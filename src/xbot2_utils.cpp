#include "xbot2_utils.hpp"

using namespace Xbot2Utils;

//************* IqRosGetter *************//

IqRosGetter::IqRosGetter(double dt,
                         bool verbose,
                         double mov_avrg_cutoff_freq)
    :_dt{dt}, _verbose{verbose}, _mov_avrg_cutoff_freq{mov_avrg_cutoff_freq}
{

}

void IqRosGetter::init_vars()
{
    _n_active_jnts = _jnt_names.size();

    _iq_out_fb = Eigen::VectorXd::Zero(_n_active_jnts);
    _timestamps = Eigen::VectorXd::Zero(_n_active_jnts);

    _mov_avrg_filter = MovAvrgFilt(_n_active_jnts, _dt, _mov_avrg_cutoff_freq);

}

void IqRosGetter::get_last_iq_out(Eigen::VectorXd& iq_out_fb)
{
    if(_vars_were_initialized)
    { // only assign output if _iq_out_fb has a nonzero dimension
      // (guaranteed after init_vars is called)
        iq_out_fb = _iq_out_fb;
    }
}

void IqRosGetter::get_last_iq_out_filt(Eigen::VectorXd& iq_out_fb_filt)
{
    if(_vars_were_initialized)
    { // only assign output if _iq_out_fb has a nonzero dimension
      // (guaranteed after init_vars is called)

        _mov_avrg_filter.add_sample(_iq_out_fb);
        _mov_avrg_filter.get(iq_out_fb_filt);
    }
}

void IqRosGetter::get_last_iq_out_stamps(Eigen::VectorXd& timestamps)
{
    timestamps = _timestamps;
}

void IqRosGetter::get_time_reference(double& t_ref)
{
    t_ref = _time_ref;
}

void IqRosGetter::set_jnt_names(std::vector<std::string> jnt_names)
{

    _jnt_names = jnt_names;

    _set_jnt_names_from_ros = false;

    _were_jnt_names_set = true;

}

void IqRosGetter::get_jnt_names(std::vector<std::string>& jnt_names)
{

    jnt_names = _jnt_names;

}

bool IqRosGetter::is_iq_out_topic_active()
{

    if(!_is_first_aux_sig)
    { // we have received at least one aux signal --> which means that at this point we know for sure
      // the joint names associated with each message

        return true;
    }
    else
    {
        return false;
    }

}

void IqRosGetter::on_js_signal_received(const xbot_msgs::JointState& js_sig)
{

    if(_set_jnt_names_from_ros && !_were_jnt_names_set)
    {
        _jnt_names = js_sig.name;

        _were_jnt_names_set = true; // to signal that we now have read the joint names
    }

    if (_verbose)
    {
        fprintf( stderr, "\n js message received \n");
    }

}

void IqRosGetter::on_aux_signal_received(const xbot_msgs::CustomState& aux_sig)
{

    if (_verbose)
    {
        fprintf( stderr, "\n aux message received \n");
    }

    if(_is_first_aux_sig)
    {
        init_vars();

        _vars_were_initialized = true;
    }

    auto remapped_aux_tuple = aux_mapper(aux_sig); // remapping aux types

    std::vector<double> ordered_vals = std::get<1>(remapped_aux_tuple); // ordered as _jnt_names
    double* ptr = &ordered_vals[0];

    _iq_out_fb = Eigen::Map<Eigen::VectorXd>(ptr, ordered_vals.size());

}

template <typename T, typename t_v >
int IqRosGetter::find_index(std::vector<T> input_v, t_v value)
{
    /**
    Finds the index of a value in an array.

    @param input_v Input vector.
    @param value Value to be searched within the input vector.

    @return  The index of the element (-1 if not present).
    */

    auto it = find(input_v.begin(), input_v.end(), value);

    // If element was found
    if (it != input_v.end())
    {
        // calculating the index
        int index = it - input_v.begin();
        return index;
    }
    else
    {
        // The element is not present in the vector
        return -1;
    }
}

template <typename T>
std::vector<int> IqRosGetter::map_indices(std::vector<T> input_v1, std::vector<T> input_v2)
{
    /**
    Maps the elements of the first input array to the second.

    @param input_v1 First input vector.
    @param input_v2 Second input vector.

    @return  A vector of indices. indices_map[i] contains the index where the i-th element of input_v1 is present in input_v2.

    */

    int v1_size = input_v1.size(); //initialize output
    std::vector<int> indices_map(v1_size, -1);

    // Here insert a check (size of input_v1 must equal size of input_v2)
    for (int i = 0; i < input_v1.size(); i++)
    {
        indices_map[i] = find_index(input_v2, input_v1[i] );
    }
    return indices_map;
}

int IqRosGetter::get_aux_type_code(std::string msg_type)
{
    /**
    Computes the code ID associated with a given message type and saves this code to the .mat file,
    so that other programs can interpret the signal-type information. To avoid the use of strings, a unique _aux_code suffix is employed.

    @param msg_type The message type name (std::string).
    @return The message code to be associated with the message type name
    */

    int msg_code;

    if (_aux_msg_type_map.find(msg_type) == _aux_msg_type_map.end()) // msg_type not present
    {
        _aux_types_encode_number++; // increment last signal code number by 1
        _aux_msg_type_map.insert({msg_type, _aux_types_encode_number}); // put the code in the dictionary
    }

    msg_code = _aux_msg_type_map.at(msg_type); // read the key associated with the msg_type from the aux code map

    return msg_code;
}

std::tuple<std::vector<int>, std::vector<double>> IqRosGetter::aux_mapper(const xbot_msgs::CustomState& aux_sig)
{

    /**
    Loops through the chains and their joints and, based on the received message, assigns the IDs associated with each message type.

    @param msg Input message
    @return A vector of signal IDs with dimension (number of chains)*(number of joints)
    */

    int n_jnts = aux_sig.name.size(); // number of joints

    std::vector<double> msg_value_remapped(n_jnts, -1.0); // output vector for msg values
    std::vector<int> msg_type_remapped(n_jnts, -1.0); // output vector for msg types

    if (_is_first_aux_sig)
    { // this runs only the first time an aux message is received
      // (joint mapping is assumed to be constant throughout the life of this
      // object)

        _indices = map_indices(_jnt_names, aux_sig.name);

        _is_first_aux_sig = false; // mapping not needed anymore

        _time_ref = aux_sig.header.stamp.toSec(); // this time will be used as a reference
    }

    for (int i = 0; i < n_jnts; i++) // mapping
    {
        if (aux_sig.type[i].find(_iq_sig_basename) != std::string::npos)
        { // we are reading the iq aux type

            int encoded_type = get_aux_type_code(aux_sig.type[i]);
            msg_type_remapped[_indices[i]] = encoded_type;
            msg_value_remapped[_indices[i]] = aux_sig.value[i];

            _timestamps[_indices[i]] = aux_sig.header.stamp.toSec() - _time_ref;// getting timestamp

        }

    }

    return make_tuple(msg_type_remapped, msg_value_remapped);
}
