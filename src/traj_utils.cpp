#include <math.h> 

#include <iostream>
#include <fstream>

#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "awesome_utils/traj_utils.hpp"

using namespace TrajUtils;

//***************************** PeisekahTrans *****************************//

PeisekahTrans::PeisekahTrans(){};

PeisekahTrans::PeisekahTrans(Eigen::VectorXd start_point, Eigen::VectorXd end_point, double exec_time, double dt)
: _start_point{start_point}, _end_point{end_point}, _exec_time{exec_time}, _dt{dt}
{

    check_input_dim();

    compute_traj();

}

Eigen::VectorXd PeisekahTrans::eval_at(int node)
{

    return _traj.col(node);

}

double PeisekahTrans::get_exec_time()
{
    return _exec_time;
}

double PeisekahTrans::get_traj_dt()
{
    return _dt;
}

double PeisekahTrans::get_n_nodes()
{
    return _n_nodes;
}

double PeisekahTrans::get_n_dim()
{
    return _n_dim;
}

void PeisekahTrans::check_input_dim()
{   
    int start_size = _start_point.size(); 
    int end_size = _end_point.size();
    if ( start_size != end_size)
    { 

        std::string exception = std::string("The starting point and end point dimensions do not match: ") + 
                                std::to_string(start_size) + std::string(" VS ") +
                                std::to_string(end_size) + std::string("\n");
                                
        throw std::invalid_argument(exception);
        
    }

    _n_dim = start_size;
}

void PeisekahTrans::rate_adapter()
{
    double dt_des = _dt;

    int n_int_aux = floor(_exec_time / dt_des);

    double dt1 = _exec_time / n_int_aux;
    double dt2 = _exec_time / (n_int_aux + 1);
    double abs_diff1 = abs(_exec_time - dt_des * (n_int_aux));
    double abs_diff2 = abs(_exec_time - dt_des * (n_int_aux + 1));

    if (abs_diff1 < abs_diff2)
    { 
        _n_nodes = n_int_aux + 1;
        _dt = dt1;
    }
    else
    {
        _n_nodes = n_int_aux + 2;  
        _dt = dt2;
    }

}

double PeisekahTrans::compute_peisekah_val(double phase, double start_point, double end_point)
{   

    double common_part_traj = (126.0 * pow(phase, 5) - 420.0 * pow(phase, 6) + 
                            540.0 * pow(phase, 7) - 315.0 * pow(phase, 8) + 
                            70.0 * pow(phase, 9));

    auto value = start_point + (end_point - start_point) *  common_part_traj;

    return value;

}
Eigen::VectorXd PeisekahTrans::compute_peisekah_vect_val(double phase, Eigen::MatrixXd start_point,  Eigen::MatrixXd end_point)
{   
    int n_dim_start =  start_point.rows() >= start_point.cols() ? start_point.rows(): start_point.cols();
    bool column_wise_strt =  start_point.rows() >= start_point.cols() ? true: false;

    int n_dim_trgt =  end_point.rows() >= end_point.cols() ? end_point.rows(): end_point.cols();
    bool column_wise_trgt =  end_point.rows() >= end_point.cols() ? true: false;

    if (n_dim_start != n_dim_trgt)
    {
        std::string exception = std::string("compute_peisekah_vect_val: dimension mismatch in the provided points: ") + 
                                std::string("start point dim: ") + std::to_string(n_dim_start) + std::string(", ") +
                                std::string("final point dim: ") + std::to_string(n_dim_trgt) + std::string(".") ;

        throw std::invalid_argument(exception);
    }

    if (!(column_wise_strt && column_wise_trgt))
    {
        std::string exception = std::string("compute_peisekah_vect_val: input points have to be both column-major or row-major!");

        throw std::invalid_argument(exception);
    }



    Eigen::VectorXd peisekah_sample(n_dim_start);

    for (int k = 0; k < n_dim_start; k++)
    { 
        if (column_wise_strt)
        {
            peisekah_sample(k) = compute_peisekah_val(phase, start_point(k, 0), end_point(k, 0));
        }
        else
        {
            peisekah_sample(k) = compute_peisekah_val(phase, start_point(0, k), end_point(0, k));
        }
        
    }

    return peisekah_sample;
}

void PeisekahTrans::compute_traj()
{
    rate_adapter(); // adapt _dt and _n_nodes to match the required _exec_time

    _traj = Eigen::MatrixXd::Zero(_n_dim, _n_nodes);

    for (int k = 0; k < _n_dim; k++)
    { // loop through joints (rows)

        for (int i = 0; i < _n_nodes; i++)
        { // loop through samples (columns)

            double phase = (double) i/(_n_nodes - 1);
            
            _traj(k, i) = compute_peisekah_val(phase, _start_point(k), _end_point(k)); 

        }

    }

}

Eigen::MatrixXd PeisekahTrans::get_traj()
{

    return _traj;
    
}

//***************************** TrajLinInterp *****************************//

TrajLinInterp::TrajLinInterp(){};

TrajLinInterp::TrajLinInterp(Eigen::VectorXd sample_time, Eigen::MatrixXd input_traj, int interp_dir, double time_check_tol)
: _sample_times{sample_time}, _traj{input_traj}, _interp_dir{interp_dir}, _time_check_tol{time_check_tol}
{
    
    check_dim_match(); // at this point, it is guaranteed that dimensions along the interpolation axis
    // (which now is for sure valid) match

    _n_dims = (_interp_dir == 1) ? _traj.rows(): _traj.cols();

}

Eigen::MatrixXd TrajLinInterp::eval_at(Eigen::VectorXd interp_times)
{

    // check interp_times is within _sample_times
    check_time_vector(interp_times);

    int n_int_samples = interp_times.size();
    Eigen::MatrixXd inter_mat = Eigen::MatrixXd::Zero(_n_dims, n_int_samples);

    for (int i = 0; i < _n_dims; i++)
    {

        if (_interp_dir == 1)
        { //  row-wise

            inter_mat.row(i) = interp_1d(interp_times, i);

        }
        else
        { // col-wise

            inter_mat.col(i) = interp_1d(interp_times, i);

        }

    }

    return inter_mat;

}

    
void TrajLinInterp::check_time_vector(Eigen::VectorXd interp_times)
{
    int n_int_samples = interp_times.size();

    if ( (interp_times(0) <  _sample_times(0) - _time_check_tol) ||
        (interp_times(n_int_samples - 1) > _sample_times(_n_samples - 1) + _time_check_tol ) )
    { // checking interp_times is actually within _sample_times (with a tolerance given by _time_check_tol)  

        std::string exception = std::string("check_time_vector: The provided interpolation array ") + 
                                std::string("[") + std::to_string(interp_times(0)) + std::string(", ") +
                                std::to_string(interp_times(n_int_samples - 1)) + std::string("]") +
                                std::string(" is outside the trajectory time array ") +
                                std::string("[") + std::to_string(_sample_times(0)) + std::string(", ") +
                                std::to_string(_sample_times(_n_samples - 1)) + std::string("]\n");

        throw std::invalid_argument(exception);
        
    }

}

int TrajLinInterp::get_traj_n_samples()
{

    if (_interp_dir == 0)
    { // colum-wise
        int sample_number = _traj.rows();

        return sample_number; 
    }
    else
    {
        if (_interp_dir == 1)
        { // row-wise
            int sample_number = _traj.cols();
        
            return sample_number; 
        }
        else
        {
            std::string exception = "get_traj_n_samples: Invalid interpolation direction " +
                                    std::to_string(_interp_dir) + "provided. " +
                                    "Allowed values are 0(row-wise) and 1(colum-wise)\n";

            _interp_dir = -1;

            throw std::invalid_argument(exception);

        }
    }                

}

void TrajLinInterp::check_dim_match()
{
    _n_samples = get_traj_n_samples();

    if ( !(_n_samples == _sample_times.rows()) || (_n_samples == _sample_times.cols()) )
    { // checking that time vector and trajectory dimensions match 

        std::string exception = std::string("check_dim_match: Time array and trajectory dimension do not match. ") + 
                                std::string("They are, respectively,") + 
                                std::string(" (") + std::to_string(_sample_times.rows()) + std::string(", ") + 
                                std::to_string(_sample_times.cols()) + std::string(")") +
                                std::string("and") + 
                                std::string(" (") + std::to_string(_traj.rows()) + ", " + std::to_string(_traj.cols()) +
                                std::string(")\n");

        throw std::invalid_argument(exception);

    }
}
            
Eigen::VectorXd TrajLinInterp::interp_1d(Eigen::VectorXd& interp_times, int dim_index)
{   
    int n_int_samples = interp_times.size();

    Eigen::VectorXd interp_vect = Eigen::VectorXd::Zero(n_int_samples);
    
    for (int i = 0; i < n_int_samples; i++)
    {
        int first_indx = -1;
        int second_indx = -1;

        get_closest_points(interp_times(i), first_indx, second_indx);

        double interval_dt = _sample_times(second_indx) - _sample_times(first_indx);

        double t_norm = interp_times(i) - _sample_times(first_indx);

        interp_vect(i) = interp_0d(t_norm, interval_dt, first_indx, second_indx, dim_index);
    }

    return interp_vect;
}

double TrajLinInterp::interp_0d(double t_norm, double interval_dt, int first_indx, int second_indx, int dim_index)
{   
    double interp_val;

    if (_interp_dir == 1)
    { // row-wise

        interp_val = _traj(dim_index, first_indx) + t_norm / interval_dt * 
                    (_traj(dim_index, second_indx) - _traj(dim_index, first_indx));

    } // col-wise
    else{

        interp_val = _traj(first_indx, dim_index) + t_norm / interval_dt * 
                    (_traj(second_indx, dim_index) - _traj(first_indx, dim_index));

    }
    
    return interp_val;
}

void TrajLinInterp::get_closest_points(double inter_time, int& first_indx, int& second_indx)
{
    
    Eigen::VectorXd diff_array = Eigen::VectorXd::Zero(_n_samples); 

    int closest_sampl_idx = -1;
    double previous_pos_diff = 1000000000000000000.0; // auxiliary value, initialized to very high value

    for (int i = 0; i < _n_samples; i++)
    {
        diff_array(i) = inter_time - _sample_times(i); 

        if (diff_array(i) > 0 && diff_array(i) < previous_pos_diff)
        {
            previous_pos_diff = diff_array(i);

            closest_sampl_idx = i;
        }
    }

    if (closest_sampl_idx == -1) // all negative differences --> we are before the first time sample
    {
        first_indx = 0;
        second_indx = 1;
    }

    if (closest_sampl_idx == (_n_samples - 1 )) // we are past the last time sample
    {
        first_indx = _n_samples - 2;
        second_indx = _n_samples - 1;
    }

    if(closest_sampl_idx != -1 && closest_sampl_idx != (_n_samples - 1 ))
    {
        first_indx = closest_sampl_idx;
        second_indx = closest_sampl_idx + 1;
    }

}


//***************************** TrajLoader *****************************//

TrajLoader::TrajLoader(){};

TrajLoader::TrajLoader(std::string data_path, bool column_major, double resample_err_tol, bool load_from_csv)
:_data_path{data_path}, _column_major_order{column_major}, _resample_err_tol{resample_err_tol}, _load_from_csv{load_from_csv}
{

    if (!_load_from_csv)
    {
        load_data_from_mat(data_path);
    }
    else
    {
        throw std::invalid_argument(std::string("Reading data from CSV files is not supported anymore. Use .mat files instead!\n"));
    }
    
    check_loaded_data_dims();

    _n_nodes = get_n_samples(_q_p);
    _n_jnts = get_n_jnts(_q_p);

    _sample_times = Eigen::VectorXd::Zero(_n_nodes);
    for (int i = 0; i < (_n_nodes - 1); i++)
    {
        _sample_times(i + 1) = _sample_times(i) + _dt_opt(0, i);
    }
    _exec_time = _sample_times(_n_nodes - 1) - _sample_times(0);

    int interp_dir = (_column_major_order) ? 1 : 0;

    opt_traj.emplace(_q_p_name, TrajLinInterp(_sample_times, _q_p, interp_dir));
    opt_traj.emplace(_q_p_dot_name, TrajLinInterp(_sample_times, _q_p_dot, interp_dir));
    opt_traj.emplace(_efforts_name,
                    TrajLinInterp(_sample_times, _tau, interp_dir));   
    opt_traj.emplace(_f_cont_name,
                    TrajLinInterp(_sample_times, _f_cont, interp_dir));  
    
}

Eigen::MatrixXd TrajLoader::read_data_from_csv(std::string data_path)
{ // keep it public so it can also be used outside class instances
    using namespace Eigen;
    using namespace std;

    vector<double> matrixEntries;

    // in this object we store the data from the matrix
    ifstream matrixDataFile(data_path);

    // this variable is used to store the row of the matrix that contains commas 
    string matrixRowString;

    // this variable is used to store the matrix entry;
    string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;


    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and
    // store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

        while (getline(matrixRowStringStream, matrixEntry, ' ')) // here we read pieces of the stream matrixRowStringStream
        // until every comma, and store the resulting character into the matrixEntry
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector
            // storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }

    // here we convet the vector variable into the matrix and return the resulting object, 
    // note that matrixEntries.data() is the pointer to the first memory location at which the entries of
    // the vector matrixEntries are stored;
    return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>
            (matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}

int TrajLoader::get_n_jnts()
{
    return _n_jnts;
}

int TrajLoader::get_n_nodes()
{
    return _n_nodes;
}

Eigen::VectorXd TrajLoader::get_sample_times()
{
    return _sample_times;
}

void TrajLoader::get_loaded_traj(Eigen::MatrixXd& q_p, Eigen::MatrixXd& q_p_dot, Eigen::MatrixXd& tau,
                                Eigen::MatrixXd& dt_opt, Eigen::MatrixXd& f_cont)
{
    q_p = _q_p;
    q_p_dot = _q_p_dot;
    tau = _tau;
    dt_opt = _dt_opt;
    f_cont = _f_cont;
}

void TrajLoader::resample(double res_dt, Eigen::MatrixXd& q_p_res, Eigen::MatrixXd& q_p_dot_res,
                        Eigen::MatrixXd& tau_res, Eigen::MatrixXd& f_cont_res)
{

    Eigen::VectorXd times = compute_res_times(res_dt);

    double n_res_nodes = times.size();

    q_p_res =  opt_traj[_q_p_name].eval_at(times);

    q_p_dot_res = opt_traj[_q_p_dot_name].eval_at(times);

    tau_res =  opt_traj[_efforts_name].eval_at(times.head(n_res_nodes - 1)); // tau is resampled excluding the last instant of time

    tau_res.conservativeResize(tau_res.rows(), tau_res.cols() + 1);
    tau_res.col(tau_res.cols() - 1) = Eigen::VectorXd::Zero(_n_jnts); // to be able to potentially send the whole trajectory concurrently
    // // a dummy null control input is added on the last sample time

    f_cont_res =  opt_traj[_f_cont_name].eval_at(times.head(n_res_nodes - 1)); // tau is resampled excluding the last instant of time

    f_cont_res.conservativeResize(f_cont_res.rows(), f_cont_res.cols() + 1);
    f_cont_res.col(f_cont_res.cols() - 1) = Eigen::VectorXd::Zero(_n_jnts); // to be able to potentially send the whole trajectory concurrently
    // // a dummy null control input is added on the last sample time

    // updating number of nodes
    _n_nodes = get_n_samples(q_p_res);

}

std::string TrajLoader::get_file_extension(std::string file)
{
    std::vector<std::string> token_list;
    boost::split(token_list, file, [](char c){return c == '.';});
    
    if(token_list.size() > 1)
    {
        return token_list.back();
    }
    
    return "";
}

int TrajLoader::get_n_jnts(Eigen::MatrixXd& mat)
{
    if(_column_major_order)
    {
        return mat.rows();
    }
    else
    {
        return mat.cols();
    }
}

int TrajLoader::get_n_samples(Eigen::MatrixXd& mat)
{
    if(_column_major_order)
    {
        return mat.cols();
    }
    else
    {
        return mat.rows();
    }
}

int TrajLoader::get_takeoff_index(double epsi)
{
    int takeoff_index = -1;

    for (int n = 0; n < (_f_cont.cols()); n++)
    {
        if (_f_cont(2, n) < epsi)
        {
            takeoff_index = n;
            break;
        }
    }

    return takeoff_index;
}

double TrajLoader::get_exec_time()
{
    return _exec_time;
}

void TrajLoader::get_opt_dt(double& dt)
{
    dt = _dt_opt(0);
}

void TrajLoader::check_loaded_data_dims()
{   
    // By convention, the tau value at the last node is zero.
    // This way the q_p, q_p_dot and tau all have the same number of samples
    if ( !( ( get_n_jnts(_q_p) == get_n_jnts(_q_p_dot) ) && (get_n_jnts(_q_p_dot) == get_n_jnts(_tau)) ) )
    { // check number of joints consistency between data

        throw std::invalid_argument(std::string("check_loaded_data_dims: The number of rows ") + 
                                    std::string("(i.e. joints) of the loaded data does not match!\n"));

    }

    if ( !( (get_n_samples(_q_p) == get_n_samples(_q_p_dot)) && (get_n_samples(_q_p_dot) == (get_n_samples(_tau))) 
        && (get_n_samples(_tau) == (get_n_samples(_f_cont)))) )
    { // check cols (torque matrix)

        throw std::invalid_argument(std::string("check_loaded_data_dims: ") +
                                    std::string("The number of columns (i.e. samples) of the loaded data does not match!\n" )+
                                    std::string("q_p samples: ") + std::to_string(get_n_samples(_q_p)) + std::string("\n") + 
                                    std::string("q_p_dot samples: ") + std::to_string(get_n_samples(_q_p_dot)) + std::string("\n") +
                                    std::string("tau samples: ") + std::to_string(get_n_samples(_tau)) + std::string("\n") + 
                                    std::string("f contact samples: ") + std::to_string(get_n_samples(_f_cont)) + std::string("\n"));

    }

    if ((get_n_samples(_tau) - 1) != _dt_opt.size()) // note that size() returns the product 
    // cols * rows (ok for _dt_opt because it is a 1-D matrix)
    {
        int a = get_n_samples(_tau);

        throw std::invalid_argument(std::string("check_loaded_data_dims:") +
                                    std::string("The size of the loaded dt vector does not match the other data!.\n"));
    }

}

// void TrajLoader::load_data_from_csv(std::string data_path)
// {

//     std::string q_p_path = data_path + _q_p_name + std::string(".csv");
//     std::string q_p_dot_path = data_path + _q_p_dot_name + std::string(".csv");
//     std::string tau_path = data_path + _efforts_name + std::string(".csv");
//     std::string dt_path = data_path + _dt_name + std::string(".csv");

//     _q_p = read_data_from_csv(q_p_path);
//     _q_p_dot = read_data_from_csv(q_p_dot_path);
//     _tau = read_data_from_csv(tau_path);

//     Eigen::MatrixXd dt_opt_aux1 = read_data_from_csv(dt_path);
//     Eigen::Map<Eigen::VectorXd> dt_opt_aux2(dt_opt_aux1.data(), dt_opt_aux1.size()); // converting to vector
//     _dt_opt = dt_opt_aux2;

//     if (_q_p.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find q_p at ") + q_p_path);
//     }
//     if (_q_p_dot.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find q_p_dot at ") + q_p_dot_path);
//     }
//     if (_tau.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find tau at ") + tau_path);
//     }
//     if (_dt_opt.size() == 0)
//     { // reading failed    
//         throw std::runtime_error(std::string("load_data_from_csv: Failed to find dt_opt at ") + dt_path);
//     }

// }

void TrajLoader::load_data_from_mat(std::string math_path)
{
    
    // throw std::invalid_argument(std::string("load_data_from_mat: 
    // Reading from mat databases is not supported yet!! \n")); // to be removed upon new MatLogger2 merge

    XBot::MatLogger2::Options opts;
    opts.load_file_from_path = true; // enable reading
    auto _load_logger = XBot::MatLogger2::MakeLogger(math_path, opts);

    int slices; // not needed, used just to call the method properly 
    bool q_p_read_ok = _load_logger->readvar(_q_p_name, _q_p, slices);
    bool q_p_dot_read_ok = _load_logger->readvar(_q_p_dot_name, _q_p_dot, slices);
    bool tau_read_ok = _load_logger->readvar(_efforts_name, _tau, slices);
    bool f_cont_read_ok = _load_logger->readvar(_f_cont_name, _f_cont, slices);

    _tau.conservativeResize(_tau.rows(), _tau.cols()+1); // appending a vector of zero torques for the last sample
    // (input always null at the last trajectory node)
    _tau.col(_tau.cols() - 1) = Eigen::VectorXd::Zero(_tau.rows());

    _f_cont.conservativeResize(_f_cont.rows(), _f_cont.cols()+1); // appending a vector of zero torques for the last sample
    // (input always null at the last trajectory node)
    _f_cont.col(_f_cont.cols() - 1) = Eigen::VectorXd::Zero(_f_cont.rows());
    
    bool dt_read_ok = _load_logger->readvar(_dt_name, _dt_opt, slices); // here fix _dt_opt (should change to MatrixXd)

    if (!q_p_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find q_p from mat database at ") + math_path);
    }
    if (!q_p_dot_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find q_p_dot from mat database at ") + math_path);
    }
    if (!tau_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find tau from mat database at ") + math_path);
    }
    if (!f_cont_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find f_contact from mat database at ") + math_path);
    }
    if (!dt_read_ok)
    { // reading failed    
        throw std::runtime_error(std::string("Failed to find dt_opt from mat database at ") + math_path);
    }

    _load_logger.reset();

}

Eigen::VectorXd TrajLoader::compute_res_times(double dt_res)
{
    // in case a _exec_time / dt_res has not zero remainder, the resulting resampled trajectory will be replayed by the plugin
    // with a different execution time w.r.t. the expected one. The error is _exec_time - n_nodes * dt_plugin 
    
    int n_nodes = round(_exec_time / dt_res) + 1; // if _exec_time is exactly divisible by dt_res, 
    // round returns the right number of nodes. If not, the number which allows the smallest deviation
    // from the nominal execution time
    
    double exec_time_res_error = _exec_time - (n_nodes - 1) * dt_res;

    if (abs(exec_time_res_error) > _resample_err_tol)
    { // the resulting execution error is beyond the set threshold -> throw error

        std::string error = std::string("compute_res_times: The error on the execution time (" + std::to_string(_exec_time)) +
                            std::string(" s) resulting from resampling at \n") + 
                            std::to_string(dt_res) + std::string( "s is ") + 
                            std::to_string(exec_time_res_error) + std::string("s, which in absolute value greater") +
                            std::string("than the allowed one: i.e. ") +
                            std::to_string(_resample_err_tol) + std::string("s,\n");

        throw std::invalid_argument(error);
    }

    Eigen::VectorXd times_res = Eigen::VectorXd::Zero(n_nodes);
    for (int i = 0; i < (n_nodes - 1); i++)
    {
        times_res(i + 1) = times_res(i) + dt_res;
    }

    return times_res;
    
}
