#include "awesome_utils/model_interface.hpp"

using namespace ModelInterface;

Model::Model()
{

}

Model::Model(std::string urdf_path)
{
    try{
        pinocchio::urdf::buildModel(urdf_path, _pin_model);
    }
    catch (int err_code){

        _pin_model_init_ok = false;
    }

    _pin_model_init_ok = true;

    _pin_data = pinocchio::Data(_pin_model);

    _q = pinocchio::randomConfiguration(_pin_model);

    _n_jnts = (pinocchio::JointIndex)_pin_model.njoints;
    _jnt_names = std::vector<std::string>(_n_jnts);

    _nq = _q.size();
    _nv = _q_dot.size();

    for(pinocchio::JointIndex joint_id = 0; joint_id < _n_jnts; ++joint_id)
    {

        _jnt_names[joint_id] = _pin_model.names[joint_id];

    }
}

bool Model::was_model_init_ok()
{
    return _pin_model_init_ok;
}
