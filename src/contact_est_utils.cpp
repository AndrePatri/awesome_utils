#include "awesome_utils/contact_est_utils.hpp"

#include <math.h>

using namespace ContactEstUtils;

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt)
    : _model_ptr{model_ptr}, _dt{data_dt}
{

}

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt, double bandwidth)
    : _model_ptr{model_ptr}, _dt{data_dt}, _bandwidth{bandwidth}
{
    if(!_model_ptr->was_model_init_ok())
    {
        std::string exception = std::string("ContactEstUtils::MomentumBasedFObs::MomentumBasedFObs(): \n") +
                                std::string("The provided model object is not initialized properly!!\n");

        throw std::invalid_argument(exception);
    }

    _nv = _model_ptr->get_nv();

    _k = - _bandwidth * log(1 - _BW_red_factor);

    _K = MatrixXd::Identity(_nv, _nv);

    // numerical intergrators

    _integrator_tau = NumInt(_nv, _dt, _dt); // integrates only in a dt
    _integrator_g = NumInt(_nv, _dt, _dt);
    _integrator_C_T_q_dot = NumInt(_nv, _dt, _dt);

}
