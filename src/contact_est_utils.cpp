#include "awesome_utils/contact_est_utils.hpp"

#include <math.h>

using namespace ContactEstUtils;

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt)
    : _model_ptr{model_ptr}, _dt{data_dt}
{

}

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt, double bandwidth,
                                     double lambda, bool regularize_f)
    : _model_ptr{model_ptr}, _dt{data_dt}, _bandwidth{bandwidth}, _lambda{lambda},
      _regularize_f{regularize_f}
{
    if(!_model_ptr->was_model_init_ok())
    {
        std::string exception = std::string("ContactEstUtils::MomentumBasedFObs::MomentumBasedFObs(): \n") +
                                std::string("The provided model object is not initialized properly!!\n");

        throw std::invalid_argument(exception);
    }

    _nv = _model_ptr->get_nv();

    _k = - abs(_bandwidth) * log(1 - _BW_red_factor);

    _K = _k * MatrixXd::Identity(_nv, _nv); // diagonal matrix of
    // constant values

    // state transition matrices for the integration of the observer
    // dynamic (constant, so initialzied here)
    _Skp1 = MatrixXd::Identity(_nv, _nv) + _dt/2.0 * _K;
    _Skp1_inv = _Skp1.inverse();
    _Sk = MatrixXd::Identity(_nv, _nv) - _dt/2.0 * _K;

    _integrator = NumInt(_nv, _dt, _dt); // numerical integrator

    _tau_c_k = VectorXd::Zero(_nv);

    _p_km1 = VectorXd::Zero(_nv);

    _A = MatrixXd::Zero(_nv + 6, 6);
    _b = VectorXd::Zero(_nv + 6);

    _I_lambda = MatrixXd::Identity(6, 6);
    _b_lambda = VectorXd::Zero(_I_lambda.rows());

    // A regularization block can be assigned once and for all
    _A.block(_nv, 0, _I_lambda.rows(), _I_lambda.cols()) = std::sqrt(_lambda) * _I_lambda; // adding regularization

    _w_c = VectorXd::Zero(6);
    _w_c_reg = VectorXd::Zero(6);

}

void MomentumBasedFObs::update(std::string contact_framename)
{
    VectorXd v, tau, g, p;
    MatrixXd C, J_c;
    VectorXd to_be_integrated, integral,
             tau_c_kp1;

    _model_ptr->get_C(C);
    _model_ptr->get_v(v);
    _model_ptr->get_g(g);
    _model_ptr->get_p(p);
    _model_ptr->get_tau(tau);

    to_be_integrated = tau - g + C.transpose() * v;
    _integrator.add_sample(to_be_integrated);
    _integrator.get(integral);

    tau_c_kp1 = _Skp1_inv * ( _Sk * _tau_c_k +
                              (p - _p_km1) -
                              integral);

    _tau_c_k = tau_c_kp1; // update current estimate

    // retrieve the contact jacobian at the prescribed link (from v to v link wrt world frame)
    _model_ptr->get_jac(contact_framename,
                        Model::ReferenceFrame::LOCAL_WORLD_ALIGNED,
                        J_c);

    // solving QP for retrieving the force + wrench estimate
    // basically solving J.T * f_c = tau_c (but with some regularization)
    _A.block(0, 0, _nv, _A.cols()) = J_c.transpose();
    _b.segment(0, _nv) = _tau_c_k;
    _b_lambda = _w_c_reg; // the regularization is done around the regularization vector
    _b.segment(_nv, _I_lambda.rows()) = std::sqrt(_lambda) * _b_lambda;

    // exploiting Eigen builtin method for regression problems
    _w_c = _A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(_b);

    _p_km1 = p; // assigning joint-space momentum for the next update call
    if(_regularize_f)
    {
        _w_c_reg = _w_c; // will use previous solution to regularize the new solution
        // instead of using always a constant value
    }
}

void MomentumBasedFObs::get_tau_obs(VectorXd& tau_c)
{
    tau_c = _tau_c_k;
}

void MomentumBasedFObs::get_f_est(VectorXd& f_c)
{
    f_c = _w_c.segment(0, 3);
}

void MomentumBasedFObs::get_t_est(VectorXd& t_c)
{
    t_c = _w_c.segment(3, 3);
}
