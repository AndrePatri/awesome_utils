#include "awesome_utils/contact_est_utils.hpp"

#include <math.h>

using namespace ContactEstUtils;

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt)
    : _model_ptr{model_ptr}, _dt{data_dt}
{

}

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt, double bandwidth,
                                     double lambda, bool regularize_f, std::vector<int> selector)
    : _model_ptr{model_ptr}, _dt{data_dt}, _bandwidth{bandwidth}, _lambda{lambda},
      _regularize_f{regularize_f}, _selector{selector}
{
    if(!_model_ptr->was_model_init_ok())
    {
        std::string exception = std::string("ContactEstUtils::MomentumBasedFObs::MomentumBasedFObs(): \n") +
                                std::string("The provided model object is not initialized properly!!\n");

        throw std::invalid_argument(exception);
    }

    // let's make sure selector is valid
    _selector.erase(std::remove_if(
        _selector.begin(), _selector.end(),
        [](const int& x) {
            return x > 5;
        }), _selector.end()); // removing elements above 5 (a wrench can have a max. of 6 components)
    _selector.erase( unique( _selector.begin(), _selector.end() ), _selector.end() ); // erasing duplicates
    std::sort(_selector.begin(), _selector.end()); // sorting selector in growing order

    // _Jt_selector is used to make 0 the columns of J^T corresponding to
    // components of the contact wrench which are not to be stimated.
    // For the same reason, also the corresponding elements of w_c_reg
    // are set to 0 (we want the components of the wrench which are not to be
    // estimated to be 0).
    for (auto i = _selector.begin(); i != _selector.end(); ++i)
    { // we go through each element of the provided selector and
      // remove the indeces for the Jt_selector

        int index = i - _selector.begin();

        auto it = std::find(_Jt_selector.begin(), _Jt_selector.end(), _selector[index]);

        if ( it != _Jt_selector.end()) { // element found

            _Jt_selector.erase(it); // removing it

        }
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
    MatrixXd J_c, J_c_transp;
    compute_tau_c(); // computing observed residual joint efforts

    // retrieve the contact jacobian at the prescribed link (from v to vel. of link wrt world frame)
    _model_ptr->get_jac(contact_framename,
                        Model::ReferenceFrame::LOCAL_WORLD_ALIGNED,
                        J_c);

    J_c_transp = J_c.transpose();
    apply_selector(J_c_transp); // will set to 0 columns of the axis we don't want to
    // estimate.
    // Those components value will not influence the QP

    // solving QP for retrieving the force + wrench estimate
    // basically solving J.T * f_c = tau_c (but with some regularization)
    _A.block(0, 0, _nv, _A.cols()) = J_c_transp;
    _b.segment(0, _nv) = _tau_c_k;
    _b_lambda = _w_c_reg; // the regularization is done around the regularization vector
    _b.segment(_nv, _I_lambda.rows()) = std::sqrt(_lambda) * _b_lambda;

    // exploiting Eigen builtin method for regression problems
    _w_c = _A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(_b);

    if(_regularize_f)
    {
        _w_c_reg = _w_c; // will use previous solution to regularize the new solution
        // instead of using always a constant value
        apply_selector(_w_c_reg); // will set to 0 the elements corresponding to the axes we don't need
        // to estimate --> this, in conjunction to apply_selector(J_c_transp), will ensure that
        // unestimated components will always converge to 0
        // (this is not needed if !_regularize_f, since the default value for _w_c_reg is a vect. of 0s)
    }
}

void MomentumBasedFObs::apply_selector(VectorXd& vector)
{
    if (vector.size() != 6)
    {
        std::string exception = std::string("ContactEstUtils::MomentumBasedFObs::apply_selector(): \n") +
                                std::string("The input vector must be of dimension 6!!\n");

        throw std::invalid_argument(exception);
    }

    for (int i = 0; i < _Jt_selector.size(); i++)
    {
        vector(_Jt_selector[i]) = 0.0;
    }
}

void MomentumBasedFObs::apply_selector(MatrixXd& matrix)
{

    if (matrix.cols() != 6)
    {
        std::string exception = std::string("ContactEstUtils::MomentumBasedFObs::apply_selector(): \n") +
                                std::string("The input matrix must have exactly 6 columns!!\n");

        throw std::invalid_argument(exception);
    }

    for (int i = 0; i < _Jt_selector.size(); i++)
    {
        matrix.block(0, _Jt_selector[i], matrix.rows(), 1) = VectorXd::Zero(matrix.rows());
    }

}

void MomentumBasedFObs::get_tau_obs(VectorXd& tau_c)
{
    tau_c = _tau_c_k;
}

void MomentumBasedFObs::get_w_est(VectorXd& w_c)
{
    w_c = _w_c;
}

void MomentumBasedFObs::get_f_est(VectorXd& f_c)
{
    f_c = _w_c.segment(0, 3);
}

void MomentumBasedFObs::get_t_est(VectorXd& t_c)
{
    t_c = _w_c.segment(3, 3);
}

void MomentumBasedFObs::compute_tau_c()
{

    // getting quantities from the model
    VectorXd v, tau, g, p;
    MatrixXd C;
    VectorXd to_be_integrated, integral,
             tau_c_kp1;

    _model_ptr->get_C(C);
    _model_ptr->get_v(v);
    _model_ptr->get_g(g);
    _model_ptr->get_p(p);
    _model_ptr->get_tau(tau);

    // computing equation (6)  --> see header for more info on this
    to_be_integrated = g - C.transpose() * v - tau;
    _integrator.add_sample(to_be_integrated);
    _integrator.get(integral);

    tau_c_kp1 = _Skp1_inv * ( _Sk * _tau_c_k +
                              _K * ((p - _p_km1) + integral) );

    _tau_c_k = tau_c_kp1; // update current estimate

    _p_km1 = p; // assigning joint-space momentum for the next update call


}
