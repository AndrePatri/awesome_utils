#include "awesome_utils/contact_est_utils.hpp"

#include <math.h>

using namespace ContactEstUtils;

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt,
                                     std::vector<std::string> contact_framenames,
                                     double bandwidth,
                                     double lambda, bool regularize_f, std::vector<int> selector)
    : _model_ptr{model_ptr}, _dt{data_dt},
      _contact_framenames{contact_framenames},
      _bandwidth{bandwidth},
      _regularize_f{regularize_f}, _selector{selector}
{
    if(!_model_ptr->was_model_init_ok())
    {
        std::string exception = std::string("ContactEstUtils::MomentumBasedFObs::MomentumBasedFObs(): \n") +
                                std::string("The provided model object is not initialized properly!!\n");

        throw std::invalid_argument(exception);
    }

    for (int i = 0; i < _lambda.size(); i++)
    {
        _lambda(i) = lambda; // in case a scalar lambda is provided,
        // we assign the same value to all
    }

    process_selector();

    process_contactnames();

    _nv = _model_ptr->get_nv();

    compute_bandwidth();

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

    _v = VectorXd::Zero(_nv);
    _tau = VectorXd::Zero(_nv);
    _g = VectorXd::Zero(_nv);
    _p = VectorXd::Zero(_nv);
    _to_be_integrated = VectorXd::Zero(_nv);
    _integral = VectorXd::Zero(_nv);
    _C = MatrixXd::Zero(_nv, _nv);

    _A = MatrixXd::Zero(_nv + _lambda.size() * _nc, _lambda.size() * _nc);
    _b = VectorXd::Zero(_nv + _lambda.size());

    _A_lambda = MatrixXd::Zero(_lambda.size() * _nc, _lambda.size() * _nc);
    _b_lambda = VectorXd::Zero(_A_lambda.rows());
    for (int i = 0; i < _nc; i++)
    { // we assign _lambda on the diagonal of A_lambda (replicating it
      // for the number of contacts)
        for (int j = 0; j < _lambda.size(); j++)
        {
            _A_lambda(_lambda.size() * i + j,
                      _lambda.size() * i + j) = std::sqrt(_lambda(i));
        }

    }

    // A regularization block can be assigned once and for all
    _A.block(_nv, 0, _A_lambda.rows(), _A_lambda.cols()) = _A_lambda; // adding regularization

    _W = VectorXd::Zero(_lambda.size() * _nc);
    _W_reg = VectorXd::Zero(_lambda.size() * _nc);

    _J_c_tot = MatrixXd::Zero(_lambda.size() * _nc, _nv);

}

MomentumBasedFObs::MomentumBasedFObs(Model::Ptr model_ptr, double data_dt,
                  std::vector<std::string> contact_framenames,
                  double bandwidth,
                  Reg6D lambda, bool regularize_f,
                  std::vector<int> selector)
    : _model_ptr{model_ptr}, _dt{data_dt},
      _contact_framenames{contact_framenames},
      _bandwidth{bandwidth},
      _lambda{lambda},
      _regularize_f{regularize_f}, _selector{selector}
{

    process_selector();

    process_contactnames();

    _nv = _model_ptr->get_nv();

    compute_bandwidth();

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

    _A = MatrixXd::Zero(_nv + _lambda.size() * _nc, _lambda.size() * _nc);
    _b = VectorXd::Zero(_nv + _lambda.size());

    _A_lambda = MatrixXd::Zero(_lambda.size() * _nc, _lambda.size() * _nc);
    _b_lambda = VectorXd::Zero(_A_lambda.rows());
    for (int i = 0; i < _nc; i++)
    { // we assign _lambda on the diagonal of A_lambda (replicating it
      // for the number of contacts)
        for (int j = 0; j < _lambda.size(); j++)
        {
            _A_lambda(_lambda.size() * i + j,
                      _lambda.size() * i + j) = std::sqrt(_lambda(j));
        }

    }

    // A regularization block can be assigned once and for all
    _A.block(_nv, 0, _A_lambda.rows(), _A_lambda.cols()) = _A_lambda; // adding regularization

    _W = VectorXd::Zero(_lambda.size() * _nc);
    _W_reg = VectorXd::Zero(_lambda.size() * _nc);

    _J_c_tot = MatrixXd::Zero(_lambda.size() * _nc, _nv);

}

void MomentumBasedFObs::compute_bandwidth()
{
    _k = 2 * _PI * abs(_bandwidth);

}

void MomentumBasedFObs::get_contact_framenames(std::vector<std::string> names)
{
    names = _valid_contact_names;
}

void MomentumBasedFObs::get_contact_indeces(std::vector<int> indeces)
{
    indeces = _contact_indeces;
}

void MomentumBasedFObs::process_contactnames()
{
    int framename_code = 0; // contact indeces are assigned incrementally
    for (int i = 0; i < _contact_framenames.size(); i++)
    {

        if(_model_ptr->frame_exists(_contact_framenames[i]))
        {   // it's valid, we insert it

            _contact_map[_contact_framenames[i]] = framename_code;
            _valid_contact_names.push_back(_contact_framenames[i]);
            _contact_indeces.push_back(framename_code);
            _active_contacts.push_back(true); // by default, all contacts
            // are active
            framename_code = framename_code + 1;
        }
    }

    _nc = _valid_contact_names.size();

    if (_nc == 0)
    {
        std::string exception = std::string("MomentumBasedFObs::process_contactnames(): no valid frame names where provided \"");

        throw std::invalid_argument(exception);
    }

}

void MomentumBasedFObs::process_selector()
{
    // let's make sure selector is valid
    _selector.erase(std::remove_if(
        _selector.begin(), _selector.end(),
        [](const int& x) {
            return x > 5;
        }), _selector.end()); // removing elements above 5 (a wrench can have a max. of 6 components)
    _selector.erase( unique( _selector.begin(), _selector.end() ), _selector.end() ); // erasing duplicates
    std::sort(_selector.begin(), _selector.end()); // sorting selector in growing order

    // _deselector is used to make 0 the columns of J^T corresponding to
    // components of the contact wrench which are not to be stimated.
    // For the same reason, also the corresponding elements of w_c_reg
    // are set to 0 (we want the components of the wrench which are not to be
    // estimated to be 0).
    for (auto i = _selector.begin(); i != _selector.end(); ++i)
    { // we go through each element of the provided selector and
      // remove the indeces for the Jt_selector

        int index = i - _selector.begin();

        auto it = std::find(_deselector.begin(), _deselector.end(), _selector[index]);

        if ( it != _deselector.end()) { // element found

            _deselector.erase(it); // removing it

        }
    }

}

void MomentumBasedFObs::set_active_contacts(std::vector<int> contact_indeces)
{

    std::vector<int> filtered_indeces = contact_indeces;

    filtered_indeces.erase(std::remove_if(
        filtered_indeces.begin(), filtered_indeces.end(),
        [this](const int& x) {
            return x > _contact_indeces.back();
        }), filtered_indeces.end()); // removing elements above 5 (a wrench can have a max. of 6 components)
    filtered_indeces.erase( unique( filtered_indeces.begin(), filtered_indeces.end() ), filtered_indeces.end() ); // erasing duplicates
    std::sort(filtered_indeces.begin(), filtered_indeces.end()); // sorting selector in growing order

    for (int i = 0; i < filtered_indeces.size(); i++)
    {
        _active_contacts[filtered_indeces[i]] = true; // activating contact
    }

}

void MomentumBasedFObs::get_contact_jacobians()
{
    for (int i = 0; i < _nc; i++)
    {

        if(_active_contacts[i])
        { // i-th contact is active -> we get the jacobian
          // and then set to zero the columns specified by the _selector

            _model_ptr->get_jac(_valid_contact_names[i],
                               _J_buffer);

            apply_component_selector(_J_buffer); // will set to zero the rows corresponding
            // to the components of the wrench which are NOT to be estimated

        }
        else
        { // the contact is not active
            _J_buffer = MatrixXd::Zero(_lambda.size(), _nv);
        }

        // assign to the total jacobian of the contact
        _J_c_tot.block(i * _lambda.size(), 0, _lambda.size(), _nv) = _J_buffer;
    }
}

void MomentumBasedFObs::assign_regression_matrices()
{
    _A.block(0, 0, _nv, _A.cols()) = _J_c_tot.transpose();
    _b.segment(0, _nv) = _tau_c_k;
    _b_lambda = _A_lambda * _W_reg; // the regularization is done around the regularization vector
    _b.segment(_nv, _A_lambda.rows()) = _b_lambda;
}

void MomentumBasedFObs::update()
{

    compute_tau_c(); // computing observed residual joint efforts

    // retrieve the contact jacobian at the prescribed link (from v to vel. of link wrt world frame)
    get_contact_jacobians();

    assign_regression_matrices(); // we assign matrix A and vector b for the least-squares regression problem A W = b

    // exploiting Eigen builtin method for regression problems
    _W = _A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(_b);

    if(_regularize_f)
    {
        _W_reg = _W; // will use previous solution to regularize the new solution
        // instead of using always a constant value
        for (int i = 0; i < _nc; i++)
        {
            _w_buff = _W_reg.segment(i * _lambda.size(), _lambda.size());

            apply_component_selector(_w_buff); // will set to 0 the elements corresponding to the axes we don't need
            // to estimate --> this, in conjunction to apply_selector(J_c), will ensure that
            // unestimated components will always converge to 0 thanks to the regularization term
            // (this is not needed if !_regularize_f, since the default value for _W_reg is a vect. of 0s)

            _W_reg.segment(i * _lambda.size(), _lambda.size()) = _w_buff; // re-assign filtered segment
        }

    }
}

void MomentumBasedFObs::apply_component_selector(Model::Wrench& vector)
{

    for (int i = 0; i < _deselector.size(); i++)
    {
        vector(_deselector[i]) = 0.0;
    }

}

void MomentumBasedFObs::apply_component_selector(Model::SpatialJac& J)
{

    for (int i = 0; i < _deselector.size(); i++)
    {
        J.block(_deselector[i], 0,  1, J.cols()) = MatrixXd::Zero(1, J.cols()); // set row to zero
    }

}

void MomentumBasedFObs::get_tau_obs(VectorXd& tau_c)
{
    tau_c = _tau_c_k;
}

void MomentumBasedFObs::get_w(Eigen::VectorXd& w)
{
    w = _W;
}

void MomentumBasedFObs::get_w_est_at(std::string contact_framename, Model::Wrench& w_c)
{
    w_c = _W.segment(_contact_map[contact_framename] * _lambda.size(), _lambda.size());
}

void MomentumBasedFObs::get_w_est_at(int contact_index, Model::Wrench& w_c)
{
    w_c = _W.segment(contact_index * _lambda.size(), _lambda.size());
}

void MomentumBasedFObs::get_f_est_at(std::string contact_framename, Model::Force3D& f_c)
{
    f_c = _W.segment(_contact_map[contact_framename] * _lambda.size(), 3);
}

void MomentumBasedFObs::get_f_est_at(int contact_index, Model::Force3D& f_c)
{
    f_c = _W.segment(contact_index * _lambda.size(), 3);
}

void MomentumBasedFObs::get_t_est_at(std::string contact_framename, Model::Torque3D& t_c)
{
    t_c = _W.segment(3 + _contact_map[contact_framename] * _lambda.size(), 3);
}

void MomentumBasedFObs::get_t_est_at(int contact_index, Model::Torque3D& t_c)
{
    t_c = _W.segment(3 + contact_index * _lambda.size(), 3);
}

void MomentumBasedFObs::compute_tau_c()
{

    _model_ptr->get_C(_C);
    _model_ptr->get_v(_v);
    _model_ptr->get_g(_g);
    _model_ptr->get_p(_p);
    _model_ptr->get_tau(_tau); // gets the latest model tau the user has set (e.g. from
        // a measurement)

    // computing equation (6)  --> see header for more info on this
    _to_be_integrated = _g - _C.transpose() * _v - _tau;
    _integrator.add_sample(_to_be_integrated);
    _integrator.get(_integral);

    _tau_c_k = _Skp1_inv * ( _Sk * _tau_c_k +
                             _K * ((_p - _p_km1) + _integral) ); // update current estimate

    _p_km1 = _p; // assigning joint-space momentum for the next update call

}
