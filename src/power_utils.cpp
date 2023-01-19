#include<power_utils.hpp>
#include<typedefs.hpp>

using namespace PowerUtils;

RegEnergy::RegEnergy(IqRosGetter::Ptr iq_meas,
                    IqEstimator::Ptr iq_est,
                    Eigen::VectorXd R,
                    Eigen::VectorXd L_leak, Eigen::VectorXd L_m,
                    double dt,
                    bool use_iq_meas)
    :_iq_meas{iq_meas},
     _iq_est{iq_est},
     _R{R},
     _L_leak{L_leak}, _L_m{L_m},
     _dt{dt},
     _use_iq_meas{use_iq_meas}
{
    _n_jnts = iq_est->get_n_jnts();
    iq_est->get_Kt(_Kt);

    // we check for dimension consistency between R, L_leak, Lm, K_t and the number of joints
    int err = 0;
    if(_R.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(_L_leak.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(_L_m.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (err != 0)
    {
        std::string exception = std::string("RegEnergy::RegEnergy(): dimension mismatch of input data with error code: ") +
                                std::to_string(err) + std::string("\n");

        throw std::invalid_argument(exception);
    }

    _R_q = _R;

    _L_q = _L_leak + 3/2 * _L_m;

    _iq_0 = Eigen::VectorXd::Zero(_n_jnts);
    _iq_k = Eigen::VectorXd::Zero(_n_jnts);

    _ek = Eigen::VectorXd::Zero(_n_jnts);
    _pk = Eigen::VectorXd::Zero(_n_jnts);

    _pk_joule = Eigen::VectorXd::Zero(_n_jnts);
    _pk_mech = Eigen::VectorXd::Zero(_n_jnts);
    _pk_indct_est = Eigen::VectorXd::Zero(_n_jnts);

    _ek_joule = Eigen::VectorXd::Zero(_n_jnts);
    _ek_mech = Eigen::VectorXd::Zero(_n_jnts);
    _ek_indct = Eigen::VectorXd::Zero(_n_jnts);

    _omega_r = Eigen::VectorXd::Zero(_n_jnts);

    _num_diff_iq = NumDiff(_n_jnts, _dt);

    _num_int_joule = NumInt(_n_jnts, _dt, _dt);
    _num_int_mech = NumInt(_n_jnts, _dt, _dt);

    _mov_filter = MovAvrgFilt(_n_jnts, _dt, _filter_cutoff_freq);

}

void RegEnergy::set_e0(double e0)
{
    if(_is_first_update)
    {

        _e0 = e0;

    }
}

void RegEnergy::set_omega_r(Eigen::VectorXd omega_r)
{

    int err = 0;
    if(omega_r.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (err != 0)
    {
        std::string exception = std::string("RegEnergy::set_omega_r(): dimension mismatch of input data!\n");

        throw std::invalid_argument(exception);
    }

    if(_use_iq_meas)
    {
        _omega_r = omega_r;

    }else{ // we do nothing

        std::cout << utils_defs::Colors::kYellow << "RegEnergy::set_omega_r(): omega_r won't be set since use_iq_meas==true" << utils_defs::Colors::kEndl << std::endl;
    }

}

void RegEnergy::update()
{

    if(_is_first_update)
    { // set initial state
        if(!_use_iq_meas)
        {

            _iq_est->get_iq(_iq_0); // gets latest iq estimate sample from the estimator

        }else{

            _iq_meas->get_last_iq_out(_iq_0); // gets latest iq estimate sample from the actual measurements

        }

        _is_first_update = false;
    }

    if(!_use_iq_meas)
    {

        _iq_est->get_iq(_iq_k);
        _iq_est->get_omega_r(_omega_r);


    }else{

        _iq_meas->get_last_iq_out(_iq_k);
    }

    compute(); // updates energy and power values

}

void RegEnergy::compute()
{
    _num_diff_iq.add_sample(_iq_k);
    _num_diff_iq.dot(_iq_dot_est);

    compute_power(); // we compute all the terms of the power balance and then we integrate them
    // to get the energy

    _num_int_joule.add_sample(_pk_joule);
    _num_int_joule.get(_ek_joule);
    _num_int_mech.add_sample(_pk_mech);
    _num_int_joule.get(_ek_mech);

    compute_energy();
}

void RegEnergy::compute_power()
{
    _pk_joule = 3.0/2.0 * _R_q.array() * (_iq_k.array().pow(2));

    _pk_mech = (_Kt.array() * _iq_k.array()) * _omega_r.array();

    _pk_indct_est = 3.0/2.0 * _L_q.array() * (_iq_k.array() * _iq_dot_est.array());
}

void RegEnergy::compute_energy()
{
    _ek_indct = 3.0/4.0 * _L_q.array() * (_iq_k.array().pow(2) - _iq_0.array().pow(2));

    _ek = _e0 + _ek_joule.array() + _ek_mech.array() + _ek_indct.array();

    _ek_tot = _ek.sum();

}


