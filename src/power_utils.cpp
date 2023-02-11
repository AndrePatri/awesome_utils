#include<power_utils.hpp>

#include<typedefs.hpp>

using namespace PowerUtils;

RegEnergy::RegEnergy(IqRosGetter::Ptr iq_meas,
                    IqEstimator::Ptr iq_est,
                    Eigen::VectorXd R,
                    Eigen::VectorXd L_leak, Eigen::VectorXd L_m,
                    double dt,
                    bool use_iq_meas,
                    bool dump_data2mat,
                    std::string dump_path)
    :_iq_meas{iq_meas},
     _iq_est{iq_est},
     _R{R},
     _L_leak{L_leak}, _L_m{L_m},
     _dt{dt},
     _use_iq_meas{use_iq_meas},
     _dump_data2mat{dump_data2mat},
     _dump_path{dump_path}
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

    _aux_1p1_vect_energy =  Eigen::VectorXd::Zero(1);
    _aux_1p1_vect_power =  Eigen::VectorXd::Zero(1);

    _iq_0 = Eigen::VectorXd::Zero(_n_jnts);
    _iq_k = Eigen::VectorXd::Zero(_n_jnts);

    _ek = Eigen::VectorXd::Zero(_n_jnts);
    _pk = Eigen::VectorXd::Zero(_n_jnts);
    _e_recov = Eigen::VectorXd::Zero(_n_jnts);

    _pk_joule = Eigen::VectorXd::Zero(_n_jnts);
    _pk_mech = Eigen::VectorXd::Zero(_n_jnts);
    _pk_indct_est = Eigen::VectorXd::Zero(_n_jnts);

    _ek_joule = Eigen::VectorXd::Zero(_n_jnts);
    _ek_mech = Eigen::VectorXd::Zero(_n_jnts);
    _ek_indct = Eigen::VectorXd::Zero(_n_jnts);

    _omega_r = Eigen::VectorXd::Zero(_n_jnts);

    _num_diff_iq = NumDiff(_n_jnts, _dt);

    _num_int_joule = NumIntRt(_n_jnts, _dt);
    _num_int_mech = NumIntRt(_n_jnts, _dt);
    _recov_energy_integrators = std::vector<std::unique_ptr<NumIntRt>>(_n_jnts);

    for(int i = 0; i < _n_jnts; i++)
    {
        _recov_energy_integrators[i] = std::make_unique<NumIntRt>(1, _dt); // integrator for 1D variable
    }

    _mov_filter = MovAvrgFilt(_n_jnts, _dt, _filter_cutoff_freq);

    if(_dump_data2mat)
    {
        //  Initializing logger
        MatLogger2::Options opt;
        opt.default_buffer_size = _matlogger_buffer_size; // set default buffer size
        opt.enable_compression = true; // enable ZLIB compression

        std::string _dump_fullpath = _dump_path + std::string("/reg_energy");
        _logger = MatLogger2::MakeLogger(_dump_fullpath, opt); // date-time automatically appended

        _logger->set_buffer_mode(XBot::VariableBuffer::Mode::circular_buffer);

        _logger->add("dt", _dt);

        _logger->add("L_leak", _L_leak);
        _logger->add("L_m", _L_m);
        _logger->add("R", _R);
        _logger->add("L_q", _L_q);
        _logger->add("R_q", _R_q);
        _logger->add("Kt", _Kt);

        _logger->create("iq_k", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("iq_dot_est", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("omega_r", _n_jnts, 1, _matlogger_buffer_size);

        _logger->create("pk_joule", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("pk_mech", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("pk_indct_est", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("ek_joule", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("ek_mech", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("ek_indct", _n_jnts, 1, _matlogger_buffer_size);

        _logger->create("ek", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("pk", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("e_recov", _n_jnts, 1, _matlogger_buffer_size);

        _logger->create("ek_tot", 1, 1, _matlogger_buffer_size);
        _logger->create("pk_tot", 1, 1, _matlogger_buffer_size);
        _logger->create("e_recov_tot", 1, 1, _matlogger_buffer_size);

    }

}

RegEnergy::~RegEnergy()
{
    if(_dump_data2mat)
    {
        _logger.reset();
    }
}

void RegEnergy::set_e0(double& e0)
{
    if(_is_first_update)
    {

        _e0 = e0;

        if(_dump_data2mat)
        {
            _logger->add("e0", _e0);
        }

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

        if(_dump_data2mat)
        {
            _logger->add("iq_0", _iq_0);
        }

        _is_first_update = false;
    }

    if(!_use_iq_meas)
    {

        _iq_est->get_iq(_iq_k);
        _iq_est->get_omega_r(_omega_r);


    }else{

        if(_use_filt_iq_meas)
        {

            _iq_meas->get_last_iq_out_filt(_iq_k);

        }
        else{

            _iq_meas->get_last_iq_out(_iq_k);

        }

    }

    compute(); // updates energy and power values

}

void RegEnergy::compute()
{
    _num_diff_iq.add_sample(_iq_k);
    _num_diff_iq.dot(_iq_dot_est);

    compute_power(); // we compute all the terms of the power balance and then we integrate them
    // to get the energy

    evaluate_recov_energy(); // we get the regenerated energy

    _num_int_joule.add_sample(_pk_joule);
    _num_int_joule.get(_ek_joule);

    _num_int_mech.add_sample(_pk_mech);
    _num_int_mech.get(_ek_mech);

    compute_energy();
}

void RegEnergy::evaluate_recov_energy()
{
    // we integrated the total power if it is positive (i.e. towards the power source --> regenerative) and
    // only if we have triggered the evaluation of the recovered energy

    for(int i = 0; i < _n_jnts; i++)
    {

        if(_start_rec_energy_monitor && !_stop_rec_energy_monitor && _pk(i)>= 0)
        { // power from actuator i is flowing towards the power source

            _aux_1p1_vect_energy(0) = _pk(i); // we use this auxiliary vars to be able to pass by reference
            // a 1x1 eigen vector to the integration
            _recov_energy_integrators[i]->add_sample(_aux_1p1_vect_energy);
            _recov_energy_integrators[i]->get(_aux_1p1_vect_energy);

            _e_recov(i) = _aux_1p1_vect_energy(0);
        }

    }

    _e_recov_tot = _e_recov.sum();

}

void RegEnergy::compute_power()
{
    _pk_joule = 3.0/2.0 * _R_q.array() * (_iq_k.array().pow(2));

    _pk_mech = (_Kt.array() * _iq_k.array()) * _omega_r.array();

    _pk_indct_est = 3.0/2.0 * _L_q.array() * (_iq_k.array() * _iq_dot_est.array());

    _pk = - _pk_joule - _pk_indct_est - _pk_mech;

    _pk_tot = _pk.sum();
}

void RegEnergy::compute_energy()
{
    _ek_indct = 3.0/4.0 * _L_q.array() * (_iq_k.array().pow(2) - _iq_0.array().pow(2));

    _ek = _e0 - _ek_joule.array() -_ek_mech.array() - _ek_indct.array(); // using array to allow sum with a scalar (e0)

    _ek_tot = _ek.sum();

}

void RegEnergy::get_current_e_recov(Eigen::VectorXd& e_recov)
{
    e_recov = _e_recov;
}

void RegEnergy::get_current_e_recov(double& e_recov_tot)
{
    e_recov_tot = _e_recov_tot;
}

void RegEnergy::get(Eigen::VectorXd& ek)
{
    ek = _ek;
}

void RegEnergy::get(Eigen::VectorXd& ek, Eigen::VectorXd& pk)
{
    ek = _ek;
    pk = _pk;
}

void RegEnergy::get_p_terms(Eigen::VectorXd& pk_joule,
                 Eigen::VectorXd& pk_mech,
                 Eigen::VectorXd& pk_indct_est)
{
    pk_joule = _pk_joule;
    pk_mech = _pk_mech;
    pk_indct_est = _pk_indct_est;
}

void RegEnergy::get_e_terms(Eigen::VectorXd& ek_joule,
                 Eigen::VectorXd& ek_mech,
                 Eigen::VectorXd& ek_indct)
{
    ek_joule = _ek_joule;
    ek_mech = _ek_mech;
    ek_indct = _ek_indct;
}

void RegEnergy::get_p(double& p)
{
    p = _pk_tot;
}

void RegEnergy::get_e(double& e)
{
    e = _ek_tot;
}

void RegEnergy::set_log_buffsize(double size)
{
    _matlogger_buffer_size = abs(size);
}

void RegEnergy::add2log()
{
    if(_dump_data2mat)
    {
        _logger->add("iq_k", _iq_k);
        _logger->add("iq_dot_est", _iq_dot_est);
        _logger->add("omega_r", _omega_r);

        _logger->add("pk_joule", _pk_joule);
        _logger->add("pk_mech", _pk_mech);
        _logger->add("pk_indct_est", _pk_indct_est);
        _logger->add("ek_joule", _ek_joule);
        _logger->add("ek_mech", _ek_mech);
        _logger->add("ek_indct", _ek_indct);

        _logger->add("ek", _ek);
        _logger->add("pk", _pk);
        _logger->add("e_recov", _e_recov);

        _logger->add("ek_tot", _ek_tot);
        _logger->add("pk_tot", _pk_tot);
        _logger->add("e_recov_tot", _e_recov_tot);

    }

}

void RegEnergy::use_filt_iq_meas(bool filter_it)
{
    _use_filt_iq_meas =  filter_it;
}

void RegEnergy::enable_rec_energy_monitoring()
{
    _start_rec_energy_monitor = true;
    _stop_rec_energy_monitor = false;

}

void RegEnergy::disable_rec_energy_monitoring()
{
    _start_rec_energy_monitor = false;
    _stop_rec_energy_monitor = true;

}

void RegEnergy::reset_rec_energy()
{
    _e_recov.setZero();

    _e_recov_tot = 0.0;

    for(int i = 0; i < _n_jnts; i++)
    {
        _recov_energy_integrators[i]->reset(); // we reset each integrator
    }

}








