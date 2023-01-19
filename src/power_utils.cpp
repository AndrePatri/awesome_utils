#include<power_utils.hpp>

using namespace PowerUtils;

RegEnergy::RegEnergy(IqRosGetter::Ptr iq_meas,
                    IqEstimator::Ptr iq_est,
                    Eigen::VectorXd R,
                    Eigen::VectorXd L_leak, Eigen::VectorXd L_m,
                    double dt)
    :_iq_meas{iq_meas},
     _iq_est{iq_est},
     _R{R},
     _L_leak{L_leak}, _L_m{L_m},
     _dt{dt}
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
    _e_0 = Eigen::VectorXd::Zero(_n_jnts);
    _iq_k = Eigen::VectorXd::Zero(_n_jnts);

    _num_diff = NumDiff(_n_jnts, _dt);
    _num_int = NumInt(_n_jnts, _dt, _dt);

    _mov_filter = MovAvrgFilt(_n_jnts, _dt, _filter_cutoff_freq);

}

void RegEnergy::start(bool use_iq_meas)
{
    _use_iq_meas = _use_iq_meas;

    if(!_use_iq_meas)
    {

        _iq_est->get_iq(_iq_0); // gets latest iq estimate sample from the estimator

    }
    else{

        _iq_meas->get_last_iq_out(_iq_0); // gets latest iq estimate sample from the actual measurements

    }

    _was_start_called = true;
}
