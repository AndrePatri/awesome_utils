#include "include/awesome_utils/calib_utils.hpp"
#include "include/awesome_utils/sign_proc_utils.hpp"

using std::tanh;
using namespace std::chrono;

using namespace CalibUtils;
using namespace SignProcUtils;

//************* IqEstimator *************//

IqEstimator::IqEstimator(Eigen::VectorXd K_t,
                         Eigen::VectorXd K_d0, Eigen::VectorXd K_d1,
                         Eigen::VectorXd rot_MoI,
                         Eigen::VectorXd red_ratio,
                         int alpha,
                         double q_dot_3sigma,
                         bool dump_data2mat,
                         std::string dump_path)
  :_K_t{K_t},
   _K_d0{K_d0}, _K_d1{K_d1},
   _rot_MoI{rot_MoI},
   _red_ratio{red_ratio},
   _alpha{alpha},
   _q_dot_3sigma{q_dot_3sigma},
   _dump_data2mat{dump_data2mat},
   _dump_path{dump_path}
{

    _n_jnts = _K_t.size();

    _q_dot = Eigen::VectorXd::Zero(_n_jnts);
    _q_ddot = Eigen::VectorXd::Zero(_n_jnts);
    _tau = Eigen::VectorXd::Zero(_n_jnts);
    _tau_friction_linkside = Eigen::VectorXd::Zero(_n_jnts);
    _tau_friction_rotorside = Eigen::VectorXd::Zero(_n_jnts);

    _iq_est = Eigen::VectorXd::Zero(_n_jnts);

    // initialize the smooth sign function
    _smooth_sign = SmoooothSign(_q_dot_3sigma, _alpha, _use_thresholded_sign);

    if(_dump_data2mat)
    {
        //  Initializing logger
        MatLogger2::Options opt;
        opt.default_buffer_size = _matlogger_buffer_size; // set default buffer size
        opt.enable_compression = true; // enable ZLIB compression

        std::string _dump_fullpath = _dump_path + std::string("/iq_est");
        _logger = MatLogger2::MakeLogger(_dump_fullpath, opt); // date-time automatically appended

        _logger->set_buffer_mode(XBot::VariableBuffer::Mode::circular_buffer);

        _logger->add("q_dot_3sigma", _q_dot_3sigma);
        _logger->add("q_dot_3sigma", _q_dot_3sigma);

        _logger->create("Kt", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("K_d0", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("K_d1", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("rot_MoI", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("red_ratio", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("iq_est", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("tau_friction_linkside", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("tau_friction_rotorside", _n_jnts, 1, _matlogger_buffer_size);

        _logger->create("q_dot", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("q_ddot", _n_jnts, 1, _matlogger_buffer_size);
        _logger->create("tau", _n_jnts, 1, _matlogger_buffer_size);

    }
}

IqEstimator::~IqEstimator()
{
    if(_dump_data2mat)
    {
        _logger.reset();
    }
}

void IqEstimator::get_iq_estimate(std::vector<float>& iq_est)
{
    int err = 0;
    if(iq_est.size() != _n_jnts)
    {
        err = err + 1;
    }

    if (err != 0)
    {
        std::string exception = std::string("IqEstimator::get_iq_estimate(): dimension mismatch of input data \n");

        throw std::invalid_argument(exception);
    }

    compute_iq_estimates();

    iq_est = std::vector<float>(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
      iq_est[i] = _iq_est[i]; // mapping Eigen Vector to std::vector(brutal way)
    }

}

void IqEstimator::get_iq_estimate(Eigen::VectorXd& iq_est)
{

    int err = 0;
    if(iq_est.size() != _n_jnts)
    {
        err = err + 1;
    }

    if (err != 0)
    {
        std::string exception = std::string("IqEstimator::get_iq_estimate(): dimension mismatch of input data \n");

        throw std::invalid_argument(exception);
    }

    compute_iq_estimates();

    iq_est = _iq_est;

}

void IqEstimator::get_iq_estimate(std::vector<float>& iq_est,
                                  Eigen::VectorXd& K_d0, Eigen::VectorXd& K_d1,
                                  Eigen::VectorXd& rot_MoI, Eigen::VectorXd& K_t)
{
    int err = 0;
    if(K_d0.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(K_d1.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(rot_MoI.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(K_t.size() != _n_jnts)
    {
        err = err + 1;
    }

    if (err != 0)
    {
        std::string exception = std::string("IqEstimator::compute_iq_estimates(): dimension mismatch of input data -> \n") +
                                std::string("K_d0 length: ") + std::to_string(K_d0.size()) + std::string("\n") +
                                std::string("K_d1 length: ") + std::to_string(K_d1.size()) + std::string("\n") +
                                std::string("rot_MoI length: ") + std::to_string(rot_MoI.size()) + std::string("\n") +
                                std::string("K_t length: ") + std::to_string(K_t.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        throw std::invalid_argument(exception);
    }

    // update K_d0 and K_d1 with the provided values
    _K_d0 = K_d0;
    _K_d1 = K_d1;
    _rot_MoI = rot_MoI;
    _K_t = K_t;

    compute_iq_estimates(); // compute iq estimate with the runtime set gains

    iq_est = std::vector<float>(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
       iq_est[i] = _iq_est[i]; // mapping Eigen Vector to std::vector (brutal way)
    }

}

void IqEstimator::get_iq_estimate(Eigen::VectorXd& iq_est,
                                  Eigen::VectorXd& K_d0, Eigen::VectorXd& K_d1,
                                  Eigen::VectorXd& rot_MoI, Eigen::VectorXd& K_t)
{
    int err = 0;
    if(K_d0.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(K_d1.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(rot_MoI.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(K_t.size() != _n_jnts)
    {
        err = err + 1;
    }

    if (err != 0)
    {
        std::string exception = std::string("IqEstimator::compute_iq_estimates(): dimension mismatch of input data -> \n") +
                                std::string("K_d0 length: ") + std::to_string(K_d0.size()) + std::string("\n") +
                                std::string("K_d1 length: ") + std::to_string(K_d1.size()) + std::string("\n") +
                                std::string("rot_MoI length: ") + std::to_string(rot_MoI.size()) + std::string("\n") +
                                std::string("K_t length: ") + std::to_string(K_t.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        throw std::invalid_argument(exception);
    }

    // update K_d0 and K_d1 with the provided values
    _K_d0 = K_d0;
    _K_d1 = K_d1;
    _rot_MoI = rot_MoI;
    _K_t = K_t;

    compute_iq_estimates(); // compute iq estimate with the runtime set gains

    iq_est = _iq_est;

}

void IqEstimator::get_iq(Eigen::VectorXd& iq_est)
{

    iq_est = _iq_est;

}

void IqEstimator::update()
{

    compute_iq_estimates(); // compute iq estimate with the runtime set gains

}

void IqEstimator::update(Eigen::VectorXd& K_d0, Eigen::VectorXd& K_d1,
                         Eigen::VectorXd& rot_MoI, Eigen::VectorXd& K_t)
{
    int err = 0;
    if(K_d0.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(K_d1.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(rot_MoI.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(K_t.size() != _n_jnts)
    {
        err = err + 1;
    }

    if (err != 0)
    {
        std::string exception = std::string("IqEstimator::compute_iq_estimates(): dimension mismatch of input data -> \n") +
                                std::string("K_d0 length: ") + std::to_string(K_d0.size()) + std::string("\n") +
                                std::string("K_d1 length: ") + std::to_string(K_d1.size()) + std::string("\n") +
                                std::string("rot_MoI length: ") + std::to_string(rot_MoI.size()) + std::string("\n") +
                                std::string("K_t length: ") + std::to_string(K_t.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        throw std::invalid_argument(exception);
    }

    // update K_d0 and K_d1 with the provided values
    _K_d0 = K_d0;
    _K_d1 = K_d1;
    _rot_MoI = rot_MoI;
    _K_t = K_t;

    compute_iq_estimates(); // compute iq estimate with the runtime set gains

}

void IqEstimator::get_tau_link(Eigen::VectorXd& tau)
{
    tau = _tau;
}

void IqEstimator::get_tau_friction(Eigen::VectorXd& tau_friction)
{
    tau_friction = _tau_friction_linkside;
}

void IqEstimator::get_q_ddot(Eigen::VectorXd& q_ddot)
{
    q_ddot = _q_ddot;
}

void IqEstimator::compute_iq_estimates()
{

    for (int i = 0; i < _n_jnts; i++)
    {

        double static_friction_effort_linkside = - _K_d0(i) * _smooth_sign.sign(_q_dot(i));
        double dynamic_friction_effort_linkside  = - _K_d1(i) * _q_dot(i);

        _tau_friction_linkside(i) = static_friction_effort_linkside + dynamic_friction_effort_linkside;

        _tau_friction_rotorside(i) = _tau_friction_linkside(i) * _red_ratio(i);

        double total_torque_on_motor = - _tau(i) * _red_ratio(i) + _tau_friction_rotorside(i);

        double motor_omega_dot = _q_ddot(i) / _red_ratio(i);

        double required_motor_torque = _rot_MoI(i) * motor_omega_dot - total_torque_on_motor;

        _iq_est(i) = required_motor_torque / _K_t(i);
    }

}

void IqEstimator::set_current_state(Eigen::VectorXd& q_dot, Eigen::VectorXd& q_ddot, Eigen::VectorXd& tau)
{
  int err = 0;
  if(q_dot.size() != _n_jnts)
  {
      err = err + 1;
  }
  if(q_ddot.size() != _n_jnts)
  {
      err = err + 1;
  }
  if(tau.size() != _n_jnts)
  {
      err = err + 1;
  }

  if (err != 0)
  {
      std::string exception = std::string("IqEstimator::set_current_state(): dimension mismatch of input data -> \n") +
                              std::string("q_dot length: ") + std::to_string(q_dot.size()) + std::string("\n") +
                              std::string("q_ddot length: ") + std::to_string(q_ddot.size()) + std::string("\n") +
                              std::string("tau length: ") + std::to_string(tau.size()) + std::string("\n") +
                              std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

      throw std::invalid_argument(exception);
  }
  else
  {
      _q_dot = q_dot;
      _q_ddot = q_ddot;
      _tau = tau;

  }

}

int IqEstimator::get_n_jnts()
{
    return _n_jnts;
}

void IqEstimator::get_Kt(Eigen::VectorXd& Kt)
{
    Kt = _K_t;
}

void IqEstimator::get_rot_MoI(Eigen::VectorXd& rot_MoI)
{
    rot_MoI = _rot_MoI;
}

void IqEstimator::get_red_ratio(Eigen::VectorXd& red_ratio)
{
    red_ratio = _red_ratio;
}

void IqEstimator::get_omega_r(Eigen::VectorXd& omega_r)
{
    omega_r = _q_dot.array() * _red_ratio.array().inverse();
}


void IqEstimator::set_log_buffsize(double size)
{
    _matlogger_buffer_size = abs(size);
}

void IqEstimator::add2log()
{
    if(_dump_data2mat)
    {
         _logger->add("Kt", _K_t);
         _logger->add("K_d0", _K_d0);
         _logger->add("K_d1", _K_d1);
         _logger->add("rot_MoI", _rot_MoI);
         _logger->add("red_ratio", _red_ratio);
         _logger->add("iq_est", _iq_est);
         _logger->add("tau_friction_linkside", _tau_friction_linkside);
         _logger->add("tau_friction_rotorside", _tau_friction_rotorside);
         _logger->add("q_dot", _q_dot);
         _logger->add("q_ddot", _q_ddot);
         _logger->add("tau", _tau);

    }
}

//************* IqCalib *************//

IqCalib::IqCalib()
{

}

IqCalib::IqCalib(int window_size,
                 Eigen::VectorXd K_t,
                 Eigen::VectorXd rot_MoI,
                 Eigen::VectorXd red_ratio,
                 Eigen::VectorXd ig_Kd0,
                 Eigen::VectorXd ig_Kd1,
                 int alpha,
                 double q_dot_3sigma,
                 double lambda,
                 bool verbose)
  :_window_size{window_size},
    _K_t{K_t}, _rot_MoI{rot_MoI},
    _red_ratio{red_ratio},
    _alpha{alpha},
    _q_dot_3sigma{q_dot_3sigma},
    _verbose{verbose},
    _lambda{lambda},
    _ig_Kd0{ig_Kd0},
    _ig_Kd1{ig_Kd1}
{

  // the linear regression problem (for a single joint) is written as
  // A * Kd = tau_friction_measured
  // where A is obtained as [alpha_d0, alpha_d1]
  // and tau_friction_measured is the "measurement" of the fictitious additional
  // friction torque
  // The (unconstrained) optimization problem to be solved is
  // min_{Kd} ||(A * Kd - tau_friction_measured)||^2
  // which is solved by Kd_opt = A_+ * tau_friction_measured
  // This least squared problem can be easily solved employing the builtin
  // utilities of Eigen library (@ https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html)

  int err = 0;
//  int err2 = 0;
  int err3 = 0;

  // collecting errors in input dimensions(if any)
  _n_jnts = K_t.size();
  if(rot_MoI.size() != _n_jnts)
  {
      err = err + 1;
  }
  if(red_ratio.size() != _n_jnts)
  {
      err = err + 1;
  }

//  if(_lb_Kd.size() != _n_opt_vars)
//  {
//      err2 = err2 + 1;
//  }
//  if(_ub_Kd.size() != _n_opt_vars)
//  {
//      err2 = err2 + 1;
//  }

  if(_ig_Kd0.size() != _n_jnts)
  {
      err3 = err3 + 1;
  }
  if(_ig_Kd1.size() != _n_jnts)
  {
      err3 = err3 + 1;
  }

  // throwing err. in case an error occurred
  if (err != 0)
  {
      std::string exception = std::string("IqCalib::IqCalib(): dimension mismatch in one or more of the input data -> \n") +
                              std::string("K_t length: ") + std::to_string(_K_t.size()) + std::string("\n") +
                              std::string("rot_MoI length: ") + std::to_string(_rot_MoI.size()) + std::string("\n") +
                              std::string("red_ratio length: ") + std::to_string(_red_ratio.size()) + std::string("\n") +
                              std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

      throw std::invalid_argument(exception);
  }
//  if (err2 != 0)
//  {
//      std::string exception = std::string("IqCalib::IqCalib(): dimension mismatch in one or more of the input data -> \n") +
//                              std::string("lb_Kd length: ") + std::to_string(_lb_Kd.size()) + std::string("\n") +
//                              std::string("ub_Kd length: ") + std::to_string(_ub_Kd.size()) + std::string("\n") +
//                              std::string("which do not match the required length of: ") + std::to_string(_n_opt_vars);

//      throw std::invalid_argument(exception);
//  }
  if (err3 != 0)
  {
      std::string exception = std::string("IqCalib::IqCalib(): dimension mismatch in one or more of the input data -> \n") +
                              std::string("ig_Kd0 length: ") + std::to_string(_ig_Kd0.size()) + std::string("\n") +
                              std::string("ig_Kd1 length: ") + std::to_string(_ig_Kd1.size()) + std::string("\n") +
                              std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

      throw std::invalid_argument(exception);
  }
  else
  {
    _Alpha = Eigen::MatrixXd::Zero(_window_size * _n_jnts, _n_opt_vars);

    _I_lambda = Eigen::MatrixXd::Identity(_n_opt_vars, _n_opt_vars);
    _b_lambda = Eigen::VectorXd::Zero(_I_lambda.rows());

    _A = Eigen::MatrixXd::Zero(_window_size + _I_lambda.rows(), _n_opt_vars);
    _b = Eigen::VectorXd::Zero(_window_size + _I_lambda.rows());

    _tau_total = Eigen::VectorXd::Zero(_n_jnts);
    _tau_friction = Eigen::VectorXd::Zero(_window_size * _n_jnts);

    _Kd0 = Eigen::VectorXd::Zero(_n_jnts);
    _Kd1 = Eigen::VectorXd::Zero(_n_jnts);

    _ig_Kd = Eigen::VectorXd::Zero(_n_opt_vars);
    _lb_Kd = Eigen::VectorXd::Zero(_n_opt_vars);
    _ub_Kd = Eigen::VectorXd::Zero(_n_opt_vars);

    _alpha_d0 = Eigen::VectorXd::Zero(_window_size * _n_jnts);
    _alpha_d1 = Eigen::VectorXd::Zero(_window_size * _n_jnts);

    _q_dot = Eigen::VectorXd::Zero(_n_jnts);
    _q_ddot = Eigen::VectorXd::Zero(_n_jnts);
    _iq = Eigen::VectorXd::Zero(_n_jnts);
    _tau = Eigen::VectorXd::Zero(_n_jnts);

    _sol_time = Eigen::VectorXd::Zero(_n_jnts);

    // initialize the smooth sign function
    _smooth_sign = SignProcUtils::SmoooothSign(_q_dot_3sigma, _alpha, _use_thresholded_sign);


  }

}

void IqCalib::shift_data(Eigen::VectorXd& data,
                         bool towards_back)
{
    // shifting data towards the back of the qeue
    // (starting from the penultimate sample) of EACH JOINT
    // (recall data from each joint is attached AFTER the data of the previous one)

    int last_sample_index = _window_size - 1; // index of last sample WITHIN each
    // joint

    for (int jnt = 0; jnt < _n_jnts; jnt++)
    { // shifting data for each joint
        if (towards_back)
        {
            int last_sample_index = _window_size - 1;

            for (int i = last_sample_index - 1; i >= 0; i--)
            {
                data(i + 1 + _window_size * jnt) = data(i + _window_size * jnt);
            }
        }
        else
        {
            for (int i = 1; i <= last_sample_index; i++)
            {
                data(i - 1 + _window_size * jnt) = data(i + _window_size * jnt);
            }
        }
    }

}

void IqCalib::shift_data(Eigen::MatrixXd& data,
                         bool rowwise,
                         bool towards_back)
{

    // shifting data towards the back of the qeue
    // (starting from the penultimate sample)
    // by default, shift rows towards bottom

    int last_sample_index = _window_size - 1;

    for (int jnt = 0; jnt < _n_jnts; jnt++)
    { // shifting data for each joint

        if (towards_back)
        {
            for (int i = last_sample_index - 1; i >= 0; i--)
            {

                if (rowwise)
                {
                    data.block(i + 1 + _window_size * jnt, 0, 1, data.cols()) = data.block(i + _window_size * jnt, 0, 1, data.cols());
                }
                else
                {
                    data.block(0, i + 1 + _window_size * jnt, data.rows(), 1) = data.block(0, i + _window_size * jnt, data.rows(), 1);
                }
            }
        }
        else
        {
            for (int i = 1; i <= last_sample_index; i++)
            {

                if (rowwise)
                {
                    data.block(i - 1 + _window_size * jnt, 0, 1, data.cols()) = data.block(i + _window_size * jnt, 0, 1, data.cols());
                }
                else
                {
                    data.block(0, i - 1 + _window_size * jnt, data.rows(), 1) = data.block(0, i + _window_size * jnt, data.rows(), 1);
                }
            }
        }
    }
}

void IqCalib::solve_iq_cal_QP(int jnt_index)
{
    // jnt_index -> 0-based indexing

    // extracting data of joint jnt_index
    _A.block(0, 0, _window_size, _A.cols()) = _Alpha.block(_window_size * jnt_index, 0, _window_size, _Alpha.cols());
    _A.block(_window_size, 0, _I_lambda.rows(), _I_lambda.cols()) = std::sqrt(_lambda) * _I_lambda; // adding regularization

    _b.segment(0, _window_size) = _tau_friction.segment(_window_size * jnt_index, _window_size);
    _b_lambda = _ig_Kd; // the regularization is done around _ig_Kd
    _b.segment(_window_size, _I_lambda.rows()) = std::sqrt(_lambda) * _b_lambda;

    _sol_start = high_resolution_clock::now(); // profiling solution time

    // solving unconstrained linear regression problem with Eigen builtin method
    Eigen::VectorXd opt_Kd = _A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(_b);

    _sol_stop = high_resolution_clock::now(); // profiling solution time

    _sol_time(jnt_index) = duration_cast<milliseconds>(_sol_stop - _sol_start).count();

    _Kd0(jnt_index) = opt_Kd(0);
    _Kd1(jnt_index) = opt_Kd(1);
}

void IqCalib::set_ig(Eigen::VectorXd& ig_Kd0,
                     Eigen::VectorXd& ig_Kd1)
{
    // check on ig dimensions
    int err = 0;
    if (ig_Kd0.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (ig_Kd1.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (err == 0)
    { // all ok --> set igs

        _ig_Kd0 = ig_Kd0;
        _ig_Kd1 = ig_Kd1;

    }
    else
    {
        std::string exception = std::string("IqCalib::get_current_optimal_Kd(): dimension mismatch in one or more of the input data -> \n") +
                                std::string("ig_Kd0 length: ") + std::to_string(ig_Kd0.size()) + std::string("\n") +
                                std::string("ig_Kd1 length: ") + std::to_string(ig_Kd1.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        throw std::invalid_argument(exception);
    }

}

void IqCalib::get_sol_millis(Eigen::VectorXd& millis)
{
    millis = _sol_time;
}

void IqCalib::get_current_optimal_Kd(Eigen::VectorXd& Kd0_opt,
                                     Eigen::VectorXd& Kd1_opt)

{
    // just in case the user provides uninitialized vector
    // or vector of wrong dimension

    Kd0_opt = Eigen::VectorXd::Zero(_n_jnts);
    Kd1_opt = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        _ig_Kd << _ig_Kd0(i), _ig_Kd1(i); // using the last set ig for Kd
        solve_iq_cal_QP(i);
    }

    // assign output
    Kd0_opt = _Kd0;
    Kd1_opt = _Kd1;

}

void IqCalib::add_sample(Eigen::VectorXd& q_dot,
                         Eigen::VectorXd& q_ddot,
                         Eigen::VectorXd& iq,
                         Eigen::VectorXd& tau)
{

    _q_dot = q_dot;
    _q_ddot = q_ddot;
    _iq = iq;
    _tau = tau;

    // shift data one sample behind
    shift_data(_alpha_d0);
    shift_data(_alpha_d1);
    shift_data(_tau_friction);

    // assign current sample values
    compute_alphad0();
    compute_alphad1();
    assemble_Alpha();
    compute_tau_friction();

}


void IqCalib::compute_alphad0()
{

    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_d0(i * _window_size) =  - _smooth_sign.sign(_q_dot(i)); // assign last sample
    }
}

void IqCalib::compute_alphad1()
{
    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_d1(i * _window_size) = - _q_dot(i); // assign last sample
    }
}

void IqCalib::assemble_Alpha()
{
    _Alpha.block(0, 0, _Alpha.rows(), 1) = _alpha_d0;
    _Alpha.block(0, 1, _Alpha.rows(), 1) = _alpha_d1;
}

void IqCalib::compute_tau_friction()
{
    for (int i = 0; i < _n_jnts; i++)
    {

        _tau_total(i) = 1.0 / _red_ratio(i) *
              ( _rot_MoI(i) * _q_ddot(i) / _red_ratio(i) - _K_t(i) * _iq(i)); // total torque acting
        // on the motor rotor estimated using the iq measurement and the estimate on the motor axis acceleration.
        // The difference between this component and the torque measured on the link side gives the cumulative unmodeled
        // effort on the rotor caused by dissipative actions present between the link and the rotor itself. Ideally,
        // this difference should be zero.
        // To model this dissipative effects, we use a simple static friction + dynamic friction model:
        // tau_friction = Kd0 * sign(q_dot) * Kd1 * q_dot

        _tau_friction(i * _window_size) = _tau_total(i) + _tau(i); // assign to last sample
    }
}

void IqCalib::get_current_tau_total(Eigen::VectorXd& tau_total)
{
    // just in case the user provides uninitialized vector
    // or vector of wrong dimension

    tau_total = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        tau_total(i) = _tau_total(i);
    }

}

void IqCalib::get_current_tau_friction(Eigen::VectorXd& tau_friction)
{
    // just in case the user provides uninitialized vector
    // or vector of wrong dimension

    tau_friction = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        tau_friction(i) = _tau_friction(i * _window_size);
    }

}

void IqCalib::get_current_alpha(Eigen::VectorXd& alpha_d0, Eigen::VectorXd& alpha_d1)
{

    // just in case the user provides uninitialized vector
    // or vector of wrong dimension

    alpha_d0 = Eigen::VectorXd::Zero(_n_jnts);
    alpha_d1 = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        alpha_d0(i) = _alpha_d0(i * _window_size);
        alpha_d1(i) = _alpha_d1(i * _window_size);
    }

}

//************* RotDynCal *************//

RotDynCal::RotDynCal()
{

}

RotDynCal::RotDynCal(int window_size,
                 Eigen::VectorXd red_ratio,
                 Eigen::VectorXd ig_Kt,
                 Eigen::VectorXd ig_rot_MoI,
                 Eigen::VectorXd ig_Kd0,
                 Eigen::VectorXd ig_Kd1,
                 double lambda,
                 int alpha,
                 double q_dot_3sigma,
                 bool verbose)
  :_window_size{window_size},
    _red_ratio{red_ratio},
    _ig_Kt{ig_Kt},
    _ig_rot_MoI{ig_rot_MoI},
    _ig_Kd0{ig_Kd0},
    _ig_Kd1{ig_Kd1},
    _alpha{alpha},
    _q_dot_3sigma{q_dot_3sigma},
    _verbose{verbose}
{

    // the linear regression problem (for a single joint) is written as
    // A * Kd = tau_friction_measured
    // where A is obtained as [alpha_d0, alpha_d1]
    // and tau_friction_measured is the "measurement" of the fictitious additional
    // friction torque
    // The (unconstrained) optimization problem to be solved is
    // min_{Kd} ||(A * Kd - tau_friction_measured)||^2
    // which is solved by Kd_opt = A_+ * tau_friction_measured
    // This least squared problem can be easily solved employing the builtin
    // utilities of Eigen library (@ https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html)

    _n_jnts = red_ratio.size();

    check_dims(); // checks input data dimensions; if not ok, throws error and stops

    init_vars();

    for (int i = 0; i < _n_opt_vars; i++)
    { // regularization vector
        _lambda(i) = lambda;
    }

}

RotDynCal::RotDynCal(int window_size,
                 Eigen::VectorXd red_ratio,
                 Eigen::VectorXd ig_Kt,
                 Eigen::VectorXd ig_rot_MoI,
                 Eigen::VectorXd ig_Kd0,
                 Eigen::VectorXd ig_Kd1,
                 Eigen::VectorXd lambda,
                 int alpha,
                 double q_dot_3sigma,
                 bool verbose)
  :_window_size{window_size},
    _red_ratio{red_ratio},
    _ig_Kt{ig_Kt},
    _ig_rot_MoI{ig_rot_MoI},
    _ig_Kd0{ig_Kd0},
    _ig_Kd1{ig_Kd1},
    _alpha{alpha},
    _q_dot_3sigma{q_dot_3sigma},
    _verbose{verbose}
{

    // the linear regression problem (for a single joint) is written as
    // A * Kd = tau_friction_measured
    // where A is obtained as [alpha_d0, alpha_d1]
    // and tau_friction_measured is the "measurement" of the fictitious additional
    // friction torque
    // The (unconstrained) optimization problem to be solved is
    // min_{Kd} ||(A * Kd - tau_friction_measured)||^2
    // which is solved by Kd_opt = A_+ * tau_friction_measured
    // This least squared problem can be easily solved employing the builtin
    // utilities of Eigen library (@ https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html)

    _n_jnts = red_ratio.size();

    check_dims(); // checks input data dimensions; if not ok, throws error and stops

    init_vars();

    set_lambda(lambda); // sets regularization vector

}

bool RotDynCal::check_dims()
{

    int err = 0;

    if(_ig_Kt.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(_ig_rot_MoI.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(_ig_Kd0.size() != _n_jnts)
    {
        err = err + 1;
    }
    if(_ig_Kd1.size() != _n_jnts)
    {
        err = err + 1;
    }

    if (err != 0)
    {
        std::string exception = std::string("IqCalib::IqCalib(): dimension mismatch in one or more of the input data -> \n") +
                                std::string("ig_Kt length: ") + std::to_string(_ig_Kt.size()) + std::string("\n") +
                                std::string("ig_rot_MoI length: ") + std::to_string(_ig_rot_MoI.size()) + std::string("\n") +
                                std::string("ig_Kd0 length: ") + std::to_string(_ig_Kd0.size()) + std::string("\n") +
                                std::string("ig_Kd1 length: ") + std::to_string(_ig_Kd1.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        throw std::invalid_argument(exception);
    }

    return true;

}

void RotDynCal::init_vars()
{

    _I_lambda = Eigen::MatrixXd::Identity(_n_opt_vars, _n_opt_vars);
    _Lambda_reg = _I_lambda;
    _b_lambda = Eigen::VectorXd::Zero(_I_lambda.rows());
    _lambda = Eigen::VectorXd::Zero(_I_lambda.rows());
    _lambda_des = _lambda;

    _A = Eigen::MatrixXd::Zero(_window_size + _I_lambda.rows(), _n_opt_vars);
    _b = Eigen::VectorXd::Zero(_window_size + _I_lambda.rows());

    _sol = Eigen::VectorXd::Zero(_n_opt_vars);

    _ig = Eigen::VectorXd::Zero(_n_opt_vars);
    _lb = Eigen::VectorXd::Zero(_n_opt_vars);
    _ub = Eigen::VectorXd::Zero(_n_opt_vars);

    _K_t = Eigen::VectorXd::Zero(_n_jnts);
    _rot_MoI = Eigen::VectorXd::Zero(_n_jnts);
    _Kd0 = Eigen::VectorXd::Zero(_n_jnts);
    _Kd1 = Eigen::VectorXd::Zero(_n_jnts);

    _Alpha = Eigen::MatrixXd::Zero(_window_size * _n_jnts, _n_opt_vars);
    _alpha_d0 = Eigen::VectorXd::Zero(_window_size * _n_jnts);
    _alpha_d1 = Eigen::VectorXd::Zero(_window_size * _n_jnts);
    _alpha_kt = Eigen::VectorXd::Zero(_window_size * _n_jnts);
    _alpha_inertial = Eigen::VectorXd::Zero(_window_size * _n_jnts);
    _alpha_tlink = Eigen::VectorXd::Zero(_window_size * _n_jnts);

    _q_dot = Eigen::VectorXd::Zero(_n_jnts);
    _q_ddot = Eigen::VectorXd::Zero(_n_jnts);
    _iq = Eigen::VectorXd::Zero(_n_jnts);
    _tau = Eigen::VectorXd::Zero(_n_jnts);

    _tau_friction = Eigen::VectorXd::Zero(_n_jnts);
    _tau_mot = Eigen::VectorXd::Zero(_n_jnts);
    _tau_inertial = Eigen::VectorXd::Zero(_n_jnts);
    _tau_lm= Eigen::VectorXd::Zero(_n_jnts);

    _sol_time = Eigen::VectorXd::Zero(_n_jnts);

    _sol_mask = std::vector<bool>(_n_opt_vars);

    // initialize the smooth sign function
    _smooth_sign = SignProcUtils::SmoooothSign(_q_dot_3sigma, _alpha, _use_thresholded_sign);

}

void RotDynCal::add_sample(Eigen::VectorXd& q_dot,
                         Eigen::VectorXd& q_ddot,
                         Eigen::VectorXd& iq,
                         Eigen::VectorXd& tau)
{

    // assigning state
    _q_dot = q_dot;
    _q_ddot = q_ddot;
    _iq = iq;
    _tau = tau;

    // shift old data one sample behind
    shift_data(_alpha_d0);
    shift_data(_alpha_d1);
    shift_data(_alpha_kt);
    shift_data(_alpha_inertial);
    shift_data(_alpha_tlink);

    // adding new data
    compute_alphad0();
    compute_alphad1();
    compute_alpha_inertial();
    compute_alpha_tlink();
    compute_alpha_kt();

    assemble_Alpha();

}

void RotDynCal::shift_data(Eigen::VectorXd& data,
                         bool towards_back)
{
    // shifting data towards the back of the qeue
    // (starting from the penultimate sample) of EACH JOINT
    // (recall data from each joint is attached AFTER the data of the previous one)

    int last_sample_index = _window_size - 1; // index of last sample WITHIN each
    // joint

    for (int jnt = 0; jnt < _n_jnts; jnt++)
    { // shifting data for each joint
        if (towards_back)
        {
            int last_sample_index = _window_size - 1;

            for (int i = last_sample_index - 1; i >= 0; i--)
            {
                data(i + 1 + _window_size * jnt) = data(i + _window_size * jnt);
            }
        }
        else
        {
            for (int i = 1; i <= last_sample_index; i++)
            {
                data(i - 1 + _window_size * jnt) = data(i + _window_size * jnt);
            }
        }
    }

}

void RotDynCal::shift_data(Eigen::MatrixXd& data,
                         bool rowwise,
                         bool towards_back)
{

    // shifting data towards the back of the qeue
    // (starting from the penultimate sample)
    // by default, shift rows towards bottom

    int last_sample_index = _window_size - 1;

    for (int jnt = 0; jnt < _n_jnts; jnt++)
    { // shifting data for each joint

        if (towards_back)
        {
            for (int i = last_sample_index - 1; i >= 0; i--)
            {

                if (rowwise)
                {
                    data.block(i + 1 + _window_size * jnt, 0, 1, data.cols()) = data.block(i + _window_size * jnt, 0, 1, data.cols());
                }
                else
                {
                    data.block(0, i + 1 + _window_size * jnt, data.rows(), 1) = data.block(0, i + _window_size * jnt, data.rows(), 1);
                }
            }
        }
        else
        {
            for (int i = 1; i <= last_sample_index; i++)
            {

                if (rowwise)
                {
                    data.block(i - 1 + _window_size * jnt, 0, 1, data.cols()) = data.block(i + _window_size * jnt, 0, 1, data.cols());
                }
                else
                {
                    data.block(0, i - 1 + _window_size * jnt, data.rows(), 1) = data.block(0, i + _window_size * jnt, data.rows(), 1);
                }
            }
        }
    }
}

void RotDynCal::apply_solution_mask(int jnt_index)
{
    for(int i = 0; i < _sol_mask.size(); i ++)
    {
        if(!_sol_mask[i])
        { // we don't want to estimate this parameter for joint "jnt_index" --> we set an high reg.
          // around the last provided setpoint

            _lambda(i) = _very_high_regularization;
        }
        if(_sol_mask[i])
        { // we use the last lambda setpoint
            _lambda(i) = _lambda_des(i);
        }
    }
}

void RotDynCal::solve_iq_cal_QP(int jnt_index)
{

    apply_solution_mask(jnt_index); // modifies _lambda, so that inactive parameters will be promoted to converge
    // to the provided initial guesses (which in this case are to be interpreted as nominal values)

    // extracting data of joint jnt_index
    _A.block(0, 0, _window_size, _A.cols()) = _Alpha.block(_window_size * jnt_index, 0, _window_size, _Alpha.cols());

    for(int i = 0; i < _I_lambda.cols(); i ++)
    {
        _Lambda_reg(i, i) = std::sqrt(_lambda(i));
    }
    _A.block(_window_size, 0, _I_lambda.rows(), _I_lambda.cols()) = _Lambda_reg; // adding regularization

    _b.segment(0, _window_size) = - _alpha_tlink.segment(_window_size * jnt_index, _window_size);
    _b_lambda = _ig; // the regularization is done around _ig

    _b.segment(_window_size, _I_lambda.rows()) = _lambda.array().sqrt() * _b_lambda.array();

    _sol_start = high_resolution_clock::now(); // profiling solution time

    // solving unconstrained linear regression problem with Eigen builtin method
    _sol = _A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(_b);

    _sol_stop = high_resolution_clock::now(); // profiling solution time

    _sol_time(jnt_index) = duration_cast<milliseconds>(_sol_stop - _sol_start).count();

    _K_t(jnt_index) = _sol(0);
    _Kd0(jnt_index) = _sol(1);
    _Kd1(jnt_index) = _sol(2);
    _rot_MoI(jnt_index) = _sol(3);

}

void RotDynCal::set_lambda(Eigen::VectorXd& lambda)
{

    if(lambda.size() != _lambda_des.size())
    {
        std::string warning = std::string("RotDynCal::set_lambda(): dimension mismatch -> \n") +
                                std::string("provided lambda size: ") + std::to_string(lambda.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_lambda_des.size());

        std::cout << Colors::kYellow << warning <<  Colors::kEndl << std::endl;
    }
    else
    {
        _lambda_des = lambda;
    }

}

void RotDynCal::set_solution_mask(std::vector<bool>& mask)
{

    if(mask.size() != _sol_mask.size())
    {
        std::string warning = std::string("RotDynCal::set_solution_mask(): dimension mismatch -> \n") +
                                std::string("provided mask size: ") + std::to_string(mask.size()) + std::string("\n") +
                                std::string("which do not match the required size of: ") + std::to_string(_sol_mask.size());

        std::cout << Colors::kYellow << warning <<  Colors::kEndl << std::endl;
    }
    else
    {
        _sol_mask = mask;
    }

}

void RotDynCal::set_ig_Kd0(Eigen::VectorXd& ig_Kd0)
{
    // check on ig dimensions
    int err = 0;
    if (ig_Kd0.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (err == 0)
    { // all ok --> set igs

        _ig_Kd0 = ig_Kd0;

    }
    else
    {
        std::string warning = std::string("RotDynCal::set_ig_Kd0(): dimension mismatch -> \n") +
                                std::string("ig_Kd0 length: ") + std::to_string(ig_Kd0.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        std::cout << Colors::kYellow << warning <<  Colors::kEndl << std::endl;
    }

}

void RotDynCal::set_ig_Kd1(Eigen::VectorXd& ig_Kd1)
{
    // check on ig dimensions
    int err = 0;
    if (ig_Kd1.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (err == 0)
    { // all ok --> set igs

        _ig_Kd1 = ig_Kd1;

    }
    else
    {
        std::string warning = std::string("RotDynCal::set_ig_Kd1(): dimension mismatch -> \n") +
                                std::string("ig_Kd1 length: ") + std::to_string(ig_Kd1.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        std::cout << Colors::kYellow << warning <<  Colors::kEndl << std::endl;
    }

}

void RotDynCal::set_ig_Kt(Eigen::VectorXd& ig_Kt)
{
    // check on ig dimensions
    int err = 0;
    if (ig_Kt.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (err == 0)
    { // all ok --> set igs

        _ig_Kt = ig_Kt;

    }
    else
    {
        std::string warning = std::string("RotDynCal::set_ig_Kt(): dimension mismatch -> \n") +
                                std::string("ig_Kt length: ") + std::to_string(ig_Kt.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        std::cout << Colors::kYellow << warning <<  Colors::kEndl << std::endl;
    }

}

void RotDynCal::set_ig_MoI(Eigen::VectorXd& ig_rot_MoI)
{
    // check on ig dimensions
    int err = 0;
    if (ig_rot_MoI.size() != _n_jnts)
    {
        err = err + 1;
    }
    if (err == 0)
    { // all ok --> set igs

        _ig_rot_MoI = ig_rot_MoI;

    }
    else
    {
        std::string warning = std::string("RotDynCal::set_ig_Kt(): dimension mismatch -> \n") +
                                std::string("ig_rot_MoI length: ") + std::to_string(_ig_rot_MoI.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        std::cout << Colors::kYellow << warning <<  Colors::kEndl << std::endl;
    }

}

void RotDynCal::solve()
{
    for (int i = 0; i < _n_jnts; i++)
    {
        _ig << _ig_Kt(i), _ig_Kd0(i), _ig_Kd1(i), _ig_rot_MoI(i); // assigning to the ig vector the
        // latest set i.g. values (either from the user of the defaults)

        solve_iq_cal_QP(i);
    }
}

void RotDynCal::get_sol_millis(Eigen::VectorXd& millis)
{
    millis = Eigen::VectorXd::Zero(_n_jnts);

    millis = _sol_time;
}

void RotDynCal::get_opt_Kd0(Eigen::VectorXd& Kd0_opt)

{
    Kd0_opt = Eigen::VectorXd::Zero(_n_jnts);

    Kd0_opt = _Kd0;

}

void RotDynCal::get_opt_Kd1(Eigen::VectorXd& Kd1_opt)

{
    Kd1_opt = Eigen::VectorXd::Zero(_n_jnts);

    Kd1_opt = _Kd1;

}

void RotDynCal::get_opt_Kt(Eigen::VectorXd& Kt)

{
    Kt = Eigen::VectorXd::Zero(_n_jnts);

    Kt = _K_t;

}

void RotDynCal::get_opt_rot_MoI(Eigen::VectorXd& rot_MoI)

{
    rot_MoI = Eigen::VectorXd::Zero(_n_jnts);

    rot_MoI = _rot_MoI;

}

void RotDynCal::compute_alphad0()
{

    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_d0(i * _window_size) =  - _smooth_sign.sign(_q_dot(i));
    }
}

void RotDynCal::compute_alphad1()
{
    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_d1(i * _window_size) = - _q_dot(i);
    }
}

void RotDynCal::compute_alpha_kt()
{
    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_kt(i * _window_size) = _iq(i);
    }
}

void RotDynCal::compute_alpha_inertial()
{
    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_inertial(i * _window_size) = - _q_ddot(i);
    }
}

void RotDynCal::compute_alpha_tlink()
{
    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_tlink(i * _window_size) = _tau(i) * _red_ratio(i);
    }
}

void RotDynCal::assemble_Alpha()
{
    // ordering: [Kt, Kd0, Kd1, rot_MoI]

    _Alpha.block(0, 0, _Alpha.rows(), 1) = _alpha_kt;
    _Alpha.block(0, 1, _Alpha.rows(), 1) = _alpha_d0;
    _Alpha.block(0, 2, _Alpha.rows(), 1) = _alpha_d1;
    _Alpha.block(0, 3, _Alpha.rows(), 1) = _alpha_inertial;
}

void RotDynCal::get_tau_friction(Eigen::VectorXd& tau_friction)
{

    tau_friction = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        tau_friction(i) = - _Kd0(i) * _smooth_sign.sign(_q_dot(i)) - _Kd1(i) * _q_dot(i);
    }

}

void RotDynCal::get_tau_motor(Eigen::VectorXd& tau_mot)
{

    tau_mot = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        tau_mot(i) = _K_t(i) * _iq(i);
    }

}

void RotDynCal::get_tau_inertial(Eigen::VectorXd& tau_inertial)
{

    tau_inertial = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        tau_inertial(i) = - _rot_MoI(i) * _q_ddot(i);
    }

}

void RotDynCal::get_alpha_d(Eigen::VectorXd& alpha_d0, Eigen::VectorXd& alpha_d1)
{

    // just in case the user provides uninitialized vector
    // or vector of wrong dimension

    alpha_d0 = Eigen::VectorXd::Zero(_n_jnts);
    alpha_d1 = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        alpha_d0(i) = _alpha_d0(i * _window_size);
        alpha_d1(i) = _alpha_d1(i * _window_size);
    }

}

void RotDynCal::get_alpha_inertial(Eigen::VectorXd& alpha_inertial)
{
    alpha_inertial = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        alpha_inertial(i) = _alpha_inertial(i * _window_size);
    }
}

void RotDynCal::get_alpha_kt(Eigen::VectorXd& alpha_kt)
{
    alpha_kt = Eigen::VectorXd::Zero(_n_jnts);

    for (int i = 0; i < _n_jnts; i++)
    {
        alpha_kt(i) = _alpha_kt(i * _window_size);
    }
}

void RotDynCal::get_cal_mask(std::vector<bool>& cal_mask)
{
    cal_mask = std::vector<bool>(_sol_mask.size());

    cal_mask = _sol_mask;

}

void RotDynCal::get_lambda(Eigen::VectorXd& lambda)
{
    lambda = Eigen::VectorXd::Zero(_lambda.size());

    lambda = _lambda;
}

void RotDynCal::get_lambda_des(Eigen::VectorXd& lambda_des)
{
    lambda_des = Eigen::VectorXd::Zero(_lambda.size());

    lambda_des = _lambda_des;
}

void RotDynCal::get_ig_Kd0(Eigen::VectorXd& ig_Kd0)
{
    ig_Kd0 = Eigen::VectorXd::Zero(_n_jnts);

    ig_Kd0 = _ig_Kd0;

}

void RotDynCal::get_ig_Kd1(Eigen::VectorXd& ig_Kd1)
{
    ig_Kd1 = Eigen::VectorXd::Zero(_n_jnts);

    ig_Kd1 = _ig_Kd1;
}

void RotDynCal::get_ig_Kt(Eigen::VectorXd& ig_Kt)
{
    ig_Kt = Eigen::VectorXd::Zero(_n_jnts);

    ig_Kt = _ig_Kt;
}

void RotDynCal::get_ig_MoI(Eigen::VectorXd& ig_rot_MoI)
{
    ig_rot_MoI = Eigen::VectorXd::Zero(_n_jnts);

    ig_rot_MoI = _ig_rot_MoI;
}
