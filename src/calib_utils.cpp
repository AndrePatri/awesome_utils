#include "awesome_utils/calib_utils.hpp"
#include "awesome_utils/sign_proc_utils.hpp"

using std::tanh;
using namespace std::chrono;

using namespace CalibUtils;
using namespace SignProcUtils;

//************* IqEstimator *************//

IqEstimator::IqEstimator(Eigen::VectorXd K_t,
                         Eigen::VectorXd K_d0, Eigen::VectorXd K_d1,
                         Eigen::VectorXd rot_MoI,
                         Eigen::VectorXd red_ratio,
                         double tanh_coeff,
                         double q_dot_3sigma)
  :_K_t{K_t},
   _K_d0{K_d0}, _K_d1{K_d1},
   _rot_MoI{rot_MoI},
   _red_ratio{red_ratio},
   _tanh_coeff{tanh_coeff},
   _q_dot_3sigma{q_dot_3sigma}
{

    _n_jnts = _K_t.size();

    _q_dot = Eigen::VectorXd::Zero(_n_jnts);
    _q_ddot = Eigen::VectorXd::Zero(_n_jnts);
    _tau = Eigen::VectorXd::Zero(_n_jnts);
    _tau_friction = Eigen::VectorXd::Zero(_n_jnts);

    _iq_est = Eigen::VectorXd::Zero(_n_jnts);

    // initialize the awesome sign function
    _sign_with_memory = SignProcUtils::SignWithMem(_q_dot_3sigma, _tanh_coeff);
}

IqEstimator::IqEstimator()
{

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
    }

    compute_iq_estimates();

    iq_est = _iq_est;

}

void IqEstimator::get_iq_estimate(std::vector<float>& iq_est,
                                  Eigen::VectorXd K_d0, Eigen::VectorXd K_d1)
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

    if (err != 0)
    {
        std::string exception = std::string("IqEstimator::compute_iq_estimates(): dimension mismatch of input data -> \n") +
                                std::string("K_d0 length: ") + std::to_string(K_d0.size()) + std::string("\n") +
                                std::string("K_d1 length: ") + std::to_string(K_d1.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        throw std::invalid_argument(exception);
    }
    else
    {
        // update K_d0 and K_d1 with the provided values
        _K_d0 = K_d0;
        _K_d1 = K_d1;

        compute_iq_estimates(); // compute iq estimate with the runtime set gains

        iq_est = std::vector<float>(_n_jnts);

        for (int i = 0; i < _n_jnts; i++)
        {
           iq_est[i] = _iq_est[i]; // mapping Eigen Vector to std::vector (brutal way)
        }

    }

}

void IqEstimator::get_iq_estimate(Eigen::VectorXd& iq_est,
                                  Eigen::VectorXd K_d0, Eigen::VectorXd K_d1)
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

    if (err != 0)
    {
        std::string exception = std::string("IqEstimator::compute_iq_estimates(): dimension mismatch of input data -> \n") +
                                std::string("K_d0 length: ") + std::to_string(K_d0.size()) + std::string("\n") +
                                std::string("K_d1 length: ") + std::to_string(K_d1.size()) + std::string("\n") +
                                std::string("which do not match the required length of: ") + std::to_string(_n_jnts);

        throw std::invalid_argument(exception);
    }
    else
    {
        // update K_d0 and K_d1 with the provided values
        _K_d0 = K_d0;
        _K_d1 = K_d1;

        compute_iq_estimates(); // compute iq estimate with the runtime set gains

        iq_est = _iq_est;

    }

}

void IqEstimator::get_tau_link(Eigen::VectorXd& tau)
{
    tau = _tau;
}

void IqEstimator::get_tau_friction(Eigen::VectorXd& tau_friction)
{
    tau_friction = _tau_friction;
}

void IqEstimator::get_q_ddot(Eigen::VectorXd& q_ddot)
{
    q_ddot = _q_ddot;
}

void IqEstimator::compute_iq_estimates()
{

    for (int i = 0; i < _n_jnts; i++)
    {

        double static_friction_effort = _K_d0(i) * _sign_with_memory.sign(_q_dot(i));
        double dynamic_friction_effort = _K_d1(i) * _q_dot(i);

        _tau_friction(i) = static_friction_effort + dynamic_friction_effort;

        double total_torque_on_motor = _tau(i) + _tau_friction(i);

        double motor_omega_dot = _q_ddot(i) / _red_ratio(i);

        double required_motor_torque = _rot_MoI(i) * motor_omega_dot + total_torque_on_motor * _red_ratio(i);

        _iq_est(i) = required_motor_torque / _K_t(i);
    }

}

void IqEstimator::set_current_state(Eigen::VectorXd q_dot, Eigen::VectorXd q_ddot, Eigen::VectorXd tau)
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
                 double tanh_coeff,
                 double q_dot_3sigma,
                 double lambda,
                 bool verbose)
  :_window_size{window_size},
    _K_t{K_t}, _rot_MoI{rot_MoI},
    _red_ratio{red_ratio},
    _tanh_coeff{tanh_coeff},
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

    // initialize the awesome sign-with-memory function
    _sign_with_memory = SignProcUtils::SignWithMem(_q_dot_3sigma, _tanh_coeff);


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

void IqCalib::set_ig(Eigen::VectorXd ig_Kd0,
                     Eigen::VectorXd ig_Kd1)
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

void IqCalib::add_sample(Eigen::VectorXd q_dot,
                         Eigen::VectorXd q_ddot,
                         Eigen::VectorXd iq,
                         Eigen::VectorXd tau)
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
        _alpha_d0(i * _window_size) =  (double) _sign_with_memory.sign(_q_dot(i)); // assign last sample
    }
}

void IqCalib::compute_alphad1()
{
    for (int i = 0; i < _n_jnts; i++)
    {
        _alpha_d1(i * _window_size) = _q_dot(i); // assign last sample
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
              ( - _rot_MoI(i) * _q_ddot(i) / _red_ratio(i) + _K_t(i) * _iq(i)); // total torque acting
        // on the motor rotor estimated using the iq measurement and the estimate on the motor axis acceleration.
        // The difference between this component and the torque measured on the link side gives the cumulative unmodeled
        // effort on the rotor caused by dissipative actions present between the link and the rotor itself. Ideally,
        // this difference should be zero.
        // To model this dissipative effects, we use a simple static friction + dynamic friction model:
        // tau_friction = Kd0 * sign(q_dot) * Kd1 * q_dot

        _tau_friction(i * _window_size) = _tau_total(i) - _tau(i); // assign to last sample
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

