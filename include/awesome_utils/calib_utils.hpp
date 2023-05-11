#ifndef CALIB_UTILS_H
#define CALIB_UTILS_H

#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include <map>
#include <vector>

#include <algorithm>
#include <chrono>
#include <memory>

#include "sign_proc_utils.hpp"
#include "typedefs.hpp"

#include <matlogger2/matlogger2.h>

using namespace SignProcUtils;
using namespace XBot;
using namespace utils_defs;

namespace CalibUtils{

    //************* Notes on rotor-side dynamics *************//

    /// The dynamics of the rotor is:
    ///
    /// (1) tau_m + tau_r = J_r * q_m_ddot
    ///
    /// where
    ///
    /// (2) tau_m = K_t * i_q
    ///
    /// and tau_r is the resultant torque on the rotor side coming from the reducer
    /// (or directly from the link, if there's no transmission).
    /// tau_r can be simply modeled as
    ///
    /// (3) tau_r = tau_lm + tau_f
    ///
    /// where
    /// tau_lm is the link-side measured torque, reflected to the rotor-side and
    /// tau_f is the disturbance torque on the rotor, which can be attributed to
    /// frictional effects.
    ///
    /// (4) tau_lm = - tau_l * eta
    ///
    /// where tau_l is the measured link-side torque and
    /// eta is the reduction ratio of the tranmission
    ///
    /// Specifically, we employ a simple Coulomb model for the tau_f:
    ///
    /// (5) tau_f = - ( k_d0 * sign(q_m_dot) + K_d1 * q_m_dot)

    //************* first order friction observer *************//

    /// It is trivial, for example, to build a first order friction observer exploiting the model (1).
    /// Specifically, one can obtain an estimate for tau_f as
    ///
    /// tau_f = - (tau_m + tau_lm - J_r * q_m_ddot)
    ///
    /// where q_m_ddot can be estimated numerically from measurements of q_m_dot
    /// The dynamics of the observer is
    ///
    /// y_dot = K * (tau_f - y)
    /// -> (5) y_dot = K * (- K_t * i_q + tau_l * eta + J_r * q_m_ddot) - y)
    ///
    /// where K is a positive definite matrix (usually a diagonal matrix)
    /// We integrate numerically (5) over a control interval to obtain
    ///
    /// y_k - y_km1 = K * ( int_{t_km1}^{t_k}{- K_t * i_q + tau_l * eta} * dt + J_r * (q_m_dot_k - q_m_dot_km1) - int_{t_km1}^{t_k}{ y } * dt )
    ///
    /// Let us approximate the integral of y over dt as (y_k + y_km1)/2 * h
    /// where h = t_k - t_km1 (i.e. trapezoidal integration).
    ///
    /// --> (I + K * h/2.0) * y_k = (I - K * h/2.0) * y_km1 + b_k
    ///
    /// where
    ///
    /// b_k = K * ( int_{t_km1}^{t_k}{-K_t * i_q + tau_l * eta} * dt + J_r * (q_m_dot_k - q_m_dot_km1) )
    ///
    /// which means that the update law for the observer is given by
    ///
    /// (6) y_k = (I + K * h/2.0)^{-1} * (I - K * h/2.0) * y_km1 + (I + K * h/2.0)^{-1} * b_k
    ///
    /// (6) is a discrete time-invariant system of the type
    ///
    /// y_k = A * y_km1 + B * u_k
    ///
    /// If we approximate (6) with its continous equivalent (5) we can compute the bandwidth of the
    /// observer considering that the filter is a first order system. -->
    ///
    /// 2 * PI * BW = K
    ///
    /// (7) BW = K / (2.0 * PI)
    ///
    /// Additional remarks:
    /// - the friction observer is a first order filter and, as such, introduced
    ///   a lag in the estimate of the friction torque. If this torque is used to compensate
    ///   friction, the lag can introduce instabilities in the closed loop system.
    /// - The BW computed with (7) is an approximation of the actual one, which should
    ///   be computed using the discrete dynamics (6)
    ///
    ///
    ///
    /// \brief Class to computed a friction-compensated estimate for the quadrature
    /// current of a classical three-phase BLDC actuator, without the need to
    /// employ the current measurement.
    ///
    /// The computation of the model requires a series of actuator parameters:
    /// - K_t --> motor torque constant
    /// - K_d0, K_d1 --> estimates for the static and dynamic friction coefficients (Coulomb friction model)
    /// - rot_MoI --> rotor axial moment of inertia
    /// - red_ratio --> actuator reduction ratio
    /// - tanh_coeff, q_dot_3sigma --> paramters needed by the awesome SignWithMem class

    //************* iq model calibration-related stuff *************//

    class IqEstimator
    {

      public:

        typedef std::weak_ptr<IqEstimator> WeakPtr;
        typedef std::shared_ptr<IqEstimator> Ptr;
        typedef std::unique_ptr<IqEstimator> UniquePtr;

        IqEstimator() = default;

        ~IqEstimator();

        IqEstimator(Eigen::VectorXd K_t,
                    Eigen::VectorXd K_d0, Eigen::VectorXd K_d1,
                    Eigen::VectorXd rot_MoI,
                    Eigen::VectorXd red_ratio,
                    int alpha = 10,
                    double q_dot_3sigma = 0.001,
                    bool dump_data2mat = false,
                    std::string dump_path = "/tmp");

        void set_current_state(Eigen::VectorXd& q_dot, Eigen::VectorXd& q_ddot, Eigen::VectorXd& tau); // link side state

        void get_iq_estimate(std::vector<float>& iq_est); // updates + gets estimate
        void get_iq_estimate(Eigen::VectorXd& iq_est); // updates + gets estimate
        void get_iq_estimate(std::vector<float>& iq_est,
                             Eigen::VectorXd& K_d0, Eigen::VectorXd& K_d1,
                             Eigen::VectorXd& rot_MoI, Eigen::VectorXd& K_t); // updates + gets estimate
        void get_iq_estimate(Eigen::VectorXd& iq_est,
                             Eigen::VectorXd& K_d0, Eigen::VectorXd& K_d1,
                             Eigen::VectorXd& rot_MoI, Eigen::VectorXd& K_t); // updates + gets estimate

        void update(Eigen::VectorXd& K_d0, Eigen::VectorXd& K_d1,
                    Eigen::VectorXd& rot_MoI, Eigen::VectorXd& K_t); // only updates
        void update(); // only updates

        void get_iq(Eigen::VectorXd& iq_est); // gets the current iq

        void get_tau_link(Eigen::VectorXd& tau);
        void get_tau_friction(Eigen::VectorXd& tau_friction);
        void get_q_ddot(Eigen::VectorXd& q_ddot);

        int get_n_jnts();

        void get_Kt(Eigen::VectorXd& Kt);

        void get_rot_MoI(Eigen::VectorXd& rot_MoI);

        void get_red_ratio(Eigen::VectorXd& red_ratio);

        void get_omega_r(Eigen::VectorXd& omega_r);

        void set_log_buffsize(double size);

        void add2log();

      private:

        Eigen::VectorXd _K_t, _K_d0, _K_d1, _rot_MoI, _red_ratio,
                        _iq_est,
                        _tau_friction_linkside,
                        _tau_friction_rotorside;

        Eigen::VectorXd _q_dot, _q_ddot, _tau;

        int _alpha = 5;
        double _q_dot_3sigma = 0.001; // max amplitute of the noise contained in the velocity signal
        // (basically equal to 3 * sigma, where sigma is the standard deviation of the noise)

        int _n_jnts;

        bool _use_thresholded_sign = true,
            _dump_data2mat = false;

        std::string _dump_path = "\tmp";

        double _matlogger_buffer_size = 1e5;

        void compute_iq_estimates();

        SmoooothSign _smooth_sign;

        MatLogger2::Ptr _logger;

    };

    /// \brief DEPRECATED!!!
    /// Class to calibrate the approximate model of quadrature current
    /// using measurements from a real experiment (suitable for real-time implementation).
    /// Currently, we only optimize for the static and dynamic Coulomb friction coefficients
    /// K_d0 and K_d1

    class IqCalib
    {

      public:

        IqCalib();

        IqCalib(int window_size,
                Eigen::VectorXd K_t,
                Eigen::VectorXd rot_MoI,
                Eigen::VectorXd red_ratio,
                Eigen::VectorXd ig_Kd0,
                Eigen::VectorXd ig_Kd1,
                int alpha = 10,
                double q_dot_3sigma = 0.001,
                double lambda = 2.0,
                bool verbose = false
                );

        void add_sample(Eigen::VectorXd& q_dot,
                   Eigen::VectorXd& q_ddot,
                   Eigen::VectorXd& iq,
                   Eigen::VectorXd& tau); // q_dot, q_ddot, etc... are link side values

        void set_ig(Eigen::VectorXd& ig_Kd0,
                    Eigen::VectorXd& ig_Kd1);

        void get_current_optimal_Kd(Eigen::VectorXd& Kd0_opt,
                               Eigen::VectorXd& Kd1_opt);

        void get_current_tau_total(Eigen::VectorXd& tau_total);
        void get_current_tau_friction(Eigen::VectorXd& tau_friction);
        void get_current_alpha(Eigen::VectorXd& alpha_d0, Eigen::VectorXd& alpha_d1);

        void get_sol_millis(Eigen::VectorXd& millis);

      private:

        bool _verbose = false; // whether to print info messages

        int _window_size = 300; // number of samples which will be retained
                                // and used to solve the calibration problem

        int _n_jnts = - 1; // dimension of the input signal ( = number of joints
                           // on which calibration is run)

        int _n_opt_vars = 2;

        int _alpha = 5; // handtuned coefficient used to approximate the
                            // ideal sign() function with a C^{inf} hyperbolic tangent function.
                            // The higher tanh_coeff, the steeper the transition from -1 to 1 is.
        double _q_dot_3sigma = 0.001; // max amplitute of the noise contained in the velocity signal
        // (basically equal to 3 * sigma, where sigma is the standard deviation of the noise)

        double _lambda = 2.0; // regularization gain for the least square problem

        bool _use_thresholded_sign = true;

        std::chrono::time_point<std::chrono::high_resolution_clock> _sol_start, _sol_stop;

        Eigen::VectorXd _sol_time; // current solution time of the regression problem,
        // (one solution time for each joint)

        Eigen::VectorXd _alpha_d0, _alpha_d1; // we choose to model the tau_friction
                                              // (choice dictated by the observations
                                              // on the measured mismatch between tau_total and tau)
                                              // as
                                              // tau_friction = Kd0 * sign(q_dot) + Kd1 * q_dot
                                              // which is a friction torque made of a static component (Kd0 * sign(q_dot))
                                              // and a dynamic component ( Kd1 * q_dot)

        Eigen::MatrixXd _I_lambda; // regularization identity matrix for
        // the least square calibration problem
        Eigen::VectorXd _b_lambda; // regularization vector(normally a vector of 0s)

        Eigen::MatrixXd _Alpha; // least square problem TOTAL matrix (for all joints)
                                // _Alpha_i * _Kd_i = _tau_friction_i, where i is the i-th joint
                                // _Alpha is obtained stacking up [_alpha_d0, _alpha_d1]
        Eigen::MatrixXd _A; // least square problem actual matrix (single joint + regularization)

        Eigen::VectorXd _tau_friction; // tau_friction is equal to the model mismatch
                                        // between tau_total and tau(measured torque). Ideally, these two quantities
                                        // coincide
        Eigen::VectorXd _b; // least square problem actual meas. vector (single joint + regularization)

        Eigen::VectorXd _tau_total; // total torque on the rotor (reported to
                                    // the link-side) which is computed (with measurements) as
                                    // tau_total = 1/red_ratio * (rot_MoI * q_ddot / red_ratio - K_t * iq)

        Eigen::VectorXd _Kd0, _Kd1; // current optimical iq model calibration coefficient
                                    // Kd0 -> static friction component
                                    // Kd1 -> dynamic friction component
        Eigen::VectorXd _lb_Kd, _ub_Kd; // upper and lower bound for the coefficients of the friction torque model
        Eigen::VectorXd _ig_Kd,
                        _ig_Kd0, _ig_Kd1; // initial guess for the optimal params to be used in the QP
        // (the QP is regularized not around 0 but around _ig_Kd)

        Eigen::VectorXd _q_dot, _q_ddot, _iq, _tau; // measurements necessary for the computation
                                                    // of the calibration coefficients

        Eigen::VectorXd _K_t,
                        _rot_MoI,
                        _red_ratio; // actuator paramters (supposed to be perfectly known in
                                    // advance)

        SmoooothSign _smooth_sign;

        void shift_data(Eigen::VectorXd& data,
                        bool towards_back = true); // shift vector data towards the
                                                   // back of the data qeue
        void shift_data(Eigen::MatrixXd& data,
                        bool rowwise = true,
                        bool towards_back = true); // shift vector data towards the
                                                   // back of the data qeue

        void solve_iq_cal_QP(int jnt_index); // solve the calibration QP for a single joint

        void compute_alphad0();
        void compute_alphad1();
        void assemble_Alpha();

        void compute_tau_friction();

    };

    /// \brief
    /// Class to calibrate the rotor-side dynamics of a BLDC actuator, employing a
    /// simple Couloumb-like friction model.
    /// The rotor-side eq. of motion is
    /// Kt * iq - [Kd0 * sign(q_dot_l) + Kd1 * q_dot_l] * eta - I_r * q_ddot_l / eta = - tau_l * eta
    ///
    /// where
    ///
    /// Kt is torque constant of the motor (to be estimated)
    /// iq is the quadrature current (measured)
    /// Kd0 is the static friction coeff. (to be estimated)
    /// Kd1 is the dynamic friction coeff. (to be estimated)
    /// eta is the reduction ratio of the actuator (0 < eta <= 1; e.g. 1/50)
    /// I_r is the rotational inertia of the motor along it's axis of rotation (to be estimated)
    /// tau_l is the link-side torque (measured)
    /// q_dot_l is the link-side velocity (measured)
    /// q_ddot_r is the rotor acceleration
    ///
    /// Given measurements of iq, q_dot_l, tau_l and an estimate/measuremtn of q_ddot_l, we can
    /// set up a simple regression prb over an horizon of data and estimate Kt, Kd0, Kd1, I_r or a subset
    /// of this set of paramters.
    /// By running the optimization in real-time with fresh new data at each solution sample we perform a
    /// Moving Horizon Estimation (MHE) of the rotor dynamics.
    ///
    class RotDynCal
    {

      public:

        typedef std::weak_ptr<RotDynCal> WeakPtr;
        typedef std::shared_ptr<RotDynCal> Ptr;
        typedef std::unique_ptr<RotDynCal> UniquePtr;

        RotDynCal();

        RotDynCal(int window_size,
                Eigen::VectorXd red_ratio,
                Eigen::VectorXd ig_Kt,
                Eigen::VectorXd ig_rot_MoI,
                Eigen::VectorXd ig_Kd0,
                Eigen::VectorXd ig_Kd1,
                double lambda = 2.0,
                int alpha = 10,
                double q_dot_3sigma = 0.001,
                bool verbose = false,
                bool compute_regr_err = false
                );

        RotDynCal(int window_size,
                Eigen::VectorXd red_ratio,
                Eigen::VectorXd ig_Kt,
                Eigen::VectorXd ig_rot_MoI,
                Eigen::VectorXd ig_Kd0,
                Eigen::VectorXd ig_Kd1,
                Eigen::VectorXd lambda,
                int alpha = 10,
                double q_dot_3sigma = 0.001,
                bool verbose = false,
                bool compute_regr_err = false
                );

        void add_sample(Eigen::VectorXd& q_dot,
                   Eigen::VectorXd& q_ddot,
                   Eigen::VectorXd& iq,
                   Eigen::VectorXd& tau); // q_dot, q_ddot, etc... are link side values

        void set_ig_Kd0(Eigen::VectorXd& ig_Kd0);
        void set_ig_Kd1(Eigen::VectorXd& ig_Kd1);
        void set_ig_Kt(Eigen::VectorXd& ig_Kt);
        void set_ig_MoI(Eigen::VectorXd& ig_rot_MoI);

        void set_lambda(Eigen::VectorXd& lambda); // sets regularization [Kt, Kd0, Kd1, rot_MoI]

        void set_lambda_high(Eigen::VectorXd& lambda_high); // high regularization aroung ig for inactive parameters

        void set_solution_mask(std::vector<bool>& mask); // mask to select which
        // rotor dynamics paramter/s to calibrate. Inactive paramters are assigned
        // a big regularization value around the ig and hence assume values arbitrarily
        // close to ig. To select which value the paraemter should converge to,
        // first set the associated ig.
        // ordering: [Kt, Kd0, Kd1, rot_MoI]

        void solve();

        void get_regr_error(Eigen::MatrixXd& regr_err);

        void get_opt_Kd0(Eigen::VectorXd& Kd0_opt);
        void get_opt_Kd1(Eigen::VectorXd& Kd1_opt);
        void get_opt_Kt(Eigen::VectorXd& Kt);
        void get_opt_rot_MoI(Eigen::VectorXd& rot_MoI);

        void get_tau_friction(Eigen::VectorXd& tau_friction);
        void get_Alpha(Eigen::MatrixXd& Alpha);
        void get_alpha_d(Eigen::VectorXd& alpha_d0,
                               Eigen::VectorXd& alpha_d1);
        void get_alpha_d0(Eigen::VectorXd& alpha_d0);
        void get_alpha_d1(Eigen::VectorXd& alpha_d1);
        void get_alpha_inertial(Eigen::VectorXd& alpha_inertial);
        void get_alpha_kt(Eigen::VectorXd& alpha_kt);

        void get_tau_motor(Eigen::VectorXd& tau_mot);
        void get_tau_inertial(Eigen::VectorXd& tau_inertial);

        void get_sol_millis(Eigen::VectorXd& millis);

        void get_cal_mask(std::vector<bool>& cal_mask);

        void get_lambda(Eigen::VectorXd& lambda);
        void get_lambda_des(Eigen::VectorXd& lambda_des);
        void get_lambda_high(Eigen::VectorXd& lambda_high);

        void get_ig_Kd0(Eigen::VectorXd& ig_Kd0);
        void get_ig_Kd1(Eigen::VectorXd& ig_Kd1);
        void get_ig_Kt(Eigen::VectorXd& ig_Kt);
        void get_ig_MoI(Eigen::VectorXd& ig_rot_MoI);

        void reset_window();
        bool is_window_full();

      private:

        bool _verbose = false; // whether to print info messages

        int _window_size = 300; // number of samples which will be retained
                                // and used to solve the calibration problem

        int _window_fill_counter = 0; // counts how many fresh samples have been added
        // to the window data (1 <= _window_fill_counter <= _windows_size)
        bool _is_window_full = false; // whether the window is full or not

        bool _was_solve_called = false; // whether the solve method was called after
        // having added new data

        bool _compute_regr_err = false; // whether to also compute the regression error

        int _n_jnts = - 1; // dimension of the input signal ( = number of joints
                           // on which calibration is run)

        int _n_opt_vars = 4; // we optimize for [Kt, rot_MoI, Kd0, Kd1]

        int _alpha = 5; // handtuned coefficient used to approximate the
                            // ideal sign() function with a C^{inf} hyperbolic tangent function.
                            // The higher tanh_coeff, the steeper the transition from -1 to 1 is.
        double _q_dot_3sigma = 0.001; // max amplitute of the noise contained in the velocity signal
        // (basically equal to 3 * sigma, where sigma is the standard deviation of the noise)

        double _very_high_regularization_kt = 1e6;
        double _very_high_regularization_rot_MoI = 1e10;
        double _very_high_regularization_kd0 = 1e6;
        double _very_high_regularization_kd1 = 1e6;

        Eigen::VectorXd _lambda_high,
                        _lambda, _lambda_des; // regularization gains for the least square problem

        bool _use_thresholded_sign = true;

        std::vector<bool> _sol_mask;

        std::chrono::time_point<std::chrono::high_resolution_clock> _sol_start, _sol_stop;

        Eigen::VectorXd _sol_time; // current solution time of the regression problem,
        // (one solution time for each joint)

        Eigen::VectorXd _sol; // calibration solution (for a single joint)

        Eigen::VectorXd _alpha_kt, _alpha_inertial, _alpha_d0, _alpha_d1, _alpha_tlink;

        Eigen::MatrixXd _I_lambda; // regularization identity matrix for
        // the least square calibration problem
        Eigen::MatrixXd _Lambda_reg;

        Eigen::VectorXd _b_lambda; // regularization vector

        Eigen::MatrixXd _Alpha; // least square problem TOTAL matrix (for all joints)
                                // _Alpha_i * _Kd_i = _tau_friction_i, where i is the i-th joint
                                // _Alpha is obtained stacking up [_alpha_d0, _alpha_d1]
        Eigen::MatrixXd _A; // least square problem actual matrix (single joint + regularization)

        Eigen::VectorXd _tau_mot; // motor torque on the rotor produced by the magn. field
        Eigen::VectorXd _tau_inertial; // inertial torque on the rotor
        Eigen::VectorXd _tau_friction; // estimated friction torque, seen at the rotor
        Eigen::VectorXd _tau_lm; // measured tau at the link, reported on the rotor

        Eigen::VectorXd _b; // least square problem actual meas. vector (single joint + regularization)


        Eigen::VectorXd _K_t,
                        _rot_MoI,
                        _Kd0, _Kd1; // current optimical iq model calibration coefficient
                                    // K_t -> torque constant of the motor
                                    // rot_MoI -> axial moment of inertia of the rotor
                                    // Kd0 -> static friction component
                                    // Kd1 -> dynamic friction component

        Eigen::VectorXd _lb, _ub, _ig; // upper and lower bound for the coefficients of the friction torque model
        // and full i.g. for the QP

        Eigen::VectorXd _ig_Kd0, _ig_Kd1, _ig_Kt, _ig_rot_MoI; // initial guess for the optimal params to be used in the QP

        Eigen::VectorXd _q_dot, _q_ddot, _iq, _tau; // measurements necessary for the computation
                                                    // of the calibration coefficients

        Eigen::VectorXd _red_ratio; // reduction ratio (0 < red_ratio <= 1)

        Eigen::MatrixXd _regr_error; // (n_jnts x window_size) --> regression errors with current
        // solution
        Eigen::VectorXd _prediction_aux;

        SmoooothSign _smooth_sign; // used to approximate the sign function (used by the Coulomb-like friction estimation)

        void shift_data(); // shift vector data towards the
                                                   // back of the data qeue

        void solve_mhe(int jnt_index); // solve the calibration QP for a single joint

        void compute_alphad0();
        void compute_alphad1();
        void compute_alpha_kt();
        void compute_alpha_inertial();
        void compute_alpha_tlink();

        void assemble_Alpha();

        void apply_solution_mask(int jnt_index);

        bool check_dims();
        void init_vars();

    };

}

#endif
