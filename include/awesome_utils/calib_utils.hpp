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

#include "sign_proc_utils.hpp"

using namespace SignProcUtils;

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

        IqEstimator(Eigen::VectorXd K_t,
                    Eigen::VectorXd K_d0, Eigen::VectorXd K_d1,
                    Eigen::VectorXd rot_MoI,
                    Eigen::VectorXd red_ratio,
                    int alpha = 10,
                    double q_dot_3sigma = 0.001);

        IqEstimator();

        IqEstimator(Eigen::VectorXd K_t);

        void set_current_state(Eigen::VectorXd q_dot, Eigen::VectorXd q_ddot, Eigen::VectorXd tau);

        void get_iq_estimate(std::vector<float>& iq_est);
        void get_iq_estimate(Eigen::VectorXd& iq_est);
        void get_iq_estimate(std::vector<float>& iq_est,
                             Eigen::VectorXd K_d0, Eigen::VectorXd K_d1);
        void get_iq_estimate(Eigen::VectorXd& iq_est,
                             Eigen::VectorXd K_d0, Eigen::VectorXd K_d1);

        void get_tau_link(Eigen::VectorXd& tau);
        void get_tau_friction(Eigen::VectorXd& tau_friction);
        void get_q_ddot(Eigen::VectorXd& q_ddot);

      private:

        Eigen::VectorXd _K_t, _K_d0, _K_d1, _rot_MoI, _red_ratio,
                        _iq_est, _tau_l,
                        _tau_friction_linkside,
                        _tau_friction_rotorside;

        Eigen::VectorXd _q_dot, _q_ddot, _tau;

        int _alpha = 5;
        double _q_dot_3sigma = 0.001; // max amplitute of the noise contained in the velocity signal
        // (basically equal to 3 * sigma, where sigma is the standard deviation of the noise)

        int _n_jnts;

        bool _use_thresholded_sign = true;

        void compute_iq_estimates();

        SmoooothSign _smooth_sign;

    };

    /// \brief Class to calibrate the approximate model of quadrature current
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

        void add_sample(Eigen::VectorXd q_dot,
                   Eigen::VectorXd q_ddot,
                   Eigen::VectorXd iq,
                   Eigen::VectorXd tau);

        void set_ig(Eigen::VectorXd ig_Kd0,
                    Eigen::VectorXd ig_Kd1);

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

}

#endif
