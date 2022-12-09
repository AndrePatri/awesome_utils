#ifndef SIGN_PROC_UTILS_H
#define SIGN_PROC_UTILS_H

#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include <map>
#include <vector>

namespace SignProcUtils{

    /// \brief Class to perform numerical differentiation of a
    /// constant rate signal.
    class NumDiff
    {
      public:

        NumDiff();

        NumDiff(int n_jnts, double dt, int order = 1);

        void add_sample(Eigen::VectorXd sample);

        void dot(Eigen::VectorXd& sample_dot); // get an estimate of derivative
        // of the current sample

      private:

        // coefficient from n-th order num. differentiation
        // ordered from left to right, starting from the current sample

        // depending on _order, the estimate for the derivative of the samples
        // at the current time is given by the following backward differentiation equation:
        // v_dot(k) = { sum_{i=0}^{_order + 1} [ v(k -i) * k_n(_order) ) / _dt

        Eigen::VectorXd _k_n1, _k_n2, _k_n3,
                        _k_n4, _k_n5, _k_n6;

        Eigen::VectorXd _k_n; // assigned on the basis of the selected order of the estimate

        Eigen::MatrixXd _window_data;

        double _dt = 0.0; // assuming evenly spaced samples

        int _n_jnts, _order;
    };

    /// \brief Class to perform numerical integration of a
    /// constant rate signal over a fixed time horizon.
    class NumInt
    {
      public:

        NumInt();

        NumInt(int n_jnts, double dt, double T_horizon);

        void add_sample(Eigen::VectorXd sample);

        void get(Eigen::VectorXd& _int_sample); // get an estimate of integral
        // of the variable along the specified T_horizon (from 0 to T_horizon)
        // (assuming constant dt)

      private:

        Eigen::MatrixXd _window_data;

        double _dt = 0.0; // assuming evenly spaced samples
        double _T_horizon = 0.0;

        int _n_jnts = -1, _n_intervals = -1, _n_samples = -1;

    };

    /// \brief Class to perform moving-average filtering of a signal s(t).
    ///
    /// To get an approximate expression for this filter bandwidth (valid for large window size)
    /// consider the continous version of the filter:
    ///
    /// (0.0) y(t) = -1/T_w * int_0^{T_w} s(t - tau) d_tau
    ///
    /// say s = A * sin(omega * t)
    ///
    /// the result of the integral is
    ///
    /// y(t) = S/(omega * T_w) * sqrt(2 * (1 - cos(omega * T_w))) * sin(omega * t + phi)
    ///
    /// the bandwidth is defined as y(t) / (A * sin(omega * t)) = 1/sqrt(2)
    ///
    /// by manipulating the expression one gets
    ///
    /// sin^2(omega * T_W/2.0) = 1/2.0 * (omega * T_w/2.0)^2
    ///
    /// considering that omega * T_w is positive
    /// the equation to be solved for is sin x = x/sqrt(2.0)
    ///
    /// which shows that necessarily 0 <= x <= sqrt(2.0)
    /// This equation is trivial to solve numerically. The result is
    ///
    /// x \approx k:= 1.3915576345446
    ///
    /// This means that an approximate formula for the bandwidth is
    ///
    /// (1) f_bw \approx k/(PI * dt * (window_size - 1))
    ///
    /// or, equivalently, for the window size
    ///
    /// (2)  window_size \approx 1 + k/(PI * f_bw * dt)
    ///
    /// A less approximate expression for (1) and (2) can be found by working
    /// directly with the discrete version of (0.0).
    /// See https://dsp.stackexchange.com/questions/9966/what-is-the-cut-off-frequency-of-a-moving-average-filter
    ///
    /// (1.1) f_bw \approx sqrt( (_magic_number / (_samples_dt * _cutoff_freq))^2 + 1)
    ///
    /// (1.2) window_size \approx (k/PI) /dt * 1.0 / sqrt(_window_size^2 - 1)

    class MovAvrgFilt
    {
        public:

          MovAvrgFilt();

          MovAvrgFilt(int n_jnts, double dt, int window_size = 10);
          MovAvrgFilt(int n_jnts, double dt, double cutoff_freq = 15);

          void add_sample(Eigen::VectorXd sample);

          void get(Eigen::VectorXd& filt_sample);

          void get_cutoff_freq(double& cutoff_f);
          void get_window_size(int& window_size);

        private:

          Eigen::MatrixXd _window_data;

          double _k = 1.3915576345446;
          double _PI =  2 * std::acos(0);
          double _magic_number = _k / _PI; // reference
          // @ https://dsp.stackexchange.com/questions/9966/what-is-the-cut-off-frequency-of-a-moving-average-filter
          int _n_jnts, _window_size;

          double _cutoff_freq = -1.0;

          double _samples_dt = -1.0;

        bool _is_first_run = true;

    };

    /// \brief Simple implementation of a sign with memory.
    /// Instead of using directly the sign of the signal to detect sign changes,
    /// we use a continous approximation of the sign function (tanh function)
    /// which has two benefits -> one: it is naturally normalized ([-1, 1])
    /// and is continous and differentiable in 0. On the contrary, using directly a
    /// normalized version of the signal would cause numerical instability close to the origin.
    ///
    /// The sign change is triggered when the
    /// continous approximation of the sign function goes above or below
    /// the inteval [ -tanh(_tanh_coeff * signal_3sigma), tanh(_tanh_coeff * signal_3sigma)]
    ///
    /// where signal_3sigma is 3 * standard deviation of the signal noise.
    ///
    /// This way we are pretty much sure that the sign implementation will be immune to
    /// instability due to noise (assuming some characteristics of the noise are known).
    ///
    /// This mechanism introduces a lower bound on the sensitivity of the sign function.
    ///
    /// With an ideal and clean signal, one would set signal_3sigma to 0 and have infinite
    /// sensitivity.

    class SignWithMem
    {
        public:

            SignWithMem();

            SignWithMem(double signal_3sigma,
                        double tanh_coeff);

            int sign(double value);

        private:

            double _tanh_thresh = 0.05;

            double _tanh_coeff = 10.0; // this is not necessary
            // it's just a scaling to make the value of the threshold prettier

            double _signal_3sigma = 0.001; // 3 * standard deviation
            // of the noise present in the signal of which it's necessary to compute the sign.

            int _sign, _previous_sign = 1;

            double _value = 0.0;
            double _approx_sign = 0.0; // continous approximation of the sign
            // function

            void sign_with_memory();
            double approx_sign(double value);

    };

}

#endif // SIGN_PROC_UTILS_H
