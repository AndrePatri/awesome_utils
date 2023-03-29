#ifndef SIGN_PROC_UTILS_H
#define SIGN_PROC_UTILS_H

#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include <map>
#include <vector>
#include <memory>

namespace SignProcUtils{

    /// \brief Class to perform numerical differentiation of a
    /// constant rate signal.
    class NumDiff
    {
      public:

        typedef std::weak_ptr<NumDiff> WeakPtr;
        typedef std::shared_ptr<NumDiff> Ptr;
        typedef std::unique_ptr<NumDiff> UniquePtr;

        NumDiff();

        NumDiff(int n_jnts, double dt, int order = 1);

        void add_sample(Eigen::VectorXd& sample);

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

        typedef std::weak_ptr<NumInt> WeakPtr;
        typedef std::shared_ptr<NumInt> Ptr;
        typedef std::unique_ptr<NumInt> UniquePtr;

        NumInt();

        NumInt(int n_jnts, double dt, double T_horizon);

        void add_sample(Eigen::VectorXd& sample);

        void get(Eigen::VectorXd& _int_sample); // get an estimate of integral
        // of the variable along the specified T_horizon (from 0 to T_horizon)
        // (assuming constant dt)

        void reset();

      private:

        Eigen::MatrixXd _window_data;

        double _dt = 0.0; // assuming evenly spaced samples
        double _T_horizon = 0.0;

        int _n_jnts = -1, _n_intervals = -1, _n_samples = -1;

    };

    /// \brief Class to perform numerical real time integration of a
    /// constant rate signal (exploits NumInt over a window of two samples).
    /// The two most recent samples are integrated and added to the current
    /// value of the integral.
    class NumIntRt
    {
      public:

        typedef std::weak_ptr<NumIntRt> WeakPtr;
        typedef std::shared_ptr<NumIntRt> Ptr;
        typedef std::unique_ptr<NumIntRt> UniquePtr;

        NumIntRt();

        NumIntRt(int n_jnts, double dt);

        void add_sample(Eigen::VectorXd& sample);

        void get(Eigen::VectorXd& _int_sample); // get an estimate of integral
        // of the variable along the specified T_horizon (from 0 to T_horizon)
        // (assuming constant dt)

        void reset(); // resets the current value of the integral to 0 and
        // also removes all the previously added samples

      private:

        Eigen::VectorXd _int_k, _int_km1; // integral value
        // at current sample and at the previous sample time

        int _n_jnts = -1;

        NumInt _num_int;

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
    ///
    /// As a side note, the exact equations to be solved for obtaining the bandwidth would be
    ///
    /// { window_size * 1/sqrt(2) = sum_{j = 0}^{N-1} cos(omega * dt * j)
    /// {
    /// {
    /// {sum_{j = 0}^{N-1} sin(omega * dt * j) = 0
    ///
    /// which is solved by the trascendental equation
    ///
    /// N * 1/sqrt(2) * sin(omega * h / 2.0)= sin(omega * dt * window_size /2)

    class MovAvrgFilt
    {
        public:

          typedef std::weak_ptr<MovAvrgFilt> WeakPtr;
          typedef std::shared_ptr<MovAvrgFilt> Ptr;
          typedef std::unique_ptr<MovAvrgFilt> UniquePtr;

          MovAvrgFilt();

          MovAvrgFilt(int n_jnts, double dt, int window_size = 10);
          MovAvrgFilt(int n_jnts, double dt, double cutoff_freq = 15);

          void add_sample(Eigen::VectorXd& sample);

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

           typedef std::weak_ptr<SignWithMem> WeakPtr;
           typedef std::shared_ptr<SignWithMem> Ptr;
           typedef std::unique_ptr<SignWithMem> UniquePtr;

            SignWithMem();

            SignWithMem(double signal_3sigma,
                        double tanh_coeff);

            int sign(double& value);

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
            double approx_sign(double& value);

    };

    /// \brief An implementation of continous sign function approximation.
    /// instead of computing
    ///
    /// sign(value)
    ///
    /// or sign(value/abs(value))
    ///
    /// we compute the sign as
    ///
    /// sign = tanh(k * value)
    ///
    /// How to give k a "smart" value?
    /// we approximate the derivative of the continous sign function at the origin as
    ///
    /// k/cosh^2(k * value)|_0 = k
    ///
    /// We want tanh(k * value) to be approximately 1 after a user defined threshold,
    /// which should depend on the uncertainty we have of the signal (i.e. noise).
    /// To do that, we use the linear approximation of the function and set it to 1
    /// after a predefined input threshold:
    ///
    /// tanh(k * signal_3sigma * alpha) \approx 1
    ///
    /// or, equivalently,
    ///
    /// tanh(k * signal_3sigma * alpha) = beta
    ///
    /// where
    /// - beta is a value indicating how much the sign function will be close to 1
    ///   after signal_3sigma * alpha (e.g. 0.95 a.k.a. 95%)
    /// - signal_3sigma is 3 * standard deviation of the signal noise
    /// - alpha indicates after how many signal_3sigma the sign output should be close to 1
    ///   and should be [1, +inf] where higher values produce a slower transition
    ///
    /// -> (1) k = atanh(beta) / (signal_3sigma * alpha)

    class SmoooothSign
    {
        public:

            typedef std::weak_ptr<SmoooothSign> WeakPtr;
            typedef std::shared_ptr<SmoooothSign> Ptr;
            typedef std::unique_ptr<SmoooothSign> UniquePtr;

            SmoooothSign();

            SmoooothSign(double signal_3sigma,
                           int alpha,
                           double beta,
                           bool use_threshold = false);

            SmoooothSign(double signal_3sigma,
                           int alpha,
                           bool use_threshold = false);

            double sign(double& value);

        private:

            double _signal_3sigma = 1e-8; // 3 * standard deviation
            // of the noise present in the signal of which it's necessary to compute the sign.

            int _alpha = 5;

            double _threshold = 1e-8;

            bool _use_threshold = false;

            double _k = 1.0, _beta = 0.95;

            double _sign = 0.0;

            double _value = 0.0;

            double smooooth_sign(double& value);
    };


}

#endif // SIGN_PROC_UTILS_H
