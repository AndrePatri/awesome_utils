#ifndef SIGN_PROC_UTILS_H
#define SIGN_PROC_UTILS_H

#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include <map>
#include <vector>

namespace SignProcUtils{

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

          double _magic_number = 0.442946470689452340308369; // reference
          // @ https://dsp.stackexchange.com/questions/9966/what-is-the-cut-off-frequency-of-a-moving-average-filter

          int _n_jnts, _window_size;

          double _cutoff_freq = -1.0;

          double _samples_dt = -1.0;

        bool _is_first_run = true;

    };

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
