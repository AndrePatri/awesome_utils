#ifndef POWER_UTILS_HPP
#define POWER_UTILS_HPP

#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include <map>
#include <vector>

#include <algorithm>
#include <chrono>

#include<calib_utils.hpp>
#include<xbot2_utils.hpp>

using namespace CalibUtils;
using namespace Xbot2Utils;
using namespace SignProcUtils;

namespace PowerUtils{

    /// \brief Class estimate the energy flow towards the power bus
    /// when using FOC actuators
    class RegEnergy
    {
        public:

            RegEnergy() = default;

            RegEnergy(IqRosGetter::Ptr iq_meas,
                      IqEstimator::Ptr iq_est,
                      Eigen::VectorXd R,
                      Eigen::VectorXd L_leak, Eigen::VectorXd L_m,
                      double dt);

            void start(bool use_iq_meas = false); // start "monitoring" the energy flow on the bus

            void get(double& e_current); // getter for the total (joint-wise sum) energy flowgin towards the power bus
            void get(Eigen::VectorXd& e_jnt); // getter for the energy flowign towards the power bus from each joint

        private:

            int _n_jnts;

            bool _was_start_called = false,
                 _use_iq_meas = false;

            double _dt; // [s]
            double _filter_cutoff_freq = 15.0; // [Hz]

            Eigen::VectorXd _L_leak, _L_m, _R,
                            _L_q, _R_q,
                            _Kt;

            Eigen::VectorXd _iq_k, _iq_0,
                            _e_0;

            IqRosGetter::Ptr _iq_meas;
            IqEstimator::Ptr _iq_est;

            NumDiff _num_diff;
            NumInt _num_int;

            MovAvrgFilt _mov_filter;

            void update(bool use_iq_meas = false);


    };

}

#endif // POWER_UTILS_HPP


