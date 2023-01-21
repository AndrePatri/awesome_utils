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
#include <memory>

#include <matlogger2/matlogger2.h>

#include "calib_utils.hpp"
#include "../xbot2_utils/xbot2_utils.hpp"

using namespace CalibUtils;
using namespace Xbot2Utils;
using namespace SignProcUtils;

namespace PowerUtils{

    /// \brief Class estimate the energy flow towards the power bus
    /// when using FOC actuators
    class RegEnergy
    {
        public:

            typedef std::weak_ptr<RegEnergy> WeakPtr;
            typedef std::shared_ptr<RegEnergy> Ptr;
            typedef std::unique_ptr<RegEnergy> UniquePtr;

            RegEnergy() = default;

            ~RegEnergy();

            RegEnergy(IqRosGetter::Ptr iq_meas,
                      IqEstimator::Ptr iq_est,
                      Eigen::VectorXd R,
                      Eigen::VectorXd L_leak, Eigen::VectorXd L_m,
                      double dt,
                      bool use_iq_meas = false,
                      bool dump_data2mat = false,
                      std::string dump_path = "/tmp");

            void set_e0(double e0); // set initial energy level

            void set_omega_r(Eigen::VectorXd omega_r); // set current rotor velocity externally (only callable if use_iq_meas == true)

            void update();

            double get_p(); // getter for the total (joint-wise sum) energy flowgin towards the power bus
            double get_e(); // getter for the total (joint-wise sum) power

            void get(Eigen::VectorXd& ek);

            void get(Eigen::VectorXd& ek, Eigen::VectorXd& pk);

            void get_p_terms(Eigen::VectorXd _pk_joule,
                             Eigen::VectorXd _pk_mech,
                             Eigen::VectorXd _pk_indct_est);
            void get_e_terms(Eigen::VectorXd _ek_joule,
                             Eigen::VectorXd _ek_mech,
                             Eigen::VectorXd _ek_indct);

            void set_log_buffsize(double size);

            void add2log();

        private:

            int _n_jnts;

            bool _is_first_update = true,
                 _use_iq_meas = false,
                 _dump_data2mat = false;

            std::string _dump_path = "\tmp";

            double _matlogger_buffer_size = 1e5;

            double _dt; // [s]
            double _filter_cutoff_freq = 15.0; // [Hz]
            double _e0 = 0.0; // initial energy level
            double _ek_tot, _pk_tot; // total energy towards bus

            Eigen::VectorXd _L_leak, _L_m, _R,
                            _L_q, _R_q,
                            _Kt;

            Eigen::VectorXd _iq_k, _iq_0,
                            _iq_dot_est,
                            _omega_r;

            Eigen::VectorXd _pk_joule, _pk_mech, _pk_indct_est;
            Eigen::VectorXd _ek_joule, _ek_mech, _ek_indct;

            Eigen::VectorXd _ek, _pk; // joint-wise

            IqRosGetter::Ptr _iq_meas;
            IqEstimator::Ptr _iq_est;

            NumDiff _num_diff_iq;

            NumIntRt _num_int_joule, _num_int_mech;

            MovAvrgFilt _mov_filter;

            MatLogger2::Ptr _logger;

            void compute();

            void compute_power();

            void compute_energy();
    };

}

#endif // POWER_UTILS_HPP


