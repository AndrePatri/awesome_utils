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

#if defined(WITH_XBOT2)
#include "../xbot2_utils/xbot2_utils.hpp"
using namespace Xbot2Utils;
using namespace XBot;
#endif

using namespace CalibUtils;

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

            #if defined(WITH_XBOT2)

            RegEnergy(IqOutRosGetter::Ptr iq_meas,
                      IqEstimator::Ptr iq_est,
                      Eigen::VectorXd R,
                      Eigen::VectorXd L_leak, Eigen::VectorXd L_m,
                      double bus_p_leak, // constant (assumed so) leak power due to bus losses
                      double dt,
                      bool use_iq_meas = false,
                      bool dump_data2mat = false,
                      std::string dump_path = "/tmp");

            #else

            RegEnergy(IqEstimator::Ptr iq_est,
                      Eigen::VectorXd R,
                      Eigen::VectorXd L_leak, Eigen::VectorXd L_m,
                      double bus_p_leak, // constant (assumed so) leak power due to bus losses
                      double dt,
                      bool dump_data2mat = false,
                      std::string dump_path = "/tmp");

            #endif

            void set_e0(double& e0); // set initial energy level

            void set_omega_r(Eigen::VectorXd& omega_r); // set current rotor velocity externally (only callable if use_iq_meas == true)

            void update();

            void get_p(double& p); // getter for the total (joint-wise sum) energy flowing towards the power bus
            void get_e(double& e); // getter for the total (joint-wise sum) power

            void get(Eigen::VectorXd& ek);

            void get(Eigen::VectorXd& ek,
                     Eigen::VectorXd& pk);

            void get_p_terms(Eigen::VectorXd& _pk_joule,
                             Eigen::VectorXd& _pk_mech,
                             Eigen::VectorXd& _pk_indct_est);
            void get_e_terms(Eigen::VectorXd& _ek_joule,
                             Eigen::VectorXd& _ek_mech,
                             Eigen::VectorXd& _ek_indct);

            void get_current_e_recov(Eigen::VectorXd& e_recov);
            void get_current_e_recov(double& e_recov_tot);

            void enable_rec_energy_monitoring();
            void disable_rec_energy_monitoring();
            void reset_rec_energy();

            void set_log_buffsize(double size);

            #if defined(WITH_XBOT2)
            void use_filt_iq_meas(bool filter_it = true);
            #endif

            void add2log();

        private:

            int _n_jnts;

            bool _is_first_update = true,
                 _dump_data2mat = false,
                 _start_rec_energy_monitor = false,
                 _stop_rec_energy_monitor = false;

            bool _use_iq_meas = false,
                _use_filt_iq_meas = true;

            std::string _dump_path = "\tmp";

            double _matlogger_buffer_size = 1e5;

            double _bus_p_leak = 0.0;

            double _dt; // [s]
            double _filter_cutoff_freq = 15.0; // [Hz]
            double _e0 = 0.0; // initial energy level
            double _ek_tot, _pk_tot; // total energy towards bus
            double _e_recov_tot = 0.0; // total amount of energy recovered

            Eigen::VectorXd _aux_1p1_vect_power,
                            _aux_1p1_vect_energy; // auxiliary vars

            Eigen::VectorXd _L_leak, _L_m, _R,
                            _L_q, _R_q,
                            _Kt;

            Eigen::VectorXd _iq_k, _iq_0,
                            _iq_dot_est,
                            _omega_r;

            Eigen::VectorXd _pk_joule, _pk_mech, _pk_indct_est;
            Eigen::VectorXd _ek_joule, _ek_mech, _ek_indct;

            Eigen::VectorXd _ek, _pk,
                            _e_recov; // joint-wise

            Eigen::VectorXd _dummy_eigen_scalar; // aux variable

            #if defined(WITH_XBOT2)
            IqRosGetter::Ptr _iq_meas;
            #endif

            IqEstimator::Ptr _iq_est;

            NumDiff _num_diff_iq;

            NumIntRt _num_int_joule, _num_int_mech,
                    _num_int_tot_pow, _num_int_tot_pow_reg;

            std::vector<std::unique_ptr<NumIntRt>> _recov_energy_integrators; // one for each joint
            // since we a joint could be in regeneration mode while others may not

            MovAvrgFilt _mov_filter;

            MatLogger2::Ptr _logger;

            void compute();

            void compute_power();

            void compute_energy();

            void evaluate_recov_energy();
    };

}

#endif // POWER_UTILS_HPP


