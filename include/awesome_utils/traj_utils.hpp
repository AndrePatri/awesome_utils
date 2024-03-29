#ifndef PLUGIN_UTILS_H
#define PLUGIN_UTILS_H

#include <math.h> 
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include <vector>
#include <map>

#include <matlogger2/matlogger2.h>

#include "sign_proc_utils.hpp"

namespace TrajUtils{

    class PeisekahTrans
    {
        public:

            PeisekahTrans(); // default constructor

            PeisekahTrans(Eigen::VectorXd start_point, Eigen::VectorXd end_point, double exec_time, double dt);
            
            Eigen::VectorXd eval_at(int node);

            void compute_peisekah_val(const double& phase, const double& start_point, const double& end_point,
                                        double& val);

            void compute_peisekah_val_dot(const double& phase, const double& start_point, const double& end_point, const double& t_exec,
                                                 double& val_dot);

            void compute_peisekah_vect_val(const double& phase, const Eigen::VectorXd& start_point, const Eigen::VectorXd& end_point,
                                           Eigen::VectorXd& val);

            void compute_peisekah_vect_val_dot(const double& phase, const Eigen::VectorXd& start_point, const Eigen::VectorXd& end_point, const double& t_exec,
                                                    Eigen::VectorXd& val_dot);
            double get_exec_time();
            double get_traj_dt();
            double get_n_nodes();
            double get_n_dim();
            void get_traj(Eigen::MatrixXd& traj);

        private:

            Eigen::VectorXd _start_point, _end_point;
            Eigen::MatrixXd _traj;
            Eigen::VectorXd _current_sample;

            double _exec_time, _dt;
            double _common_part_traj, _common_part_traj_dot; // preallocation
            int _n_nodes, _n_dim;

            void check_input_dim();
            void rate_adapter();
            void compute_traj();

    };

    class TrajLinInterp
    {
        public:

            TrajLinInterp(); // default constructor

            TrajLinInterp(Eigen::VectorXd sample_time, Eigen::MatrixXd input_traj, int interp_dir = 1,
                        double time_check_tol = 0.01);

            Eigen::MatrixXd eval_at(Eigen::VectorXd interp_times);

        private: 

            double _time_check_tol; // parameter used by check_time_vector

            Eigen::MatrixXd _traj; 
            Eigen::VectorXd _sample_times;

            int _n_samples, _n_dims;

            int _interp_dir;

            void check_time_vector(Eigen::VectorXd interp_times);
            int get_traj_n_samples();
            void check_dim_match();
            Eigen::VectorXd interp_1d(Eigen::VectorXd& interp_times, int dim_index);
            double interp_0d(double t_norm, double interval_dt, int first_indx, int second_indx, int dim_index);
            void get_closest_points(double inter_time, int& first_indx, int& second_indx);
    };

    class TrajLoader
    {
        public:

            TrajLoader(); // default constructor

            TrajLoader(std::string data_path, bool column_major = true, double resample_err_tol = 0.0001, bool load_from_csv = false);

            Eigen::MatrixXd read_data_from_csv(std::string data_path);

            int get_n_jnts();
            int get_n_nodes();
            int get_takeoff_index(double epsi = 1e-15);

            Eigen::VectorXd get_sample_times();
            Eigen::VectorXd compute_res_times(double dt_res);
            double get_exec_time();

            void get_loaded_traj(Eigen::MatrixXd& q_p, Eigen::MatrixXd& q_p_dot, Eigen::MatrixXd& tau, Eigen::MatrixXd& dt_opt, Eigen::MatrixXd& f_cont_res, Eigen::MatrixXd& iq);
            
            void resample(double res_dt, Eigen::MatrixXd& q_p_res, Eigen::MatrixXd& q_p_dot_res,
                        Eigen::MatrixXd& tau_res, Eigen::MatrixXd& f_cont_res, Eigen::MatrixXd& iq_res);

            void get_opt_dt(double& dt);

        private:

            std::string _data_path; // path to directory where data is stored. If a file name and mat extension is provided, data il loaded from a .mat database, otherwise from CSV file
            bool _column_major_order; // if true, the n. joints is given by the rows of the data and the number of samples by the columns, viceversa otherwise
            double _resample_err_tol; // acceptable resample execution time tolerance below which the resampled trajectory is considered valid
            bool _load_from_csv;

            std::string _q_p_name = "q_p"; // these names have to match the ones of the loaded data
            std::string _q_p_dot_name = "q_p_dot";
            std::string _efforts_name = "tau";
            std::string _f_cont_name = "f_contact";
            std::string _dt_name = "dt_opt";
            std::string _iq_name = "i_q";

            Eigen::MatrixXd _q_p;
            Eigen::MatrixXd _q_p_dot;
            Eigen::MatrixXd _tau;
            Eigen::MatrixXd _f_cont;
            Eigen::MatrixXd _iq;

            Eigen::MatrixXd _dt_opt;
            Eigen::VectorXd _sample_times;

            double _exec_time;

            int _n_nodes, _n_jnts;

            std::map<std::string, TrajLinInterp> opt_traj; 

            std::string get_file_extension(std::string file);

            int get_n_jnts(Eigen::MatrixXd& mat);

            int get_n_samples(Eigen::MatrixXd& mat);

            void check_loaded_data_dims();

            void load_data_from_csv(std::string data_path);

            void load_data_from_mat(std::string math_path);
    };

    class SweepCos
    {
        /// \brief Implements a simple sweep cosinusoidal trajectory:
        ///
        /// sweep(t) = q_bar + (q_ub - q_lb)/2.0 * cos(phi(t)) = q_bar + (q_ub - q_lb)/2.0 * cos(int_0^{t}{omega}\,dt)
        ///
        /// sweep_dot(t) = - (q_ub - q_lb)/2.0 * sin(phi(t)) * phi_dot = - (q_ub - q_lb)/2.0 * sin(int_0^{t}{omega}\,dt) * omega
        ///
        /// sweep_ddot(t) = - (q_ub - q_lb)/2.0 * (cos(phi(t)) * phi_dot^2 + sin(phi(t)) * phi_ddot)
        ///
        public:

            SweepCos();

            SweepCos(double& omega0, double& omegaf, double& T_omega,
                     double& q_lb, double& q_ub,
                     double& dt);

            void eval_at(double& time, double& val, double& val_dot, double& val_ddot);

            void get_stuff(double& phase_omega, double& _ramp_up, double& omega_k, double& time_ref);

            void reset();

        private:

            bool _ramp_up = true, _not_revered_yet = true;

            double _omega0 = 0.0, _omegaf = 0.0, _T_omega = 0.0, _phase_omega = 0.0,
                   _q_lb = 0.0, _q_ub = 0.0, _q_bar = 0.0,
                   _time_ref = 0.0, _ramp_time = 0.0,

                   _dt = -1.0;

            double _omega_k_int = 0.0, _omega_k = 0.0, _omega_dot_k = 0.0,  _omega_ddot_k = 0.0;

            SignProcUtils::NumIntRt _num_int_omega;

            Eigen::VectorXd _aux_vect;

            PeisekahTrans _peisekah_utils;

    };
}

#endif
