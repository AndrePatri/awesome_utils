#ifndef XBOT2_UTILS_HPP
#define XBOT2_UTILS_HPP

#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>

#include <xbot2/ros/ros_support.h>
#include <xbot_msgs/CustomState.h>
#include <xbot_msgs/JointState.h>
#include <xbot2/system/clock.h>

#if defined(EC_XBOT2_CLIENT_FOUND)
    #include <ec_xbot2/joint_ec.h>
#endif

#include <map>
#include <vector>
#include <memory>

#include "../awesome_utils/sign_proc_utils.hpp"

using namespace SignProcUtils;
using namespace XBot::chrono;

namespace Xbot2Utils{

    using namespace XBot;

    class IqRosGetter
    { // **** DEPRECATED ****//

      public:

        typedef std::weak_ptr<IqRosGetter> WeakPtr;
        typedef std::shared_ptr<IqRosGetter> Ptr;
        typedef std::unique_ptr<IqRosGetter> UniquePtr;

        IqRosGetter() = default;

        IqRosGetter(double dt,
                    bool verbose = false,
                    double mov_avrg_cutoff_freq = 15.0);

        void on_aux_signal_received(const xbot_msgs::CustomState& aux_sig);

        void on_js_signal_received(const xbot_msgs::JointState& js_sig);

        void get_last_iq_out(Eigen::VectorXd& iq_out_fb);

        void get_last_iq_out_filt(Eigen::VectorXd& iq_out_fb_filt);

        void get_last_iq_out_stamps(Eigen::VectorXd& timestamps);

        void get_time_reference(double& t_ref);

        void set_jnt_names(std::vector<std::string> jnt_names);

        bool is_iq_out_topic_active();

        void get_jnt_names(std::vector<std::string>& jnt_names);

      private:

        bool _is_first_aux_sig = true,
             _vars_were_initialized = false,
             _were_jnt_names_set = false, _set_jnt_names_from_ros = true,
             _verbose = false;

        double _time_ref = 0.0;

        double _dt;

        double _mov_avrg_cutoff_freq = 15.0;

        std::string _iq_sig_basename = "iq_out_fb";

        int _aux_types_encode_number = 0; // global counter used to assign incremental values to aux signal types

        int _n_active_jnts = 0;

        std::vector<int> _indices; // for holding the joint mapping

        std::map<std::string, int> _aux_msg_type_map;
        std::vector<std::string> _jnt_names; // auxiliary vector where the joint names (as visible in the js message) are saved.

        Eigen::VectorXd _iq_out_fb, _timestamps;

        MovAvrgFilt _mov_avrg_filter;

        void init_vars();

        template <typename T, typename t_v >
        int find_index(std::vector<T> input_v, t_v value);

        template <typename T>
        std::vector<int> map_indices(std::vector<T> input_v1, std::vector<T> input_v2);

        int get_aux_type_code(std::string msg_type);

        std::tuple<std::vector<int>, std::vector<double>> aux_mapper(const xbot_msgs::CustomState& aux_sig);

    };

    class IqOutRosGetter
    {

      public:

        typedef std::weak_ptr<IqRosGetter> WeakPtr;
        typedef std::shared_ptr<IqRosGetter> Ptr;
        typedef std::unique_ptr<IqRosGetter> UniquePtr;

        IqOutRosGetter() = default;

        IqOutRosGetter(std::vector<std::string> jnt_names,
                       double dt,
                       double mov_avrg_cutoff_freq = 15.0,
                       bool verbose = false);

        void on_aux_signal_received_ros(const xbot_msgs::CustomState& aux_sig); // gets iq measurements from
        // the ROS interface (high latencies are expected)

        #if defined(EC_XBOT2_CLIENT_FOUND)
        void on_aux_signal_received(const XBot::Hal::JointEcAux& aux_sig); // uses internal xbot2 topics
        // for minimum latency
        #endif

        void get_last_iq_out(Eigen::VectorXd& iq_out_fb); // get latest read iq measurement

        void get_last_iq_out_filt(Eigen::VectorXd& iq_out_fb_filt); // gets latest (filtered) iq measurement

        void get_last_iq_out_stamps(Eigen::VectorXd& timestamps); // gets timestamps of latest iq measurement

        void get_time_reference(double& t_ref); // gets the time reference w.r.t. the timestamps are returned

        bool is_joint_mapping_done(); // check whether the mapping between the aux signal and the user-set
        // desired joints was performed (which is true after the first aux signal was received)

        void get_jnt_names(std::vector<std::string>& jnt_names); // returns the joints names given by the user upon
        // class initialization

      private:

        bool _jnt_mapping_done = false,
             _vars_were_initialized = false,
             _verbose = false;

        double _time_ref = 0.0;

        double _dt;

        double _mov_avrg_cutoff_freq = 15.0;

        std::string _iq_out_sig_basename = "iq_out_fb"; // name associated with the
        // iq measurement

        std::string _exception;

        int _aux_types_encode_number = 0; // global counter used to assign incremental values to aux signal types

        int _n_jnts_req = -1, _n_jnts_aux_sig = -1;

        std::vector<int> _indices; // for holding the joint mapping between the user-provided joint names and
        // the ones present in the aux signal

        std::map<std::string, int> _aux_msg_type_map; // map between the message type name and its unique ID
        std::vector<std::string> _jnt_names; // desired joints, to be provided by the user upon initialization

        std::vector<double> _msg_value_remapped; // holds the reordered message values
        std::vector<int> _msg_type_remapped; // for each message value, holds the ID associated with the value

        Eigen::VectorXd _iq_out_fb; // vector of ordered (as _jnt_names) iq measurements
        Eigen::VectorXd _timestamps; // vector of ordered (as _jnt_names) iq timestamps

        MovAvrgFilt _mov_avrg_filter; // internal filter for the measurements

        void init_vars();

        template <typename T, typename t_v >
        int find_index(t_v value, std::vector<T> input_v);

        template <typename T>
        std::vector<int> map_indices(std::vector<T> input_v1, std::vector<T> input_v2);

        int get_aux_type_code(std::string msg_type);

        void aux_mapper_ros(const xbot_msgs::CustomState& aux_sig);

        #if defined(EC_XBOT2_CLIENT_FOUND)
        void aux_mapper(const XBot::Hal::JointEcAux& aux_sig);
        #endif

    };

}

#endif // XBOT2_UTILS_HPP
