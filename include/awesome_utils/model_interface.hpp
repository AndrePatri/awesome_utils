#ifndef MODEL_INTERFACE_HPP
#define MODEL_INTERFACE_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <string>
#include <vector>

#include <pinocchio/algorithm/aba-derivatives.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/center-of-mass-derivatives.hpp>
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/centroidal-derivatives.hpp>
#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/contact-dynamics.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/energy.hpp>
#include <pinocchio/algorithm/frames-derivatives.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/geometry.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/model.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/regressor.hpp>

#include <pinocchio/parsers/urdf.hpp>

#include "typedefs.hpp"
#include "sign_proc_utils.hpp"

using namespace Eigen;
using namespace SignProcUtils;

namespace ModelInterface
{

    /// \brief Simple model interface wrapper class
    /// around the awesome Pinocchio's model implementation.

    class Model
    {
        public:

            enum ReferenceFrame
            {
                WORLD = 0, //This is spatial in world frame
                LOCAL = 1, //This is spatial in local frame
                LOCAL_WORLD_ALIGNED = 2 //This is classical in world frame
            };

            typedef std::weak_ptr<Model> WeakPtr;
            typedef std::shared_ptr<Model> Ptr;
            typedef std::unique_ptr<Model> UniquePtr;

            Model();

            Model(std::string _urdf_path, bool add_floating_jnt = false);

            bool was_model_init_ok();

            void update();
            void update_forward_kin();
            void update_frames_forward_kin(); // updates the joint placements and spatial
            // velocities according to the current joint configuration and velocity
            void update_B();
            void update_B_inv();
            void update_C();
            void update_b();
            void update_g();
            void update_tau_rnea();
            void update_p();

            void set_neutral();
            void set_random();
            void set_q(VectorXd& q);
            void set_v(VectorXd& v);
            void set_a(VectorXd& a);
            void set_tau(VectorXd& tau);

            void get_state(VectorXd& q, VectorXd& v, VectorXd& a,
                           VectorXd& tau);

            void get_q(VectorXd& q);
            void get_v(VectorXd& v);
            void get_a(VectorXd& a);

            void get_tau(VectorXd& tau);

            void get_B(MatrixXd& B); // joint-space inertia matrix
            void get_B_inv(MatrixXd& B_inv); // inverse of joint-space inertia matrix

            void get_C(MatrixXd& C); // Coriolis-matrix
            void get_g(VectorXd& g);

            void get_p(VectorXd& p); // generalized momentum (B * v)
            void get_b(VectorXd& b); // bias vector (C * v)
            void get_jac(std::string& frame_name,
                         utils_defs::SpatialJac& J,
                         ReferenceFrame ref = ReferenceFrame::LOCAL_WORLD_ALIGNED);
            void get_jac_dot(std::string& frame_name,
                         utils_defs::SpatialJacDot& J_dot,
                         ReferenceFrame ref = ReferenceFrame::LOCAL_WORLD_ALIGNED); // time derivative of Jacobian

            void get_frame_pose(std::string& frame_name,
                                utils_defs::PosVec3D& position, utils_defs::RotMat3D& rotation);
            void get_frame_pose(std::string& frame_name,
                                Affine3d& pose);
            void get_frame_vel(std::string& frame_name,
                               utils_defs::Twist& vel,
                               ReferenceFrame ref = ReferenceFrame::LOCAL_WORLD_ALIGNED);
            void get_frame_vel(std::string& frame_name,
                               utils_defs::LinVel& lin_vel, utils_defs::AngVel& omega,
                               ReferenceFrame ref = ReferenceFrame::LOCAL_WORLD_ALIGNED);

            void get_robot_mass(double& mass);

            int get_nq();
            int get_nv();
            int get_jnt_number();

            bool frame_exists(std::string& framename);

            std::string get_urdf_path();

            void get_jnt_lim(VectorXd& q_min, VectorXd& q_max);

            std::vector<std::string> get_jnt_names();

        private:

            bool _pin_model_init_ok = false;

            bool _is_floating_base = false;

            std::string _urdf_path;

            int _nq = 0, _nv = 0, _n_jnts = 0;

            double mass = -1.0;

            VectorXd _q_min, _q_max;

            VectorXd _q_min_aux, _q_max_aux; // these are only used internally
            // to avoid using local variables and when setting random configurations
            // for floating base robots

            VectorXd _effort_limits, _vel_limits;

            VectorXd _q, _v, _a;

            MatrixXd _B, _B_inv, _C;
            VectorXd _g, _tau,
                    _p, _b;

            std::vector<std::string> _jnt_names;

            std::vector<int> _nqs, _nvs;

            double _mass = -1.0;

            pinocchio::Model _pin_model;
            pinocchio::Data _pin_data;

            void update_all();
            void B(); // joint-space inertia matrix
            void B_inv(); // inverse of joint-space inertia matrix
            void C(); // Coriolis
            void g();
            void p(); // joint space momentum of the system
            void b(); // bias forces (C * v)

            void jacobian(std::string& frame_name, Model::ReferenceFrame ref,
                          utils_defs::SpatialJac& J);

            void jacobian_dot(std::string& frame_name, Model::ReferenceFrame ref,
                          utils_defs::SpatialJacDot& J_dot);

            void rnea(); // The Recursive Newton-Euler algorithm.
            // It computes the inverse dynamics, aka the joint torques \\
            // according to the current state of the system and the desired
            // joint accelerations.
            void aba(); // The Articulated-Body algorithm.
            // It computes the forward dynamics, aka the joint accelerations
            // given the current state and actuation of the model.

            void crba(); // Computes the upper triangular part of the joint space inertia matrix M by
            // using the Composite Rigid Body Algorithm

            void centroidal_dyn();

            void kin_energy();

            void pot_energy();

    };
}

#endif // MODEL_INTERFACE_HPP
