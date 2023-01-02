#ifndef CARTESIAN_IMP_UTILS_HPP
#define CARTESIAN_IMP_UTILS_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "typedefs.hpp"
#include "model_interface.hpp"

using namespace ModelInterface;
using namespace Eigen;

namespace CartesianImpUtils
{

    /**
    * @brief Class to represent a cartesian task.
    *
    * A cartesian task is specified by three reference quantities:
    *
    * - The task(positional), which is represented by a desired orientation and position
    *
    * - The derivative of the task, represented by a 6D vector made of translational vel and
    *   and rotational one
    *
    * - The double derivative of the task, which is also a 6D vector representing translational
    *   and angular acceleration
    *
    * The class holds a methods for computing the difference of the reference task wrt
    * to a given task (e.g. measured). The error for the velocity and acceleration parts are
    * computed by simply making the difference of the inputs, while for the rototranslational
    * representation of the task we make use of the logarithmic map, so that the resulting
    * error vector is a 6D vector, where the first three elements represent the translational
    * error, while the remaining represent the orientation error along the coordinated axis of the
    * referene orientation frame.
    */

    class CartesianTask
    {
        public:

            typedef std::weak_ptr<CartesianTask> WeakPtr;
            typedef std::shared_ptr<CartesianTask> Ptr;
            typedef std::unique_ptr<CartesianTask> UniquePtr;

            struct CartTask
            {
                utils_defs::RotMat3D rot = MatrixXd::Identity(3, 3); // orientation part of the cartesian task

                utils_defs::PosVec3D pos = VectorXd::Zero(3);
            };

            typedef Matrix<double, 6, 1> CartTaskDot;
            typedef Matrix<double, 6, 1> CartTaskDdot;

            typedef Matrix<double, 6, 1> CartTaskErr;
            typedef Matrix<double, 6, 1> CartTaskDotErr;
            typedef Matrix<double, 6, 1> CartTaskDdotErr;

            CartesianTask();

            void update_ref(CartTask chi_ref,
                        CartTaskDot chi_dot_ref,
                        CartTaskDdot chi_ddot_ref);

            void update_ref(utils_defs::PosVec3D pos_ref, utils_defs::RotMat3D rot_ref,
                        CartTaskDot chi_dot_ref,
                        CartTaskDdot chi_ddot_ref);

            void update_meas(CartTask chi_meas,
                        CartTaskDot chi_dot_meas,
                        CartTaskDdot chi_ddot_meas);

            void update_meas(utils_defs::PosVec3D pos_meas, utils_defs::RotMat3D rot_meas,
                        CartTaskDot chi_dot_meas,
                        CartTaskDdot chi_ddot_meas);

            // compute task errors between last set reference task and the input
            CartTaskErr task_err();
            CartTaskDotErr task_dot_err();
            CartTaskDdotErr task_ddot_err();
            CartTaskDdot task_ddot_ref();

            void get_task(CartTask& task_ref, CartTask& task_meas);
            void get_task_dot(CartTaskDot& task_dot_ref, CartTaskDot& task_dot_meas);
            void get_task_ddot(CartTaskDdot& task_ddot_ref, CartTaskDdot& task_ddot_meas);

        private:

            CartTask _chi_ref, _chi_meas;
            CartTaskDot _chi_dot_ref, _chi_dot_meas;
            CartTaskDdot _chi_ddot_ref, _chi_ddot_meas;

            void set_chi_ref(CartTask chi_ref);
            void set_chi_ref(utils_defs::PosVec3D pos_ref, utils_defs::RotMat3D rot_ref);
            void set_chi_dot_ref(CartTaskDot chi_dot_ref);
            void set_chi_ddot_ref(CartTaskDdot chi_ddot_ref);

            void set_chi_meas(CartTask chi_meas);
            void set_chi_meas(utils_defs::PosVec3D pos_meas, utils_defs::RotMat3D rot_meas);
            void set_chi_dot_meas(CartTaskDot chi_dot_meas);
            void set_chi_ddot_meas(CartTaskDdot chi_ddot_meas);

    };

    /**
    * @brief Class to perform model-based cartesian impedance control at
    * a specific robot frame.
    *
    * Let's look at the derivation of the control law.
    *
    * Given a task Chi, it's dynamics is given by
    *
    * (1.0) Chi_ddot = J * a + J_dot * v
    *
    * (Chi_dot = J_dot * v)
    *
    * where a can be computed with the forward dynamics.
    *
    * Specifically,
    *
    *   B(q) * a + C(q, v) * v + g(q) = tau_cmd + tau_d + tau_c + tau_i
    * ->B(q) * a + C(q, v) * q_dot + g(q) = tau_cmd + tau_d + sum_j^{nc} {J_j^T * f_j} + J^T * f_i
    *
    * where f_i is the interaction force at the cartesian control frame,
    * expressed in the frame coherent with the jacobian J.
    *
    * tau_d is a "disturbance" torque vector which, for example,
    * can include friction torques
    *
    * tau_c is the sum of the contact forces reflected at the joint level (f_i excluded)
    *
    * Given the invertibility of B, the forward dynamics is computed with
    *
    * (1.1) a = B_inv * (J^T * f + tau_cmd + tau_d - C * v - g)
    *
    * Substitution of (1.1) into (1.0) gives, after a little rearranging
    *
    *   Chi_ddot = J * B_inv * (J^T * f + tau_cmd + tau_d - C * v - g) + J_dot * v
    *-->Chi_ddot = (J * B_inv * J^T) * f + J * B_inv * (tau_cmd + tau_d - C * v - g) + J_dot * v
    *-->Chi_ddot = (J * B_inv * J^T) * f + J * B_inv * (tau_cmd + tau_d + tau_c - g) + (J_dot - J * B_inv * C) * v
    *
    * Let define the cartesian inertia matrix as
    *
    * (1.2) Lambda:= (J * B_inv * J^T)^{-1}
    *
    * The dynamics of the task Chi can hence be written as
    *
    *    Lambda * Chi_ddot = f_i + Lambda * J * B_inv * (tau_cmd + tau_d - g) + Lambda * (J_dot - J * B_inv * C) * v
    *--> Lambda * Chi_ddot = f_i + Lambda * J * B_inv * (tau_cmd + tau_d - g) - Lambda * (J * B_inv * C - J_dot) * v
    *
    *--> Lambda * Chi_ddot + Lambda * (J * B_inv * C - J_dot) * v + Lambda * J * B_inv * g = f_i + Lambda * J * B_inv * (tau_cmd + tau_d)
    *
    * we define the right weighted pseudo inverse of J as
    *
    * (1.4) J_rps_w := Lambda * J * B_inv
    *
    * The final version of the task dynamics is
    *
    * (1.5) Lambda * Chi_ddot + Lambda * (J * B_inv * C - J_dot) * v + J_rps_w * g = f_i + J_rps_w * (tau_cmd + tau_d + tau_c)
    *
    * Suppose we choose tau_cmd as
    *
    * tau_cmd = - tau_d - tau_c + tau^ =  - tau_d - tau_c + J^T * f^
    *
    * (1.5) becomes
    *
    * (1.5.1) Lambda * Chi_ddot + Lambda * (J * B_inv * C - J_dot) * v + J_rps_w * g =
    * = f_i + f_star
    *
    * We'd like to impose the MSD-like error dynamics given by
    *
    * (1.6) Lambda_ref * Chi_err_ddot + K_D^ref * Chi_err_dot + K_P^ref * Chi_err = f_i
    *
    * We extract Chi_ddot from (1.6), plug it in (1.5.1) and extract f^ :
    *
    * (1.7) f^ = (Lambda * Lambda_ref^{-1} - I) * f_i + Lambda * Chi_ddot_ref + h(q, v) + J_rps_w * g +
    *           - Lambda * Lambda_ref^{-1} * (K_D_ref * Chi_err_dot + K_P^ref * Chi_err)
    *
    * where h(q, v) := (J_rps_w * C - Lambda * J_dot) * v
    *
    * In particular, if we choose Lambda_ref = Lambda the (1.7) simplifies to
    *
    * (1.7.1) f^ = Lambda * Chi_ddot_ref + h(q, v) + J_rps_w * g +
    *              - (K_D_ref * Chi_err_dot + K_P^ref * Chi_err)
    *
    * which has the advantage of not needing a measurement of f_i
    *
    * Under the hypothesis that Lambda is diagonally dominant and the gain matrices
    * K_D and K_P are chosen to be diagonal, then we have approximately a 6 decoupled oscillators.
    * We can simply choose the gains to that K_P(j, j) = (K_D(j, j))^2 / (4 * Lambda(j, j))
    *
    * For a more in depth analysis refer to
    *
    * ## F. Angelini et al., "Online Optimal Impedance Planning for Legged Robots," 2019 IEEE/RSJ
    * International Conference on Intelligent Robots and Systems (IROS), 2019, pp. 6028-6035,
    * doi: 10.1109/IROS40897.2019.8967696. ##
    *
    */

    class CartesianImpController
    {
        public:

            typedef std::weak_ptr<CartesianImpController> WeakPtr;
            typedef std::shared_ptr<CartesianImpController> Ptr;
            typedef std::unique_ptr<CartesianImpController> UniquePtr;

            CartesianImpController();

            CartesianImpController(Model::Ptr model_ptr,
                                   CartesianTask::Ptr cart_task);
            CartesianImpController(Model::Ptr model_ptr,
                                   CartesianTask::Ptr cart_task,
                                   std::string cart_cntrl_framename);

            void update(VectorXd tau_d,
                        VectorXd tau_c);

            void update(std::string cart_cntrl_framename,
                        VectorXd tau_d,
                        VectorXd tau_c);

            void set_cart_impedance(utils_defs::CartStiffMat stifness_mat,
                                    utils_defs::CartDampMat damping_mat);
            void set_cart_impedance(utils_defs::CartStiffVect stifness_vect,
                                    utils_defs::CartDampVect damping_vect);

            void set_cart_impedance(utils_defs::CartStiffVect stifness_vect); // damping is
            // computed to give an (approximately) critically damped response

            VectorXd h();
            utils_defs::Wrench f_star();
            utils_defs::JacRightPseudoInv J_rps_w();
            utils_defs::CartInertiaMat Lambda();
            utils_defs::CartInertiaMat Lambda_inv();

            VectorXd tau_cmd(); // cartesian impedance control torque

        private:

            Model::Ptr _model_ptr;

            int _nq = 0, _nv = 0;

            std::string _cart_cntrl_framename = "";

            bool _was_cntrl_framename_set = false;

            bool _auto_critical_damp = true; // whether to employ a
            // simple approximation of the critical damping gains

            CartesianTask::Ptr _cart_task;

            utils_defs::CartStiffMat _cart_stiff;
            utils_defs::CartDampMat _cart_damp;

            utils_defs::CartStiffVect _cart_stiff_vect;
            utils_defs::CartDampVect _cart_damp_vect;

            utils_defs::CartInertiaMat _Lambda,
                           _Lambda_inv; // cartesian impedance matrix

            utils_defs::SpatialJacDot _J_dot;
            utils_defs::SpatialJac _J;

            utils_defs::JacRightPseudoInv _J_rps_w;

            utils_defs::Wrench _f_star;

            VectorXd _q, _v;
            MatrixXd _B_inv,
                     _C; // inverse of generalized inertia matrix and bias Matrix
            VectorXd _g; // joint-space gravitational vector

            utils_defs::CartVect _h; // auxiliary vector

            VectorXd _tau_d; // disturbance torques on the joints (e.g. friction torques)
            VectorXd _tau_c; // contact torques
            VectorXd _tau_cmd; // cartesian impedance control torque

            void map_impedance_vect2mat();

            void compute_quantities();

            void compute_critically_damped_gains();

            void compute_J_rps_w();

            void compute_lambda_inv(); // exposed here for possible
            // external usage
            void compute_lambda();

            void compute_tau_cmd();

            void compute_f_star();

            void compute_h();

            void update_internal_states(VectorXd tau_d, VectorXd tau_c);

    };
}

#endif // CARTESIAN_IMP_UTILS_HPP

