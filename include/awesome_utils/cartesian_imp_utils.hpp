#ifndef CARTESIAN_IMP_UTILS_HPP
#define CARTESIAN_IMP_UTILS_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "include/awesome_utils/typedefs.hpp"
#include "include/awesome_utils/model_interface.hpp"

using namespace ModelInterface;

namespace CartesianImpUtils
{

    class CartesianTask
    {
        public:

            typedef std::weak_ptr<CartesianTask> WeakPtr;
            typedef std::shared_ptr<CartesianTask> Ptr;
            typedef std::unique_ptr<CartesianTask> UniquePtr;

            struct CartTask
            {
                utils_defs::RotMat3D rot; // orientation part of the cartesian task

                utils_defs::PosVec3D pos;
            };

            typedef Matrix<double, 6, 1> CartTaskDot;
            typedef Matrix<double, 6, 1> CartTaskDdot;

            typedef Matrix<double, 6, 1> CartTaskErr;
            typedef Matrix<double, 6, 1> CartTaskDotErr;
            typedef Matrix<double, 6, 1> CartTaskDdotErr;

            CartesianTask();

            void update(CartTask chi_ref,
                        CartTaskDot chi_dot_ref,
                        CartTaskDdot chi_ddot_ref);

            void update(utils_defs::PosVec3D pos_ref, utils_defs::RotMat3D rot_ref,
                        CartTaskDot chi_dot_ref,
                        CartTaskDdot chi_ddot_ref);

            // compute task error between last set reference task and the input
            CartTaskErr compute_task_err(utils_defs::PosVec3D pos, utils_defs::RotMat3D rot);
            CartTaskErr compute_task_err(CartTask cart_task);

            CartTaskDotErr compute_task_dot_err(CartTaskDot cart_task_dot);

            CartTaskDdotErr compute_task_ddot_err(CartTaskDdot cart_task_ddot);

        private:

            CartTask _chi_ref;
            CartTaskDot _chi_dot_ref;
            CartTaskDdot _chi_ddot_ref;

            void set_chi_ref(CartTask chi_ref);
            void set_chi_ref(utils_defs::PosVec3D pos_ref, utils_defs::RotMat3D rot_ref);
            void set_chi_dot_ref(CartTaskDot chi_dot_ref);
            void set_chi_ddot_ref(CartTaskDdot chi_ddot_ref);

    };

    /**
    * @brief Class to perform cartesian impedance control.
    *
    */

    class CartesianImpController
    {
        public:

            typedef std::weak_ptr<CartesianImpController> WeakPtr;
            typedef std::shared_ptr<CartesianImpController> Ptr;
            typedef std::unique_ptr<CartesianImpController> UniquePtr;

            CartesianImpController();

            CartesianImpController(Model::Ptr model_ptr);

            void update();

        private:

            int _nq = 0, _nv = 0;

    };
}

#endif // CARTESIAN_IMP_UTILS_HPP

