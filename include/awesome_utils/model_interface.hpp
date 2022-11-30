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

namespace ModelInterface
{
    using namespace Eigen;

    class Model
    {
        public:

            Model();

            Model(std::string _urdf_path);

            bool was_model_init_ok();

        private:

            bool _pin_model_init_ok = false;

            std::string _urdf_path;

            int _nq = 0, _nv = 0, _n_jnts = 0;

            VectorXd _q_min, _q_max;

            VectorXd _q, _q_dot, _q_ddot;

            std::vector<std::string> _jnt_names;

            double _mass = -1.0;

            pinocchio::Model _pin_model;
            pinocchio::Data _pin_data;

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

            void ccrba();
            void fk();
            void CoM();
            void CoM_dot();
            void jac();
            void kin_energy();
            void pot_energy();

    };
}

#endif // MODEL_INTERFACE_HPP
