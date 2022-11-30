#ifndef MODEL_INTERFACE_HPP
#define MODEL_INTERFACE_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <string>
#include <vector>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/geometry.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"

namespace ModelInterface
{
    using namespace Eigen;

    class Model
    {
        public:

            Model();

            Model(std::string _urdf_path);

        private:

            std::string _urdf_path;

            int _nq, _nv;

            VectorXd _q_min, _q_max;

            std::vector<std::string> _jnt_names;

            double _mass;

            pinocchio::Model _pin_model;

            void rnea();
            void centroidal_dyn();
            void crba();
            void ccrba();
            void fk();
            void CoM();
            void jac();
            void kin_energy();
            void pot_energy();
            void aba();

    };
}

#endif // MODEL_INTERFACE_HPP
