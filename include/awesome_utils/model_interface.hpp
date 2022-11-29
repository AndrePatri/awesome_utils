#ifndef MODEL_INTERFACE_HPP
#define MODEL_INTERFACE_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <string>

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/geometry.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"

namespace ModelInterface
{
    class Model
    {
        public:

            Model();

        private:

            std::string _urdf_path;


    };
}

#endif // MODEL_INTERFACE_HPP
