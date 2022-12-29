#ifndef CARTESIAN_IMP_UTILS_HPP
#define CARTESIAN_IMP_UTILS_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "model_interface.hpp"

using namespace ModelInterface;

namespace CartesianImpUtils
{

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

