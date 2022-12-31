#ifndef ORIENTATION_UTILS_HPP
#define ORIENTATION_UTILS_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cmath>

#include "typedefs.hpp"

using namespace Eigen;

namespace RotErr
{

    typedef Matrix<double, 3, 1> RotErr3d;

    /**
    * @brief Class to compute an orientation error between two frames
    * using the logarithmic map (i.e. a map between a rotation matrix and
    * a 3D vector).
    * This formulation has some benefits on different orientation error
    * formulation. In particular, it decouples the error along the 3 coordinate
    * axis, which is very useful when performing, e.g. cartesian impedance control,
    * or when performing TO with pose constraints.
    *
    * The logarithmic map has the following properties:
    * - It is an isomorphism, meaning that it is a one-to-one
    *   correspondence between rotations and their corresponding
    *   3-dimensional vector representations.
    *   This means that every rotation has a unique 3-dimensional
    *   vector representation, and vice versa.
    * - It is smooth, meaning that it is a continuous function.
    *   This means that small rotations are represented by small
    *   3-dimensional vectors, and large rotations are represented
    *   by large 3-dimensional vectors.
    * - It is not a linear function.
    *
    */

    RotErr3d LogMap(utils_defs::RotMat3D R, utils_defs::RotMat3D R_ref);

}

#endif // ORIENTATION_UTILS_HPP
