#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

namespace utils_defs
{
    typedef Matrix<double, 3, 3> Mat3D;

    typedef Matrix<double, 3, 3> RotMat3D;
    typedef Affine3d Affine3D;
    typedef Matrix<double, 3, 1> PosVec3D;
    typedef Matrix<double, 6, 1> Twist;
    typedef Matrix<double, 3, 1> LinVel;
    typedef Matrix<double, 3, 1> AngVel;
    typedef Matrix<double, 6, 1> Wrench;
    typedef Matrix<double, 3, 1> Force3D;
    typedef Matrix<double, 3, 1> Torque3D;

    typedef Matrix<double, 6, -1> SpatialJac;
    typedef Matrix<double, -1, 6> SpatialJacT;
    typedef Matrix<double, 6, -1> SpatialJacDot;

    typedef Matrix<double, -1, 6> JacRightPseudoInv;
    typedef Matrix<double, 6, -1> JacLeftPseudoInv;

    typedef Matrix<double, 6, 1> CartVect;
    typedef Matrix<double, 6, 6> CartInertiaMat;
    typedef Matrix<double, 6, 6> CartStiffMat;
    typedef Matrix<double, 6, 6> CartDampMat;
    typedef Matrix<double, 6, 1> CartStiffVect;
    typedef Matrix<double, 6, 1> CartDampVect;

}

#endif // TYPEDEFS_HPP
