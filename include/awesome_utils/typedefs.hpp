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

    namespace Colors {
    const std::string kRed = "\033[0;31m";
    const std::string kBoldRed = "\033[1;31m";
    const std::string kGreen = "\033[0;32m";
    const std::string kBoldGreen = "\033[1;32m";
    const std::string kYellow = "\033[0;33m";
    const std::string kBoldYellow = "\033[1;33m";
    const std::string kBlue = "\033[0;34m";
    const std::string kBoldBlue = "\033[1;34m";
    const std::string kMagenta = "\033[0;35m";
    const std::string kBoldMagenta = "\033[1;35m";
    const std::string kCyan = "\033[0;36m";
    const std::string kBoldCyan = "\033[1;36m";
    const std::string kEndl = "\033[0m";
    };

}

#endif // TYPEDEFS_HPP
