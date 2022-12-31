#include "orientation_utils.hpp"

RotErr::RotErr3d RotErr::LogMap(utils_defs::RotMat3D R, utils_defs::RotMat3D R_ref)
{
    utils_defs::RotMat3D R_err = R_ref.transpose() * R; // orientation or actual frame w.r.t. target frame

    utils_defs::Mat3D S = (R_err - R_err.transpose())/ 2.0;  // logarithmic map

    RotErr3d r_err;

    r_err(0) = S(2, 1);
    r_err(1) = S(0, 2);
    r_err(2) = S(1, 0);

    r_err *= 1/std::sqrt(1 + R_err.trace());

    return r_err;

};
