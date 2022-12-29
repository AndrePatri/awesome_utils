#include "include/awesome_utils/orientation_utils.hpp"

using namespace RotErr;

RotErr3D LogMap(Mat3D R, Mat3D R_ref)
{
    Mat3D R_err = R_ref.transpose() * R; // orientation or actual frame w.r.t. target frame

    Mat3D S = (R_err - R_err.transpose())/ 2.0;  // logarithmic map

    RotErr3D r_err;

    r_err(0) = S(2, 1);
    r_err(1) = S(0, 2);
    r_err(2) = S(1, 0);

    r_err *= 1/std::sqrt(1 + R_err.trace());

    return r_err;

};
