// Copyright (C) 2023  Andrea Patrizi (AndrePatri, andreapatrizi1b6e6@gmail.com)
// 
// This file is part of awesome_utils and distributed under the General Public License version 2 license.
// 
// awesome_utils is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// awesome_utils is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with awesome_utils.  If not, see <http://www.gnu.org/licenses/>.
// 
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
