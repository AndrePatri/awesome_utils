# Copyright (C) 2023  Andrea Patrizi (AndrePatri, andreapatrizi1b6e6@gmail.com)
# 
# This file is part of awesome_utils and distributed under the General Public License version 2 license.
# 
# awesome_utils is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# 
# awesome_utils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with awesome_utils.  If not, see <http://www.gnu.org/licenses/>.
# 
#!/usr/bin/env python3

import awesome_utils.awesome_pyutils as utils
import numpy as np
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Test file for the python bindings for the Awesome Utils cpp library')

    parser.add_argument('--urdf_path', '-p', type = str, help = 'path to test URDF')

    args = parser.parse_args()

    urdf_path = args.urdf_path

    add_floating_joint = False

    frame_name = "tip1"

    model = utils.Model(urdf_path, add_floating_joint)

    mass = model.get_robot_mass()

    nq = model.get_nq()

    nv = model.get_nv()

    n_jnts = model.get_jnt_number()

    q = np.zeros(nq)

    v = np.zeros(nv)

    model.set_q(q)

    model.set_v(v)

    model.update()

    print("\n \\\****  Debug print ****\\\ \n")
    print("Joint names: ", model.get_jnt_names())
    print("Model mass: " + str(mass) + " Kg")
    print("nq: " + str(nq))
    print("nv: " + str(nv))
    print("n_jnts: " + str(n_jnts))
    print("p: " + str(model.get_p()))
    print("b: " + str(model.get_b()))
    print("B: " + str(model.get_B()))
    print("C: " + str(model.get_C()))
    print("g: " + str(model.get_g()))

    print("J: " + str(model.get_J(frame_name, utils.ReferenceFrame.LOCAL_WORLD_ALIGNED)))
    print("J_dot: " + str(model.get_J_dot(frame_name, utils.ReferenceFrame.LOCAL_WORLD_ALIGNED)))


