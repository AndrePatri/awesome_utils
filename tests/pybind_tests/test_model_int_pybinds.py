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

	model = utils.Model(urdf_path, add_floating_joint)

	mass = model.get_robot_mass()

	nq = model.get_nq()
	nv = model.get_nv()
	n_jnts = model.get_jnt_number()

	q = np.zeros(nq)
	v = np.zeros(nv)
	a = np.zeros(nv)
	tau = np.zeros(nv)

	model.set_q(q)
	model.set_v(v)
	model.set_a(a)

