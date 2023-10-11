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
#include <gtest/gtest.h>
#include <chrono>

#include "tests/generated/cmake_config.h"

#include "model_interface.hpp"
#include "cartesian_imp_utils.hpp"
#include "typedefs.hpp"

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace ModelInterface;
using namespace Eigen;
using namespace CartesianImpUtils;

std::string urdf_path;
std::string urdf_path_floating;

IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

namespace
{
    template <typename Func>
    double measure_sec(Func f)
    {
        auto tic = std::chrono::high_resolution_clock::now();

        f();

        auto toc = std::chrono::high_resolution_clock::now();

        return std::chrono::duration_cast<std::chrono::nanoseconds>(toc-tic).count()*1e-9;
    }
}

class TestCartImp: public ::testing::Test {


    protected:

         TestCartImp(){

         }

         virtual ~TestCartImp() {
         }

         virtual void SetUp() {

         }

         virtual void TearDown() {
         }


};

TEST_F(TestCartImp, test_cart_cntrl)
{
    std::string tip_framename = "tip1";
    std::string base_link_frame_name = "base_link";
    std::string test_rig_frame_name = "test_rig";

    CartesianTask::CartTask chi_ref, chi_meas;
    CartesianTask::CartTaskDot chi_dot_ref, chi_dot_meas;
    CartesianTask::CartTaskDdot chi_ddot_ref, chi_ddot_meas;

    Model::Ptr model_ptr(new Model(urdf_path_floating));
    CartesianTask::Ptr task_ptr(new CartesianTask());
    CartesianImpController::Ptr cart_controller(new CartesianImpController(model_ptr,
                                                                           task_ptr,
                                                                           tip_framename));
    chi_ref.pos << 0.0, 0.0, 0.0;
    chi_meas.pos << 1.0, 1.0, 1.0;
    chi_dot_ref << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    chi_dot_meas << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    chi_ddot_ref << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    chi_ddot_meas << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    task_ptr->update_ref(chi_ref, chi_dot_ref, chi_ddot_ref);
    task_ptr->update_meas(chi_meas, chi_dot_meas, chi_ddot_meas);

    VectorXd q, v, a, tau;
    model_ptr->get_state(q, v, a, tau);
    Eigen::VectorXd base_position = Eigen::VectorXd::Zero(3);
    Eigen::VectorXd base_orient = Eigen::VectorXd::Zero(4);
    base_position << 0.0, 0.0, 0.0; // z height doesn't matter
    base_orient << 0.0, 0.0, 0.0, 1.0; // qx, qy, qz, w --> pinocchio follows Eigen's convention
    q << base_position, base_orient, 1.0, -1.0; // leg vertical with horizontal base link
    model_ptr->set_q(q);
    model_ptr->set_v(v);
    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    VectorXd tau_d = VectorXd::Zero(tau.size());
    VectorXd tau_c = VectorXd::Zero(tau.size());

    utils_defs::CartStiffVect stiffness_vect;
    stiffness_vect << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    cart_controller->set_cart_impedance(stiffness_vect);

    cart_controller->update(tau_d, tau_c);

    utils_defs::Wrench f_star = cart_controller->f_star();
    utils_defs::JacRightPseudoInv J_rpi = cart_controller->J_rps_w();
    utils_defs::CartInertiaMat Lambda_inv = cart_controller->Lambda_inv();;
    utils_defs::CartInertiaMat Lambda = cart_controller->Lambda();
    VectorXd tau_cmd = cart_controller->tau_cmd();

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;

    std::cout << "** f_star: \n" << f_star.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J_rpi: \n" << J_rpi.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Lambda_inv: \n" << Lambda_inv.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Lambda: \n" << Lambda.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tau_cmd: \n" << tau_cmd.format(CleanFmt) << "\n " << std::endl;

}

TEST_F(TestCartImp, test_impact_severity_ratio)
{
    std::string tip_framename = "tip1";

    CartesianTask::CartTask chi_ref, chi_meas;
    CartesianTask::CartTaskDot chi_dot_ref, chi_dot_meas;
    CartesianTask::CartTaskDdot chi_ddot_ref, chi_ddot_meas;

    Model::Ptr model_ptr(new Model(urdf_path));
    CartesianTask::Ptr task_ptr(new CartesianTask());
    CartesianImpController::Ptr cart_controller(new CartesianImpController(model_ptr,
                                                                           task_ptr,
                                                                           tip_framename));
    chi_ref.pos << 0.0, 0.0, 0.0;
    chi_meas.pos << 1.0, 1.0, 1.0;
    chi_dot_ref << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    chi_dot_meas << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    chi_ddot_ref << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    chi_ddot_meas << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    task_ptr->update_ref(chi_ref, chi_dot_ref, chi_ddot_ref);
    task_ptr->update_meas(chi_meas, chi_dot_meas, chi_ddot_meas);

    VectorXd q, v, a, tau;
    model_ptr->set_neutral();
    model_ptr->get_state(q, v, a, tau);
    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    VectorXd tau_d = VectorXd::Zero(tau.size());
    VectorXd tau_c = VectorXd::Zero(tau.size());

    utils_defs::CartStiffVect stiffness_vect;
    stiffness_vect << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    cart_controller->set_cart_impedance(stiffness_vect);

    cart_controller->update(tau_d, tau_c);

    utils_defs::JacRightPseudoInv J_rpi = cart_controller->J_rps_w();
    utils_defs::CartInertiaMat Lambda_inv = cart_controller->Lambda_inv();

    double k_n = 0.0; // normal restitution coefficient
    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    double rho = 0.0;

    int n_conf_test = 19;
    VectorXd q_low = VectorXd::Zero(q.size());
    VectorXd q_high = VectorXd::Zero(q.size());

    double pi = 3.1416;
    q_low <<0, -0.5 * pi, -0.5 * pi;
    q_high <<0, 0.5 * pi, 0.5* pi;
    MatrixXd Q = MatrixXd::Zero(q_low.size(), n_conf_test);
    VectorXd Rho = VectorXd::Zero(n_conf_test);
    MatrixXd LAMBDA = MatrixXd::Zero(Lambda_inv.rows(), Lambda_inv.cols() * n_conf_test);
    VectorXd interval = 1/((double)(n_conf_test - 1)) * (q_high - q_low);

    double mass = 0.0;
    model_ptr->get_robot_mass(mass);
    for(int i = 0; i < n_conf_test; i++)
    {
      Q.block(0, i, Q.rows(), 1) = q_low + i * interval;

      model_ptr->set_q(Q.block(0, i, Q.rows(), 1));
      model_ptr->update();

      cart_controller->update(tau_d, tau_c);

      Lambda_inv = cart_controller->Lambda_inv();

      LAMBDA.block(0, i * (Lambda_inv.cols()), Lambda_inv.rows(), Lambda_inv.cols()) = Lambda_inv;

      Rho(i) =  (1 + k_n)/ abs(Lambda_inv(2, 2));

    }

    std::cout << "** J_rpi: \n" << J_rpi.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Q: \n" << Q.format(CleanFmt) << "\n " << std::endl;

    for(int i = 0; i < n_conf_test; i++)
    {
        std::cout << "** Lambda_inv "  << i << " \n"<< LAMBDA.block(0, i * (Lambda_inv.cols()), Lambda_inv.rows(), Lambda_inv.cols()).format(CleanFmt) << "\n " << std::endl;
    }

    std::cout << "** Rho: \n" << Rho.format(CleanFmt) << "\n " << std::endl;

    std::cout << "** 1/robot_mass: \n" << 1.0/mass << "\n " << std::endl;

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string urdf_name = "awesome_leg_floating";
    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_fullpath = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
    urdf_path = URDF_PATH + std::string("awesome_leg") + std::string(".urdf");
    urdf_path_floating = URDF_PATH + std::string("awesome_leg_floating") + std::string(".urdf");

    return RUN_ALL_TESTS();
}

