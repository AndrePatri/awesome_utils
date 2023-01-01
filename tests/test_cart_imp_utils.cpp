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

class TestModelInterface: public ::testing::Test {


    protected:

         TestModelInterface(){

         }

         virtual ~TestModelInterface() {
         }

         virtual void SetUp() {

         }

         virtual void TearDown() {
         }


};

TEST_F(TestModelInterface, test_cart_cntrl)
{
    std::string tip_framename = "tip1";
    std::string base_link_frame_name = "base_link";
    std::string test_rig_frame_name = "test_rig";

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
    model_ptr->get_state(q, v, a, tau);
    q << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, -1.0;
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

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string urdf_name = "awesome_leg_floating";
    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_fullpath = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
    urdf_path = urdf_fullpath;

    return RUN_ALL_TESTS();
}
