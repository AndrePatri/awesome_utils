#include <gtest/gtest.h>
#include <chrono>

#include "tests/generated/cmake_config.h"

#include "include/awesome_utils/model_interface.hpp"
#include "include/awesome_utils/cartesian_imp_utils.hpp"

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace ModelInterface;

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
    Model::Ptr model_ptr(new Model(urdf_path));

    std::string tip_framename = "tip1";
    std::string base_link_frame_name = "base_link";
    std::string test_rig_frame_name = "test_rig";

    double mass = -1.0;
    Eigen::VectorXd q, v, a, tau,
                    g, p, b;
    Eigen::MatrixXd B, C;
    utils_defs::SpatialJac J;
    utils_defs::SpatialJac J_dot;
    utils_defs::PosVec3D position;
    utils_defs::RotMat3D rotation;
    utils_defs::Twist vel;
    utils_defs::Affine3D pose;

    model_ptr->get_state(q, v, a, tau);
    model_ptr->set_q(q);
    model_ptr->set_v(v);
    model_ptr->set_a(a);
    model_ptr->set_tau(tau);

    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string urdf_name = "awesome_leg";
    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_fullpath = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
    urdf_path = urdf_fullpath;

    return RUN_ALL_TESTS();
}
