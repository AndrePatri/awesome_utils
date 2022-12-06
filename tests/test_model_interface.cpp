#include <gtest/gtest.h>
#include <chrono>

#include "tests/generated/cmake_config.h"

#include "awesome_utils/model_interface.hpp"

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


TEST_F(TestModelInterface, load_model)
{
    Model::Ptr model_ptr(new Model(urdf_path));

    ASSERT_TRUE(model_ptr->was_model_init_ok());

    auto jnt_names = model_ptr->get_jnt_names();

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Joint names: **"  << std::endl;
    for (std::string i: jnt_names)
        std::cout << "--> " << i << std::endl;
    std::cout << "** joint number: " <<model_ptr->get_jnt_number() << "\n " << std::endl;
    std::cout << "** nq: " << model_ptr->get_nq() << std::endl;
    std::cout << "** nv: " << model_ptr->get_nv() << std::endl;

}

TEST_F(TestModelInterface, compute_quantities)
{
    Model::Ptr model_ptr(new Model(urdf_path));

    std::string tip_framename = "tip1";
    std::string base_link_frame_name = "base_link";
    std::string test_rig_frame_name = "test_rig";

    Eigen::VectorXd q, v, a, tau,
                    g, p, b;
    Eigen::MatrixXd B, C;
    Model::SpatialJac J;
    Model::PosVec3D position;
    Model::RotMat3D rotation;
    Model::Twist vel;
    Model::Affine3D pose;

    model_ptr->get_state(q, v, a, tau);
    model_ptr->set_q(q);
    model_ptr->set_v(v);
    model_ptr->set_a(a);
    model_ptr->set_tau(tau);

    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    model_ptr->get_B(B);
    model_ptr->get_C(C);
    model_ptr->get_g(g);
    model_ptr->get_b(b);
    model_ptr->get_p(p);
    model_ptr->get_jac(tip_framename,
                  J,
                  Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);

    model_ptr->get_frame_pose(tip_framename,
                              position, rotation);

    model_ptr->get_frame_pose(test_rig_frame_name,
                              pose);

    model_ptr->get_frame_vel(tip_framename,
                             vel);

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** B: \n" << B.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** C: \n" << C.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** g: \n" << g.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** b: \n" << b.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** p: \n" << p.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J (q_dot -> " << tip_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** frame position: \n" << position.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** frame rotation matrix: \n" << rotation.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** frame position from Affine3D: \n" << pose.translation().format(CleanFmt) << "\n " << std::endl;
    std::cout << "** frame orientation from Affine3D: \n" << pose.rotation().format(CleanFmt) << "\n " << std::endl;

    std::cout << "** frame generalized velocity: \n" << vel.format(CleanFmt) << "\n " << std::endl;

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string urdf_name = "awesome_leg";
    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_fullpath = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
    urdf_path = urdf_fullpath;

    return RUN_ALL_TESTS();
}
