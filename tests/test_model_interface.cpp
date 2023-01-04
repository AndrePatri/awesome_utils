#include <gtest/gtest.h>
#include <chrono>

#include "tests/generated/cmake_config.h"

#include "model_interface.hpp"

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace ModelInterface;

std::string urdf_path_fixed;
std::string urdf_path_floating;

std::string urdf_path_anymal;
std::string urdf_path_quadruped;

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

TEST_F(TestModelInterface, load_fixed_base_leg)
{
    Model::Ptr model_ptr(new Model(urdf_path_fixed));

    ASSERT_TRUE(model_ptr->was_model_init_ok());

    auto jnt_names = model_ptr->get_jnt_names();

    std::cout << "** FIXED BASE LEG DEBUG PRINTS**\n"  << std::endl;
    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Joint names: **"  << std::endl;
    for (std::string i: jnt_names)
        std::cout << "--> " << i << std::endl;
    std::cout << "** joint number: " <<model_ptr->get_jnt_number() << "\n " << std::endl;
    std::cout << "** nq: " << model_ptr->get_nq() << std::endl;
    std::cout << "** nv: " << model_ptr->get_nv() << std::endl;

    std::string tip_framename = "tip1";
    std::string base_link_frame_name = "base_link";
    std::string test_rig_frame_name = "test_rig";

    double mass = -1.0;
    Eigen::VectorXd q, v, a, tau,
                    g, p, b;
    Eigen::MatrixXd B, B_inv, C;
    utils_defs::SpatialJac J;
    utils_defs::SpatialJac J_dot;
    utils_defs::PosVec3D position;
    utils_defs::RotMat3D rotation;
    utils_defs::Twist vel;
    utils_defs::Affine3D pose;

    model_ptr->set_neutral(); // sets q to a neutral configuration vector
    model_ptr->get_state(q, v, a, tau);
//    model_ptr->set_q(q);
//    model_ptr->set_v(v);
//    model_ptr->set_a(a);
//    model_ptr->set_tau(tau);

    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    model_ptr->get_robot_mass(mass);
    model_ptr->get_B(B);
    model_ptr->get_B_inv(B_inv);
    model_ptr->get_C(C);
    model_ptr->get_g(g);
    model_ptr->get_b(b);
    model_ptr->get_p(p);
    model_ptr->get_jac(tip_framename,
                  J,
                  Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);
    model_ptr->get_jac_dot(tip_framename, J_dot, Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);

    model_ptr->get_frame_pose(tip_framename,
                              position, rotation);

    model_ptr->get_frame_pose(test_rig_frame_name,
                              pose);

    model_ptr->get_frame_vel(tip_framename,
                             vel);

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Robot mass: \n" << mass << "\n " << std::endl;
    std::cout << "** q: \n" << q.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** v: \n" << v.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** a: \n" << a.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tau: \n" << tau.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B: \n" << B.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B_inv: \n" << B_inv.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** C: \n" << C.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** g: \n" << g.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** b: \n" << b.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** p: \n" << p.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J (q_dot -> " << tip_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J_dot (q_dot -> " << tip_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J_dot.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame position: \n" << position.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame rotation matrix: \n" << rotation.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame position from Affine3D: \n" << pose.translation().format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame orientation from Affine3D: \n" << pose.rotation().format(CleanFmt) << "\n " << std::endl;

    std::cout << "** tip frame generalized velocity: \n" << vel.format(CleanFmt) << "\n " << std::endl;

}

TEST_F(TestModelInterface, load_floating_base_leg)
{
    Model::Ptr model_ptr(new Model(urdf_path_floating));

    ASSERT_TRUE(model_ptr->was_model_init_ok());

    auto jnt_names = model_ptr->get_jnt_names();

    std::cout << "** FLOATING BASE LEG DEBUG PRINTS**\n"  << std::endl;
    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Joint names: **"  << std::endl;
    for (std::string i: jnt_names)
        std::cout << "--> " << i << std::endl;
    std::cout << "** joint number: " <<model_ptr->get_jnt_number() << "\n " << std::endl;
    std::cout << "** nq: " << model_ptr->get_nq() << std::endl;
    std::cout << "** nv: " << model_ptr->get_nv() << std::endl;

    std::string tip_framename = "tip1";
    std::string base_link_framename = "base_link";

    double mass = -1.0;
    Eigen::VectorXd q, v, a, tau,
                    g, p, b;
    Eigen::MatrixXd B, B_inv, C;
    utils_defs::SpatialJac J_tip, J_base;
    utils_defs::SpatialJac J_dot_tip, J_dot_base;
    utils_defs::PosVec3D tip_position, base_position;
    utils_defs::RotMat3D tip_rotation, base_rotation;
    utils_defs::Twist tip_vel, base_vel;
    utils_defs::Affine3D tip_pose, base_pose;

    model_ptr->set_neutral(); // sets q to a neutral configuration vector
    model_ptr->get_state(q, v, a, tau);
//    model_ptr->set_q(q);
//    model_ptr->set_v(v);
//    model_ptr->set_a(a);
//    model_ptr->set_tau(tau);

    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    model_ptr->get_robot_mass(mass);
    model_ptr->get_B(B);
    model_ptr->get_B_inv(B_inv);
    model_ptr->get_C(C);
    model_ptr->get_g(g);
    model_ptr->get_b(b);
    model_ptr->get_p(p);
    model_ptr->get_jac(tip_framename,
                  J_tip,
                  Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);
    model_ptr->get_jac_dot(tip_framename, J_dot_tip, Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);
    model_ptr->get_jac(base_link_framename,
                  J_base,
                  Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);
    model_ptr->get_jac_dot(base_link_framename, J_dot_base, Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);

    model_ptr->get_frame_pose(tip_framename,
                              tip_position, tip_rotation);
    model_ptr->get_frame_pose(base_link_framename,
                              base_position, base_rotation);

    model_ptr->get_frame_vel(tip_framename,
                             tip_vel);
    model_ptr->get_frame_vel(base_link_framename,
                             base_vel);

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Robot mass: \n" << mass << "\n " << std::endl;
    std::cout << "** q: \n" << q.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** v: \n" << v.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** a: \n" << a.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tau: \n" << tau.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B: \n" << B.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B_inv: \n" << B_inv.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** C: \n" << C.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** g: \n" << g.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** b: \n" << b.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** p: \n" << p.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J (q_dot -> " << tip_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J_tip.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J_dot (q_dot -> " << tip_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J_dot_tip.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame position: \n" << tip_position.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame rotation matrix: \n" << tip_rotation.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame position from Affine3D: \n" << tip_pose.translation().format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tip frame orientation from Affine3D: \n" << tip_pose.rotation().format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J (q_dot -> " << base_link_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J_base.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J_dot (q_dot -> " << base_link_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J_dot_base.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** base frame position: \n" << base_position.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** base frame rotation matrix: \n" << base_rotation.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** base frame position from Affine3D: \n" << base_pose.translation().format(CleanFmt) << "\n " << std::endl;
    std::cout << "** base frame orientation from Affine3D: \n" << base_pose.rotation().format(CleanFmt) << "\n " << std::endl;
    std::cout << "** frame generalized velocity: \n" << base_vel.format(CleanFmt) << "\n " << std::endl;

}

TEST_F(TestModelInterface, load_floating_base_anymal)
{
    Model::Ptr model_ptr(new Model(urdf_path_anymal));

    ASSERT_TRUE(model_ptr->was_model_init_ok());

    auto jnt_names = model_ptr->get_jnt_names();

    std::cout << "** FLOATING BASE ANYMAL DEBUG PRINTS**\n"  << std::endl;
    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Joint names: **"  << std::endl;
    for (std::string i: jnt_names)
        std::cout << "--> " << i << std::endl;
    std::cout << "** joint number: " <<model_ptr->get_jnt_number() << "\n " << std::endl;
    std::cout << "** nq: " << model_ptr->get_nq() << std::endl;
    std::cout << "** nv: " << model_ptr->get_nv() << std::endl;

    double mass = -1.0;
    Eigen::VectorXd q, v, a, tau,
                    g, p, b;
    Eigen::MatrixXd B, B_inv, C;

    model_ptr->set_neutral(); // sets q to a neutral configuration vector
    model_ptr->get_state(q, v, a, tau);

    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    model_ptr->get_robot_mass(mass);
    model_ptr->get_B(B);
    model_ptr->get_B_inv(B_inv);
    model_ptr->get_C(C);
    model_ptr->get_g(g);
    model_ptr->get_b(b);
    model_ptr->get_p(p);

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Robot mass: \n" << mass << "\n " << std::endl;
    std::cout << "** q: \n" << q.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** v: \n" << v.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** a: \n" << a.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tau: \n" << tau.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B: \n" << B.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B_inv: \n" << B_inv.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** C: \n" << C.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** g: \n" << g.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** b: \n" << b.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** p: \n" << p.format(CleanFmt) << "\n " << std::endl;

}

TEST_F(TestModelInterface, load_floating_base_quadruped)
{
    Model::Ptr model_ptr(new Model(urdf_path_quadruped));

    ASSERT_TRUE(model_ptr->was_model_init_ok());

    auto jnt_names = model_ptr->get_jnt_names();

    std::cout << "** FLOATING BASE QUADRUPED DEBUG PRINTS**\n"  << std::endl;
    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Joint names: **"  << std::endl;
    for (std::string i: jnt_names)
        std::cout << "--> " << i << std::endl;
    std::cout << "** joint number: " <<model_ptr->get_jnt_number() << "\n " << std::endl;
    std::cout << "** nq: " << model_ptr->get_nq() << std::endl;
    std::cout << "** nv: " << model_ptr->get_nv() << std::endl;

    double mass = -1.0;
    Eigen::VectorXd q, v, a, tau,
                    g, p, b;
    Eigen::MatrixXd B, B_inv, C;

    model_ptr->set_random(); // sets q to a neutral configuration vector
    model_ptr->get_state(q, v, a, tau);

    model_ptr->update(); // computes all terms of the dynamics
    // and updates the forward kinematis

    model_ptr->get_robot_mass(mass);
    model_ptr->get_B(B);
    model_ptr->get_B_inv(B_inv);
    model_ptr->get_C(C);
    model_ptr->get_g(g);
    model_ptr->get_b(b);
    model_ptr->get_p(p);

    std::cout << "\nLoaded URDF at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** Robot mass: \n" << mass << "\n " << std::endl;
    std::cout << "** q: \n" << q.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** v: \n" << v.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** a: \n" << a.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tau: \n" << tau.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B: \n" << B.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** B_inv: \n" << B_inv.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** C: \n" << C.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** g: \n" << g.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** b: \n" << b.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** p: \n" << p.format(CleanFmt) << "\n " << std::endl;

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string urdf_name = "awesome_leg";
    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_fullpath_fixed = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
    const std::string urdf_fullpath_floating = (argc<=1) ? URDF_PATH + urdf_name + std::string("_floating.urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");

    urdf_path_fixed = urdf_fullpath_fixed;
    urdf_path_floating = urdf_fullpath_floating;

    urdf_path_anymal = URDF_PATH + std::string("anymal") + std::string(".urdf");
    urdf_path_quadruped = URDF_PATH + std::string("quadruped") + std::string(".urdf");

    return RUN_ALL_TESTS();
}
