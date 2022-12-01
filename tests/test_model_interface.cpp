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
    Model model = Model(urdf_path);

    ASSERT_TRUE(model.was_model_init_ok());

    auto jnt_names = model.get_jnt_names();

    std::cout << "\nLoaded URDF at: "<< model.get_urdf_path() << "\n " << std::endl;
    std::cout << "** Joint names: **"  << std::endl;
    for (std::string i: jnt_names)
        std::cout << "--> " << i << std::endl;
    std::cout << "** joint number: " << model.get_jnt_number() << "\n " << std::endl;
    std::cout << "** nq: " << model.get_nq() << std::endl;
    std::cout << "** nv: " << model.get_nv() << std::endl;

}

TEST_F(TestModelInterface, compute_quantities)
{
    Model model = Model(urdf_path);

    std::string tip_framename = "tip1";

    Eigen::VectorXd q, v, a, tau,
                    g, p, b;
    Eigen::MatrixXd B, C, J;

    model.get_state(q, v, a, tau);

    model.update(q, v, tau, a); // computes all terms

    model.get_B(B);
    model.get_C(C);
    model.get_g(g);
    model.get_b(b);
    model.get_p(p);
    model.get_jac(tip_framename, Model::ReferenceFrame::LOCAL_WORLD_ALIGNED,
                  J);

    std::cout << "\nLoaded URDF at: "<< model.get_urdf_path() << "\n " << std::endl;
    std::cout << "** B: \n" << B.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** C: \n" << C.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** g: \n" << g.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** b: \n" << b.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** p: \n" << p.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** J (q_dot -> " << tip_framename << " - LOCAL_WORLD_ALIGNED) :\n " << J.format(CleanFmt) << "\n " << std::endl;

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string urdf_name = "awesome_leg";
    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_fullpath = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
    urdf_path = urdf_fullpath;

    return RUN_ALL_TESTS();
}
