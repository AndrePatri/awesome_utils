#include <gtest/gtest.h>
#include <chrono>

#include "tests/generated/cmake_config.h"

#include "awesome_utils/model_interface.hpp"

#include <iostream>

using namespace ModelInterface;

std::string urdf_path;

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

    std::cout << "Loaded URDF at: "<< model.get_urdf_path() << "\n " << std::endl;
    std::cout << "** Joint names: **" << "\n " << std::endl;
    for (std::string i: jnt_names)
        std::cout << "--> " << i << "\n";
    std::cout << "** Joint number: **" << model.get_jnt_number() << "\n " << std::endl;
    std::cout << "** nq: **" << model.get_nq() << std::endl;
    std::cout << "** nv: **" << model.get_nv() << std::endl;

}

TEST_F(TestModelInterface, compute_quantities)
{
    Model model = Model(urdf_path);

    ASSERT_TRUE(model.was_model_init_ok());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    std::string urdf_name = "awesome_leg";
    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_fullpath = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
    urdf_path = urdf_fullpath;

    return RUN_ALL_TESTS();
}
