#include <gtest/gtest.h>
#include <chrono>

#include "tests/generated/cmake_config.h"

#include "xbot2_utils.hpp"

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Xbot2Utils;
using namespace Eigen;

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

class TestIqOutGetter: public ::testing::Test {


    protected:

         TestIqOutGetter(){

         }

         virtual ~TestIqOutGetter() {
         }

         virtual void SetUp() {

         }

         virtual void TearDown() {
         }


};

TEST_F(TestIqOutGetter, test_read_from_ros)
{

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
