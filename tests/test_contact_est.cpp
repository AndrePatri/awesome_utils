#include <gtest/gtest.h>
#include <chrono>

#include "tests/generated/cmake_config.h"

#include "model_interface.hpp"
#include "contact_est_utils.hpp"
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace ModelInterface;
using namespace ContactEstUtils;

std::string urdf_path_fixed;
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

class TestContactEst: public ::testing::Test {


    protected:

         TestContactEst(){

         }

         virtual ~TestContactEst() {
         }

         virtual void SetUp() {

         }

         virtual void TearDown() {
         }


};

TEST_F(TestContactEst, compute_quantities_fixed_base)
{
    Eigen::VectorXd q, v, a, tau, p,
                    tau_c;

    utils_defs::Force3D f_c;
    utils_defs::Torque3D t_c;
    utils_defs::SpatialJac J;

    Eigen::MatrixXd B, C;

    Model::Ptr model_ptr(new Model(urdf_path_fixed));

    model_ptr->get_state(q, v, a, tau);
    Eigen::VectorXd base_position = Eigen::VectorXd::Zero(3);
    Eigen::VectorXd base_orient = Eigen::VectorXd::Zero(4);
    q << 0.0, 0.0, 0.0; // leg vertical with horizontal base link
    v << 0.0, 0.0, 0.0; // static
    a << 0.0, 0.0, 0.0; // static
    tau << 0.0, 0.0, 0.0;
    model_ptr->set_q(q);
    model_ptr->set_v(v);
    model_ptr->set_a(a);
    model_ptr->set_tau(tau);

    model_ptr->update(); // computes all terms

    double dt = 0.001; // control loop dt
    double BW = 100.0; // observer bandwitdth [Hz]

    MomentumBasedFObs::Reg6D lambda = Eigen::VectorXd::Zero(6);
    lambda << 0.001, 0.001, 0.001, 0.01, 0.01, 0.01;

    std::vector<int> selector{0, 1, 2}; // only force

    std::vector<std::string> contacts{"tip1"};

    bool regularize_w = true;

    MomentumBasedFObs::Ptr f_obs_ptr(new MomentumBasedFObs(model_ptr, dt,
                                                           contacts,
                                                           BW,
                                                           lambda, regularize_w,
                                                           selector));

    f_obs_ptr->update(); // compute estimates using the current state in model_ptr

    f_obs_ptr->get_tau_obs(tau_c);
    f_obs_ptr->get_f_est_at(contacts[0], f_c);
    f_obs_ptr->get_t_est_at(contacts[0], t_c);

    std::cout << "** FIXED BASE DEBUG PRINTS**\n"  << std::endl;
    std::cout << "\nURDF loaded at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** tau_c: \n" << tau_c << "\n " << std::endl;
    std::cout << "** f_c: \n" << f_c << "\n " << std::endl;
    std::cout << "** w_c: \n" << t_c << "\n " << std::endl;

}

TEST_F(TestContactEst, compute_quantities_floating_base)
{
    Eigen::VectorXd q, v, a, tau, p,
                    tau_c;

    utils_defs::Force3D f_c;
    utils_defs::Torque3D t_c;
    utils_defs::SpatialJac J;

    Eigen::MatrixXd B, C;

    Model::Ptr model_ptr(new Model(urdf_path_floating));

    model_ptr->set_neutral();
    model_ptr->get_state(q, v, a, tau);

    model_ptr->update(); // computes all terms

    double dt = 0.001; // control loop dt
    double BW = 100.0; // observer bandwitdth [Hz]

    MomentumBasedFObs::Reg6D lambda = Eigen::VectorXd::Zero(6);
    lambda << 0.001, 0.001, 0.001, 0.01, 0.01, 0.01;

    std::vector<int> selector{0, 1, 2}; // only force

    std::vector<std::string> contacts{"tip1"};

    bool regularize_w = true;

    MomentumBasedFObs::Ptr f_obs_ptr(new MomentumBasedFObs(model_ptr, dt,
                                                           contacts,
                                                           BW,
                                                           lambda, regularize_w,
                                                           selector));

    f_obs_ptr->update(); // compute estimates using the current state in model_ptr

    f_obs_ptr->get_tau_obs(tau_c);
    f_obs_ptr->get_f_est_at(contacts[0], f_c);
    f_obs_ptr->get_t_est_at(contacts[0], t_c);

    std::cout << "** FLOATING BASE DEBUG PRINTS**\n"  << std::endl;
    std::cout << "\nURDF loaded at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** tau_c: \n" << tau_c << "\n " << std::endl;
    std::cout << "** f_c: \n" << f_c << "\n " << std::endl;
    std::cout << "** w_c: \n" << t_c << "\n " << std::endl;

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  std::string urdf_name = "awesome_leg";
  // You should change here to set up your own URDF file or just pass it as an argument of this example.
  const std::string urdf_fullpath_fixed = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
  const std::string urdf_fullpath_floating = (argc<=1) ? URDF_PATH + urdf_name + std::string("_floating.urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");

  urdf_path_fixed = urdf_fullpath_fixed;
  urdf_path_floating = urdf_fullpath_floating;

}
