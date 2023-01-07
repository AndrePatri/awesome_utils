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
std::string urdf_anymal;
std::string urdf_quadruped;

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

//TEST_F(TestContactEst, compute_quantities_fixed_base)
//{
//    Eigen::VectorXd q, v, a, tau, p,
//                    tau_c;

//    utils_defs::Force3D f_c;
//    utils_defs::Torque3D t_c;
//    utils_defs::SpatialJac J;

//    Eigen::MatrixXd B, C;

//    Model::Ptr model_ptr(new Model(urdf_path_fixed));

//    model_ptr->get_state(q, v, a, tau);
//    Eigen::VectorXd base_position = Eigen::VectorXd::Zero(3);
//    Eigen::VectorXd base_orient = Eigen::VectorXd::Zero(4);
//    q << 0.0, 0.0, 0.0; // leg vertical with horizontal base link
//    v << 0.0, 0.0, 0.0; // static
//    model_ptr->set_q(q);
//    model_ptr->set_v(v);

//    model_ptr->update(); // computes all terms

//    double dt = 0.0001; // control loop dt
//    double BW = 10000.0; // observer bandwitdth [Hz]

//    MomentumBasedFObs::Reg6D lambda = Eigen::VectorXd::Zero(6);
//    lambda << 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001;

//    std::vector<int> selector{0, 1, 2}; // only force

//    std::vector<std::string> contacts{"tip1"};

//    bool regularize_w = true;

//    MomentumBasedFObs::Ptr f_obs_ptr(new MomentumBasedFObs(model_ptr, dt,
//                                                           contacts,
//                                                           BW,
//                                                           lambda, regularize_w,
//                                                           selector));

//    f_obs_ptr->update(); // compute estimates using the current state in model_ptr

//    f_obs_ptr->get_tau_obs(tau_c);
//    f_obs_ptr->get_f_est_at(contacts[0], f_c);
//    f_obs_ptr->get_t_est_at(contacts[0], t_c);

//    std::cout << "** FIXED BASE DEBUG PRINTS**\n"  << std::endl;
//    std::cout << "\nURDF loaded at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
//    std::cout << "** tau_c: \n" << tau_c << "\n " << std::endl;
//    std::cout << "** f_c: \n" << f_c << "\n " << std::endl;
//    std::cout << "** w_c: \n" << t_c << "\n " << std::endl;

//}

//TEST_F(TestContactEst, compute_quantities_floating_base)
//{
//    Eigen::VectorXd q, v, a, tau, p,
//                    tau_c, g, CT_v,
//                    tau_c_raw_static;

//    utils_defs::Force3D f_c;
//    utils_defs::Torque3D t_c;
//    utils_defs::SpatialJac J;

//    Eigen::MatrixXd B, C;

//    Model::Ptr model_ptr(new Model(urdf_path_floating));

//    model_ptr->set_neutral();
//    model_ptr->get_state(q, v, a, tau);

//    model_ptr->update(); // computes all terms

//    model_ptr->get_g(g);
//    model_ptr->get_tau(tau);
//    tau_c_raw_static = g - tau; // raw disturbance torques (not filtered

//    double dt = 0.0001; // control loop dt
//    double BW = 10000.0; // observer bandwitdth [Hz]

//    MomentumBasedFObs::Reg6D lambda = Eigen::VectorXd::Zero(6);
//    lambda << 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001;

//    std::vector<int> selector{0, 1, 2}; // only force

//    std::vector<std::string> contacts{"tip1"};

//    bool regularize_w = true;

//    MomentumBasedFObs::Ptr f_obs_ptr(new MomentumBasedFObs(model_ptr, dt,
//                                                           contacts,
//                                                           BW,
//                                                           lambda, regularize_w,
//                                                           selector));

//    f_obs_ptr->update(); // compute estimates using the current state in model_ptr

//    f_obs_ptr->get_tau_obs(tau_c);
//    f_obs_ptr->get_f_est_at(contacts[0], f_c);
//    f_obs_ptr->get_t_est_at(contacts[0], t_c);

//    std::cout << "** FLOATING BASE DEBUG PRINTS**\n"  << std::endl;
//    std::cout << "\nURDF loaded at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
//    std::cout << "** tau_c: \n" << tau_c << "\n " << std::endl;
//    std::cout << "** tau_c_raw: \n" << tau_c_raw_static << "\n " << std::endl;
//    std::cout << "** f_c: \n" << f_c << "\n " << std::endl;
//    std::cout << "** w_c: \n" << t_c << "\n " << std::endl;

//}

//TEST_F(TestContactEst, test_contact_est_anymal)
//{
//    Eigen::VectorXd q, v, a, tau, p,
//                    tau_c, g, CT_v,
//                    tau_c_raw_static;

//    utils_defs::Force3D fc_lf, fc_lh, fc_rf, fc_rh;
//    utils_defs::Torque3D tc_lf, tc_lh, tc_rf, tc_rh;

//    Model::Ptr model_ptr(new Model(urdf_anymal));

//    model_ptr->set_neutral();
//    model_ptr->get_state(q, v, a, tau);

//    model_ptr->update(); // computes all terms

//    model_ptr->get_g(g);
//    model_ptr->get_tau(tau);
//    tau_c_raw_static = g - tau; // raw disturbance torques (not filtered

//    int nv = model_ptr->get_nv();

//    double dt = 0.0001; // control loop dt
//    double BW = 10000.0; // observer bandwitdth [Hz]

//    int rollout_number = 15; // number of samples for testing steady state convergence
//    // of the observer over multiple sample intervals

//    Eigen::MatrixXd Fc_lf, Fc_lh, Fc_rf, Fc_rh;
//    Eigen::MatrixXd Tc_lf, Tc_lh, Tc_rf, Tc_rh;

//    Eigen::MatrixXd Tau_c, Tau_c_raw;

//    Fc_lf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Fc_lh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Fc_rf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Fc_rh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Tc_lf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Tc_lh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Tc_rf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Tc_rh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
//    Tau_c = Eigen::MatrixXd::Zero(nv, rollout_number);
//    Tau_c_raw = Eigen::MatrixXd::Zero(nv, rollout_number);

//    MomentumBasedFObs::Reg6D lambda = Eigen::VectorXd::Zero(6);
//    lambda << 0.000000000000001,
//              0.000000000000001,
//              0.000000000000001,
//              0.000000000000001,
//              0.000000000000001,
//              0.000000000000001;

//    std::vector<int> selector{0, 1, 2}; // only force

//    std::vector<std::string> contacts{"LF_FOOT", "LH_FOOT", "RF_FOOT", "RH_FOOT"};

//    bool regularize_w = true;

//    MomentumBasedFObs::Ptr f_obs_ptr(new MomentumBasedFObs(model_ptr, dt,
//                                                           contacts,
//                                                           BW,
//                                                           lambda, regularize_w,
//                                                           selector));

//    for (int i = 0; i < rollout_number; i++)
//    {
//        // we leave the state of the model fixed, so we are looking
//        // at the homogenous part of the observer dynamics.

//        f_obs_ptr->update();

//        f_obs_ptr->get_tau_obs(tau_c);
//        f_obs_ptr->get_f_est_at(contacts[0], fc_lf);
//        f_obs_ptr->get_f_est_at(contacts[1], fc_lh);
//        f_obs_ptr->get_f_est_at(contacts[2], fc_rf);
//        f_obs_ptr->get_f_est_at(contacts[3], fc_rh);
//        f_obs_ptr->get_t_est_at(contacts[0], tc_lf);
//        f_obs_ptr->get_t_est_at(contacts[1], tc_lh);
//        f_obs_ptr->get_t_est_at(contacts[2], tc_rf);
//        f_obs_ptr->get_t_est_at(contacts[3], tc_rh);

//        Fc_lf.block(0, i, Fc_lf.rows(), 1) = fc_lf;
//        Fc_lh.block(0, i, Fc_lf.rows(), 1) = fc_lh;
//        Fc_rf.block(0, i, Fc_lf.rows(), 1) = fc_rf;
//        Fc_rh.block(0, i, Fc_lf.rows(), 1) = fc_rh;
//        Tc_lf.block(0, i, Fc_lf.rows(), 1) = tc_lf;
//        Tc_lh.block(0, i, Fc_lf.rows(), 1) = tc_lh;
//        Tc_rf.block(0, i, Fc_lf.rows(), 1) = tc_rf;
//        Tc_rh.block(0, i, Fc_lf.rows(), 1) = tc_rh;

//        Tau_c.block(0, i, Tau_c.rows(), 1) = tau_c;
//        Tau_c_raw.block(0, i, Tau_c.rows(), 1) = tau_c_raw_static;

//    }

//    std::cout << "** ANYMAL DEBUG PRINTS**\n"  << std::endl;
//    std::cout << "\nURDF loaded at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
//    std::cout << "** tau_c: \n" << Tau_c << "\n " << std::endl;
//    std::cout << "** tau_c_raw: \n" << Tau_c_raw << "\n " << std::endl;
//    std::cout << "** Fc_lf: \n" << Fc_lf << "\n " << std::endl;
//    std::cout << "** Fc_lh: \n" << Fc_lh << "\n " << std::endl;
//    std::cout << "** Fc_rf: \n" << Fc_rf << "\n " << std::endl;
//    std::cout << "** Fc_rh: \n" << Fc_rh << "\n " << std::endl;
//    std::cout << "** Tc_lf: \n" << Tc_lf << "\n " << std::endl;
//    std::cout << "** Tc_lh: \n" << Tc_lh << "\n " << std::endl;
//    std::cout << "** Tc_rf: \n" << Tc_rf << "\n " << std::endl;
//    std::cout << "** Tc_rh: \n" << Tc_rh << "\n " << std::endl;

//}

TEST_F(TestContactEst, test_contact_est_quadruped)
{
    Eigen::VectorXd q, v, a, tau, p,
                    tau_c, g, CT_v,
                    tau_c_raw_static;

    utils_defs::Force3D fc_lf, fc_lh, fc_rf, fc_rh;
    utils_defs::Torque3D tc_lf, tc_lh, tc_rf, tc_rh;

    Model::Ptr model_ptr(new Model(urdf_quadruped));

    model_ptr->set_neutral();
    model_ptr->get_state(q, v, a, tau);

    model_ptr->update(); // computes all terms

    model_ptr->get_g(g);
    model_ptr->get_tau(tau);
    tau_c_raw_static = g - tau; // raw disturbance torques (not filtered

    int nv = model_ptr->get_nv();

    double dt = 0.00001; // control loop dt
    double BW = 10000.0; // observer bandwitdth [Hz]

    int rollout_number = 15; // number of samples for testing steady state convergence
    // of the observer over multiple sample intervals

    Eigen::MatrixXd Fc_lf, Fc_lh, Fc_rf, Fc_rh;
    Eigen::MatrixXd Tc_lf, Tc_lh, Tc_rf, Tc_rh;

    Eigen::MatrixXd Tau_c, Tau_c_raw;

    Fc_lf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Fc_lh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Fc_rf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Fc_rh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Tc_lf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Tc_lh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Tc_rf = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Tc_rh = Eigen::MatrixXd::Zero(fc_lf.size(), rollout_number);
    Tau_c = Eigen::MatrixXd::Zero(nv, rollout_number);
    Tau_c_raw = Eigen::MatrixXd::Zero(nv, rollout_number);

    Eigen::MatrixXd J_c_tot;

    double lambda = 1e-6;

    std::vector<int> selector{0, 1, 2}; // only force

    std::vector<std::string> contacts{"ball_1", "ball_2", "ball_3", "ball_4"};

    bool regularize_delta_w = true;

    MomentumBasedFObs::Ptr f_obs_ptr(new MomentumBasedFObs(model_ptr, dt,
                                                           contacts,
                                                           BW,
                                                           lambda, regularize_delta_w,
                                                           selector));

    Eigen::MatrixXd A_lambda = Eigen::MatrixXd::Zero(2 * fc_lf.size() * contacts.size(), 2 * fc_lf.size() * contacts.size());
    Eigen::MatrixXd B_lambda = Eigen::MatrixXd::Zero(2 * fc_lf.size() * contacts.size(), rollout_number);
    Eigen::VectorXd b_lambda = Eigen::VectorXd::Zero(2 * fc_lf.size() * contacts.size());


    for (int i = 0; i < rollout_number; i++)
    {
        // we leave the state of the model fixed, so we are looking
        // at the homogenous part of the observer dynamics.

        f_obs_ptr->update();

        f_obs_ptr->get_tau_obs(tau_c);
        f_obs_ptr->get_f_est_at(contacts[0], fc_lf);
        f_obs_ptr->get_f_est_at(contacts[1], fc_lh);
        f_obs_ptr->get_f_est_at(contacts[2], fc_rf);
        f_obs_ptr->get_f_est_at(contacts[3], fc_rh);
        f_obs_ptr->get_t_est_at(contacts[0], tc_lf);
        f_obs_ptr->get_t_est_at(contacts[1], tc_lh);
        f_obs_ptr->get_t_est_at(contacts[2], tc_rf);
        f_obs_ptr->get_t_est_at(contacts[3], tc_rh);

        Fc_lf.block(0, i, Fc_lf.rows(), 1) = fc_lf;
        Fc_lh.block(0, i, Fc_lf.rows(), 1) = fc_lh;
        Fc_rf.block(0, i, Fc_lf.rows(), 1) = fc_rf;
        Fc_rh.block(0, i, Fc_lf.rows(), 1) = fc_rh;
        Tc_lf.block(0, i, Fc_lf.rows(), 1) = tc_lf;
        Tc_lh.block(0, i, Fc_lf.rows(), 1) = tc_lh;
        Tc_rf.block(0, i, Fc_lf.rows(), 1) = tc_rf;
        Tc_rh.block(0, i, Fc_lf.rows(), 1) = tc_rh;

        Tau_c.block(0, i, Tau_c.rows(), 1) = tau_c;
        Tau_c_raw.block(0, i, Tau_c.rows(), 1) = tau_c_raw_static;

        f_obs_ptr->get_reg_matrices(A_lambda, b_lambda);

        B_lambda.block(0, i, A_lambda.rows(), 1) = b_lambda;

    }

    f_obs_ptr->get_J_c_tot(J_c_tot);

    std::cout << "** QUADRUPED DEBUG PRINTS**\n"  << std::endl;
    std::cout << "\nURDF loaded at: "<< model_ptr->get_urdf_path() << "\n " << std::endl;
    std::cout << "** tau_c: \n" << Tau_c.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** tau_c_raw: \n" << Tau_c_raw.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Fc_lf: \n" << Fc_lf.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Fc_lh: \n" << Fc_lh.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Fc_rf: \n" << Fc_rf.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Fc_rh: \n" << Fc_rh.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Tc_lf: \n" << Tc_lf.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Tc_lh: \n" << Tc_lh.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Tc_rf: \n" << Tc_rf.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Tc_rh: \n" << Tc_rh.format(CleanFmt) << "\n " << std::endl;
    std::cout << "** Jc_tot transpose: \n" << J_c_tot.transpose().format(CleanFmt)<< "\n " << std::endl;
    std::cout << "** A_lambda: \n" << A_lambda.format(CleanFmt)<< "\n " << std::endl;
    std::cout << "** B_lambda: \n" << B_lambda.format(CleanFmt)<< "\n " << std::endl;


}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  std::string urdf_name = "awesome_leg";
  // You should change here to set up your own URDF file or just pass it as an argument of this example.
  const std::string urdf_fullpath_fixed = (argc<=1) ? URDF_PATH + urdf_name + std::string(".urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");
  const std::string urdf_fullpath_floating = (argc<=1) ? URDF_PATH + urdf_name + std::string("_floating.urdf"): URDF_PATH + std::string(argv[1]) + std::string(".urdf");

  urdf_path_fixed = urdf_fullpath_fixed;
  urdf_path_floating = urdf_fullpath_floating;
  urdf_anymal = URDF_PATH + std::string("anymal.urdf");
  urdf_quadruped = URDF_PATH + std::string("quadruped.urdf");

  return RUN_ALL_TESTS();
}
