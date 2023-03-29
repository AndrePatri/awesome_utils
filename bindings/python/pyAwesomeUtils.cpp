#if defined(WITH_MODEL_INTERFACE)
#include "include/awesome_utils/model_interface.hpp"
#include "include/awesome_utils/contact_est_utils.hpp"
#endif

#include "include/awesome_utils/sign_proc_utils.hpp"
#include "include/awesome_utils/traj_utils.hpp"
#include "include/awesome_utils/calib_utils.hpp"

#if defined(WITH_XBOT2)
#include "include/awesome_utils/power_utils.hpp"
#endif

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

//************* Signal processing utilities bindings *************//

using namespace SignProcUtils;

// NumDiff
auto construct_num_diff = [](int n_jnts, double dt, int order)
{

    NumDiff::Ptr num_diff;

    num_diff.reset(new NumDiff(n_jnts, dt, order));

    return num_diff;

};
auto dot = [](NumDiff& self)
{
    Eigen::VectorXd sample_dot;

    self.dot(sample_dot);

    return sample_dot;
};

// NumInt
auto construct_num_int = [](int n_jnts, double dt, double T_horizon)
{

    NumInt::Ptr num_int;

    num_int.reset(new NumInt(n_jnts, dt, T_horizon));

    return num_int;

};
namespace num_int {
auto get(NumInt& self)
{
    Eigen::VectorXd int_sample;

    self.get(int_sample);

    return int_sample;
}
}
// NumIntRt
auto construct_num_int_rt = [](int n_jnts, double dt)
{

    NumIntRt::Ptr num_int_rt;

    num_int_rt.reset(new NumIntRt(n_jnts, dt));

    return num_int_rt;

};
namespace num_int_rt {
auto get(NumIntRt& self)
{
    Eigen::VectorXd int_sample;

    self.get(int_sample);

    return int_sample;
}
}

// MovAvrgFilt
auto construct_mov_avrg_filt = [](int n_jnts, double dt, int window_size)
{

    MovAvrgFilt::Ptr mov_avrg_filt;

    mov_avrg_filt.reset(new MovAvrgFilt(n_jnts, dt, window_size));

    return mov_avrg_filt;

};
auto construct_mov_avrg_filt2 = [](int n_jnts, double dt, double cutoff_freq)
{

    MovAvrgFilt::Ptr mov_avrg_filt;

    mov_avrg_filt.reset(new MovAvrgFilt(n_jnts, dt, cutoff_freq));

    return mov_avrg_filt;

};
namespace mov_avrg_filt {
auto get(MovAvrgFilt& self)
{
    Eigen::VectorXd filt_sample;

    self.get(filt_sample);

    return filt_sample;
}
}

// SignWithMem
auto construct_sign_with_mem = [](double signal_3sigma,
                                  double tanh_coeff)
{

    SignWithMem::Ptr sign_with_mem;

    sign_with_mem.reset(new SignWithMem(signal_3sigma, tanh_coeff));

    return sign_with_mem;

};

// SmoooothSign
auto construct_smooth_sign = [](double signal_3sigma,
                                    int alpha,
                                    double beta,
                                    bool use_threshold)
{

    SmoooothSign::Ptr smooth_sign;

    smooth_sign.reset(new SmoooothSign(signal_3sigma,
                                       alpha,
                                       beta,
                                       use_threshold));

    return smooth_sign;

};
auto construct_smooth_sign2 = [](double signal_3sigma,
                                    int alpha,
                                    bool use_threshold)
{

    SmoooothSign::Ptr smooth_sign;

    smooth_sign.reset(new SmoooothSign(signal_3sigma,
                                       alpha,
                                       use_threshold));

    return smooth_sign;

};

//************* Calibration utilities bindings *************//

using namespace CalibUtils;

// IqEstimator

auto construct_iq_estimator = [](Eigen::VectorXd K_t,
                                Eigen::VectorXd K_d0, Eigen::VectorXd K_d1,
                                Eigen::VectorXd rot_MoI,
                                Eigen::VectorXd red_ratio,
                                int alpha = 10,
                                double q_dot_3sigma = 0.001,
                                bool dump_data2mat = false,
                                std::string dump_path = "/tmp")
{
    IqEstimator::Ptr iq_est;

    iq_est.reset(new IqEstimator(K_t,
                                 K_d0, K_d1,
                                 rot_MoI,
                                 red_ratio,
                                 alpha,
                                 q_dot_3sigma,
                                 dump_data2mat,
                                 dump_path));

    return iq_est;
};

namespace iq_estimator{

    auto get_iq_estimate(IqEstimator& self)
    {

        Eigen::VectorXd iq_est;

        self.get_iq_estimate(iq_est);

        return iq_est;
    }

    auto get_iq_meas(IqEstimator& self)
    {

        Eigen::VectorXd iq_meas;

        self.get_iq(iq_meas);

        return iq_meas;
    }

};

// RotDynCal
auto construct_rot_dyn_cal = [](int window_size,
        Eigen::VectorXd red_ratio,
        Eigen::VectorXd ig_Kt,
        Eigen::VectorXd ig_rot_MoI,
        Eigen::VectorXd ig_Kd0,
        Eigen::VectorXd ig_Kd1,
        double lambda,
        int alpha,
        double q_dot_3sigma,
        bool verbose)
{

    RotDynCal::Ptr rot_cal_ptr;

    rot_cal_ptr.reset(new RotDynCal(window_size,
                                    red_ratio,
                                    ig_Kt,
                                    ig_rot_MoI,
                                    ig_Kd0,
                                    ig_Kd1,
                                    lambda,
                                    alpha,
                                    q_dot_3sigma,
                                    verbose));

    return rot_cal_ptr;
};

auto construct_rot_dyn_cal2 = [](int window_size,
        Eigen::VectorXd red_ratio,
        Eigen::VectorXd ig_Kt,
        Eigen::VectorXd ig_rot_MoI,
        Eigen::VectorXd ig_Kd0,
        Eigen::VectorXd ig_Kd1,
        int alpha,
        double q_dot_3sigma,
        bool verbose)
{

    RotDynCal::Ptr rot_cal_ptr;

    rot_cal_ptr.reset(new RotDynCal(window_size,
                                    red_ratio,
                                    ig_Kt,
                                    ig_rot_MoI,
                                    ig_Kd0,
                                    ig_Kd1,
                                    alpha,
                                    q_dot_3sigma,
                                    verbose));

    return rot_cal_ptr;
};

namespace rot_dyn_cal{

    auto get_opt_Kd0(RotDynCal& self)
    {
        Eigen::VectorXd Kd0;

        self.get_opt_Kd0(Kd0);

        return Kd0;
    }

    auto get_opt_Kd1(RotDynCal& self)
    {
        Eigen::VectorXd Kd1;

        self.get_opt_Kd1(Kd1);

        return Kd1;
    }

    auto get_opt_Kt(RotDynCal& self)
    {
        Eigen::VectorXd Kt;

        self.get_opt_Kt(Kt);

        return Kt;
    }

    auto get_opt_rot_MoI(RotDynCal& self)
    {
        Eigen::VectorXd rot_MoI;

        self.get_opt_rot_MoI(rot_MoI);

        return rot_MoI;
    }

    auto get_tau_friction(RotDynCal& self)
    {
        Eigen::VectorXd tau_friction;

        self.get_tau_friction(tau_friction);

        return tau_friction;
    }

    auto get_alpha_d(RotDynCal& self)
    {
        Eigen::VectorXd alpha_d0;
        Eigen::VectorXd alpha_d1;

        self.get_alpha_d(alpha_d0, alpha_d1);

        auto output = std::make_tuple(alpha_d0, alpha_d1);

        return output;
    }

    auto get_alpha_inertial(RotDynCal& self)
    {
        Eigen::VectorXd alpha_inertial;

        self.get_alpha_inertial(alpha_inertial);

        return alpha_inertial;
    }

    auto get_alpha_kt(RotDynCal& self)
    {
        Eigen::VectorXd alpha_kt;

        self.get_alpha_kt(alpha_kt);

        return alpha_kt;
    }

    auto get_tau_motor(RotDynCal& self)
    {
        Eigen::VectorXd tau_motor;

        self.get_tau_motor(tau_motor);

        return tau_motor;
    }

    auto get_tau_inertial(RotDynCal& self)
    {
        Eigen::VectorXd tau_inertial;

        self.get_tau_inertial(tau_inertial);

        return tau_inertial;
    }

    auto get_sol_millis(RotDynCal& self)
    {
        Eigen::VectorXd sol_millis;

        self.get_sol_millis(sol_millis);

        return sol_millis;
    }

    auto get_cal_mask(RotDynCal& self)
    {
        std::vector<bool> cal_mask;

        self.get_cal_mask(cal_mask);

        return cal_mask;
    }

    auto get_lambda(RotDynCal& self)
    {
        Eigen::VectorXd lambda;

        self.get_lambda(lambda);

        return lambda;
    }

    auto get_lambda_des(RotDynCal& self)
    {
        Eigen::VectorXd lambda;

        self.get_lambda_des(lambda);

        return lambda;
    }

    auto get_lambda_high(RotDynCal& self)
    {
        Eigen::VectorXd lambda;

        self.get_lambda_high(lambda);

        return lambda;
    }

    auto get_ig_Kd0(RotDynCal& self)
    {
        Eigen::VectorXd Kd0;

        self.get_ig_Kd0(Kd0);

        return Kd0;
    }

    auto get_ig_Kd1(RotDynCal& self)
    {
        Eigen::VectorXd kd1;

        self.get_ig_Kd1(kd1);

        return kd1;
    }

    auto get_ig_Kt(RotDynCal& self)
    {
        Eigen::VectorXd Kt;

        self.get_ig_Kt(Kt);

        return Kt;
    }

    auto get_ig_MoI(RotDynCal& self)
    {
        Eigen::VectorXd MoI;

        self.get_ig_MoI(MoI);

        return MoI;
    }

};

#if defined(WITH_XBOT2)

//************* Power utils bindings *************//

using namespace PowerUtils;

#endif

#if defined(WITH_MODEL_INTERFACE)

//************* Model interface bindings *************//

using namespace ModelInterface;

auto construct_model_interface = [](std::string urdf_path, bool add_floating_jnt)
{

    Model::Ptr model_ptr;

    model_ptr.reset(new Model(urdf_path, add_floating_jnt));

    return model_ptr;
};

auto get_robot_mass = [](Model& self)
{

    double mass = -1.0;

    self.get_robot_mass(mass);

    return mass;
};

auto get_q = [](Model& self)
{

    Eigen::VectorXd q;

    self.get_q(q);

    return q;
};

auto get_v = [](Model& self)
{

    Eigen::VectorXd v;

    self.get_v(v);

    return v;

};

auto get_a = [](Model& self)
{

    Eigen::VectorXd a;

    self.get_a(a);

    return a;

};

auto get_tau = [](Model& self)
{

    Eigen::VectorXd tau;

    self.get_tau(tau);

    return tau;

};

auto get_p = [](Model& self)
{

    Eigen::VectorXd p;

    self.get_p(p);

    return p;

};

auto get_b = [](Model& self)
{

    Eigen::VectorXd b;

    self.get_b(b);

    return b;

};

auto get_g = [](Model& self)
{

    Eigen::VectorXd g;

    self.get_g(g);

    return g;

};

auto get_B = [](Model& self)
{

    Eigen::MatrixXd B;

    self.get_B(B);

    return B;

};

auto get_C = [](Model& self)
{

    Eigen::MatrixXd C;

    self.get_C(C);

    return C;

};

auto get_J = [](Model& self,
                std::string frame_name,
                Model::ReferenceFrame ref_frame = Model::ReferenceFrame::LOCAL_WORLD_ALIGNED)
{

    utils_defs::SpatialJac J;

    self.get_jac(frame_name, J, ref_frame);

    return J;

};

auto get_J_dot = [](Model& self,
                    std::string frame_name,
                    Model::ReferenceFrame ref_frame = Model::ReferenceFrame::LOCAL_WORLD_ALIGNED)
{

    utils_defs::SpatialJac J_dot;

    self.get_jac_dot(frame_name, J_dot, ref_frame);

    return J_dot;

};

#endif

PYBIND11_MODULE(awesome_pyutils, m) {

    // signal processing utilities

    py::class_<NumDiff, std::shared_ptr<NumDiff>>(m, "NumDiff")
            .def(py::init(construct_num_diff),
                 py::arg("n_jnts"),
                 py::arg("dt"),
                 py::arg("order") = 1
                 )
            .def("add_sample", &NumDiff::add_sample,
                 py::arg("sample"))

            .def("dot", dot)
            ;

    py::class_<NumInt, std::shared_ptr<NumInt>>(m, "NumInt")
            .def(py::init(construct_num_int),
                 py::arg("n_jnts"),
                 py::arg("dt"),
                 py::arg("T_horizon")
                 )

            .def("add_sample", &NumInt::add_sample,
                 py::arg("sample"))

            .def("get", num_int::get)

            ;

    py::class_<NumIntRt, std::shared_ptr<NumIntRt>>(m, "NumIntRt")
            .def(py::init(construct_num_int_rt),
                 py::arg("n_jnts"),
                 py::arg("dt")
                 )

            .def("add_sample", &NumIntRt::add_sample,
                 py::arg("sample"))

            .def("get", num_int_rt::get)

            ;

    py::class_<MovAvrgFilt, std::shared_ptr<MovAvrgFilt>>(m, "MovAvrgFilt")
            .def(py::init(construct_mov_avrg_filt),
                 py::arg("n_jnts"),
                 py::arg("dt"),
                 py::arg("window_size")
                 )
            .def(py::init(construct_mov_avrg_filt2),
                 py::arg("n_jnts"),
                 py::arg("dt"),
                 py::arg("cutoff_freq")
                 )

            .def("add_sample", &MovAvrgFilt::add_sample,
                 py::arg("sample"))

            .def("get", mov_avrg_filt::get)

            ;

    py::class_<SignWithMem, std::shared_ptr<SignWithMem>>(m, "SignWithMem")
            .def(py::init(construct_sign_with_mem),
                 py::arg("signal_3sigma"),
                 py::arg("tanh_coeff")
                 )

            .def("sign", &SignWithMem::sign,
                 py::arg("value"))

            ;

    py::class_<SmoooothSign, std::shared_ptr<SmoooothSign>>(m, "SmoooothSign")
            .def(py::init(construct_smooth_sign),
                 py::arg("signal_3sigma"),
                 py::arg("alpha"),
                 py::arg("beta"),
                 py::arg("use_threshold")
                 )
            .def(py::init(construct_smooth_sign2),
                 py::arg("signal_3sigma"),
                 py::arg("alpha"),
                 py::arg("use_threshold")
                 )

            .def("sign", &SmoooothSign::sign,
                 py::arg("value"))

            ;

    // calib. utilities

    py::class_<RotDynCal, std::shared_ptr<RotDynCal>>(m, "RotDynCal")
            .def(py::init(construct_rot_dyn_cal),
                 py::arg("window_size"),
                 py::arg("red_ratio"),
                 py::arg("ig_Kt"),
                 py::arg("ig_rot_MoI"),
                 py::arg("ig_Kd0"),
                 py::arg("ig_Kd1"),
                 py::arg("lambda") = 2.0,
                 py::arg("alpha") = 10,
                 py::arg("q_dot_3sigma") = 0.001,
                 py::arg("verbose") = false)
            .def(py::init(construct_rot_dyn_cal2),
                 py::arg("window_size"),
                 py::arg("red_ratio"),
                 py::arg("ig_Kt"),
                 py::arg("ig_rot_MoI"),
                 py::arg("ig_Kd0"),
                 py::arg("ig_Kd1"),
                 py::arg("alpha") = 10,
                 py::arg("q_dot_3sigma") = 0.001,
                 py::arg("verbose") = false)

            .def("add_sample", &RotDynCal::add_sample,
                 py::arg("q_dot"),
                 py::arg("q_ddot"),
                 py::arg("iq"),
                 py::arg("tau"))

            .def("set_ig_Kd0", &RotDynCal::set_ig_Kd0,
                 py::arg("ig_Kd0"))
            .def("set_ig_Kd1", &RotDynCal::set_ig_Kd1,
                 py::arg("ig_Kd1"))
            .def("set_ig_Kt", &RotDynCal::set_ig_Kt,
                 py::arg("ig_Kt"))
            .def("set_ig_MoI", &RotDynCal::set_ig_MoI,
                 py::arg("ig_rot_MoI"))

            .def("set_lambda", &RotDynCal::set_lambda,
                 py::arg("lambda"))

            .def("set_lambda_high", &RotDynCal::set_lambda_high,
                 py::arg("lambda_high"))

            .def("set_solution_mask", &RotDynCal::set_solution_mask,
                 py::arg("mask"))

            .def("solve", &RotDynCal::solve)

            .def("get_opt_Kd0", rot_dyn_cal::get_opt_Kd0)
            .def("get_opt_Kd1", rot_dyn_cal::get_opt_Kd1)
            .def("get_opt_Kt", rot_dyn_cal::get_opt_Kt)
            .def("get_opt_rot_MoI", rot_dyn_cal::get_opt_rot_MoI)

            .def("get_tau_friction", rot_dyn_cal::get_tau_friction)
            .def("get_alpha_d", rot_dyn_cal::get_alpha_d)
            .def("get_alpha_inertial", rot_dyn_cal::get_alpha_inertial)
            .def("get_alpha_kt", rot_dyn_cal::get_alpha_kt)

            .def("get_tau_motor", rot_dyn_cal::get_tau_motor)
            .def("get_tau_inertial", rot_dyn_cal::get_tau_inertial)

            .def("get_sol_millis", rot_dyn_cal::get_sol_millis)

            .def("get_cal_mask", rot_dyn_cal::get_cal_mask)

            .def("get_lambda", rot_dyn_cal::get_lambda)
            .def("get_lambda_des", rot_dyn_cal::get_lambda_des)
            .def("get_lambda_high", rot_dyn_cal::get_lambda_high)

            .def("get_ig_Kd0", rot_dyn_cal::get_ig_Kd0)
            .def("get_ig_Kd1", rot_dyn_cal::get_ig_Kd1)
            .def("get_ig_Kt", rot_dyn_cal::get_ig_Kt)
            .def("get_ig_MoI", rot_dyn_cal::get_ig_MoI)

            .def("reset_window", &RotDynCal::reset_window)
            .def("is_window_full", &RotDynCal::is_window_full)
            ;

    py::class_<IqEstimator, std::shared_ptr<IqEstimator>>(m, "IqEstimator")
            .def(py::init(construct_iq_estimator),
                 py::arg("K_t"),
                 py::arg("K_d0"),
                 py::arg("K_d1"),
                 py::arg("rot_MoI"),
                 py::arg("red_ratio"),
                 py::arg("alpha") = 10,
                 py::arg("q_dot_3sigma") = 0.001,
                 py::arg("dump_data2mat") = false,
                 py::arg("dump_path") = "/tmp")

            .def("set_current_state", &IqEstimator::set_current_state,
                 py::arg("q_dot"),
                 py::arg("q_ddot"),
                 py::arg("tau"))

            .def("update",
                 py::overload_cast<Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&>(&IqEstimator::update),
                 py::arg("K_d0"),
                 py::arg("K_d1"),
                 py::arg("rot_MoI"),
                 py::arg("K_t"))

            .def("update",
                 py::overload_cast<>(&IqEstimator::update))

            .def("set_log_buffsize", &IqEstimator::set_log_buffsize,
                 py::arg("size"))

            .def("add2log", &IqEstimator::add2log)

            .def("get_iq_estimate", iq_estimator::get_iq_estimate)
            .def("get_iq_meas", iq_estimator::get_iq_estimate)

            ;


    // model interface
    #if defined(WITH_MODEL_INTERFACE)

    py::enum_<Model::ReferenceFrame>(m, "ReferenceFrame", py::arithmetic())
                .value("WORLD", Model::ReferenceFrame::WORLD)
                .value("LOCAL", Model::ReferenceFrame::LOCAL)
                .value("LOCAL_WORLD_ALIGNED", Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);

    py::class_<Model, std::shared_ptr<Model>>(m, "Model")
            .def(py::init(construct_model_interface),
                 py::arg("urdf_path"),
                 py::arg("add_floating_jnt") = false)

            .def("get_robot_mass", get_robot_mass)
            .def("get_nq", &Model::get_nq)
            .def("get_nv", &Model::get_nv)
            .def("get_jnt_number", &Model::get_jnt_number)
            .def("get_jnt_names", &Model::get_jnt_names)

            .def("get_q",  get_q)
            .def("get_v",  get_v)
            .def("get_a",  get_a)
            .def("get_tau",  get_tau)
            .def("get_p",  get_p)
            .def("get_b",  get_b)
            .def("get_B",  get_B)
            .def("get_C",  get_C)
            .def("get_g",  get_g)

            .def("get_J",  get_J)
            .def("get_J_dot",  get_J_dot)

            .def("set_q",  &Model::set_q)
            .def("set_v", &Model::set_v)
            .def("set_a", &Model::set_a)
            .def("set_tau", &Model::set_tau)
            .def("update", &Model::update)

            ;
    #endif

    // power utilities
    #if defined(WITH_XBOT2)


    #endif

}



