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

//************* Calibration utilities bindings *************//

using namespace CalibUtils;

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

auto get_sol_millis = [](RotDynCal& self)
{
    Eigen::VectorXd sol_millis;

    self.get_sol_millis(sol_millis);

    return sol_millis;
};

auto get_alpha_d = [](RotDynCal& self)
{
    Eigen::VectorXd alpha_d0;
    Eigen::VectorXd alpha_d1;

    self.get_alpha_d(alpha_d0, alpha_d1);

    auto output = std::make_tuple(alpha_d0, alpha_d1);

    return output;
};

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

#if defined(WITH_XBOT2)

//************* Power utils bindings *************//

#endif

PYBIND11_MODULE(awesome_pyutils, m) {

    // signal processing utilities

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

            .def("get_opt_Kd0", &RotDynCal::get_opt_Kd0,
                 py::arg("Kd0_opt"))
            .def("get_opt_Kd1", &RotDynCal::get_opt_Kd1,
                 py::arg("Kd1_opt"))
            .def("get_opt_Kt", &RotDynCal::get_opt_Kt,
                 py::arg("Kt"))
            .def("get_opt_rot_MoI", &RotDynCal::get_opt_rot_MoI,
                 py::arg("rot_MoI"))

            .def("get_tau_friction", &RotDynCal::get_tau_friction,
                 py::arg("tau_friction"))
            .def("get_alpha_d", get_alpha_d)
            .def("get_alpha_inertial", &RotDynCal::get_alpha_inertial,
                 py::arg("alpha_inertial"))
            .def("get_alpha_kt", &RotDynCal::get_alpha_kt,
                 py::arg("alpha_kt"))

            .def("get_tau_motor", &RotDynCal::get_tau_motor,
                 py::arg("tau_mot"))
            .def("get_tau_inertial", &RotDynCal::get_tau_inertial,
                 py::arg("tau_inertial"))

            .def("get_sol_millis", get_sol_millis)

            .def("get_cal_mask", &RotDynCal::get_cal_mask,
                 py::arg("cal_mask"))

            .def("get_lambda", &RotDynCal::get_lambda,
                 py::arg("lambda"))
            .def("get_lambda_des", &RotDynCal::get_lambda_des,
                 py::arg("lambda_des"))
            .def("get_lambda_high", &RotDynCal::get_lambda_high,
                 py::arg("lambda_high"))

            .def("get_ig_Kd0", &RotDynCal::get_ig_Kd0,
                 py::arg("ig_Kd0"))
            .def("get_ig_Kd1", &RotDynCal::get_ig_Kd1,
                 py::arg("ig_Kd1"))
            .def("get_ig_Kt", &RotDynCal::get_ig_Kt,
                 py::arg("ig_Kt"))
            .def("get_ig_MoI", &RotDynCal::get_ig_MoI,
                 py::arg("ig_rot_MoI"))

            .def("reset_window", &RotDynCal::reset_window)
            .def("is_window_full", &RotDynCal::is_window_full)
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



