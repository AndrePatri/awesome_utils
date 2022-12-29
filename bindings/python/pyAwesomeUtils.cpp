#include "include/awesome_utils/calib_utils.hpp"
#include "include/awesome_utils/sign_proc_utils.hpp"
#include "include/awesome_utils/contact_est_utils.hpp"
#include "include/awesome_utils/model_interface.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace ModelInterface;

//************* Model interface bindings *************//

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

PYBIND11_MODULE(awesome_pyutils, m) {

    py::class_<Model, std::shared_ptr<Model>>(m, "Model")
            .def(py::init(construct_model_interface),
                 py::arg("urdf_path"),
                 py::arg("add_floating_jnt") = false)

            .def("get_robot_mass", get_robot_mass)
            .def("get_nq", &Model::get_nq)
            .def("get_nv", &Model::get_nv)
            .def("get_jnt_number", &Model::get_jnt_number)

            .def("get_q",  get_q)
            .def("get_v",  get_v)
            .def("get_a",  get_a)
            .def("get_tau",  get_tau)
            .def("get_p",  get_p)
            .def("get_b",  get_b)
            .def("get_B",  get_B)
            .def("get_C",  get_C)
            .def("get_g",  get_g)

            .def("set_q",  &Model::set_q)
            .def("set_v", &Model::set_v)
            .def("set_a", &Model::set_a)
            .def("set_tau", &Model::set_tau)
            .def("update", &Model::update)

            ;

}



