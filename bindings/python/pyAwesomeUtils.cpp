#include "awesome_utils/calib_utils.hpp"
#include "awesome_utils/sign_proc_utils.hpp"
#include "awesome_utils/contact_est_utils.hpp"
#include "awesome_utils/model_interface.hpp"
#include "awesome_utils/traj_utils.hpp"
#include "awesome_utils/xbot2_utils.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

auto construct_model_interface = [](std::string urdf_path, bool add_floating_jnt)
{

    Model::Ptr model_ptr;

    model_ptr.reset(new ModelInterface::Model(urdf_path, add_floating_jnt));

    return model_ptr;
};

PYBIND11_MODULE(awesome_pyutils, m) {

    py::class_<Model, std::shared_ptr<ModelInterface::Model>>(m, "Model")
            .def(py::init(construct_model_interface),
                 py::arg("urdf_path"),
                 py::arg("add_floating_jnt") = false)
            ;

}



