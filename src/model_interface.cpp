#include "awesome_utils/model_interface.hpp"

using namespace ModelInterface;

Model::Model()
{

}

Model::Model(std::string urdf_path, bool add_floating_jnt)
{
    try{

        if (add_floating_jnt)
        { // add floating joint --> useful if working with a floating base robot URDF
          // lacking a floating base root joint
            pinocchio::JointModelFreeFlyer floating_joint;
            pinocchio::urdf::buildModel(urdf_path, floating_joint, _pin_model);
        }
        else{

             pinocchio::urdf::buildModel(urdf_path, _pin_model);

        }

    }
    catch (int err_code){

        _pin_model_init_ok = false;
    }

    _pin_model_init_ok = true;
    _urdf_path = urdf_path;

    _pin_data = pinocchio::Data(_pin_model);

//    _q = pinocchio::randomConfiguration(_pin_model);

    _n_jnts = _pin_model.njoints;
    _jnt_names = _pin_model.names;

    _nq = _pin_model.nq;
    _nv = _pin_model.nv;

    _B = MatrixXd::Zero(_nv, _nv);
    _C = MatrixXd::Zero(_nv, _nv);
    _g = VectorXd::Zero(_nv);
    _p = VectorXd::Zero(_nv);
    _tau = VectorXd::Zero(_nv);
    _q = VectorXd::Zero(_nq);
    _v = VectorXd::Zero(_nv);
    _a = VectorXd::Zero(_nv);

    _q_min = _pin_model.lowerPositionLimit;
    _q_max = _pin_model.upperPositionLimit;

    _nqs = std::vector<int>(_n_jnts);
    _nvs = std::vector<int>(_n_jnts);

    _nqs = _pin_model.nqs;
    _nvs = _pin_model.nvs;

    _effort_limits = _pin_model.effortLimit;
    _vel_limits = _pin_model.velocityLimit;

}

bool Model::was_model_init_ok()
{
    return _pin_model_init_ok;
}

std::vector<std::string> Model::get_jnt_names()
{
    return _jnt_names;
}

void Model::update(VectorXd q, VectorXd v, VectorXd tau)
{
    _q = q;
    _v = v;
    _tau = tau;

    update_all();

}

void Model::update(VectorXd q, VectorXd v, VectorXd tau, VectorXd a)
{
    _q = q;
    _v = v;
    _tau = tau;
    _a = a;

    update_all();

}

void Model::update_all()
{
    // update quantities
    B();
    C();
    g();
    tau();
    p();
    b();
}

void Model::B()
{
    pinocchio::crba(_pin_model, _pin_data, _q); // only computes the upper triangular part of B
    // we copy the upper triangular half into the lower
    _pin_data.M.triangularView<Eigen::StrictlyLower>() = _pin_data.M.transpose().triangularView<Eigen::StrictlyLower>();

    _B = _pin_data.M; // we get the full joint-space inertia matrix
}

void Model::get_B(MatrixXd& B)
{
    B = _B;
}

void Model::C()
{
    pinocchio::computeCoriolisMatrix(_pin_model, _pin_data, _q, _v);

    _C = _pin_data.C;
}

void Model::get_C(MatrixXd& C)
{
    C = _C;
}

void Model::g()
{
    pinocchio::computeGeneralizedGravity(_pin_model, _pin_data, _q);
    _g = _pin_data.g;
}

void Model::get_g(VectorXd& g)
{
    g = _g;
}

void Model::tau()
{
    _tau = _pin_data.tau;
}

void Model::get_tau(VectorXd& tau)
{
    tau = _tau;
}

void Model::p()
{
    _p = _B * _v; // joint space momentum of the system
}

void Model::get_p(VectorXd& p)
{
    p = _p;
}

void Model::b()
{
    _b = _C * _v; // bias forces
}

void Model::get_b(VectorXd& b)
{
    b = _b;
}

void Model::get_state(VectorXd& q, VectorXd& v, VectorXd& a,
                      VectorXd& tau)
{
    q = _q;
    v = _v;
    a = _a;
    tau = _tau;
}

int Model::get_nq()
{
    return _nq;
}

int Model::get_nv()
{
    return _nv;
}

int Model::get_jnt_number()
{
   return _n_jnts;
}

std::string Model::get_urdf_path()
{
    return _urdf_path;
}

void Model::get_jnt_lim(VectorXd& q_min, VectorXd& q_max)
{
    q_min = _q_min;
    q_max = _q_max;
}

void Model::rnea()
{

}

void Model::aba()
{

}

void Model::crba()
{

}

void Model::centroidal_dyn()
{

}

void Model::jac()
{

}

void Model::kin_energy()
{

}

void Model::pot_energy()
{

}
