// Copyright (C) 2023  Andrea Patrizi (AndrePatri, andreapatrizi1b6e6@gmail.com)
// 
// This file is part of awesome_utils and distributed under the General Public License version 2 license.
// 
// awesome_utils is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// awesome_utils is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with awesome_utils.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "model_interface.hpp"
#include <pinocchio/spatial/se3.hpp>

using namespace ModelInterface;

Model::Model(std::string urdf_path, bool add_floating_jnt)
{
    try{

        if (add_floating_jnt)
        { // add floating joint --> useful if working with a floating base robot URDF
          // lacking a floating base root joint
            pinocchio::JointModelFreeFlyer floating_joint;
            pinocchio::urdf::buildModel(urdf_path, floating_joint, _pin_model);
        }
        else
        {

            pinocchio::urdf::buildModel(urdf_path, _pin_model);

        }

    }
    catch (int err_code){

        _pin_model_init_ok = false;
    }

    _urdf_path = urdf_path;

    _pin_data = pinocchio::Data(_pin_model);

    _n_jnts = _pin_model.njoints;
    _jnt_names = _pin_model.names;

    _nq = _pin_model.nq;
    _nv = _pin_model.nv;

    _B = MatrixXd::Zero(_nv, _nv);
    _B_inv = MatrixXd::Zero(_nv, _nv);
    _C = MatrixXd::Zero(_nv, _nv);
    _b = VectorXd::Zero(_nv);
    _g = VectorXd::Zero(_nv);
    _p = VectorXd::Zero(_nv);
    _tau = VectorXd::Zero(_nv);
    _q = VectorXd::Zero(_nq);
    _v = VectorXd::Zero(_nv);
    _a = VectorXd::Zero(_nv);

    _q_min = _pin_model.lowerPositionLimit;
    _q_max = _pin_model.upperPositionLimit;

    _q_min_aux = VectorXd::Zero(_q_min.size());
    _q_max_aux = VectorXd::Zero(_q_max.size());
    _q_min_aux << -1.0, -1.0, -1.0, _q_min.tail(_q_min.size() - 3);
    _q_max_aux << 1.0, 1.0, 1.0, _q_max.tail(_q_max.size() - 3);

    _nqs = std::vector<int>(_n_jnts);
    _nvs = std::vector<int>(_n_jnts);

    _nqs = _pin_model.nqs;
    _nvs = _pin_model.nvs;

    _effort_limits = _pin_model.effortLimit;
    _vel_limits = _pin_model.velocityLimit;

    pinocchio::computeTotalMass(_pin_model, _pin_data);
    _mass = _pin_data.mass[0];

    if (_jnt_names[0] == "universe" && _nqs[1] == 7 && _nvs[1] == 6)
    { // in case the universe joint is present, look at the second joint

        _is_floating_base = true;

    }

    if (_jnt_names[0] != "universe" && _nqs[0] == 7 && _nvs[0] == 6)
    { // in case the universe joint is not present, look directly at the first joint

        _is_floating_base = true;

    }

    _pin_model_init_ok = true;

}

bool Model::frame_exists(std::string& framename)
{

    return _pin_model.existFrame(framename);

}

void Model::get_robot_mass(double& mass)
{
    mass = _mass;
}

bool Model::was_model_init_ok()
{
    return _pin_model_init_ok;
}

std::vector<std::string> Model::get_jnt_names()
{
    return _jnt_names;
}

void Model::update()
{

    update_all();

}

void Model::set_neutral()
{
    _q = pinocchio::neutral(_pin_model);
}

void Model::set_random()
{

    if(_is_floating_base)
    {
        _q = pinocchio::randomConfiguration(_pin_model, _q_min_aux, _q_max_aux);
    }
    else
    {
        _q = pinocchio::randomConfiguration(_pin_model, _q_min, _q_max);
    }

}

void Model::set_q(VectorXd& q)
{

    _q = q;

}

void Model::set_v(VectorXd& v)
{
    _v = v;
}

void Model::set_a(VectorXd& a)
{
    _a = a;
}

void Model::set_tau(VectorXd& tau)
{
    _tau = tau;
}

void Model::update_all()
{

    update_frames_forward_kin();

    update_forward_kin();

    // update dynamics quantities
    update_B(); // joint-space inertia matrix
    update_B_inv(); // inverse of joint-space inertia matrix
    update_C(); // coriolis, centrifugal effects
    update_g(); // gravitational effects
    update_p(); // joint-space momentum (aka generalized momentum)
    update_b(); // bias vector (aka C * v)
}

void Model::update_frames_forward_kin()
{

    pinocchio::framesForwardKinematics(_pin_model, _pin_data, _q); // forward kinematics update

}

void Model::update_forward_kin()
{
    pinocchio::forwardKinematics(_pin_model,
                                 _pin_data,
                                 _q,
                                 _v);
}

void Model::update_B()
{
    B();
}

void Model::update_B_inv()
{
    B_inv();
}

void Model::B()
{
    pinocchio::crba(_pin_model, _pin_data, _q); // only computes the upper triangular part of B
    // we copy the upper triangular half into the lower
    _pin_data.M.triangularView<Eigen::StrictlyLower>() = _pin_data.M.transpose().triangularView<Eigen::StrictlyLower>();

    _B = _pin_data.M; // we get the full joint-space inertia matrix
}

void Model::B_inv()
{
    pinocchio::computeMinverse(_pin_model, _pin_data, _q); // only computes the upper triangular part
    _pin_data.Minv.triangularView<Eigen::StrictlyLower>() = _pin_data.Minv.transpose().triangularView<Eigen::StrictlyLower>();

    _B_inv = _pin_data.Minv; // we get the full joint-space inertia matrix
}

void Model::get_B(MatrixXd& B)
{
    B = _B;
}

void Model::get_B_inv(MatrixXd& B_inv)
{
    B_inv = _B_inv;
}

void Model::update_C()
{
    C();
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

void Model::update_g()
{
    g();
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

void Model::update_tau_rnea()
{
    rnea();
}

void Model::get_tau(VectorXd& tau)
{
    tau = _tau;
}

void Model::get_q(VectorXd& q)
{
    q = _q;
}

void Model::get_v(VectorXd& v)
{
    v = _v;
}

void Model::get_a(VectorXd& a)
{
    a = _a;
}

void Model::update_p()
{
    p();
}

void Model::p()
{
    _p = _B * _v; // joint space momentum of the system
}

void Model::get_p(VectorXd& p)
{
    p = _p;
}

void Model::update_b()
{
    b();
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


void Model::jacobian(std::string& frame_name, Model::ReferenceFrame ref,
                    utils_defs::SpatialJac& J)
{

    bool does_frame_exist = frame_exists(frame_name);

    if (!does_frame_exist)
    {
        std::string exception = std::string("ModelInterface::Model::jacobian(): the provided frame \"") +
                                std::string(frame_name) + std::string("\" does not exist!");

        throw std::invalid_argument(exception);
    }

    pinocchio::FrameIndex frame_idx = _pin_model.getFrameId(frame_name);

    J = MatrixXd::Zero(6, _nv);

    pinocchio::computeFrameJacobian(_pin_model, _pin_data, _q, frame_idx, pinocchio::ReferenceFrame(ref), J); // using latest set q

}

void Model::jacobian_dot(std::string& frame_name, Model::ReferenceFrame ref,
              utils_defs::SpatialJacDot& J_dot)
{

    bool does_frame_exist = frame_exists(frame_name);

    if (!does_frame_exist)
    {
        std::string exception = std::string("ModelInterface::Model::jacobian_dot(): the provided frame \"") +
                                std::string(frame_name) + std::string("\" does not exist!");

        throw std::invalid_argument(exception);
    }

    pinocchio::FrameIndex frame_idx = _pin_model.getFrameId(frame_name);

    J_dot = utils_defs::SpatialJacDot(6, _nv);

    pinocchio::computeJointJacobiansTimeVariation(_pin_model, _pin_data, _q, _v); // this is needed
    // before the call to getFrameJacobianTimeVariation

    pinocchio::getFrameJacobianTimeVariation(_pin_model,
                                             _pin_data,
                                             frame_idx,
                                             pinocchio::ReferenceFrame(ref),
                                             J_dot);

}

void Model::get_jac(std::string& frame_name,
                    utils_defs::SpatialJac& J,
                    Model::ReferenceFrame ref)
{

    jacobian(frame_name, ref, J);

}

void Model::get_jac_dot(std::string& frame_name,
             utils_defs::SpatialJacDot& J_dot,
             ReferenceFrame ref)
{

    jacobian_dot(frame_name, ref, J_dot);

}


void Model::get_frame_vel(std::string& frame_name,
                          utils_defs::Twist& vel,
                          ReferenceFrame ref)
{

    bool does_frame_exist = _pin_model.existFrame(frame_name);

    if (!does_frame_exist)
    {
        std::string exception = std::string("ModelInterface::Model::get_frame_vel(): the provided frame \"") +
                                std::string(frame_name) + std::string("\" does not exist!");

        throw std::invalid_argument(exception);
    }

    pinocchio::FrameIndex frame_idx = _pin_model.getFrameId(frame_name);

    auto full_vel = pinocchio::getFrameVelocity(_pin_model,
                     _pin_data,
                     frame_idx,
                     pinocchio::ReferenceFrame(ref));

    vel.segment(0, 3) = full_vel.linear();
    vel.segment(3, 3) = full_vel.angular();

}

void Model::get_frame_vel(std::string& frame_name,
                          utils_defs::LinVel& lin_vel, utils_defs::AngVel& omega,
                          ReferenceFrame ref)
{

    bool does_frame_exist = _pin_model.existFrame(frame_name);

    if (!does_frame_exist)
    {
        std::string exception = std::string("ModelInterface::Model::get_frame_vel(): the provided frame \"") +
                                std::string(frame_name) + std::string("\" does not exist!");

        throw std::invalid_argument(exception);
    }

    pinocchio::FrameIndex frame_idx = _pin_model.getFrameId(frame_name);

    auto full_vel = pinocchio::getFrameVelocity(_pin_model,
                     _pin_data,
                     frame_idx,
                     pinocchio::ReferenceFrame(ref));

    lin_vel = full_vel.linear();
    omega = full_vel.angular();

}

void Model::get_frame_pose(std::string& frame_name,
                    utils_defs::Affine3D& pose)
{
    bool does_frame_exist = _pin_model.existFrame(frame_name);

    if (!does_frame_exist)
    {
        std::string exception = std::string("ModelInterface::Model::get_frame_pose(): the provided frame \"") +
                                std::string(frame_name) + std::string("\" does not exist!");

        throw std::invalid_argument(exception);
    }

    pinocchio::FrameIndex frame_idx = _pin_model.getFrameId(frame_name);

    pose.translation() = _pin_data.oMf.at(frame_idx).translation();
    pose.linear() = _pin_data.oMf.at(frame_idx).rotation();

}

void Model::get_frame_pose(std::string& frame_name,
                    utils_defs::PosVec3D& position, utils_defs::RotMat3D& rotation)
{

    bool does_frame_exist = _pin_model.existFrame(frame_name);

    if (!does_frame_exist)
    {
        std::string exception = std::string("ModelInterface::Model::get_frame_pose(): the provided frame \"") +
                                std::string(frame_name) + std::string("\" does not exist!");

        throw std::invalid_argument(exception);
    }

    pinocchio::FrameIndex frame_idx = _pin_model.getFrameId(frame_name);

    position = _pin_data.oMf.at(frame_idx).translation();
    rotation = _pin_data.oMf.at(frame_idx).rotation();

}

void Model::rnea()
{

    pinocchio::rnea(_pin_model,
         _pin_data,
         _q,
         _v,
         _a);

    _tau = _pin_data.tau;

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

void Model::kin_energy()
{

}

void Model::pot_energy()
{

}
