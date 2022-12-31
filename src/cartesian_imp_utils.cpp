#include "include/awesome_utils/cartesian_imp_utils.hpp"

#include "include/awesome_utils/orientation_utils.hpp"

using namespace CartesianImpUtils;


//************* CartesianTask *************//

void CartesianTask::set_chi_ref(CartTask chi_ref)
{
    _chi_ref = chi_ref;
}

void CartesianTask::set_chi_ref(utils_defs::PosVec3D pos_ref,
                           utils_defs::RotMat3D rot_ref)
{
    _chi_ref.pos = pos_ref;
    _chi_ref.rot = rot_ref;
}

void CartesianTask::set_chi_dot_ref(CartTaskDot chi_dot_ref)
{
    _chi_dot_ref = chi_dot_ref;
}

void CartesianTask::set_chi_ddot_ref(CartTaskDdot chi_ddot_ref)
{
    _chi_ddot_ref = chi_ddot_ref;
}

void CartesianTask::set_chi_meas(CartTask chi_meas)
{
    _chi_meas = chi_meas;
}

void CartesianTask::set_chi_meas(utils_defs::PosVec3D pos_meas,
                           utils_defs::RotMat3D rot_meas)
{
    _chi_meas.pos = pos_meas;
    _chi_meas.rot = rot_meas;
}

void CartesianTask::set_chi_dot_meas(CartTaskDot chi_dot_meas)
{
    _chi_dot_meas = chi_dot_meas;
}

void CartesianTask::set_chi_ddot_meas(CartTaskDdot chi_ddot_meas)
{
    _chi_ddot_meas = chi_ddot_meas;
}

void CartesianTask::update_ref(CartTask chi_ref,
            CartTaskDot chi_dot_ref,
            CartTaskDdot chi_ddot_ref)
{
    set_chi_ref(chi_ref);
    set_chi_dot_ref(chi_dot_ref);
    set_chi_ddot_ref(chi_ddot_ref);
}

void CartesianTask::update_ref(utils_defs::PosVec3D pos_ref, utils_defs::RotMat3D rot_ref,
            CartTaskDot chi_dot_ref,
            CartTaskDdot chi_ddot_ref)
{
    set_chi_ref(pos_ref, rot_ref);
    set_chi_dot_ref(chi_dot_ref);
    set_chi_ddot_ref(chi_ddot_ref);
}

void CartesianTask::update_meas(CartTask chi_meas,
            CartTaskDot chi_dot_meas,
            CartTaskDdot chi_ddot_meas)
{
    set_chi_meas(chi_meas);
    set_chi_dot_meas(chi_dot_meas);
    set_chi_ddot_meas(chi_ddot_meas);
}

void CartesianTask::update_meas(utils_defs::PosVec3D pos_meas, utils_defs::RotMat3D rot_meas,
            CartTaskDot chi_dot_meas,
            CartTaskDdot chi_ddot_meas)
{
    set_chi_meas(pos_meas, rot_meas);
    set_chi_dot_meas(chi_dot_meas);
    set_chi_ddot_meas(chi_ddot_meas);
}

CartesianTask::CartTaskErr CartesianTask::task_err()
{

    CartTaskErr err;

    err.segment(0, 3) = _chi_meas.pos - _chi_ref.pos;

    err.segment(3, 3) = RotErr::LogMap(_chi_meas.rot, _chi_ref.rot);

    return err;

}

CartesianTask::CartTaskDotErr CartesianTask::task_dot_err()
{

    return (_chi_dot_meas - _chi_dot_ref);

}

CartesianTask::CartTaskDdotErr CartesianTask::task_ddot_err()
{

    return (_chi_ddot_meas - _chi_ddot_ref);

}

CartesianTask::CartTaskDdot CartesianTask::task_ddot_ref()
{
    return _chi_ddot_ref;
}

void CartesianTask::get_task(CartTask& task_ref, CartTask& task_meas)
{
    task_ref = _chi_ref;
    task_meas = _chi_ref;
}

void CartesianTask::get_task_dot(CartTaskDot& task_dot_ref, CartTaskDot& task_dot_meas)
{
    task_dot_ref = _chi_dot_ref;
    task_dot_meas = _chi_dot_meas;
}

void CartesianTask::get_task_ddot(CartTaskDdot& task_ddot_ref, CartTaskDdot& task_ddot_meas)
{
    task_ddot_ref = _chi_ddot_ref;
    task_ddot_meas = _chi_ddot_meas;
}

//************* CartesianImpController *************//

CartesianImpController::CartesianImpController()
{

}

CartesianImpController::CartesianImpController(Model::Ptr model_ptr,
                                               CartesianTask::Ptr cart_task,
                                               std::string cart_cntrl_framename)
{
    _model_ptr = model_ptr;

    _cart_task = cart_task;

    _nq = _model_ptr->get_nq();

    _nv = _model_ptr->get_nv();

    _cart_cntrl_framename = cart_cntrl_framename;

    _was_cntrl_framename_set = true;

    // some initializations
    _q = VectorXd::Zero(_nq);
    _v = VectorXd::Zero(_nv);
    _B_inv = MatrixXd::Zero(_nv, _nv);
    _C = MatrixXd::Zero(_nv, _nv);
    _g = VectorXd::Zero(_nv);
    _tau_d = VectorXd::Zero(_nv);
    _tau_c = VectorXd::Zero(_nv);
    _tau_cmd = VectorXd::Zero(_nv);

}

CartesianImpController::CartesianImpController(Model::Ptr model_ptr,
                                               CartesianTask::Ptr cart_task)
{
    _model_ptr = model_ptr;

    _cart_task = cart_task;

    _nq = _model_ptr->get_nq();

    _nv = _model_ptr->get_nv();

    // some initializations
    _q = VectorXd::Zero(_nq);
    _v = VectorXd::Zero(_nv);
    _B_inv = MatrixXd::Zero(_nv, _nv);
    _C = MatrixXd::Zero(_nv, _nv);
    _g = VectorXd::Zero(_nv);
    _tau_d = VectorXd::Zero(_nv);
    _tau_c = VectorXd::Zero(_nv);
    _tau_cmd = VectorXd::Zero(_nv);

}

void CartesianImpController::map_impedance_vect2mat()
{
    // diagonal
    _cart_stiff(0, 0) = _cart_stiff_vect(0);
    _cart_stiff(1, 1) = _cart_stiff_vect(1);
    _cart_stiff(2, 2) = _cart_stiff_vect(2);
    _cart_damp(0, 0) = _cart_damp_vect(0);
    _cart_damp(1, 1) = _cart_damp_vect(1);
    _cart_damp(2, 2) = _cart_damp_vect(2);
    // upper triangular part
    _cart_stiff(3, 3) = _cart_stiff_vect(3);
    _cart_stiff(4, 4) = _cart_stiff_vect(4);
    _cart_stiff(5, 5) = _cart_stiff_vect(5);
    _cart_damp(3, 3) = _cart_damp_vect(3);
    _cart_damp(4, 4) = _cart_damp_vect(4);
    _cart_damp(5, 5) = _cart_damp_vect(5);
}

void CartesianImpController::update(VectorXd tau_d, VectorXd tau_c)
{

    if (tau_d.size() != tau_c.size())
    {
        std::string exception = std::string("CartesianImpUtils::CartesianImpController::update(): tau_d size does not match tau_c");

        throw std::invalid_argument(exception);
    }

    if (tau_d.size() != _nv)
    {
        std::string exception = std::string("CartesianImpUtils::CartesianImpController::update(): tau_d and tau_c should have size equal to nv");

        throw std::invalid_argument(exception);
    }

    if (!_was_cntrl_framename_set)
    {
        std::string exception = std::string("CartesianImpUtils::CartesianImpController::update(): you have to set a control frame name!");

        throw std::invalid_argument(exception);
    }

    if(_auto_critical_damp)
    { // this is triggered depending on which overload of
      // set_cart_impedance was last called
        compute_critically_damped_gains();
    }

    update_internal_states(tau_d, tau_c);

    compute_tau_cmd();

}

void CartesianImpController::update(std::string cart_cntrl_framename,
                                    VectorXd tau_d, VectorXd tau_c)
{

    if(_auto_critical_damp)
    {// this is triggered depending on which overload of
     // set_cart_impedance was last called
        compute_critically_damped_gains();
    }

    update_internal_states(tau_d, tau_c);

    compute_tau_cmd();

}

void CartesianImpController::update_internal_states(VectorXd tau_d, VectorXd tau_c)
{
    compute_quantities();

    _tau_d = tau_d;
    _tau_c = tau_c;

    _model_ptr->get_q(_q);
    _model_ptr->get_v(_v);
}

void CartesianImpController::compute_critically_damped_gains()
{
    for (int i = 0; i < _cart_damp_vect.size(); i++)
    {
        _cart_damp_vect(i) = 2 * std::sqrt(_Lambda(i, i) * _cart_stiff_vect(i, i));
    }
}

void CartesianImpController::set_cart_impedance(utils_defs::CartStiffMat stifness_mat,
                        utils_defs::CartDampMat damping_mat)
{
    _auto_critical_damp = false; // if we call this method, we want to set the impedance
    // manually

    // we make sure inputs are symmetric
    _cart_stiff = 0.5 * (stifness_mat + stifness_mat.transpose());

    _cart_damp = 0.5 * (damping_mat + damping_mat.transpose());
}

void CartesianImpController::set_cart_impedance(utils_defs::CartStiffVect stifness_vect,
                        utils_defs::CartDampVect damping_vect)
{
    _auto_critical_damp = false; // if we call this method, we want to set the impedance
    // manually

    _cart_stiff_vect = stifness_vect;
    _cart_damp_vect = damping_vect;

    // since we assume diagonal imp. matrices
    // we reset them to 0 before adding the diagonal in case
    // theie off-diagonal elements different from 0 before

    _cart_stiff = Eigen::MatrixXd::Zero(6, 6);
    _cart_damp = Eigen::MatrixXd::Zero(6, 6);

    map_impedance_vect2mat(); // updates _cart_stiff and _cart_damp
}

void CartesianImpController::set_cart_impedance(utils_defs::CartStiffVect stifness_vect)
{
    _auto_critical_damp = true; // if we call this method, we want to set the impedance
    // automatically --> the next call to the update() method will compute approximately
    // critically damped gains

    _cart_stiff_vect = stifness_vect;

}

void CartesianImpController::compute_lambda_inv()
{
    _Lambda_inv = _J * _B_inv * _J.transpose();
}

void CartesianImpController::compute_lambda()
{
    _Lambda = _Lambda_inv.inverse();
}

void CartesianImpController::compute_J_rps_w()
{
    _J_rps_w = _B_inv * _J.transpose() * _Lambda;
}

void CartesianImpController::compute_quantities()
{
    _model_ptr->get_jac(_cart_cntrl_framename,
                        _J,
                        Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);

    _model_ptr->get_jac_dot(_cart_cntrl_framename,
                        _J_dot,
                        Model::ReferenceFrame::LOCAL_WORLD_ALIGNED);

    _model_ptr->get_B_inv(_B_inv);

    _model_ptr->get_C(_C);

    _model_ptr->get_g(_g);

    compute_lambda_inv();

    compute_lambda();

    compute_J_rps_w();

}

void CartesianImpController::compute_h()
{
    _h = (_J_rps_w.transpose() * _C - _Lambda * _J_dot) * _v;
}

void CartesianImpController::compute_f_star()
{
    compute_h();

    _f_star = _Lambda * _cart_task->task_ddot_ref() + _h + _J_rps_w.transpose() * _g -
              (_cart_stiff * _cart_task->task_err() + _cart_damp * _cart_task->task_dot_err());
}

void CartesianImpController::compute_tau_cmd()
{
    compute_f_star();

    _tau_cmd = - _tau_d - _tau_c + _J.transpose() * _f_star;
}

VectorXd CartesianImpController::h()
{
    return _h;
}

utils_defs::Wrench CartesianImpController::f_star()
{
    return _f_star;
}

utils_defs::JacRightPseudoInv CartesianImpController::J_rps_w()
{
    return _J_rps_w;
}

utils_defs::CartInertiaMat CartesianImpController::Lambda()
{
    return _Lambda;
}

utils_defs::CartInertiaMat CartesianImpController::Lambda_inv()
{
    return _Lambda_inv;
}

VectorXd CartesianImpController::tau_cmd()
{
    return _tau_cmd;
}
