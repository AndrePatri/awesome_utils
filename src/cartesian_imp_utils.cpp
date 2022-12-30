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

void CartesianTask::update(CartTask chi_ref,
            CartTaskDot chi_dot_ref,
            CartTaskDdot chi_ddot_ref)
{
    set_chi_ref(chi_ref);
    set_chi_dot_ref(chi_dot_ref);
    set_chi_ddot_ref(chi_ddot_ref);
}

void CartesianTask::update(utils_defs::PosVec3D pos_ref, utils_defs::RotMat3D rot_ref,
            CartTaskDot chi_dot_ref,
            CartTaskDdot chi_ddot_ref)
{
    set_chi_ref(pos_ref, rot_ref);
    set_chi_dot_ref(chi_dot_ref);
    set_chi_ddot_ref(chi_ddot_ref);
}

CartesianTask::CartTaskErr CartesianTask::task_err(utils_defs::PosVec3D pos, utils_defs::RotMat3D rot)
{

    CartesianTask::CartTaskErr err;

    err.segment(0, 3) = pos - _chi_ref.pos;

    err.segment(3, 3) = RotErr::LogMap(rot, _chi_ref.rot);

    return err;

}

CartesianTask::CartTaskErr CartesianTask::task_err(CartTask cart_task)
{

    CartesianTask::CartTaskErr err;

    err.segment(0, 3) = cart_task.pos - _chi_ref.pos;

    err.segment(3, 3) = RotErr::LogMap(cart_task.rot, _chi_ref.rot);

    return err;

}

CartesianTask::CartTaskDotErr CartesianTask::task_dot_err(CartTaskDot cart_task_dot)
{

    return (cart_task_dot - _chi_dot_ref);

}

CartesianTask::CartTaskDdotErr CartesianTask::task_ddot_err(CartTaskDdot cart_task_ddot)
{

    return (cart_task_ddot - _chi_ddot_ref);

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
}

CartesianImpController::CartesianImpController(Model::Ptr model_ptr,
                                               CartesianTask::Ptr cart_task)
{
    _model_ptr = model_ptr;

    _cart_task = cart_task;

    _nq = _model_ptr->get_nq();

    _nv = _model_ptr->get_nv();

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

void CartesianImpController::update()
{
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

    compute_quantities();

}

void CartesianImpController::update(std::string cart_cntrl_framename)
{

    if(_auto_critical_damp)
    {// this is triggered depending on which overload of
     // set_cart_impedance was last called
        compute_critically_damped_gains();
    }

    compute_quantities();

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
