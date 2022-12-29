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

CartesianTask::CartTaskErr CartesianTask::compute_task_err(utils_defs::PosVec3D pos, utils_defs::RotMat3D rot)
{

    CartesianTask::CartTaskErr err;

    err.segment(0, 3) = pos - _chi_ref.pos;

    err.segment(3, 3) = RotErr::LogMap(rot, _chi_ref.rot);

    return err;

}

CartesianTask::CartTaskErr CartesianTask::compute_task_err(CartTask cart_task)
{

    CartesianTask::CartTaskErr err;

    err.segment(0, 3) = cart_task.pos - _chi_ref.pos;

    err.segment(3, 3) = RotErr::LogMap(cart_task.rot, _chi_ref.rot);

    return err;

}

CartesianTask::CartTaskDotErr CartesianTask::compute_task_dot_err(CartTaskDot cart_task_dot)
{

    return (cart_task_dot - _chi_dot_ref);

}

CartesianTask::CartTaskDdotErr CartesianTask::compute_task_ddot_err(CartTaskDdot cart_task_ddot)
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

}

void CartesianImpController::set_cart_impedance(CartStiffMat stifness_mat,
                        CartDampMat damping_mat)
{
    // we make sure inputs are symmetric
    _cart_stiff = 0.5 * (stifness_mat + stifness_mat.transpose());

    _cart_damp = 0.5 * (damping_mat + damping_mat.transpose());
}

void CartesianImpController::set_cart_impedance(CartStiffVect stifness_vect,
                        CartDampVect damping_vect)
{
    _cart_stiff_vect = stifness_vect;
    _cart_damp_vect = damping_vect;

    // since we assume diagonal imp. matrices
    // we reset them to 0 before adding the diagonal in case
    // theie off-diagonal elements different from 0 before

    _cart_stiff = Eigen::MatrixXd::Zero(6, 6);
    _cart_damp = Eigen::MatrixXd::Zero(6, 6);

    map_impedance_vect2mat(); // updates _cart_stiff and _cart_damp
}

void CartesianImpController::compute_cart_inertia_mat()
{
}
