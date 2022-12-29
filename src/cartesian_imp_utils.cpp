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

void CartesianImpController::update()
{

}

