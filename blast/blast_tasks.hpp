#pragma once
#include <iostream>
#include "blast_error.hpp"
#include "manipulators/kinova_link6.h"
#include "blast_math.hpp"
namespace blast {

struct Manipulator_state {
    Array PVA;
    bool gripper_state;
    enum trajectory_generation_type { traj_gen_null = 0, traj_gen_blast = 1, traj_gen_straight = 2 };
    enum movement_type { mvt_robot = 0, mvt_gripper = 1, mvt_robot_and_gripper = 2 };
};

// Pick and place task planning
// Inputs :
//      manip       : Link6, used for forward_kinematics and jacobian functions
//      Current_manip_state  :   Manipulator_state of current position, used as first element of Manip_list output.
//      Obj_list    :   Matrix of objects to be picked up as a matrix of cartesian positions of size 7*N_objects
//                      (Last int used as dropbox number)
//      Drop_box_list : List of dropboxes as a matrix of cartesian positions of size 6*N_objects
//      hva         :   Vec3 containing height, velocity, and acceleration of intermediate point (above pick-up zone)
//                      Must all be positive. Acceleration is currently hard-coded to 0, so will have no effect.
// Outputs :
//      Manip_list  :   List of Manipulator_state which indicate fully how the robot is positionned. Going through the list
//                      two elements at a time will give a task that must be completed.
std::list<Manipulator_state> pick_and_place(Link6 &manip, const Manipulator_state &Current_manip_state, const Matrix &Obj_list, const Matrix Drop_box_list, const Vec3 hva) {
    // todo: Forward should give orientation

    // Initialization
    int N_joints = Drop_box_list.rows;
    int N_objects = Obj_list.cols;
    Array h(6);
    h[2] = hva.x;
    Array v(6);
    v[2] = hva.y;
    Array a(6);
    a[2] = hva.z;
    std::list<Manipulator_state> manip_list;
    manip_list.resize(8*N_objects + 2);
    manip_list.push_back(Current_manip_state);

    Manipulator_state manip_state;

    // Open gripper before starting
    manip_state = Current_manip_state;
    manip_state.gripper_state = 0;
    manip_state.mvt_gripper;
    manip_state.traj_gen_null;
    manip_list.push_back(manip_state);

    // The pick-and-place task has four points for each object :
    //  - Above the object
    //  - At the object
    //  - Above the bin
    //  - At the bin

    Array joint_pos1(N_joints);
    Array joint_pos2(N_joints);
    Array joint_pos3(N_joints);
    Array joint_pos4(N_joints);
    Array joint_vel1(N_joints);
    Array joint_vel2(N_joints);
    Array joint_vel3(N_joints);
    Array joint_vel4(N_joints);
    // Array joint_acc1(N_joints);
    // Array joint_acc2(N_joints);
    // Array joint_acc3(N_joints);
    // Array joint_acc4(N_joints);
    Array joint_acc(N_joints, 0);

    for (int i = 0; i < N_objects; i++) {
        // Joint positions
        auto current_obj = Obj_list.col(i);
        int drop_box_id = current_obj.back();
        current_obj.size--;
        // current_obj.back() = [];    // Destroy last element

        joint_pos2 = manip.inverse_kinematics(current_obj);
        joint_pos1 = manip.inverse_kinematics(current_obj + h, joint_pos2);

        Array current_dropbox = Drop_box_list.col(drop_box_id);
        joint_pos4 = manip.inverse_kinematics(current_dropbox);
        joint_pos3 = manip.inverse_kinematics(current_dropbox + h, joint_pos4);

        // Joint velocity
        for (int j = 0; j < N_joints; j++) {
            joint_vel2[j] = 0;
            joint_vel4[j] = 0;
        }
        joint_vel1 = pinv(manip.jacobian(joint_pos1))*v;
        joint_vel3 = pinv(manip.jacobian(joint_pos3))*v;

        // todo: If we use joint acceleration, we need to change the formula for 1 and 3 since it is not right
        // // Joint acceleration
        // for (int j = 0; j < N_joints; j++) {
        //     joint_acc2[j] = 0;
        //     joint_acc4[j] = 0;
        // }
        // joint_acc1 = pinv(jacobian(joint_pos1))*a;
        // joint_acc3 = pinv(jacobian(joint_pos3))*a;

        // Filling manip list
        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos1, -joint_vel1), -joint_acc1);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos1, -joint_vel1), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_blast;
        manip_list.push_back(manip_state);

        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos2, joint_vel2), joint_acc2);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos2, joint_vel2), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);

        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos2, joint_vel2), joint_acc2);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos2, joint_vel2), joint_acc);
        manip_state.gripper_state = 1;
        manip_state.mvt_gripper;
        manip_state.traj_gen_null;
        manip_list.push_back(manip_state);

        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos1, joint_vel1), joint_acc1);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos1, joint_vel1), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);

        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos3, -joint_vel3), -joint_acc3);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos3, -joint_vel3), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_blast;
        manip_list.push_back(manip_state);

        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos4, joint_vel4), joint_acc4);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos4, joint_vel4), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);

        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos4, joint_vel4), joint_acc4);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos4, joint_vel4), joint_acc);
        manip_state.gripper_state = 1;
        manip_state.mvt_gripper;
        manip_state.traj_gen_null;
        manip_list.push_back(manip_state);

        // manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos3, joint_vel3), joint_acc3);
        manip_state.PVA = concatenate_arrays(concatenate_arrays(joint_pos3, joint_vel3), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);
    }
    return manip_list;
}

} // namespace blast