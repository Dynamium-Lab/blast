#pragma once
#include <iostream>
#include "blast_manipulator.hpp"
#include "blast_math.hpp"
using namespace blast;

// Pick and place task planning
// Inputs : 
//      Current_pva :   Manipulator_state of current position, used as first element of Manip_list output.
//      Obj_list    :   Matrix of objects to be picked up as a matrix of joint positions of size N_objects*(N_joints + 1)
//                      (Last int used as dropbox number)
//      Drop_box_list : List of dropboxes as a matrix of joint positions of size N_objects*N_joints
//      hva         :   Vec3 containing height, velocity, and acceleration of intermediate point (above pick-up zone)
//                      Must all be positive
// Outputs :
//      Manip_list  :   List of Manipulator_state which indicate fully how the robot is positionned. Going through the list 
//                      two elements at a time will give a task that must be completed.
std::list<Manipulator_state> Link6::pick_and_place(const Manipulator_state &Current_manip_state, const Matrix &Obj_list, const Matrix Drop_box_list, const Vec3 hva) {
    // Initialization
    int N_joints = Drop_box_list.cols;
    int N_objects = Obj_list.rows;
    Array h(3, 0);
    h[2] = hva.x;
    Array v(3, 0);
    v[2] = hva.y;
    Array a(3, 0);
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
    // More tasks are defined to deal with gripper movement at points 2 and 4

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
        for (int j = 0; j < N_joints; j++) {
            joint_pos2[j] = Obj_list(i, j);
        }
        joint_pos1 = inverse_kinematics(forward_kinematics(joint_pos2)+ h, joint_pos2);

        int drop_box_id = Obj_list(i, N_joints);
        for (int j = 0; j < N_joints; j++) {
            joint_pos4 = Drop_box_list(drop_box_id, j);
        }
        joint_pos3 = inverse_kinematics(forward_kinematics(joint_pos4) + h, joint_pos4);

        // Joint velocity
        for (int j = 0; j < N_joints; j++) {
            joint_vel2[j] = 0;
            joint_vel4[j] = 0;
        }
        joint_vel1 = pinv(jacobian(joint_pos1))*v;
        joint_vel3 = pinv(jacobian(joint_pos3))*v;

        // todo: If we use joint acceleration, we need to change the formula for 1 and 3 since it is not right
        // // Joint acceleration
        // for (int j = 0; j < N_joints; j++) {
        //     joint_acc2[j] = 0;
        //     joint_acc4[j] = 0;
        // }
        // joint_acc1 = pinv(jacobian(joint_pos1))*a;
        // joint_acc3 = pinv(jacobian(joint_pos3))*a;

        // Filling manip list
        // manip_state.PVA = add_arrays(add_arrays(joint_pos1, -joint_vel1), -joint_acc1);
        manip_state.PVA = add_arrays(add_arrays(joint_pos1, -joint_vel1), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_blast;
        manip_list.push_back(manip_state);

        // manip_state.PVA = add_arrays(add_arrays(joint_pos2, joint_vel2), joint_acc2);
        manip_state.PVA = add_arrays(add_arrays(joint_pos2, joint_vel2), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);

        // manip_state.PVA = add_arrays(add_arrays(joint_pos2, joint_vel2), joint_acc2);
        manip_state.PVA = add_arrays(add_arrays(joint_pos2, joint_vel2), joint_acc);
        manip_state.gripper_state = 1;
        manip_state.mvt_gripper;
        manip_state.traj_gen_null;
        manip_list.push_back(manip_state);

        // manip_state.PVA = add_arrays(add_arrays(joint_pos1, joint_vel1), joint_acc1);
        manip_state.PVA = add_arrays(add_arrays(joint_pos1, joint_vel1), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);

        // manip_state.PVA = add_arrays(add_arrays(joint_pos3, -joint_vel3), -joint_acc3);
        manip_state.PVA = add_arrays(add_arrays(joint_pos3, -joint_vel3), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_blast;
        manip_list.push_back(manip_state);

        // manip_state.PVA = add_arrays(add_arrays(joint_pos4, joint_vel4), joint_acc4);
        manip_state.PVA = add_arrays(add_arrays(joint_pos4, joint_vel4), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);

        // manip_state.PVA = add_arrays(add_arrays(joint_pos4, joint_vel4), joint_acc4);
        manip_state.PVA = add_arrays(add_arrays(joint_pos4, joint_vel4), joint_acc);
        manip_state.gripper_state = 1;
        manip_state.mvt_gripper;
        manip_state.traj_gen_null;
        manip_list.push_back(manip_state);

        // manip_state.PVA = add_arrays(add_arrays(joint_pos3, joint_vel3), joint_acc3);
        manip_state.PVA = add_arrays(add_arrays(joint_pos3, joint_vel3), joint_acc);
        manip_state.gripper_state = 0;
        manip_state.mvt_robot;
        manip_state.traj_gen_straight;
        manip_list.push_back(manip_state);
    }
    return manip_list;
}