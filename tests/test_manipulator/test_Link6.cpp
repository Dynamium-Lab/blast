#include <blast>
using namespace blast;

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch.hpp>

Link6 manip;
auto joints = manip.joints;

TEST_CASE("Link6 Initializer") {
    CHECK(manip.pmax.size == joints);
    CHECK(manip.vmax.size == joints);
    CHECK(manip.amax.size == joints);
    CHECK(manip.tau_max.size == joints);

}

// todo:
// TEST_CASE("Link6 Set Payload") {
// }

TEST_CASE("Link6 Internal Constraints") {
    Trajectory traj(1, joints);
    auto ncon = manip.ncon(1);
    Array constraints(ncon);

    manip.dynamics(traj);
    Array torque_const_zero(joints);
    for (u32 i = 0; i < joints; i++)
        torque_const_zero[i] = (abs(manip._efforts(i, 0)) - manip.tau_max[i]) / manip.tau_max[i];

    manip.internal_constraints(traj, constraints.data);

    CHECK(constraints.size == ncon); // todo: is there a way to actually test this ?
    for (u32 i = 9; i < ncon - joints; i++)
        CHECK(constraints[i] == -1);
    for ( u32 i = ncon - joints; i < ncon; i++)
        CHECK(constraints[i] == torque_const_zero[i + joints - ncon]);

    // todo: check internal_collisions
}

// todo:
// TEST_CASE("Link6 Validate Task") {
// }

TEST_CASE("Link6", "Dynamics") {
// todo: check dyn at a specific known point
    Trajectory traj(1, joints);
    Array pos(joints);
    Array vel(joints);
    Array acc(joints);

    manip.dynamics(traj);
    auto dyn_array = manip.dynamics(pos, vel, acc);
    for (u32 i = 0; i < joints; i++)
        CHECK(manip._efforts(i, 0) == dyn_array[i]);

    for (u32 i = 0; i < 100; i++) {
        Trajectory traj_rnd(1, joints);
        auto pos_rnd = random_array(joints, PI);
        auto vel_rnd = random_array(joints, PI);
        auto acc_rnd = random_array(joints, PI);
        traj_rnd.pos = pos_rnd;
        traj_rnd.vel = vel_rnd;
        traj_rnd.acc = acc_rnd;

        manip.dynamics(traj_rnd);
        auto dyn_array = manip.dynamics(pos_rnd, vel_rnd, acc_rnd);
        for (u32 i = 0; i < joints; i++)
            CHECK(manip._efforts(i, 0) == dyn_array[i]);
    }
}

TEST_CASE("Link6 Forward Kinematics") {
    // Testing initialized at 0
    Array pos_array(joints);
    Matrix pos_matrix(joints, 1);

    Array tmp_array = {0, 0, 0, 0, 0, 0};
    Matrix tmp_matrix(tmp_array);

    CHECK(is_close(pos_array, tmp_array));
    CHECK(is_close(pos_matrix, tmp_matrix));

    // Testing for zeros
    auto f_kin_array = manip.forward_kinematics(pos_array);
    auto f_kin_matrix = manip.forward_kinematics(pos_matrix);


    CHECK(f_kin_array[0] == f_kin_matrix(0, 0));
    CHECK(f_kin_array[1] == f_kin_matrix(1, 0));
    CHECK(f_kin_array[2] == f_kin_matrix(2, 0));

    f_kin_array.size = 3; // note: don't try this at home.

    Array actual_pos(3);
    actual_pos[0] = 0.489;
    actual_pos[1] = 0.117;
    actual_pos[2] = 1.136;
    CHECK(is_close(f_kin_array, actual_pos, 0.001));

    // Testing for home position
    pos_array = deg2rad({0.0, -20.0, 110.0, 0.0, 40.0, 0.0});

    pos_matrix(0, 0) = deg2rad(0.0);
    pos_matrix(1, 0) = deg2rad(-20.0);
    pos_matrix(2, 0) = deg2rad(110.0);
    pos_matrix(3, 0) = deg2rad(0.0);
    pos_matrix(4, 0) = deg2rad(40.0);
    pos_matrix(5, 0) = deg2rad(0.0);

    f_kin_array = manip.forward_kinematics(pos_array);
    f_kin_matrix = manip.forward_kinematics(pos_matrix);

    CHECK(f_kin_array[0] == f_kin_matrix(0, 0));
    CHECK(f_kin_array[1] == f_kin_matrix(1, 0));
    CHECK(f_kin_array[2] == f_kin_matrix(2, 0));

    f_kin_array.size = 3; // note: don't try this at home.

    actual_pos[0] = 0.649;
    actual_pos[1] = 0.117;
    actual_pos[2] = 0.026; // todo: Verify?
    CHECK(is_close(f_kin_array, actual_pos, 0.001));

    // todo: check if both f_kin are the same at multiple random positions
    for (u32 i = 0; i < 100; i++) {
        auto rnd_array = random_array(joints, PI);
        Matrix rnd_matrix(rnd_array);

        auto f_kin_array_rnd = manip.forward_kinematics(rnd_array);
        auto f_kin_matrix_rnd = manip.forward_kinematics(rnd_matrix);

        CHECK(f_kin_array_rnd[0] == f_kin_matrix_rnd(0, 0));
        CHECK(f_kin_array_rnd[1] == f_kin_matrix_rnd(1, 0));
        CHECK(f_kin_array_rnd[2] == f_kin_matrix_rnd(2, 0));
    }
}

// todo:
// TEST_CASE("Link6 Jacobian"){
// }

TEST_CASE("Link6 Capsules") {
    Array pos_array(joints);
    Matrix pos_matrix(joints, 1);

    Array tmp_array = {0, 0, 0, 0, 0, 0};
    Matrix tmp_matrix(tmp_array);

    CHECK(is_close(pos_array, tmp_array));
    CHECK(is_close(pos_matrix, tmp_matrix));

    auto robot_caps_matrix = manip.robot_capsules(pos_array);
    auto robot_caps_vector = manip.robot_capsules(pos_matrix, 1);

    CHECK(robot_caps_matrix.size == (u32)robot_caps_vector.size()*(sizeof(Capsule)/sizeof(real) - 2)); // Vec3 has 4 real, with one being _pad that is not desired in this calculation

    for (u32 i = 0; i < 100; i++) {
        auto rnd_array = random_array(joints, PI);
        Matrix rnd_matrix(rnd_array);
        auto robot_caps_matrix = manip.robot_capsules(rnd_array);
        auto robot_caps_vector = manip.robot_capsules(rnd_matrix, 1);
        for (u32 j = 0; j < robot_caps_matrix.cols; j++) {
            CHECK(robot_caps_matrix(0, j) == robot_caps_vector[j].p1.x);
            CHECK(robot_caps_matrix(1, j) == robot_caps_vector[j].p1.y);
            CHECK(robot_caps_matrix(2, j) == robot_caps_vector[j].p1.z);
            CHECK(robot_caps_matrix(3, j) == robot_caps_vector[j].p2.x);
            CHECK(robot_caps_matrix(4, j) == robot_caps_vector[j].p2.y);
            CHECK(robot_caps_matrix(5, j) == robot_caps_vector[j].p2.z);
            CHECK(robot_caps_matrix(6, j) == robot_caps_vector[j].r);
        }
    }
}