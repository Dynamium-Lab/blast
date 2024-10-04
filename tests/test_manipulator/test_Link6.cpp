#include "blast.h"
using namespace blast;

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch.hpp>

Link6 manip;

TEST_CASE("Link6", "Initializer") {
    CHECK(manip.pmax.size == 6);
    CHECK(manip.vmax.size == 6);
    CHECK(manip.amax.size == 6);
    CHECK(manip.tau_max.size == 6);

}

// todo:
// TEST_CASE("Link6", "Set Payload") {

// }

// todo:
// TEST_CASE("Link6", "Internal Constraints") {

// }

// todo:
// TEST_CASE("Link6", "Validate Task") {

// }

TEST_CASE("Link6", "Ncon") {
    int points = 100;
    CHECK(manip.ncon(points) == 2700);
}

// todo:
// TEST_CASE("Link6", "Dynamics") {
// todo: check if both dyn are the same at multiple random posistions
// }

TEST_CASE("Link6", "Forward Kinematics") {
    // Testing initialized at 0
    Array pos_array(6);
    Matrix pos_matrix(6, 1);

    Array tmp_array = {0, 0, 0, 0, 0, 0};
    Matrix tmp_matrix(tmp_array);

    CHECK(pos_array == tmp_array);
    CHECK(pos_matrix == tmp_matrix);

    // Testing for zeros
    auto f_kin_array = manip.forward_kinematics(pos_array);
    auto f_kin_matrix = manip.forward_kinematics(pos_matrix);

    CHECK(f_kin_array[0] == f_kin_matrix(0, 0));
    CHECK(f_kin_array[1] == f_kin_matrix(1, 0));
    CHECK(f_kin_array[2] == f_kin_matrix(2, 0));

    Array actual_pos(3);
    actual_pos[0] = 0.477;
    actual_pos[1] = 0.128;
    actual_pos[2] = 1.145;
    CHECK(abs(f_kin_array[0] - actual_pos[0]) <= 0.001);
    CHECK(abs(f_kin_array[1] - actual_pos[1]) <= 0.001);
    CHECK(abs(f_kin_array[2] - actual_pos[2]) <= 0.001);

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

    actual_pos[0] = 0.655;
    actual_pos[1] = 0.117;
    actual_pos[2] = 0.037;
    CHECK(abs(f_kin_array[0] - actual_pos[0]) <= 0.001);
    CHECK(abs(f_kin_array[1] - actual_pos[1]) <= 0.001);
    CHECK(abs(f_kin_array[2] - actual_pos[2]) <= 0.001);

    // todo: check if both f_kin are the same at multiple random positions
    for (u32 i = 0; i < 100; i++) {
        auto rnd_array = random_array(6, PI);
        Matrix rnd_matrix(rnd_array);

        auto f_kin_array_rnd = manip.forward_kinematics(rnd_array);
        auto f_kin_matrix_rnd = manip.forward_kinematics(rnd_matrix);

        CHECK(f_kin_array_rnd[0] == f_kin_matrix_rnd(0, 0));
        CHECK(f_kin_array_rnd[1] == f_kin_matrix_rnd(1, 0));
        CHECK(f_kin_array_rnd[2] == f_kin_matrix_rnd(2, 0));
    }
}
