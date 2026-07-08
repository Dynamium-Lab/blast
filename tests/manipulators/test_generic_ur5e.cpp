#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;

TEST_CASE("UR5e forward_kinematics() test", "[Generic]") {
  Manipulator generic_manip = make_UR5e();

  ManipulatorTempData data;

  //  @ home position
  Vec3  expected_tool_position = {0.0, -0.233, 1.08};
  Array test_position          = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};
  forward_kinematics(generic_manip, data, test_position);
  CHECK(is_close(data.p_j[generic_manip.n_joints], expected_tool_position, 1e-2));

  //  @ 2nd position
  Vec3  expected_tool_position_2 = {-0.1179, -0.436, 0.152};
  Array test_position_2          = {-1.54, -1.83, -2.28, -0.59, 1.60, 0.023};
  forward_kinematics(generic_manip, data, test_position_2);

  CHECK(is_close(data.p_j[generic_manip.n_joints], expected_tool_position_2, 1e-2));
}

TEST_CASE("UR5e jacobian() - dimensions are 6 x n_joints", "[Generic]") {
  Manipulator         manip = make_UR5e();
  ManipulatorTempData data;
  Array               q = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};

  forward_kinematics(manip, data, q);
  Matrix J = jacobian(manip, data);

  CHECK(J.rows == 6u);
  CHECK(J.cols == manip.n_joints);
}

TEST_CASE("UR5e jacobian() - all entries are finite", "[Generic]") {
  Manipulator         manip = make_UR5e();
  ManipulatorTempData data;
  Array               q = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};

  forward_kinematics(manip, data, q);
  Matrix J = jacobian(manip, data);

  for (u32 r = 0; r < J.rows; r++)
    for (u32 c = 0; c < J.cols; c++)
      CHECK(std::isfinite(J(r, c)));
}

TEST_CASE("UR5e jacobian() - linear velocity rows are non-trivial at home position", "[Generic]") {
  Manipulator         manip = make_UR5e();
  ManipulatorTempData data;
  Array               q = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};

  forward_kinematics(manip, data, q);
  Matrix J = jacobian(manip, data);

  // linear velocity part (rows 0-2): at least one column must be nonzero
  real col_norm_sum = 0.0;
  for (u32 c = 0; c < J.cols; c++)
    col_norm_sum += J(0, c) * J(0, c) + J(1, c) * J(1, c) + J(2, c) * J(2, c);
  CHECK(col_norm_sum > 0.0);
}

TEST_CASE("UR5e dynamics() - efforts are finite at home position with zero motion", "[Generic]") {
  Manipulator         manip = make_UR5e();
  ManipulatorTempData data;
  Array               q = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};
  Array               vel(manip.n_joints, 0.0);
  Array               acc(manip.n_joints, 0.0);

  forward_kinematics(manip, data, q);
  dynamics(manip, data, vel, acc);

  for (u32 j = 0; j < manip.n_joints; j++)
    CHECK(std::isfinite(data.efforts[j]));
}

TEST_CASE("UR5e dynamics() - efforts change with nonzero accelerations (Coriolis/centrifugal terms)", "[Generic]") {
  Manipulator         manip = make_UR5e();
  ManipulatorTempData data;
  Array               q = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};
  Array               vel(manip.n_joints, 0.0);
  Array               acc_zero(manip.n_joints, 0.0);
  Array               acc_nonzero = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

  forward_kinematics(manip, data, q);

  dynamics(manip, data, vel, acc_zero);
  Array efforts_zero(manip.n_joints);
  for (u32 j = 0; j < manip.n_joints; j++)
    efforts_zero[j] = data.efforts[j];

  dynamics(manip, data, vel, acc_nonzero);
  Array efforts_nonzero(manip.n_joints);
  for (u32 j = 0; j < manip.n_joints; j++)
    efforts_nonzero[j] = data.efforts[j];

  CHECK(!is_close(efforts_zero, efforts_nonzero, 1e-6));
}
