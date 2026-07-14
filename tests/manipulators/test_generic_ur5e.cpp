#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;


static Manipulator make_1dof(real mass, Vec3 cog, Mat3 inertia) {
  ManipulatorLimits lim;
  lim.position_max = Array{PI};
  lim.position_min = Array{-PI};
  lim.velocity_max = Array{1.0f};

  ManipulatorKinematics kin;
  kin.joint_offsets[0]    = {0, 0, 0.5f};
  kin.joint_axes[0]       = {0, 0, 1};
  kin.static_rotations[0] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  ManipulatorDynamics dyn;
  dyn.link_masses[0]     = mass;
  dyn.cog_offsets[0]     = cog;
  dyn.inertia_tensors[0] = inertia;

  return Manipulator{1, lim, kin, &dyn};
}

Tool make_gripper() {
  Tool t;
  t.position       = {0.0, 0.0, 0.05};
  t.rotation       = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  t.mass           = 0.85;
  t.cog_offset     = {0.0, 0.0, 0.03};
  t.inertia_tensor = Mat3{
          0.001, 0.0, 0.0,
          0.0, 0.001, 0.0,
          0.0, 0.0, 0.0005};
  return t;
}

Payload make_payload() {
  Payload p;
  p.position       = {0.0, 0.0, 0.02};
  p.rotation       = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  p.mass           = 1.2;
  p.cog_offset     = {0.0, 0.0, 0.04};
  p.inertia_tensor = Mat3{
          0.0008, 0.0, 0.0,
          0.0, 0.0008, 0.0,
          0.0, 0.0, 0.0003};
  return p;
}

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

TEST_CASE("UR5e dynamics() - efforts differ with vs without gripper+payload", "[Gripper]") {
  Manipulator         manip_bare = make_UR5e();
  ManipulatorTempData data_bare;
  Array               q = {0.0, 0.0, -PI / 2, PI / 2, 0.0, 0.0};
  Array               vel(manip_bare.n_joints, 0.0);
  Array               acc(manip_bare.n_joints, 0.0);

  forward_kinematics(manip_bare, data_bare, q);
  dynamics(manip_bare, data_bare, vel, acc);
  Array efforts_bare(manip_bare.n_joints);
  for (u32 j = 0; j < manip_bare.n_joints; j++)
    efforts_bare[j] = data_bare.efforts[j];

  Manipulator manip_loaded = make_UR5e();
  manip_loaded.set_tool(make_gripper());
  manip_loaded.set_payload(make_payload());
  ManipulatorTempData data_loaded;

  forward_kinematics(manip_loaded, data_loaded, q);
  dynamics(manip_loaded, data_loaded, vel, acc);
  Array efforts_loaded(manip_loaded.n_joints);
  for (u32 j = 0; j < manip_loaded.n_joints; j++)
    efforts_loaded[j] = data_loaded.efforts[j];

  real max_loaded = std::max(std::abs(max(efforts_loaded)), std::abs(min(efforts_loaded)));
  real max_bare   = std::max(std::abs(max(efforts_bare)), std::abs(min(efforts_bare)));

  CHECK(max_loaded > max_bare);
}
