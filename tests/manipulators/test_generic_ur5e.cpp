#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;

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

  CHECK(!is_close(efforts_zero, efforts_nonzero, BLAST_EPSILON));
}

TEST_CASE("UR5e dynamics() - efforts differ with vs without tool+payload", "[Tool]") {
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

// ─── Tool + Payload tests ────────────────────────────────────────────────────

TEST_CASE("Zero-offset tool equals flange", "[Manipulator]") {
  Manipulator manip = make_UR5e();

  Tool tool;
  tool.position = {0, 0, 0};
  tool.rotation = eye();

  manip.set_tool(tool);

  ManipulatorTempData data;

  Array q = {0.2, -1.0, 0.4, -1.3, 0.2, 0.1};

  forward_kinematics(manip, data, q);

  CHECK(is_close(data.tool_position, data.p_j[manip.n_joints - 1], BLAST_EPSILON));
  CHECK(is_close(data.tool_rotation, data.rotations_mult[manip.n_joints - 1], BLAST_EPSILON));
}

TEST_CASE("Manipulator tool FK transform", "[Manipulator]") {
  Manipulator manip = make_UR5e();

  Tool tool;
  tool.position = {0.1, 0, 0};
  tool.rotation = rpy2rotation({0, 0, PI / 2});
  manip.set_tool(tool);

  ManipulatorTempData data;

  Array q = {0.2, -1.0, 0.5, -1.2, 0.3, 0.1};
  forward_kinematics(manip, data, q);

  Vec3 expected_position = data.p_j[manip.n_joints - 1] + data.rotations_mult[manip.n_joints - 1] * tool.position;
  Mat3 expected_rotation = data.rotations_mult[manip.n_joints - 1] * tool.rotation;

  CHECK(is_close(data.tool_position, expected_position, BLAST_EPSILON));
  CHECK(is_close(data.tool_rotation, expected_rotation, BLAST_EPSILON));
}

TEST_CASE("Manipulator payload FK transform", "[Manipulator]") {
  Manipulator manip = make_UR5e();

  Payload payload;
  payload.position = {0.2, 0, 0};
  manip.set_payload(payload);

  ManipulatorTempData data;

  Array q = {0.1, -1.2, 0.3, -1.1, 0.4, 0.2};
  forward_kinematics(manip, data, q);

  Vec3 expected_position = data.p_j[manip.n_joints - 1] + data.rotations_mult[manip.n_joints - 1] * payload.position;
  Mat3 expected_rotation = data.rotations_mult[manip.n_joints - 1] * payload.rotation;

  CHECK(is_close(data.payload_position, expected_position, BLAST_EPSILON));
  CHECK(is_close(data.payload_rotation, expected_rotation, BLAST_EPSILON));
}

TEST_CASE("Tool Jacobian finite difference", "[Manipulator]") {
  Manipulator manip = make_UR5e();

  Tool tool;
  tool.position             = {0.15, 0.02, 0.1};
  tool.tool_center_position = {0, 0, 0};
  manip.set_tool(tool);

  ManipulatorTempData data;

  Array q = {0.3, -1.0, 0.5, -1.5, 0.2, 0.1};
  forward_kinematics(manip, data, q);

  Matrix J = jacobian(manip, data);

  constexpr real eps = BLAST_EPSILON;

  for (int i = 0; i < manip.n_joints; i++) {
    Array qp = q;
    Array qm = q;

    qp[i] += eps;
    qm[i] -= eps;

    ManipulatorTempData dp;
    ManipulatorTempData dm;

    forward_kinematics(manip, dp, qp);
    forward_kinematics(manip, dm, qm);

    Vec3 expected_position = (dp.tool_position - dm.tool_position) / (2 * eps);
    Vec3 expected_rotation = data.rotations_mult[i] * manip.joint_axes[i];

    CHECK(is_close(expected_position, Vec3(J(0, i), J(1, i), J(2, i)), 1e-4));
    CHECK(is_close(expected_rotation, Vec3(J(3, i), J(4, i), J(5, i)), BLAST_EPSILON));
  }
}

TEST_CASE("Removing tool restores flange FK", "[Manipulator]") {
  Manipulator manip = make_UR5e();

  ManipulatorTempData data1;
  ManipulatorTempData data2;

  Array q = {0.2, -1.0, 0.3, -1.4, 0.2, 0.1};
  forward_kinematics(manip, data1, q);

  Tool tool;
  tool.position = {0.5, 0, 0};

  manip.set_tool(tool);
  manip.remove_tool();

  forward_kinematics(manip, data2, q);

  CHECK(is_close(data1.p_j[manip.n_joints], data2.p_j[manip.n_joints], BLAST_EPSILON));
}

TEST_CASE("Tool changes dynamics", "[Manipulator]") {
  Manipulator manip1 = make_UR5e();
  Manipulator manip2 = make_UR5e();

  Tool tool;

  tool.mass     = 5.0;
  tool.position = {0.2, 0, 0};
  manip2.set_tool(tool);

  ManipulatorTempData d1;
  ManipulatorTempData d2;

  Array q   = {0.1, -1.0, 0.5, -1.2, 0.3, 0.2};
  Array qd  = {0, 0, 0, 0, 0, 0};
  Array qdd = {0, 0, 0, 0, 0, 0};
  forward_kinematics(manip1, d1, q);
  forward_kinematics(manip2, d2, q);

  dynamics(manip1, d1, qd, qdd);
  dynamics(manip2, d2, qd, qdd);

  for (int i = 0; i < manip1.n_joints; i++) {
    CHECK(std::abs(d1.efforts[i]) < std::abs(d2.efforts[i]));
  }
}

TEST_CASE("Massless tool leaves dynamics unchanged", "[Manipulator]") {
  Manipulator manip1 = make_UR5e();
  Manipulator manip2 = make_UR5e();

  Tool tool;

  tool.mass           = 0;
  tool.position       = {0.3, 0, 0};
  tool.inertia_tensor = {};

  manip2.set_tool(tool);

  Array q   = {0.3, -1.0, 0.4, -1.1, 0.3, 0.1};
  Array qd  = {0.2, 0.1, -0.1, 0.2, 0.1, 0};
  Array qdd = {0.3, -0.2, 0.4, 0.1, -0.1, 0.2};

  ManipulatorTempData d1, d2;

  forward_kinematics(manip1, d1, q);
  forward_kinematics(manip2, d2, q);

  dynamics(manip1, d1, qd, qdd);
  dynamics(manip2, d2, qd, qdd);

  CHECK(is_close(d1.efforts, d2.efforts, BLAST_EPSILON));
}

TEST_CASE("Removing tool restores dynamics", "[Manipulator]") {

  Manipulator manip = make_UR5e();

  Array q   = {0.2, -1.0, 0.4, -1.2, 0.3, 0.1};
  Array qd  = {0.1, 0.2, -0.1, 0.3, 0.1, -0.2};
  Array qdd = {0.2, -0.3, 0.1, 0.4, -0.2, 0.1};

  ManipulatorTempData original;
  ManipulatorTempData restored;

  forward_kinematics(manip, original, q);
  dynamics(manip, original, qd, qdd);

  Tool tool;
  tool.mass     = 5;
  tool.position = {0.2, 0, 0};

  manip.set_tool(tool);
  manip.remove_tool();

  forward_kinematics(manip, restored, q);
  dynamics(manip, restored, qd, qdd);

  CHECK(is_close(original.efforts, restored.efforts, BLAST_EPSILON));
}

TEST_CASE("Payload changes dynamics", "[Manipulator]") {
  Manipulator manip1 = make_UR5e();
  Manipulator manip2 = make_UR5e();

  Payload payload;

  payload.mass     = 5.0;
  payload.position = {0.2, 0, 0};
  manip2.set_payload(payload);

  ManipulatorTempData d1;
  ManipulatorTempData d2;

  Array q   = {0.1, -1.0, 0.5, -1.2, 0.3, 0.2};
  Array qd  = {0, 0, 0, 0, 0, 0};
  Array qdd = {0, 0, 0, 0, 0, 0};
  forward_kinematics(manip1, d1, q);
  forward_kinematics(manip2, d2, q);

  dynamics(manip1, d1, qd, qdd);
  dynamics(manip2, d2, qd, qdd);

  for (int i = 0; i < manip1.n_joints; i++) {
    CHECK(std::abs(d1.efforts[i]) < std::abs(d2.efforts[i]));
  }
}

TEST_CASE("Removing payload restores dynamics", "[Manipulator]") {
  Manipulator manip = make_UR5e();

  Array q   = {0.2, -1.0, 0.4, -1.2, 0.3, 0.1};
  Array qd  = {0.1, 0.2, -0.1, 0.3, 0.1, -0.2};
  Array qdd = {0.2, -0.3, 0.1, 0.4, -0.2, 0.1};

  ManipulatorTempData original;
  ManipulatorTempData restored;

  forward_kinematics(manip, original, q);
  dynamics(manip, original, qd, qdd);

  Payload payload;
  payload.mass     = 5;
  payload.position = {0.2, 0, 0};

  manip.set_payload(payload);
  manip.remove_payload();

  forward_kinematics(manip, restored, q);
  dynamics(manip, restored, qd, qdd);

  CHECK(is_close(original.efforts, restored.efforts, BLAST_EPSILON));
}

TEST_CASE("Tool and payload together modify dynamics", "[Manipulator]") {
  Manipulator manip = make_UR5e();

  Tool tool;
  tool.mass     = 2.0;
  tool.position = {0.1, 0, 0};

  Payload payload;
  payload.mass     = 3.0;
  payload.position = {0.3, 0, 0};

  manip.set_tool(tool);
  manip.set_payload(payload);

  ManipulatorTempData data;

  Array q   = {0.3, -1.0, 0.4, -1.2, 0.2, 0.1};
  Array qd  = {0.2, 0.1, -0.2, 0.3, 0.1, -0.1};
  Array qdd = {0.3, -0.1, 0.2, 0.1, -0.2, 0.4};

  forward_kinematics(manip, data, q);

  CHECK_NOTHROW(dynamics(manip, data, qd, qdd));
}
