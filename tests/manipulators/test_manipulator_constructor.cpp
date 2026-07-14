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

// ─── constructor ────────────────────────────────────────────────────────────

TEST_CASE("Manipulator constructor - kinematics-only leaves dynamic fields zero-initialized", "[constructor]") {
  ManipulatorLimits lim;
  lim.position_max = Array{PI};
  lim.position_min = Array{-PI};
  lim.velocity_max = Array{1.0f};

  ManipulatorKinematics kin;
  kin.joint_offsets[0]    = {0, 0, 0.5f};
  kin.joint_axes[0]       = {0, 0, 1};
  kin.static_rotations[0] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  Manipulator m{1, lim, kin}; // no dynamics, no capsules

  CHECK(m.n_joints == 1u);
  CHECK(m._n_caps == 0);
  CHECK(m.link_masses[0] == 0.0f);
  CHECK(m._n_internal_collisions == 0);
}

TEST_CASE("Manipulator constructor - n_joints matches joint_count for standard manipulators", "[constructor]") {
  CHECK(make_Kinova_Gen3().n_joints == 7u);
  CHECK(make_UR5e().n_joints == 6u);
}

// ─── set_capsules ────────────────────────────────────────────────────────────

TEST_CASE("UR5e set_capsules() - capsule count is correct", "[constructor]") {
  Manipulator ur5e = make_UR5e();
  CHECK(ur5e._n_caps == 7);
}

TEST_CASE("UR5e set_capsules() - _n_internal_collisions counts lower-triangle pairs and base contacts", "[constructor]") {
  Manipulator ur5e = make_UR5e();
  // 8 lower-triangle collision-matrix pairs + 4 capsules that collide with the base sphere = 12
  CHECK(ur5e._n_internal_collisions == 12);
}

// ─── set_payload ─────────────────────────────────────────────────────────────

TEST_CASE("set_payload() - combined mass is sum of link and payload masses", "[constructor]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  m.has_tool    = true; // to stop warning in test
  Payload payload;
  payload.mass       = 2.0;
  payload.cog_offset = {0, 0, 0.3f};
  m.set_payload(payload);
  CHECK(is_close(m.link_masses[0], 12.0f, 1e-5f));
}

TEST_CASE("set_payload() - combined CoG is mass-weighted average", "[constructor]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  m.has_tool    = true; // to stop warning in test
  Payload payload;
  payload.mass       = 2.0;
  payload.cog_offset = {0, 0, 0.3f};
  m.set_payload(payload);
  CHECK(is_close(m.cog_offsets[0].z, 1.6f / 12.0f, 1e-5f));
  CHECK(is_close(m.cog_offsets[0].x, 0.0f, 1e-5f));
  CHECK(is_close(m.cog_offsets[0].y, 0.0f, 1e-5f));
}

TEST_CASE("set_payload() - rotated point mass is COM invariant", "[constructor]") {
  Manipulator m = make_1dof(0.0f, {0, 0, 0}, Mat3{});
  m.has_tool    = true; // to stop warning in test

  Payload p;
  p.mass       = 1.0;
  p.cog_offset = {1.0, 0.0, 0.0};

  // rotation is irrelevant for point mass inertia about itself
  p.rotation = {
          0, -1, 0,
          1, 0, 0,
          0, 0, 1};

  m.set_payload(p);

  Mat3 I = m.inertia_tensors[0];

  CHECK(is_close(I(0, 0), 0.0, 1e-9));
  CHECK(is_close(I(1, 1), 0.0, 1e-9));
  CHECK(is_close(I(2, 2), 0.0, 1e-9));
}

TEST_CASE("set_payload() - inertia updated via parallel axis theorem (diagonal terms)", "[constructor]") {
  Manipulator m = make_1dof(1.0f, {0, 0, 0}, Mat3{});
  m.has_tool    = true; // to stop warning in test

  Payload p;
  p.mass           = 2.0;
  p.cog_offset     = {0.3, 0.0, 0.0};
  p.inertia_tensor = Mat3{}; // point mass

  m.set_payload(p);

  // Ixx = m (y^2 + z^2)
  real expected_Iyy = 2.0 * 0.1 * 0.1 + 1.0 * 0.2 * 0.2;

  CHECK(is_close(m.inertia_tensors[0](0, 0), 0.0, 1e-9));
  CHECK(is_close(m.inertia_tensors[0](1, 1), expected_Iyy, 1e-9));
  CHECK(is_close(m.inertia_tensors[0](2, 2), expected_Iyy, 1e-9));
}

TEST_CASE("set_payload() - inertia updated via parallel axis theorem (off-diagonal terms)", "[constructor]") {
  Manipulator m = make_1dof(0.0f, {0, 0, 0}, Mat3{});
  m.has_tool    = true; // to stop warning in test

  // First mass (implicit link mass or dummy)
  m.link_masses[0] = 1.0;
  m.cog_offsets[0] = {0.0, 0.0, 0.0};

  Payload p;
  p.mass           = 2.0;
  p.cog_offset     = {1.0, 2.0, 3.0};
  p.inertia_tensor = Mat3{};

  m.set_payload(p);

  Mat3 I = m.inertia_tensors[0];

  real total_mass = 3.0;
  Vec3 com        = (1.0 * Vec3{0, 0, 0} + 2.0 * Vec3{1, 2, 3}) / total_mass;

  Vec3 d_link    = Vec3{0, 0, 0} - com;
  Vec3 d_payload = Vec3{1, 2, 3} - com;

  real expected_Ixy =
          1.0 * (-d_link.x * d_link.y) +
          2.0 * (-d_payload.x * d_payload.y);

  real expected_Ixz =
          1.0 * (-d_link.x * d_link.z) +
          2.0 * (-d_payload.x * d_payload.z);

  real expected_Iyz =
          1.0 * (-d_link.y * d_link.z) +
          2.0 * (-d_payload.y * d_payload.z);

  CHECK(is_close(I(0, 1), expected_Ixy, 1e-9));
  CHECK(is_close(I(0, 2), expected_Ixz, 1e-9));
  CHECK(is_close(I(1, 2), expected_Iyz, 1e-9));
}

TEST_CASE("set_payload() - mass and cog updated with pre-existing tool", "[constructor]") {
  Manipulator manip   = make_UR5e();
  int         link_id = manip.n_joints - 1;

  manip.set_tool(make_gripper());
  real mass_with_tool = manip.link_masses[link_id];

  Payload payload = make_payload();
  manip.set_payload(payload);

  CHECK(manip.has_payload);
  CHECK(manip.link_masses[link_id] == Approx(mass_with_tool + payload.mass));
}

TEST_CASE("set_payload() - inertia tensor stays symmetric and finite", "[constructor]") {
  Manipulator manip = make_UR5e();
  manip.set_tool(make_gripper());
  manip.set_payload(make_payload());

  int  link_id = manip.n_joints - 1;
  Mat3 I       = manip.inertia_tensors[link_id];

  for (u32 r = 0; r < 3; r++)
    for (u32 c = 0; c < 3; c++)
      CHECK(std::isfinite(I(r, c)));

  CHECK(I(0, 1) == Approx(I(1, 0)));
  CHECK(I(0, 2) == Approx(I(2, 0)));
  CHECK(I(1, 2) == Approx(I(2, 1)));
}

// ─── remove_payload ──────────────────────────────────────────────────────────

TEST_CASE("remove_payload() - restores tool-only mass, cog, inertia", "[constructor]") {
  Manipulator manip   = make_UR5e();
  int         link_id = manip.n_joints - 1;

  manip.set_tool(make_gripper());

  real mass_tool_only    = manip.link_masses[link_id];
  Vec3 cog_tool_only     = manip.cog_offsets[link_id];
  Mat3 inertia_tool_only = manip.inertia_tensors[link_id];

  manip.set_payload(make_payload());
  manip.remove_payload();

  CHECK(!manip.has_payload);
  CHECK(manip.link_masses[link_id] == Approx(mass_tool_only));
  CHECK(is_close(manip.cog_offsets[link_id], cog_tool_only, 1e-9));

  for (u32 r = 0; r < 3; r++)
    for (u32 c = 0; c < 3; c++)
      CHECK(manip.inertia_tensors[link_id](r, c) == Approx(inertia_tool_only(r, c)).margin(1e-9));
}

// ─── set_tool ────────────────────────────────────────────────────────────────

TEST_CASE("set_tool() - has_tool flag is set", "[tool]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  CHECK(m.has_tool == false);
  Tool tool;
  tool.mass           = 1.0f;
  tool.inertia_tensor = Mat3{};
  m.set_tool(tool);
  CHECK(m.has_tool == true);
}

TEST_CASE("set_tool() - mass and cog updated correctly", "[constructor]") {
  Manipulator manip       = make_UR5e();
  int         link_id     = manip.n_joints - 1;
  real        mass_before = manip.link_masses[link_id];
  Vec3        cog_before  = manip.cog_offsets[link_id];

  Tool gripper = make_gripper();
  manip.set_tool(gripper);

  CHECK(manip.has_tool);
  CHECK(manip.link_masses[link_id] == Approx(mass_before + gripper.mass));

  Vec3 expected_cog = (mass_before * cog_before + gripper.mass * (gripper.position + gripper.rotation * gripper.cog_offset)) / manip.link_masses[link_id];
  CHECK(is_close(manip.cog_offsets[link_id], expected_cog, 1e-9));
}

TEST_CASE("set_tool() - inertia tensor stays symmetric and finite", "[constructor]") {
  Manipulator manip = make_UR5e();
  manip.set_tool(make_gripper());

  int  link_id = manip.n_joints - 1;
  Mat3 I       = manip.inertia_tensors[link_id];

  for (u32 r = 0; r < 3; r++)
    for (u32 c = 0; c < 3; c++)
      CHECK(std::isfinite(I(r, c)));

  CHECK(I(0, 1) == Approx(I(1, 0)));
  CHECK(I(0, 2) == Approx(I(2, 0)));
  CHECK(I(1, 2) == Approx(I(2, 1)));
}

TEST_CASE("set_tool() - rotation is applied correctly", "[constructor]") {
  Manipulator m = make_1dof(0.0f, {0, 0, 0}, Mat3{});

  Tool t;
  // rotate tool 90 degrees around x-axis
  t.position       = {0, 0, 0};
  t.rotation       = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  t.inertia_tensor = {1e-4, 0, 0, 0, 1e-3, 0, 0, 0, 5e-3};
  t.cog_offset     = {0, 0, 0};
  t.mass           = 0.85;

  m.set_tool(t);

  Mat3 I_world = m.inertia_tensors[m.n_joints - 1];

  // Expected: I_world = R I R^T
  Mat3 R = t.rotation;
  // tool inertia rotated
  Mat3 I_expected = R * t.inertia_tensor * transpose(R);

  CHECK(is_close(I_world, I_expected));
}

// ─── remove_tool ─────────────────────────────────────────────────────────────

TEST_CASE("remove_tool() - restores original mass, cog, inertia (without payload)", "[constructor]") {
  Manipulator manip   = make_UR5e();
  int         link_id = manip.n_joints - 1;

  real mass_before    = manip.link_masses[link_id];
  Vec3 cog_before     = manip.cog_offsets[link_id];
  Mat3 inertia_before = manip.inertia_tensors[link_id];

  manip.set_tool(make_gripper());
  Mat3 inertia_during = manip.inertia_tensors[link_id];
  manip.remove_tool();
  Mat3 inertia_after = manip.inertia_tensors[link_id];

  CHECK(!manip.has_tool);
  CHECK(manip.link_masses[link_id] == Approx(mass_before));
  CHECK(is_close(manip.cog_offsets[link_id], cog_before, 1e-9));

  for (u32 r = 0; r < 3; r++)
    for (u32 c = 0; c < 3; c++)
      CHECK(manip.inertia_tensors[link_id](r, c) == Approx(inertia_before(r, c)).margin(1e-9));
}

TEST_CASE("remove_tool() - restores original mass, cog, inertia (with payload)", "[constructor]") {
  Manipulator manip   = make_UR5e();
  int         link_id = manip.n_joints - 1;

  real mass_before    = manip.link_masses[link_id];
  Vec3 cog_before     = manip.cog_offsets[link_id];
  Mat3 inertia_before = manip.inertia_tensors[link_id];

  manip.set_tool(make_gripper());
  manip.set_payload(make_payload());
  manip.remove_tool();

  CHECK(!manip.has_tool);
  CHECK(!manip.has_payload);
  CHECK(manip.link_masses[link_id] == Approx(mass_before));
  CHECK(is_close(manip.cog_offsets[link_id], cog_before, 1e-9));

  for (u32 r = 0; r < 3; r++)
    for (u32 c = 0; c < 3; c++)
      CHECK(manip.inertia_tensors[link_id](r, c) == Approx(inertia_before(r, c)).margin(1e-9));
}

// ─── Tool + Payload Dynamics ────────────────────────────────────────────────────────────────

TEST_CASE("set_payload() in tool and set_tool() / set_tool() and set_payload() order invariance", "[constructor]") {
  Manipulator base1 = make_UR5e();
  Manipulator base2 = make_UR5e();

  Tool    t = make_gripper();
  Payload p = make_payload();

  // pipeline A
  base1.set_tool(t);
  base1.set_payload(p);

  // pipeline B
  base2.set_tool(t);

  Mat3 I1 = base1.inertia_tensors[base1.n_joints - 1];
  Mat3 I2 = base2.inertia_tensors[base2.n_joints - 1];

  CHECK(is_close(base1.link_masses[base1.n_joints - 1], base2.link_masses[base1.n_joints - 1]));
  CHECK(is_close(base1.cog_offsets[base1.n_joints - 1], base2.cog_offsets[base1.n_joints - 1]));
  CHECK(is_close(base1.inertia_tensors[base1.n_joints - 1], base2.inertia_tensors[base1.n_joints - 1]));
}
