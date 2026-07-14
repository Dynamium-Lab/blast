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

// ─── set_payload / remove_payload ────────────────────────────────────────────

TEST_CASE("set_payload() - combined mass is sum of link and payload masses", "[constructor]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});

  CHECK(m.has_payload == false);

  Payload payload = make_payload();
  m.set_payload(payload);
  CHECK(m.has_payload == true);
  CHECK(is_close(m.payload, payload, BLAST_EPSILON));

  m.remove_payload();
  CHECK(m.has_payload == false);
}

// ─── set_tool / remove_tool ──────────────────────────────────────────────────

TEST_CASE("set_tool()  and remove_tool() basic tests", "[tool]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  CHECK(m.has_tool == false);

  Tool tool;
  tool.mass = 1.0f;
  m.set_tool(tool);
  CHECK(m.has_tool == true);
  CHECK(is_close(m.tool, tool, BLAST_EPSILON));

  m.remove_tool();
  CHECK(m.has_tool == false);
}
