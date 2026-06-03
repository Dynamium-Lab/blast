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

TEST_CASE("set_payload() - combined mass is sum of link and payload masses", "[payload]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  m.set_payload(2.0f, {0, 0, 0.3f}, Mat3{});
  CHECK(is_close(m.link_masses[0], 12.0f, 1e-5f));
}

TEST_CASE("set_payload() - combined CoG is mass-weighted average", "[payload]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  m.set_payload(2.0f, {0, 0, 0.3f}, Mat3{});
  // av_new_z = (10*0.1 + 2*0.3) / 12 = 1.6/12
  CHECK(is_close(m.cog_offsets[0].z, 1.6f / 12.0f, 1e-5f));
  CHECK(is_close(m.cog_offsets[0].x, 0.0f, 1e-5f));
  CHECK(is_close(m.cog_offsets[0].y, 0.0f, 1e-5f));
}

TEST_CASE("set_payload() - inertia updated via parallel axis theorem", "[payload]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  m.set_payload(2.0f, {0, 0, 0.3f}, Mat3{});
  // delta_av_z = 1.6/12 - 0.1 = 1/30
  // av_to_mass_z = 0.3 - 1.6/12 = 1/6
  // I_zz = 10*(1/30)^2 + 2*(1/6)^2 = 1/90 + 1/18 = 1/15
  CHECK(is_close(m.inertia_tensors[0](2, 2), 1.0f / 15.0f, 1e-5f));
  CHECK(is_close(m.inertia_tensors[0](0, 0), 0.0f, 1e-5f));
  CHECK(is_close(m.inertia_tensors[0](1, 1), 0.0f, 1e-5f));
}

// ─── add_tool ────────────────────────────────────────────────────────────────

TEST_CASE("add_tool() - has_tool flag is set", "[tool]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  CHECK(m.has_tool == false);
  Tool tool;
  tool.mass           = 1.0f;
  tool.inertia_tensor = Mat3{};
  m.add_tool(tool);
  CHECK(m.has_tool == true);
}

TEST_CASE("add_tool() - combined mass is sum of link and tool masses", "[tool]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  Tool tool;
  tool.mass           = 2.0f;
  tool.cog_offset     = {0, 0, 0.3f};
  tool.inertia_tensor = Mat3{};
  m.add_tool(tool);
  CHECK(is_close(m.link_masses[0], 12.0f, 1e-5f));
}

TEST_CASE("add_tool() - combined CoG is mass-weighted average", "[tool]") {
  Manipulator m = make_1dof(10.0f, {0, 0, 0.1f}, Mat3{});
  Tool tool;
  tool.mass           = 2.0f;
  tool.cog_offset     = {0, 0, 0.3f};
  tool.inertia_tensor = Mat3{};
  m.add_tool(tool);
  CHECK(is_close(m.cog_offsets[0].z, 1.6f / 12.0f, 1e-5f));
  CHECK(is_close(m.cog_offsets[0].x, 0.0f, 1e-5f));
  CHECK(is_close(m.cog_offsets[0].y, 0.0f, 1e-5f));
}
