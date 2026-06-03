#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>

using namespace blast;

// Returns a representative stop-to-stop task for the UR5e.
inline blast::Task make_UR5e_task() {
  blast::Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  blast::Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  return blast::Task::stop_to_stop(start, end);
}

TEST_CASE("optimize stop-to-stop default pva + tool_speed constraints active", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task); // default enables pva + tool_speed constraints
  opt.success_tolerance = 0.01f;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0u);
}

TEST_CASE("optimize stop-to-stop pva + tool_speed + torque constraints active", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task); // default enables pva + tool_speed constraints
  opt.constraints.torque = true;
  opt.success_tolerance  = 0.01f;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0u);
}

TEST_CASE("optimization stop-to-stop pva + tool_speed + torque + self_collisions active", "[Optimization") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task); // default enables pva + tool_speed constraints
  opt.constraints.torque          = true;
  opt.constraints.self_collisions = true;
  opt.success_tolerance           = 0.01f;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0u);
}

TEST_CASE("optimization stop-to-stop all constraints active", "[Optimization") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  World world;
  world.add_box(
          Vec3{0.4, 0.0, 0.6},            // centre: 40 cm in front, 60 cm high
          Vec3{0.05, 0.3, 0.3},           // half-extents: thin vertical slab
          Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1} // upright, axis-aligned
  );

  Optimization opt(robot, task); // default enables pva + tool_speed constraints
  opt.world = world;

  opt.constraints.torque              = true;
  opt.constraints.self_collisions     = true; // avoid self-contact
  opt.constraints.external_collisions = true; // avoid world obstacles

  opt.success_tolerance = 0.01f;


  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0u);
}
