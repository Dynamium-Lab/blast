#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>

using namespace blast;

// Returns a representative stop-to-stop task for the UR5e.
static const blast::Array kStart = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
static const blast::Array kEnd   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};

inline blast::Task make_UR5e_task() {
  return blast::Task::stop_to_stop(kStart, kEnd);
}


using Catch::Matchers::WithinAbs;

TEST_CASE("initialize_optimization collision constraints", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task);

  SECTION("external collisions disabled zeroes the constraint count") {
    opt.constraints.external_collisions     = false;
    opt.constraints.n_collision_constraints = 12;
    initialize_optimization(&opt);
    CHECK(opt.constraints.n_collision_constraints == 0);
  }

  SECTION("external collisions enabled leaves the count untouched") {
    opt.constraints.external_collisions     = true;
    opt.constraints.n_collision_constraints = 12;
    initialize_optimization(&opt);
    CHECK(opt.constraints.n_collision_constraints == 12);
  }
}

TEST_CASE("initialize_optimization bounds", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task);
  initialize_optimization(&opt);

  const int n = opt.x_len();
  REQUIRE(opt.lower_bounds.size == n);
  REQUIRE(opt.upper_bounds.size == n);

  for (int i = 0; i < n - 1; i++) {
    CHECK(opt.lower_bounds[i] == -HUGE_VAL);
    CHECK(opt.upper_bounds[i] == HUGE_VAL);
  }
  CHECK_THAT(opt.lower_bounds.back(), WithinAbs(0.1, BLAST_EPSILON));
  CHECK_THAT(opt.upper_bounds.back(), WithinAbs(60.0, BLAST_EPSILON));
}

TEST_CASE("initialize_optimization task angle unwrapping", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task);
  initialize_optimization(&opt); // baseline: default UR5e task must already be unwrapped

  const auto reference = opt.task;
  const int  j         = 0; // joint under test; all other joints must stay untouched

  auto check_other_joints_untouched = [&]() {
    for (int i = 0; i < opt.manip.n_joints; i++) {
      if (i == j)
        continue;
      CHECK_THAT(opt.task(i, 0), WithinAbs(reference(i, 0), BLAST_EPSILON));
      CHECK_THAT(opt.task(i, 3), WithinAbs(reference(i, 3), BLAST_EPSILON));
    }
  };

  SECTION("default UR5e task has every joint within PI") {
    for (int i = 0; i < opt.manip.n_joints; i++) {
      CHECK(std::abs(opt.task(i, 0) - opt.task(i, 3)) <= PI + BLAST_EPSILON);
    }
  }

  SECTION("distance already <= PI is unchanged") {
    opt.manip.position_min[j] = -2 * PI;
    opt.manip.position_max[j] = 2 * PI;
    opt.task(j, 0)            = -1.0;
    opt.task(j, 3)            = 1.5;
    initialize_optimization(&opt);
    CHECK_THAT(opt.task(j, 0), WithinAbs(-1.0, BLAST_EPSILON));
    CHECK_THAT(opt.task(j, 3), WithinAbs(1.5, BLAST_EPSILON));
    check_other_joints_untouched();
  }

  SECTION("start is shifted when within limits") {
    opt.manip.position_min[j] = -2 * PI;
    opt.manip.position_max[j] = 2 * PI;
    opt.task(j, 0)            = -3.0;
    opt.task(j, 3)            = 3.0;
    initialize_optimization(&opt);
    CHECK_THAT(opt.task(j, 0), WithinAbs(-3.0 + 2 * PI, BLAST_EPSILON));
    CHECK_THAT(opt.task(j, 3), WithinAbs(3.0, BLAST_EPSILON));
    check_other_joints_untouched();
  }

  SECTION("end is shifted when start shift violates limits") {
    // start + 2PI = 3.283 > 3.2, end - 2PI = -3.283 > -3.5
    opt.manip.position_min[j] = -3.5;
    opt.manip.position_max[j] = 3.2;
    opt.task(j, 0)            = -3.0;
    opt.task(j, 3)            = 3.0;
    initialize_optimization(&opt);
    CHECK_THAT(opt.task(j, 0), WithinAbs(-3.0, BLAST_EPSILON));
    CHECK_THAT(opt.task(j, 3), WithinAbs(3.0 - 2 * PI, BLAST_EPSILON));
    check_other_joints_untouched();
  }

  SECTION("neither shift fits limits: unchanged, terminates") {
    opt.manip.position_min[j] = -3.1;
    opt.manip.position_max[j] = 3.1;
    opt.task(j, 0)            = -3.0;
    opt.task(j, 3)            = 3.0;
    initialize_optimization(&opt);
    CHECK_THAT(opt.task(j, 0), WithinAbs(-3.0, BLAST_EPSILON));
    CHECK_THAT(opt.task(j, 3), WithinAbs(3.0, BLAST_EPSILON));
    check_other_joints_untouched();
  }

  SECTION("multiple wraps until distance <= PI") {
    opt.manip.position_min[j] = -20.0;
    opt.manip.position_max[j] = 20.0;
    opt.task(j, 0)            = -10.0;
    opt.task(j, 3)            = 3.0;
    initialize_optimization(&opt);
    CHECK_THAT(opt.task(j, 0), WithinAbs(-10.0 + 4 * PI, BLAST_EPSILON));
    CHECK(std::abs(opt.task(j, 0) - opt.task(j, 3)) <= PI);
    check_other_joints_untouched();
  }
}

TEST_CASE("optimize stop-to-stop default pva + tool_speed constraints active", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task); // default enables pva + tool_speed constraints
  opt.success_tolerance = 0.01f;
  // Deterministic guess: Guess::random makes `success == true` a claim about
  // the solver, not about this test.
  opt.guess.type = Guess::straight_line;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0);
}

TEST_CASE("optimize stop-to-stop pva + tool_speed + torque constraints active", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task); // default enables pva + tool_speed constraints
  opt.constraints.torque = true;
  opt.success_tolerance  = 0.01f;
  // Deterministic guess: Guess::random makes `success == true` a claim about
  // the solver, not about this test.
  opt.guess.type = Guess::straight_line;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0);
}

TEST_CASE("optimization stop-to-stop pva + tool_speed + torque + self_collisions active", "[Optimization") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task); // default enables pva + tool_speed constraints
  opt.constraints.torque          = true;
  opt.constraints.self_collisions = true;
  opt.success_tolerance           = 0.01f;
  // Deterministic guess: Guess::random makes `success == true` a claim about
  // the solver, not about this test.
  opt.guess.type = Guess::straight_line;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0);
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

  // Deterministic guess: Guess::random makes `success == true` a claim about
  // the solver, not about this test.
  opt.guess.type = Guess::straight_line;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  CHECK(result.compute_time > 0.0f);
  CHECK(result.trajectory.t.size > 0);
}
