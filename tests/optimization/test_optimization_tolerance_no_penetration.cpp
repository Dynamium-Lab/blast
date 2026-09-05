#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include "test_helper.hpp"

using namespace blast;

// Returns a representative stop-to-stop task for the UR5e. Endpoints are file-scope
// so the tests can build a deterministic guess from them.
static const blast::Array kStart = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
static const blast::Array kEnd   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};

inline blast::Task make_UR5e_task() {
  return blast::Task::stop_to_stop(kStart, kEnd);
}

// Regression test: a large success_tolerance is meant to let nlopt/native-SQP early-out
// faster, but must never let the returned trajectory violate the real (untightened) robot
// or collision limits. Before the tightening fix, success_tolerance leaked directly into
// nlopt's own feasibility tolerance AND into the post-solve success gate, so a large value
// could let the trajectory penetrate an obstacle by up to success_tolerance world-units.
TEST_CASE("large success_tolerance does not allow real constraint violation", "[Optimization]") {
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

  opt.success_tolerance = 0.05f;

  opt.guess.type = Guess::custom;                                            // deterministic; see straight_line_guess()

  opt.guess.initial_x = blast::test::straight_line_guess(opt, kStart, kEnd); // deliberately large - well above the 0.01 default
  opt.max_tries       = 10;                                                  // random initial guess; a tightened obstacle needs a few tries

  Result result = optimize(&opt);

  CHECK(result.success == true);
  REQUIRE(result.trajectory.t.size > 0);

  // Independently re-check the solved trajectory against the ORIGINAL, untightened robot/world
  // geometry (the local `robot`/`world`, never touched by optimize()) rather than trusting
  // opt's post-solve state.
  Optimization opt_check(robot, task);
  opt_check.world       = world;
  opt_check.bspline     = opt.bspline;
  opt_check.constraints = opt.constraints;
  n_con(&opt_check);

  Array constraints_points(opt_check.constraints.n_constraints);
  compute_constraints(constraints_points.data, result.x, &opt_check);

  CHECK(max(constraints_points) < 1e-6f);
}

// Regression test for the retry-loop overwrite bug: is_valid/is_valid_more/result.trajectory
// used to be recomputed every try with the `break` on first success commented out, so
// result.num_tries always equalled max_tries regardless of when a valid trajectory was found.
TEST_CASE("optimize stops at the first valid try instead of always using the last", "[Optimization]") {
  Manipulator robot = make_UR5e();
  Task        task  = make_UR5e_task();

  Optimization opt(robot, task);         // default enables pva + tool_speed constraints, no obstacles
  opt.success_tolerance = 0.01f;
  opt.guess.type        = Guess::custom; // deterministic; see straight_line_guess()
  opt.guess.initial_x   = blast::test::straight_line_guess(opt, kStart, kEnd);
  opt.max_tries         = 5;

  Result result = optimize(&opt);

  CHECK(result.success == true);
  // The regression is "num_tries always equalled max_tries". Assert THAT, not
  // "the first try always succeeds": since the accept gate became a single
  // success_tolerance (it was success_tolerance*2), a marginal first attempt can
  // legitimately fail and retry. Measured on this task, first-try failure is ~1%
  // on Linux and more on Windows, so `== 1` is a flaky assertion about the solver,
  // not about the loop.
  CHECK(result.num_tries >= 1);
  CHECK(result.num_tries < opt.max_tries);
}
