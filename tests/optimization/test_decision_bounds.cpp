#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include "test_helper.hpp"

using namespace blast;

// decision_bounds: control points within the tightened position limits widened by
// control_point_bound_margin of the range, T within [min_duration, max_duration], in every
// optimization method; a start outside the box is clamped rather than refused (NLopt returns
// NLOPT_INVALID_ARGS for it).
TEST_CASE("decision bounds are set and respected by every method", "[Optimization]") {
  Array       start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array       end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.0};
  Manipulator robot = make_UR5e();
  for (int joint = 0; joint < robot.n_joints; joint++) {
    robot.position_min[joint] = -3.1416;
    robot.position_max[joint] = 3.1416;
  }

  for (real margin: {0.0, 0.1}) {
    Optimization opt(robot, Task::stop_to_stop(start, end));
    opt.max_duration               = 3.0;
    opt.control_point_bound_margin = margin;
    initialize_optimization_with_segments(&opt);
    auto      tolerance_snapshot    = tighten_for_success_tolerance(&opt);
    const u32 n_variables           = opt.bspline.x_len(opt.task);
    const int n_free_control_points = (int) opt.bspline.n_ctrl - 6;
    Array     lower_bounds(n_variables), upper_bounds(n_variables);
    decision_bounds(&opt, lower_bounds, upper_bounds);
    CHECK(lower_bounds.back() == Approx(opt.min_duration));
    CHECK(upper_bounds.back() == Approx(3.0));
    INFO("margin " << margin);
    for (int joint = 0; joint < opt.manip.n_joints; joint++) {
      const real range = opt.manip.position_max[joint] - opt.manip.position_min[joint]; // the tightened limits
      for (int control_point = 0; control_point < n_free_control_points; control_point++) {
        CHECK(lower_bounds[joint * n_free_control_points + control_point] == Approx(opt.manip.position_min[joint] - margin * range));
        CHECK(upper_bounds[joint * n_free_control_points + control_point] == Approx(opt.manip.position_max[joint] + margin * range));
      }
    }
    restore_from_tolerance(&opt, tolerance_snapshot);
  }

  for (auto method: {OptimizationMethod::with_segments, OptimizationMethod::baseline,
                     OptimizationMethod::with_analytical_pva, OptimizationMethod::with_analytical_dynamics,
                     OptimizationMethod::broadphase, OptimizationMethod::double_broadphase}) {
    Optimization opt(robot, Task::stop_to_stop(start, end));
    opt.method                   = method;
    opt.constraints.position     = true;
    opt.constraints.velocity     = true;
    opt.constraints.acceleration = true;
    opt.max_duration             = 3.0;
    opt.guess.type               = Guess::custom;
    Array initial_x              = blast::guess_straight_line(&opt);
    initial_x.back()             = 5.0;  // T above the cap
    initial_x[0]                 = 10.0; // a control point far outside
    opt.guess.initial_x          = initial_x;
    const Result result          = optimize(&opt);
    INFO("method " << (int) method);
    CHECK(result.nlopt_exit_criteria != NLOPT_INVALID_ARGS);                 // clamped, not refused
    REQUIRE(result.x.size == initial_x.size);
    CHECK(result.x.back() <= 3.0 + 1e-3);                                    // + the millisecond round-up of the final validation
    const real bound = 3.1416 + opt.control_point_bound_margin * 2 * 3.1416; // the true limits plus the margin
    for (u32 variable = 0; variable + 1 < result.x.size; variable++) {
      CHECK(result.x[variable] <= bound);
      CHECK(result.x[variable] >= -bound);
    }
  }
}
