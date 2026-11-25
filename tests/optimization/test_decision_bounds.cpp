#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include "test_helper.hpp"

using namespace blast;

// decision_bounds: control points within the solver's position bounds (the tightened limits, widened
// past a pinned start/goal) widened by control_point_bound_margin of the tightened range, T within [min_duration, max_duration], in every
// optimization method; a start outside the box is clamped rather than refused (NLopt returns
// NLOPT_INVALID_ARGS for it).
TEST_CASE("decision bounds are set and respected by every method", "[Optimization]") {
  Array       start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array       end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3};
  Manipulator robot = make_UR5e();
  for (int joint = 0; joint < robot.n_joints; joint++) {
    robot.position_min[joint] = -(real) 3.1416;
    robot.position_max[joint] = (real) 3.1416;
  }

  for (real margin: {(real) 0, (real) 0.1}) {
    Optimization opt(robot, Task::stop_to_stop(start, end));
    opt.max_duration               = 3;
    opt.control_point_bound_margin = margin;
    initialize_optimization_with_segments(&opt);
    auto      tolerance_snapshot    = tighten_for_success_tolerance(&opt);
    const u32 n_variables           = opt.bspline.x_len(opt.task);
    const int n_free_control_points = (int) opt.bspline.n_ctrl - 6; // stop-to-stop pins 3 at each end
    Array     lower_bounds(n_variables), upper_bounds(n_variables);
    decision_bounds(&opt, lower_bounds, upper_bounds);
    CHECK(lower_bounds.back() == Approx(opt.min_duration));
    CHECK(upper_bounds.back() == Approx(3.0));
    INFO("margin " << margin);
    // Expected box computed from the test data alone, independently of the code under test: the tightened
    // limit, or halfway between a pinned endpoint outside it and the true limit; plus margin times the
    // tightened range. The endpoints here are well clear of the tightened bound, so no ulp margin is needed.
    const real true_limit = (real) 3.1416;
    const real tightened  = true_limit / (1 + opt.success_tolerance);
    int        n_widened  = 0;
    for (int joint = 0; joint < opt.manip.n_joints; joint++) {
      INFO("joint " << joint);
      const real q_hi         = std::max(start[joint], end[joint]);
      const real q_lo         = std::min(start[joint], end[joint]);
      const real expected_max = q_hi > tightened ? q_hi + (true_limit - q_hi) / 2 : tightened;
      const real expected_min = q_lo < -tightened ? q_lo - (q_lo + true_limit) / 2 : -tightened;
      n_widened += (q_hi > tightened) + (q_lo < -tightened);
      const real range = 2 * tightened;
      for (int control_point = 0; control_point < n_free_control_points; control_point++) {
        CHECK(lower_bounds[joint * n_free_control_points + control_point] == Approx(expected_min - margin * range));
        CHECK(upper_bounds[joint * n_free_control_points + control_point] == Approx(expected_max + margin * range));
      }
    }
    CHECK(n_widened == 1); // joint 5's start (-3.12883) is outside the tightened band: the widening is exercised
    restore_from_tolerance(&opt, tolerance_snapshot);
  }

  // The bounds must hold whether or not the solve succeeds, so success is not checked.
  for (auto method: {OptimizationMethod::with_segments, OptimizationMethod::baseline,
                     OptimizationMethod::with_analytical_pva, OptimizationMethod::with_analytical_dynamics,
                     OptimizationMethod::broadphase, OptimizationMethod::double_broadphase}) {
    Optimization opt(robot, Task::stop_to_stop(start, end));
    opt.method                   = method;
    opt.constraints.position     = true;
    opt.constraints.velocity     = true;
    opt.constraints.acceleration = true;
    opt.max_duration             = 3;
    opt.guess.type               = Guess::custom;
    Array initial_x              = blast::guess_straight_line(&opt);
    initial_x.back()             = 5;  // T above the cap
    initial_x[0]                 = 10; // a control point far outside
    opt.guess.initial_x          = initial_x;
    const Result result          = optimize(&opt);
    INFO("method " << (int) method);
    CHECK(result.nlopt_exit_criteria != NLOPT_INVALID_ARGS);                               // clamped, not refused
    REQUIRE(result.x.size == initial_x.size);
    CHECK(result.x.back() <= 3.0 + 1e-3);                                                  // + the millisecond round-up of the final validation
    const real bound = (real) 3.1416 + opt.control_point_bound_margin * 2 * (real) 3.1416; // the true limits plus the margin
    for (u32 variable = 0; variable + 1 < result.x.size; variable++) {
      CHECK(result.x[variable] <= bound);
      CHECK(result.x[variable] >= -bound);
    }
  }
}
