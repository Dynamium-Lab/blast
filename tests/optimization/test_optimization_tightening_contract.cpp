#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include <cmath>
#include "test_helper.hpp"

using namespace blast;

// Pins the properties tighten_for_success_tolerance() has to have.

// +/-pi joint limits: make_UR5e() uses +/-2pi, where the endpoints below sit mid-range and
// the failure cannot appear. It needs an endpoint closer to its limit than the tightening
// moves that limit.
static Manipulator UR5e_narrow_limits() {
  Manipulator m = make_UR5e();
  for (int j = 0; j < m.n_joints; j++) {
    m.position_min[j] = -3.1416;
    m.position_max[j] = 3.1416;
  }
  return m;
}

// ---------------------------------------------------------------------------
// 1. No constraint row may be violated with a zero gradient: the SQP linearizes such a row
// as c + grad.d <= 0 with grad = 0, which no step can satisfy.
// ---------------------------------------------------------------------------
TEST_CASE("tightening introduces no violated zero-gradient constraint row", "[Optimization]") {
  // 3.14 against a +/-3.1416 limit: inside by less than the tightening moves the bound.
  Array start = {3.14, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.14};

  Manipulator  robot = UR5e_narrow_limits();
  Task         task  = Task::stop_to_stop(start, end);
  Optimization opt(robot, task);
  opt.constraints.torque = true;
  opt.success_tolerance  = 0.01;
  opt.collision_buffer   = 0.001;

  REQUIRE(validate_task(&opt) == true); // legal on true geometry

  initialize_optimization_with_segments(&opt);
  n_con_with_segments(&opt);
  auto snap = tighten_for_success_tolerance(&opt);

  const u32 xlen = (u32) opt.bspline.x_len(opt.task);
  const u32 m    = (u32) opt.constraints.n_constraints;

  // A zero-gradient row does not depend on the decision variables, so one deterministic
  // guess is as strong as many.
  int offenders = 0;
  {
    Array  x = blast::test::straight_line_guess(opt, start, end);
    Array  cons(m);
    Matrix grad(xlen, m);
    for (u32 i = 0; i < xlen * m; i++)
      grad.data[i] = 0;
    constraints_and_gradients_with_segments(x, opt, cons, grad);

    for (u32 i = 0; i < m; i++) {
      if (cons[i] <= 0)
        continue;
      real gnorm = 0;
      for (u32 k = 0; k < xlen; k++)
        gnorm += grad(k, i) * grad(k, i);
      if (std::sqrt(gnorm) < 1e-12)
        offenders++; // violated AND unreachable
    }
  }
  CHECK(offenders == 0);
  restore_from_tolerance(&opt, snap);
}

// ---------------------------------------------------------------------------
// 2. A constraint within tolerance must imply the true geometry is clear:
//    c = -(d_true - buffer) * (tol/buffer) < tol  <=>  d_true > 0
// ---------------------------------------------------------------------------
TEST_CASE("a constraint within tolerance implies true clearance", "[Optimization]") {
  const real tol = 0.01;
  Array      q   = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};

  auto min_dist = [&](Manipulator m, const World& w) {
    ManipulatorTempData t;
    forward_kinematics(m, t, q);
    compute_collision_model(m, t);
    real worst = INF_REAL;
    for (int c = 0; c < m._n_caps; c++)
      for (const auto& s: w.spheres)
        worst = std::min(worst, distance(t.capsule_list[c], s));
    return worst;
  };

  for (real buffer: {0.0005, 0.001, 0.005, 0.01}) {
    int checked = 0, both_sides = 0;
    // Sweep the obstacle through contact so the implication is tested on both sides.
    for (int step = -25; step <= 25; step++) {
      const real x = 0.095 + (real) step * 0.0008;

      World w;
      w.add_sphere(Vec3{x, 0.0, 0.35}, 0.10);

      Manipulator  robot = UR5e_narrow_limits();
      Task         task  = Task::stop_to_stop(q, q);
      Optimization opt(robot, task);
      opt.world                           = w;
      opt.constraints.external_collisions = true;
      opt.success_tolerance               = tol;
      opt.collision_buffer                = buffer;

      const real d_true = min_dist(opt.manip, opt.world);
      auto       snap   = tighten_for_success_tolerance(&opt);
      const real c      = -min_dist(opt.manip, opt.world) * opt.collision_scale;
      restore_from_tolerance(&opt, snap);

      if (c < tol) // slack for float comparison at the crossing
        CHECK(d_true > -1e-9);
      checked++;
      if (d_true < 0)
        both_sides++;
    }
    CHECK(checked == 51);
    CHECK(both_sides > 0); // the sweep really did cross contact, so this is not vacuous
  }
}

// ---------------------------------------------------------------------------
// 3. restore_from_tolerance() must put the derived unit change back.
// ---------------------------------------------------------------------------
TEST_CASE("restore puts collision_scale back to 1", "[Optimization]") {
  Array        q     = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Manipulator  robot = UR5e_narrow_limits();
  Task         task  = Task::stop_to_stop(q, q);
  Optimization opt(robot, task);
  opt.success_tolerance = 0.01;
  opt.collision_buffer  = 0.001;

  auto snap = tighten_for_success_tolerance(&opt);
  CHECK(opt.collision_scale == Approx(10.0)); // tol / buffer
  restore_from_tolerance(&opt, snap);
  CHECK(opt.collision_scale == Approx(1.0));
  CHECK(opt.manip.position_max[0] == Approx(3.1416)); // position is never tightened
}
