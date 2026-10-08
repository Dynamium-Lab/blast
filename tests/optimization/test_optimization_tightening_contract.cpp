#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include <cmath>

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
    Array  x = blast::guess_straight_line(&opt);
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
  CHECK(opt.manip.position_max[0] == Approx(3.1416)); // restored to the true bound
}

// ---------------------------------------------------------------------------
// 4. Position is tightened like every other limit, but never past a task endpoint, and a
// reported success never leaves the TRUE position bounds. Untightened, the (1 + tol)
// acceptance slack let a solution exceed a joint limit by up to tol * range / 2.
// ---------------------------------------------------------------------------
TEST_CASE("position is tightened up to the endpoints and success respects the true bounds",
          "[Optimization]") {
  Array start = {3.14, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.14};

  Manipulator  robot = UR5e_narrow_limits();
  Task         task  = Task::stop_to_stop(start, end);
  Optimization opt(robot, task);
  opt.success_tolerance = 0.01;
  opt.collision_buffer  = 0.001;
  initialize_optimization_with_segments(&opt);

  auto tolerance_snapshot = tighten_for_success_tolerance(&opt);
  for (int joint = 0; joint < opt.manip.n_joints; joint++) {
    const real endpoint_max = std::max(opt.task(joint, 0), opt.task(joint, 3));
    const real endpoint_min = std::min(opt.task(joint, 0), opt.task(joint, 3));
    CHECK(opt.manip.position_max[joint] >= endpoint_max); // never cuts off an endpoint
    CHECK(opt.manip.position_min[joint] <= endpoint_min);
    CHECK(opt.manip.position_max[joint] <= (real) 3.1416);
    CHECK(opt.manip.position_min[joint] >= (real) -3.1416);
  }
  CHECK(opt.manip.position_max[2] < (real) 3.1416); // mid-range joint: really tightened
  restore_from_tolerance(&opt, tolerance_snapshot);

  Optimization run(robot, task);
  run.success_tolerance      = 0.01;
  run.collision_buffer       = 0.001;
  run.guess.type             = Guess::custom;
  run.guess.initial_x        = blast::guess_straight_line(&run);
  run.guess.initial_x.back() = 2.0;
  const Result result        = optimize(&run);
  if (result.success) {
    const Matrix& pos = result.trajectory.pos;
    for (u32 point = 0; point < pos.cols; point++)
      for (u32 joint = 0; joint < pos.rows; joint++) {
        CHECK(pos(joint, point) <= robot.position_max[joint]);
        CHECK(pos(joint, point) >= robot.position_min[joint]);
      }
  }
}

// ---------------------------------------------------------------------------
// 4. Every solve method must apply the same -d * collision_scale transform to a
// self-collision constraint. with_segments (and compute_constraints) do; broadphase
// and double_broadphase previously computed the raw, unscaled -d instead, so under
// tightening they'd accept trajectories with self-collision margins the accept gate
// never actually validated. A collision_scale of 1 (untightened) can't catch this --
// the missing factor is a no-op -- so this must run under tightening.
// ---------------------------------------------------------------------------
TEST_CASE("tightened self-collision constraint matches across with_segments/broadphase/double_broadphase", "[Optimization]") {
  Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};

  Manipulator  robot = make_UR5e();
  Task         task  = Task::stop_to_stop(start, end);
  Optimization opt(robot, task);
  opt.constraints.self_collisions = true;
  opt.success_tolerance           = 0.01;
  opt.collision_buffer            = 0.001; // collision_scale = tol/buffer = 10

  REQUIRE(validate_task(&opt) == true);    // legal on true geometry

  initialize_optimization_with_segments(&opt);
  n_con_with_segments(&opt);
  auto snap = tighten_for_success_tolerance(&opt);

  const u32 xlen = (u32) opt.bspline.x_len(opt.task);
  const u32 m    = (u32) opt.constraints.n_constraints;
  Array     x    = blast::guess_straight_line(&opt);

  Array  cons_segments(m), cons_broadphase(m), cons_double(m);
  Matrix grad_segments(xlen, m), grad_broadphase(xlen, m), grad_double(xlen, m);
  for (u32 i = 0; i < xlen * m; i++) {
    grad_segments.data[i]   = 0;
    grad_broadphase.data[i] = 0;
    grad_double.data[i]     = 0;
  }

  constraints_and_gradients_with_segments(x, opt, cons_segments, grad_segments);
  constraints_and_gradients_with_broadphase(x, opt, cons_broadphase, grad_broadphase);
  constraints_and_gradients_with_double_broadphase(x, opt, cons_double, grad_double);

  CHECK(is_close(cons_segments, cons_broadphase));
  CHECK(is_close(cons_segments, cons_double));
  CHECK(is_close(grad_segments, grad_broadphase));
  CHECK(is_close(grad_segments, grad_double));

  restore_from_tolerance(&opt, snap);
}

// ---------------------------------------------------------------------------
// 5. The collision analogue of the position test above. A task start that clears an obstacle by less than the buffer is
// acceptable (clearance > 0), but the full buffer would make its collision row violated with a zero
// gradient (the start is pinned). The first segment's rows target what the start achieves instead,
// every other row keeps the buffer, and acceptance (c < tol <=> d + buffer > 0) does not move.
// ---------------------------------------------------------------------------
// What "exactly" means for distances and rows below. Float rows carry ~1e-7 relative distance
// error times the row scale (tolerance / target, ~20 here); the bug they guard is +tolerance/2.
#if BLAST_USE_DOUBLES
constexpr real exact_tolerance = 1e-9;
#else
constexpr real exact_tolerance = 1e-5;
#endif

TEST_CASE("collision rows at a pinned endpoint target what the endpoint achieves", "[Optimization]") {
  const real tolerance = 0.01, buffer = 0.001, start_clearance = 0.0005; // start clears the sphere by half the buffer
  Array      start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array      end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.14};

  Manipulator robot        = make_UR5e();
  int         tool_capsule = robot._n_caps - 1;
  World       world;
  {
    ManipulatorTempData manip_data;
    forward_kinematics(robot, manip_data, start);
    compute_collision_model(robot, manip_data);
    const Capsule capsule       = manip_data.capsule_list[tool_capsule];
    const Vec3    axis          = (capsule.p2 - capsule.p1) / norm(capsule.p2 - capsule.p1);
    const real    sphere_radius = 0.02;
    world.add_sphere(capsule.p2 + axis * (capsule.radius + sphere_radius + start_clearance), sphere_radius); // beyond the capsule's end, on its axis
    for (int capsule_id = 0; capsule_id < robot._n_caps; capsule_id++)
      if (capsule_id != tool_capsule)
        REQUIRE(distance(manip_data.capsule_list[capsule_id], world.spheres[0]) > 2 * buffer); // only the tool capsule is close
    REQUIRE(distance(capsule, world.spheres[0]) == Approx(start_clearance).margin(exact_tolerance));
  }

  Task         task = Task::stop_to_stop(start, end);
  Optimization opt(robot, task);
  opt.world                           = world;
  opt.constraints.external_collisions = true;
  opt.success_tolerance               = tolerance;
  opt.collision_buffer                = buffer;
  REQUIRE(validate_task(&opt) == true);

  initialize_optimization_with_segments(&opt);
  n_con_with_segments(&opt);
  auto tolerance_snapshot = tighten_for_success_tolerance(&opt);

  REQUIRE(opt.endpoint_targets_active);
  CHECK(opt.endpoint_collision_target[0][tool_capsule] == Approx(start_clearance).margin(exact_tolerance)); // what the start achieves
  for (int capsule_id = 0; capsule_id < robot._n_caps; capsule_id++)
    if (capsule_id != tool_capsule)
      CHECK(opt.endpoint_collision_target[0][capsule_id] == Approx(buffer)); // far capsules keep the full buffer

  // The endpoint must satisfy its own rows. On a trajectory that holds the start (goal == start)
  // every sample IS the start, so every row of the first and last segment has to be <= 0. Without
  // the per-row target the tool capsule's row sits at +tol/2 there, and its worst sample is beside
  // the pinned start where the free control points barely reach: gradient ~1e-3 against ~4 for the
  // same row mid-trajectory. Not exactly zero, so the zero-gradient check of 1 cannot see it; the
  // SQP needs a huge step to move it, which is the round-off line-search failure (nlopt -4).
  {
    Optimization hold(robot, Task::stop_to_stop(start, start));
    hold.world                           = world;
    hold.constraints.external_collisions = true;
    hold.success_tolerance               = tolerance;
    hold.collision_buffer                = buffer;
    initialize_optimization_with_segments(&hold);
    n_con_with_segments(&hold);
    auto hold_tolerance_snapshot = tighten_for_success_tolerance(&hold);

    const u32 n_variables             = (u32) hold.bspline.x_len(hold.task);
    const u32 n_constraints           = (u32) hold.constraints.n_constraints;
    const int constraints_per_segment = hold.constraints.n_constraints_per_segment;
    const int n_segments              = (int) n_constraints / constraints_per_segment;
    Array     x                       = blast::guess_straight_line(&hold);
    x.back()                          = 2.0; // start == goal: the limit-derived duration would be 0
    Array     constraints(n_constraints);
    Matrix    gradient(n_variables, n_constraints);
    constraints_and_gradients_with_segments(x, hold, constraints, gradient);
    real worst_boundary = -INF_REAL, worst_interior = -INF_REAL;
    for (int segment = 0; segment < n_segments; segment++)
      for (int row = 0; row < constraints_per_segment; row++) {
        const real constraint = constraints[segment * constraints_per_segment + row];
        if (segment == 0 || segment == n_segments - 1)
          worst_boundary = std::max(worst_boundary, constraint);
        else
          worst_interior = std::max(worst_interior, constraint);
      }
    CHECK(worst_boundary <= exact_tolerance);
    CHECK(worst_interior > 0); // interior rows keep the full buffer: not vacuous
    restore_from_tolerance(&hold, hold_tolerance_snapshot);
  }

  // Acceptance is the untargeted row's, on both sides of distance + buffer = 0, for the targeted
  // row and for an interior one.
  const int n_segments = (int) opt.bspline.n_ctrl - (int) opt.bspline.degree;
  for (int segment: {0, n_segments / 2}) {
    const auto row = collision_row(opt, segment, n_segments, tool_capsule);
    for (real distance: {-buffer - 1e-6, -buffer + 1e-6}) {
      const real constraint = -(distance + row.shift) * row.scale;
      CHECK((constraint < tolerance) == (distance + buffer > 0));
    }
  }
  CHECK(collision_row(opt, n_segments / 2, n_segments, tool_capsule).shift == 0);

  // Point-based methods: the pinned start/goal samples take the boundary target, an interior
  // sample the full buffer. At the start's own clearance its row is satisfied (<= 0).
  {
    const u32  n_points       = opt.bspline.n_points;
    const real start_distance = start_clearance - buffer; // the start's distance in the buffered geometry
    CHECK(collision_constraint_at_point(opt, 0, tool_capsule, start_distance) <= exact_tolerance);
    CHECK(collision_constraint_at_point(opt, n_points / 2, tool_capsule, start_distance) > 0);
    for (real distance: {-buffer - 1e-6, -buffer + 1e-6})
      for (u32 point: {0u, n_points / 2, n_points - 1})
        CHECK((collision_constraint_at_point(opt, point, tool_capsule, distance) < tolerance) == (distance + buffer > 0));
  }

  restore_from_tolerance(&opt, tolerance_snapshot);
  CHECK_FALSE(opt.endpoint_targets_active);
}
