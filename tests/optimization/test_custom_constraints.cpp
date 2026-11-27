#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include <cmath>
#include <random>
#include "test_helper.hpp"

using namespace blast;

// Custom constraints in every optimization method: per-point constraints (evaluated at every
// sample, reduced to the worst sample per segment by the segment-based methods) and global
// constraints (evaluated once per trajectory). Checks the layout of their rows, their gradients
// against central finite differences, and that a solve enforces them.

#if BLAST_USE_DOUBLES
constexpr real finite_difference_step      = 1e-6;
constexpr real finite_difference_tolerance = 1e-3;
#else
constexpr real finite_difference_step      = (real) 1e-3;
constexpr real finite_difference_tolerance = (real) 1e-1;
#endif

static const Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
static const Array end   = {-0.5825, -1.0700, 1.3892, -1.3196, 0.9927, 1.14};

// ------------------------- custom constraints under test ---------------------

// Reads everything a point constraint may use: the kinematics and efforts of the state, the
// sample's velocity, acceleration, position and time. Smooth, and inactive (all values < 0).
static Array point_mix(Optimization* opt, Array* x, ManipulatorTempData* state, u32 point) {
  (void) x;
  const auto& traj     = opt->bspline.traj;
  const int   n_joints = opt->manip.n_joints;
  Array       values(4);
  values[0] = state->p_j[n_joints].z - 2;
  values[1] = (real) 0.1 * traj.vel(1, point) * traj.vel(1, point) + (real) 0.05 * traj.acc(2, point) - 5;
  values[2] = state->efforts[1] / 200 - 5;
  values[3] = (real) 0.1 * traj.t[point] + (real) 0.01 * traj.pos(0, point) - 10;
  return values;
}

// A second point constraint, to check that several are laid out one after the other.
static Array point_tool_x(Optimization* opt, Array* x, ManipulatorTempData* state, u32 point) {
  (void) x, (void) point;
  return Array(1, state->p_j[opt->manip.n_joints].x - 3);
}

// Depends on the duration directly and, through a recomputed trajectory, on every variable.
static Array global_mix(Optimization* opt, Array* x) {
  Array values(2);
  values[0] = x->back() - 50;
  opt->bspline.compute_trajectory(*x, opt->task);
  real sum = 0;
  for (u32 entry = 0; entry < opt->bspline.traj.vel.size; entry++)
    sum += opt->bspline.traj.vel.data[entry] * opt->bspline.traj.vel.data[entry];
  values[1] = (real) 1e-3 * sum - 1; // a small offset: float keeps the finite differences accurate
  return values;
}

// Binding constraints for the solves; the limits come through custom_data.
struct Limits {
  real joint_0_speed = 0;
  real min_duration  = 0;
};

static Array point_joint_0_speed(Optimization* opt, Array* x, ManipulatorTempData* state, u32 point) {
  (void) x, (void) state;
  const auto* limits = static_cast<const Limits*>(opt->custom_data);
  return Array(1, std::abs(opt->bspline.traj.vel(0, point)) / limits->joint_0_speed - 1);
}

static Array global_min_duration(Optimization* opt, Array* x) {
  const auto* limits = static_cast<const Limits*>(opt->custom_data);
  return Array(1, (limits->min_duration - x->back()) / limits->min_duration);
}

// ------------------------- helpers -------------------------------------------

static bool is_segment_based(OptimizationMethod method) {
  return method == OptimizationMethod::with_segments || method == OptimizationMethod::broadphase || method == OptimizationMethod::double_broadphase;
}

static const char* method_name(OptimizationMethod method) {
  switch (method) {
    case OptimizationMethod::with_segments: return "with_segments";
    case OptimizationMethod::baseline: return "baseline";
    case OptimizationMethod::with_analytical_pva: return "with_analytical_pva";
    case OptimizationMethod::with_analytical_dynamics: return "with_analytical_dynamics";
    case OptimizationMethod::broadphase: return "broadphase";
    case OptimizationMethod::double_broadphase: return "double_broadphase";
  }
  return "?";
}

// Every built-in constraint on (the analytical gradient fills assume all their rows exist), a few
// spheres, and the custom constraints under test.
static Optimization make_problem(OptimizationMethod method) {
  World world;
  world.add_sphere(Vec3{(real) 0.35, (real) 0.25, (real) 0.45}, (real) 0.08);
  world.add_sphere(Vec3{-(real) 0.2, (real) 0.45, (real) 0.30}, (real) 0.06);
  world.add_sphere(Vec3{(real) 0.1, -(real) 0.4, (real) 0.6}, (real) 0.07);
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, end));
  opt.method = method;
  opt.world  = world;
  create_static_bounding_volume_hierarchy(opt.world, opt.world.static_bounding_volume_hierarchy); // the broadphase methods read it
  auto& selection    = opt.constraints;
  selection.position = selection.velocity = selection.acceleration = selection.torque = selection.tool_speed = true;
  selection.self_collisions = selection.external_collisions = true;
  selection.add_point_constraint(point_mix, 4);
  selection.add_point_constraint(point_tool_x, 1);
  selection.add_global_constraint(global_mix, 2);
  if (is_segment_based(method)) {
    initialize_optimization_with_segments(&opt);
    n_con_with_segments(&opt);
  } else {
    initialize_optimization(&opt);
    n_con(&opt);
  }
  return opt;
}

static Array random_x(Optimization& opt, std::mt19937& random_engine, int trial) {
  std::normal_distribution<real> noise(0, (real) 0.3);
  Array                          x = blast::guess_straight_line(&opt);
  x.back()                         = 1 + (real) 0.2 * (real) trial;
  for (u32 variable = 0; variable + 1 < x.size; variable++)
    x[variable] += noise(random_engine);
  return x;
}

// Constraint values (and, when gradient is given, the gradient in NLopt layout grad[row * x_len + variable]) of a method.
static void evaluate(OptimizationMethod method, Optimization& opt, const Array& x, Array& constraints, real* gradient) {
  const u32 m = constraints.size, x_len = x.size;
  switch (method) {
    case OptimizationMethod::baseline: nlopt_constraints(m, constraints.data, x_len, x.data, gradient, &opt); return;
    case OptimizationMethod::with_analytical_pva: nlopt_constraints_with_analytical_pva(m, constraints.data, x_len, x.data, gradient, &opt); return;
    case OptimizationMethod::with_analytical_dynamics: nlopt_constraints_with_analytical_dynamics(m, constraints.data, x_len, x.data, gradient, &opt); return;
    case OptimizationMethod::with_segments: nlopt_constraints_with_segments(m, constraints.data, x_len, x.data, gradient, &opt); return;
    case OptimizationMethod::broadphase: nlopt_constraints_with_broadphase(m, constraints.data, x_len, x.data, gradient, &opt); return;
    case OptimizationMethod::double_broadphase: nlopt_constraints_with_double_broadphase(m, constraints.data, x_len, x.data, gradient, &opt); return;
  }
}

// Rows of the custom constraints, in constraint-vector order.
static std::vector<u32> custom_rows(const Optimization& opt, OptimizationMethod method) {
  const u32        n_point = opt.constraints.n_point_constraint_values(), n_global = opt.constraints.n_global_constraint_values();
  const u32        m       = opt.constraints.n_constraints;
  const u32        n_blocks = is_segment_based(method) ? opt.bspline.n_ctrl - opt.bspline.degree : opt.bspline.n_points;
  const u32        per_block = (m - n_global) / n_blocks;
  std::vector<u32> rows;
  for (u32 block = 0; block < n_blocks; block++)
    for (u32 row = 0; row < n_point; row++)
      rows.push_back(block * per_block + per_block - n_point + row);
  for (u32 row = 0; row < n_global; row++)
    rows.push_back(m - n_global + row);
  return rows;
}

// The point constraints at a sample, computed by the test.
static Array expected_point_values(Optimization& opt, const Array& x, u32 point) {
  ManipulatorTempData state;
  auto                pos = opt.bspline.traj.pos.col(point);
  auto                vel = opt.bspline.traj.vel.col(point);
  auto                acc = opt.bspline.traj.acc.col(point);
  forward_kinematics(opt.manip, state, pos);
  compute_collision_model(opt.manip, state);
  dynamics(opt.manip, state, vel, acc);
  Array x_copy = x;
  Array a      = point_mix(&opt, &x_copy, &state, point);
  Array b      = point_tool_x(&opt, &x_copy, &state, point);
  Array values(a.size + b.size);
  for (u32 i = 0; i < a.size; i++)
    values[i] = a[i];
  values[a.size] = b[0];
  return values;
}

static const OptimizationMethod all_methods[] = {OptimizationMethod::baseline, OptimizationMethod::with_analytical_pva, OptimizationMethod::with_analytical_dynamics,
                                                 OptimizationMethod::with_segments, OptimizationMethod::broadphase, OptimizationMethod::double_broadphase};

// ------------------------- tests ---------------------------------------------

TEST_CASE("custom constraints: registration and counts", "[Optimization][Custom]") {
  ConstraintSelection selection;
  CHECK(selection.n_point_constraint_values() == 0);
  CHECK(selection.n_global_constraint_values() == 0);
  selection.add_point_constraint(point_mix, 4);
  selection.add_point_constraint(point_tool_x, 1);
  selection.add_global_constraint(global_mix, 2);
  selection.add_constraint(global_min_duration, 1); // deprecated name of add_global_constraint
  CHECK(selection.n_point_constraint_values() == 5);
  CHECK(selection.n_global_constraint_values() == 3);
  CHECK(selection.global_constraints.size() == 2);

  for (auto method: all_methods) {
    INFO(method_name(method));
    Optimization opt     = make_problem(method);
    Optimization builtin = opt;
    builtin.constraints.point_constraints.clear();
    builtin.constraints.n_point_constraints.clear();
    builtin.constraints.global_constraints.clear();
    builtin.constraints.n_global_constraints.clear();
    if (is_segment_based(method)) {
      n_con_with_segments(&builtin);
      const u32 n_segments = opt.bspline.n_ctrl - opt.bspline.degree;
      CHECK(opt.constraints.n_constraints_per_segment == builtin.constraints.n_constraints_per_segment + 5);
      CHECK(opt.constraints.n_constraints == builtin.constraints.n_constraints + n_segments * 5 + 2);
    } else {
      n_con(&builtin);
      CHECK(opt.constraints.n_constraints == builtin.constraints.n_constraints + opt.bspline.n_points * 5 + 2);
    }
  }
}

TEST_CASE("custom constraints: values and layout in every method", "[Optimization][Custom]") {
  std::mt19937 random_engine(3);
  for (auto method: all_methods) {
    INFO(method_name(method));
    Optimization opt = make_problem(method);
    const Array  x   = random_x(opt, random_engine, 1);
    Array        constraints(opt.constraints.n_constraints);
    evaluate(method, opt, x, constraints, nullptr);

    // expected values, from the trajectory at x
    Optimization reference = opt;
    reference.bspline.compute_trajectory(x, reference.task);
    const auto rows     = custom_rows(opt, method);
    const u32  n_point  = 5;
    const u32  n_blocks = is_segment_based(method) ? opt.bspline.n_ctrl - opt.bspline.degree : opt.bspline.n_points;
    const u32  per_seg  = opt.bspline.n_points / n_blocks;
    for (u32 block = 0; block < n_blocks; block++) {
      Array expected(n_point, -INF_REAL);
      for (u32 point = block * per_seg; point < (block + 1) * per_seg; point++) { // one sample per block for the point-based methods
        const Array values = expected_point_values(reference, x, point);
        for (u32 row = 0; row < n_point; row++)
          expected[row] = std::max(expected[row], values[row]);
      }
      for (u32 row = 0; row < n_point; row++)
        CHECK(constraints[rows[block * n_point + row]] == Approx(expected[row]).margin(test::abs_tol));
    }
    Array       x_copy = x;
    const Array global = global_mix(&reference, &x_copy);
    CHECK(constraints[rows[rows.size() - 2]] == Approx(global[0]).margin(test::abs_tol));
    CHECK(constraints[rows[rows.size() - 1]] == Approx(global[1]).margin(test::abs_tol));

    // the tightening shift moves every custom row, and nothing else
    Array shifted(opt.constraints.n_constraints);
    opt.custom_constraint_shift = (real) 0.25;
    evaluate(method, opt, x, shifted, nullptr);
    opt.custom_constraint_shift = 0;
    std::vector<bool> is_custom(constraints.size, false);
    for (auto row: rows)
      is_custom[row] = true;
    int wrong = 0;
    for (u32 row = 0; row < constraints.size; row++) {
      const real expected = constraints[row] + (is_custom[row] ? (real) 0.25 : 0);
      wrong += std::abs(shifted[row] - expected) > test::approx_tol * (1 + std::abs(expected));
    }
    CHECK(wrong == 0);
  }
}

TEST_CASE("custom constraints: gradients match finite differences in every method", "[Optimization][Custom]") {
  for (auto method: all_methods) {
    INFO(method_name(method));
    Optimization       opt     = make_problem(method);
    const u32          m       = opt.constraints.n_constraints;
    const auto         rows    = custom_rows(opt, method);
    std::mt19937       random_engine(7);
    int                mismatched = 0, checked = 0;
    for (int trial = 0; trial < 4; trial++) {
      const Array x = random_x(opt, random_engine, trial);
      const u32   n = x.size;
      Array       constraints(m), gradient(m * n);
      evaluate(method, opt, x, constraints, gradient.data);
      for (u32 variable = 0; variable < n; variable++) {
        const real step    = (variable + 1 == n) ? finite_difference_step * x[variable] : finite_difference_step;
        Array      x_plus  = x;
        Array      x_minus = x;
        x_plus[variable] += step;
        x_minus[variable] -= step;
        Array constraints_plus(m), constraints_minus(m);
        evaluate(method, opt, x_plus, constraints_plus, nullptr);
        evaluate(method, opt, x_minus, constraints_minus, nullptr);
        for (auto row: rows) {
          const real finite_difference = (constraints_plus[row] - constraints_minus[row]) / (2 * step);
          // Segment rows are a max over samples: where the argmax switches inside +/-step the finite
          // difference mixes two pieces, so allow a few.
          mismatched += std::abs(gradient[row * n + variable] - finite_difference) > finite_difference_tolerance * (std::abs(finite_difference) + finite_difference_tolerance);
          checked++;
        }
      }
    }
    INFO("mismatching entries: " << mismatched << "/" << checked);
    CHECK(mismatched <= checked / 500);
  }
}

TEST_CASE("custom constraints are enforced by every method", "[Optimization][Custom]") {
  const Array near_end = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  for (auto method: all_methods) {
    INFO(method_name(method));
    auto make = [&]() {
      Optimization opt(make_UR5e(), Task::stop_to_stop(start, near_end));
      opt.method                   = method;
      opt.constraints.position     = true;
      opt.constraints.velocity     = true;
      opt.constraints.acceleration = true;
      opt.guess.type               = Guess::straight_line; // deterministic: the default random guess depends on earlier tests
      opt.max_tries                = 3;
      if (method == OptimizationMethod::with_analytical_pva) { // its gradient fill assumes every built-in row exists
        opt.constraints.torque = opt.constraints.tool_speed = opt.constraints.self_collisions = opt.constraints.external_collisions = true;
        opt.world.add_sphere(Vec3{5, 5, 5}, (real) 0.1);
      }
      return opt;
    };

    // reference: no custom constraint
    Optimization reference_opt = make();
    const Result reference     = optimize(&reference_opt);
    REQUIRE(reference.success);
    real joint_0_speed = 0;
    for (u32 point = 0; point < reference.trajectory.vel.cols; point++)
      joint_0_speed = std::max(joint_0_speed, std::abs(reference.trajectory.vel(0, point)));
    REQUIRE(joint_0_speed > 0);

    Limits limits;
    limits.joint_0_speed = (real) 0.5 * joint_0_speed;
    limits.min_duration  = (real) 1.3 * reference.x.back();

    { // point constraint
      Optimization opt = make();
      opt.custom_data  = &limits;
      opt.constraints.add_point_constraint(point_joint_0_speed, 1);
      const Result result = optimize(&opt);
      REQUIRE(result.success);
      for (u32 point = 0; point < result.trajectory.vel.cols; point++) // the validation trajectory, sampled every millisecond
        CHECK(std::abs(result.trajectory.vel(0, point)) <= limits.joint_0_speed);
      CHECK(result.x.back() > reference.x.back()); // slowing joint 0 down lengthens the time-optimal trajectory
    }

    { // global constraint
      Optimization opt = make();
      opt.custom_data  = &limits;
      opt.constraints.add_global_constraint(global_min_duration, 1);
      const Result result = optimize(&opt);
      REQUIRE(result.success);
      CHECK(result.x.back() >= limits.min_duration);
      CHECK(result.x.back() <= limits.min_duration * (real) 1.05); // the time objective pushes it to the bound
    }
  }
}
