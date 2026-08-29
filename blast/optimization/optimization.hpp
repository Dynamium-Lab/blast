#pragma once

#include <blast>
#include <utility>
#ifdef BLAST_USE_NATIVE_SQP
#include "optimization/sqp.hpp"
#else
#include "nlopt.h"
#endif
#include "optimization/constraints.hpp"
#include "optimization/initial_guess.hpp"
#include "optimization/objective.hpp"


namespace blast {

struct Result {
  bool          success       = false;
  bool          success_false = false;
  real          compute_time  = 0;
  Optimization* opt; // todo: fix this (pointer becomes out of bounds fast)
  Array         x;
  Array         x0;
  nlopt_result  nlopt_exit_criteria = NLOPT_SUCCESS;
  int           num_eval            = 0;
  int           num_tries           = 0;
  Trajectory    trajectory;

  real max_constraint_value             = 0.0;
  int  max_constraint_idx               = 0;
  real max_constraint_more_points_value = 0.0;
  int  max_constraint_more_points_idx   = 0;

  Result() = delete;

  explicit Result(Optimization* optim) {
    opt = optim;
  }
};

inline Optimization::Optimization(const Manipulator& new_manip, const Task& new_task) :
    manip(new_manip),
    bspline(new_manip.n_joints),
    task(new_task.to_matrix()),
    custom_data(nullptr) {
  // Default values
  method     = OptimizationMethod::with_segments;
  guess.type = Guess::random;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.tool_speed   = true;

  objective.time_weight = 1.0;
}

inline Optimization::Optimization(const Manipulator& new_manip, const Task& new_task, const Bspline& new_bspline) :
    manip(new_manip),
    bspline(new_bspline),
    task(new_task.to_matrix()),
    custom_data(nullptr) {
  // Default values
  method               = OptimizationMethod::with_segments;
  guess.type           = Guess::random;
  guess.n_random_shots = 100;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.tool_speed   = true;

  objective.time_weight = 1.0;
}

inline int Optimization::x_len() const {
  return (int) bspline.x_len(task);
}

inline void Optimization::set_manip(Manipulator new_manip) {
  manip = std::move(new_manip);
}

inline void Optimization::set_bspline(Bspline new_bspline) {
  bspline = std::move(new_bspline);
}

inline void Optimization::set_guess(Guess new_guess) {
  guess = std::move(new_guess);
}

inline void Optimization::set_constraints(ConstraintSelection new_constraints) {
  constraints = std::move(new_constraints);
}

inline void Optimization::set_objective(Objective new_objective) {
  objective = std::move(new_objective);
}

inline void Optimization::set_task(const Task& new_task) {
  task = new_task.to_matrix();
}

inline void Optimization::set_world(World new_world) {
  world = std::move(new_world);
}

inline void n_con(Optimization* opt) {
  const auto n_points            = opt->bspline.n_points;
  const auto n_joints            = opt->manip.n_joints;
  const auto n_constraints_basic = n_points * n_joints;
  opt->constraints.n_constraints = 0;
  if (opt->constraints.position)
    opt->constraints.n_constraints += n_constraints_basic;
  if (opt->constraints.velocity)
    opt->constraints.n_constraints += n_constraints_basic;
  if (opt->constraints.acceleration)
    opt->constraints.n_constraints += n_constraints_basic;
  if (opt->constraints.torque)
    opt->constraints.n_constraints += n_constraints_basic;

  if (opt->constraints.tool_speed)
    opt->constraints.n_constraints += n_points;
  if (opt->constraints.self_collisions)
    opt->constraints.n_constraints += n_points;
  if (opt->constraints.external_collisions)
    opt->constraints.n_constraints += n_points * opt->manip._n_caps;


  for (auto& n: opt->constraints.n_extra_constraints)
    opt->constraints.n_constraints += n;
}

inline void initialize_optimization(Optimization* opt) {
  // Constraints
  if (!opt->constraints.external_collisions) {
    opt->constraints.n_collision_constraints = 0;
  }

  // todo: swap INF for large value

  opt->lb        = Array(opt->x_len(), -HUGE_VAL);
  opt->ub        = Array(opt->x_len(), -HUGE_VAL);
  opt->lb.back() = 0.1;
  opt->ub.back() = 60.0;

  // Task
  auto task = opt->task;
  for (int i = 0; i < opt->manip.n_joints; i++) {
    while (std::abs(task(i, 0) - task(i, 3)) > PI) { // PI since the highest angle between two points on a circle is halfway through, so PI
      // Try updating start value
      real tmp_value = task(i, 0);
      if (task(i, 0) < task(i, 3)) {
        tmp_value += 2 * PI;
      } else {
        tmp_value -= 2 * PI;
      }
      if (tmp_value > opt->manip.position_min[i] && tmp_value < opt->manip.position_max[i]) {
        task(i, 0) = tmp_value;
      } else {
        // Try updating end value
        tmp_value = task(i, 3);
        if (task(i, 3) < task(i, 0)) {
          tmp_value += 2 * PI;
        } else {
          tmp_value -= 2 * PI;
        }
        if (tmp_value > opt->manip.position_min[i] && tmp_value < opt->manip.position_max[i]) {
          task(i, 0) = tmp_value;
        } else {
          break; // Nothing to be done
        }
      }
    }
  }
  opt->task = task;
}

struct ToleranceSnapshot {
  std::array<real, MAX_JOINTS>   position_min;
  std::array<real, MAX_JOINTS>   position_max;
  std::array<real, MAX_JOINTS>   velocity_max;
  std::array<real, MAX_JOINTS>   acceleration_max;
  std::array<real, MAX_JOINTS>   torque_max;
  real                           tool_speed_max;
  std::array<real, MAX_CAPSULES> capsule_radius; // [0, manip._n_caps)
  real                           base_sphere_radius;
  World                          world;
};

// Tightens the limits/collision geometry constraints are evaluated against by success_tolerance,
// so that whenever the solver's own tolerance-based early-out triggers, the original,
// untightened limits are still satisfied. Must be paired with restore_from_tolerance() before
// the Optimization is returned or reused.
//
// success_tolerance is the only bar: it is the solver's con_tol, it sets the tightening here,
// and it is the post-solve accept gate. The collision algebra then cancels exactly, for any
// tolerance and any buffer:
//
//   c = -(d_true - buffer) * (tol/buffer) < tol   <=>   d_true > 0
//
inline ToleranceSnapshot tighten_for_success_tolerance(Optimization* opt) {
  ToleranceSnapshot snap;
  snap.position_min       = opt->manip.position_min;
  snap.position_max       = opt->manip.position_max;
  snap.velocity_max       = opt->manip.velocity_max;
  snap.acceleration_max   = opt->manip.acceleration_max;
  snap.torque_max         = opt->manip.torque_max;
  snap.tool_speed_max     = opt->manip.tool_speed_max;
  snap.base_sphere_radius = opt->manip._base_sphere.radius;
  for (int i = 0; i < opt->manip._n_caps; i++)
    snap.capsule_radius[i] = opt->manip._collision_model[i].radius;
  snap.world = opt->world;

  // Snapshot taken, nothing modified, so restore_from_tolerance() is a no-op.
  if (!opt->tighten_for_tolerance) {
    opt->collision_scale = 1.0;
    return snap;
  }

  const real tol       = opt->success_tolerance;
  const real ratio_div = 1 + tol;

  // The unit change that makes the dimensionless success_tolerance land on the buffer.
  // A non-positive buffer means no unit change.
  const real buffer    = opt->collision_buffer > 0 ? opt->collision_buffer : tol;
  opt->collision_scale = tol / buffer;

  // Position is deliberately NOT tightened. Task::stop_to_stop pins the boundary control
  // points, so the position rows of the boundary segments are constant in the decision
  // variables -- gradient exactly 0. Shrinking the range can flip such a row to violated,
  // and a violated row with no gradient makes the SQP subproblem infeasible at every
  // iterate. validate_task() already checks start and goal against the true bounds.
  for (int j = 0; j < opt->manip.n_joints; j++) {
    opt->manip.velocity_max[j] /= ratio_div;
    opt->manip.acceleration_max[j] /= ratio_div;
    opt->manip.torque_max[j] /= ratio_div;
  }
  opt->manip.tool_speed_max /= ratio_div;

  // Halved so a pairwise check between two inflated objects sums to exactly the buffer.
  const real col_margin = buffer / 2;
  for (int i = 0; i < opt->manip._n_caps; i++)
    opt->manip._collision_model[i].radius += col_margin;
  opt->manip._base_sphere.radius += col_margin;

  for (auto& box: opt->world.boxes)
    box.extents += Vec3(col_margin, col_margin, col_margin);
  for (auto& sphere: opt->world.spheres)
    sphere.radius += col_margin;
  for (auto& capsule: opt->world.capsules)
    capsule.radius += col_margin;
  for (auto& dyn_box: opt->world.dynamic_boxes)
    for (auto& box: dyn_box.trajectory)
      box.extents += Vec3(col_margin, col_margin, col_margin);
  for (auto& dyn_sphere: opt->world.dynamic_spheres)
    for (auto& sphere: dyn_sphere.trajectory)
      sphere.radius += col_margin;
  for (auto& dyn_capsule: opt->world.dynamic_capsules)
    for (auto& capsule: dyn_capsule.trajectory)
      capsule.radius += col_margin;
  for (auto& door: opt->world.dynamic_doors)
    door.extents += Vec3(col_margin, col_margin, col_margin);

  return snap;
}

inline void restore_from_tolerance(Optimization* opt, const ToleranceSnapshot& snap) {
  opt->manip.position_min        = snap.position_min;
  opt->manip.position_max        = snap.position_max;
  opt->manip.velocity_max        = snap.velocity_max;
  opt->manip.acceleration_max    = snap.acceleration_max;
  opt->manip.torque_max          = snap.torque_max;
  opt->manip.tool_speed_max      = snap.tool_speed_max;
  opt->manip._base_sphere.radius = snap.base_sphere_radius;
  for (int i = 0; i < opt->manip._n_caps; i++)
    opt->manip._collision_model[i].radius = snap.capsule_radius[i];
  opt->world           = snap.world;
  opt->collision_scale = 1.0; // raw metres again outside the tightening window
}

inline Result optimize_baseline_impl(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  auto tolerance_snapshot = tighten_for_success_tolerance(opt); // restored below, next to opt->guess restore

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = opt->max_eval;
  stop.maxtime    = opt->max_time;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 30.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, opt->max_time);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, opt->max_eval);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    PROFILE_SCOPE("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("NLopt optimization");
#endif

      real f = HUGE_VAL;
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.degree, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        try_count++; // count this try before breaking, so num_tries reports tries made
        break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original
  restore_from_tolerance(opt, tolerance_snapshot);

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// ------------------------- Accelerated with segments ---------------------------

// todo: separate task and initialization
inline void initialize_optimization_with_segments(Optimization* opt) {
  // Constraints
  if (!opt->constraints.external_collisions) {
    opt->constraints.n_collision_constraints = 0;
  }

  // todo: swap INF for large value

  opt->lb        = Array(opt->x_len(), -HUGE_VAL);
  opt->ub        = Array(opt->x_len(), HUGE_VAL);
  opt->lb.back() = 0.1;
  opt->ub.back() = 60.0;

  // Task
  auto task = opt->task;
  for (int i = 0; i < opt->manip.n_joints; i++) {
    while (std::abs(task(i, 0) - task(i, 3)) > PI) { // PI since the highest angle between two points on a circle is halfway through, so PI
      // Try updating start value
      real tmp_value = task(i, 0);
      if (task(i, 0) < task(i, 3)) {
        tmp_value += 2 * PI;
      } else {
        tmp_value -= 2 * PI;
      }
      if (tmp_value > opt->manip.position_min[i] && tmp_value < opt->manip.position_max[i]) {
        task(i, 0) = tmp_value;
      } else {
        // Try updating end value
        tmp_value = task(i, 3);
        if (task(i, 3) < task(i, 0)) {
          tmp_value += 2 * PI;
        } else {
          tmp_value -= 2 * PI;
        }
        if (tmp_value > opt->manip.position_min[i] && tmp_value < opt->manip.position_max[i]) {
          task(i, 0) = tmp_value;
        } else {
          break; // Nothing to be done
        }
      }
    }
  }
  opt->task = task;
}

inline void n_con_with_segments(Optimization* opt) {
  const int n_segments = ((int) opt->bspline.n_ctrl - (int) opt->bspline.degree);

  opt->constraints.n_constraints_per_segment = 0;
  if (opt->constraints.position)
    opt->constraints.n_constraints_per_segment += opt->manip.n_joints;
  if (opt->constraints.velocity)
    opt->constraints.n_constraints_per_segment += opt->manip.n_joints;
  if (opt->constraints.acceleration)
    opt->constraints.n_constraints_per_segment += opt->manip.n_joints;
  if (opt->constraints.torque)
    opt->constraints.n_constraints_per_segment += opt->manip.n_joints;
  if (opt->constraints.tool_speed)
    opt->constraints.n_constraints_per_segment += 1;
  if (opt->constraints.self_collisions)
    opt->constraints.n_constraints_per_segment += 1;
  if (opt->constraints.external_collisions) {
    opt->constraints.n_constraints_per_segment += opt->manip._n_caps;
  }

  opt->constraints.n_constraints = n_segments * opt->constraints.n_constraints_per_segment;
}

inline Result optimize_with_segments_impl(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto T1 = get_tick_us();

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization_with_segments(opt);
  n_con_with_segments(opt);

  Result result(opt); // todo: this is expensive
  result.opt->task = opt->task;

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  auto tolerance_snapshot = tighten_for_success_tolerance(opt); // restored below, next to opt->guess restore

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = opt->max_eval;
  stop.maxtime    = opt->max_time;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 30.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_with_segments;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_with_segments, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, opt->max_time);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, opt->max_eval);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    PROFILE_SCOPE("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Initial guess");
#endif
      x         = init_guess_segments(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("NLopt optimization");
#endif

      real f = HUGE_VAL;
      // note: can we initialize grad to 0 here
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation");
#endif
      Array  constraints_points(opt->constraints.n_constraints);
      Matrix gradient;
      constraints_and_gradients_with_segments(x, *opt, constraints_points, gradient);
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.degree, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con_with_segments(&opt_val_more);
      Array  constraints_more_points(opt_val_more.constraints.n_constraints);
      Matrix gradient;
      constraints_and_gradients_with_segments(x, opt_val_more, constraints_more_points, gradient);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        try_count++; // count this try before breaking, so num_tries reports tries made
        break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original
  restore_from_tolerance(opt, tolerance_snapshot);

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// ------------------------- Accelerated functions --------------------------------

inline Result optimize_with_analytical_pva_impl(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  auto tolerance_snapshot = tighten_for_success_tolerance(opt); // restored below, next to opt->guess restore

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.0001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = opt->max_eval;
  stop.maxtime    = opt->max_time;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_with_analytical_pva;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_with_analytical_pva, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, opt->max_time);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, opt->max_eval);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    PROFILE_SCOPE("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("NLopt optimization");
#endif

      real f = HUGE_VAL;
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.degree, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        try_count++; // count this try before breaking, so num_tries reports tries made
        break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original
  restore_from_tolerance(opt, tolerance_snapshot);

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

inline Result optimize_with_analytical_dynamics_impl(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  auto tolerance_snapshot = tighten_for_success_tolerance(opt); // restored below, next to opt->guess restore

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.0001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = opt->max_eval;
  stop.maxtime    = opt->max_time;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_with_analytical_dynamics;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_with_analytical_dynamics, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, opt->max_time);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, opt->max_eval);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    PROFILE_SCOPE("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("NLopt optimization");
#endif

      real f = HUGE_VAL;
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.degree, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        try_count++; // count this try before breaking, so num_tries reports tries made
        break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original
  restore_from_tolerance(opt, tolerance_snapshot);

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}


inline Result optimize(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  switch (opt->method) {
    case OptimizationMethod::baseline:
      return optimize_baseline_impl(opt, output_steps_ms);
    case OptimizationMethod::with_analytical_pva:
      return optimize_with_analytical_pva_impl(opt, output_steps_ms);
    case OptimizationMethod::with_analytical_dynamics:
      return optimize_with_analytical_dynamics_impl(opt, output_steps_ms);
    case OptimizationMethod::with_segments:
      return optimize_with_segments_impl(opt, output_steps_ms);
  }
  return optimize_with_segments_impl(opt, output_steps_ms); // unreachable
}

} // namespace blast
