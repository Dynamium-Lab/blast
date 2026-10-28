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

  opt->lower_bounds        = Array(opt->x_len(), -HUGE_VAL);
  opt->upper_bounds        = Array(opt->x_len(), -HUGE_VAL);
  opt->lower_bounds.back() = opt->min_duration;
  opt->upper_bounds.back() = opt->max_duration;

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

// Fills lower_bounds / upper_bounds (sized x_len) for the decision vector
// [free control points, joint-major, (n_ctrl - 6) per joint | T]. Call after tighten_for_success_tolerance().
// T is bounded to [min_duration, max_duration]. For stop-to-stop tasks the control points are also
// bounded to the (tightened) position limits, widened by control_point_bound_margin of each joint's
// range; set bound_control_points = false to skip that. Other task layouts get the duration bounds only.
inline void decision_bounds(const Optimization* opt, Array& lower_bounds, Array& upper_bounds) {
  for (u32 variable = 0; variable < lower_bounds.size; variable++) {
    lower_bounds[variable] = -INF_REAL;
    upper_bounds[variable] = INF_REAL;
  }
  lower_bounds.back()             = opt->min_duration;
  upper_bounds.back()             = opt->max_duration;
  const int n_joints              = opt->manip.n_joints;
  const int n_free_control_points = (int) opt->bspline.n_ctrl - 6;
  // Non-stop-to-stop layouts have a different variable count; they get duration bounds only.
  if ((int) lower_bounds.size != n_joints * n_free_control_points + 1)
    return;
  if (!opt->bound_control_points)
    return;
  // Fallback for joints with no position limits: ±2π (screwing tasks aside, a stop-to-stop
  // trajectory never needs more than a full revolution per joint). margin is skipped for
  // unbounded joints to avoid multiplying by an infinite range.
  for (int joint = 0; joint < n_joints; joint++) {
    const real lo           = std::isfinite(opt->manip.position_min[joint]) ? opt->manip.position_min[joint] : -2 * PI;
    const real hi           = std::isfinite(opt->manip.position_max[joint]) ? opt->manip.position_max[joint] : 2 * PI;
    const bool finite_range = std::isfinite(opt->manip.position_max[joint] - opt->manip.position_min[joint]);
    const real margin       = finite_range ? opt->control_point_bound_margin * (hi - lo) : real(0);
    for (int control_point = 0; control_point < n_free_control_points; control_point++) {
      lower_bounds[joint * n_free_control_points + control_point] = lo - margin;
      upper_bounds[joint * n_free_control_points + control_point] = hi + margin;
    }
  }
}

inline void clamp_to_bounds(Array& x, const Array& lower_bounds, const Array& upper_bounds) {
  for (u32 variable = 0; variable < x.size; variable++)
    x[variable] = std::clamp(x[variable], lower_bounds[variable], upper_bounds[variable]);
}

/**
 * @brief Sets the clearance targets of the first/last segment's collision constraints to what
 *        the pinned start/goal configuration achieves, capped at `buffer`.
 *
 * Measured per capsule (external) and for self-collision, on the tightened geometry. Static
 * obstacles only. Why the endpoints need a target: docs/guides/constraints.md.
 */
inline void set_pinned_col_clearances(Optimization* opt, real buffer) {
  // Floor, as a fraction of the buffer: keeps tolerance / clearance finite for an endpoint at ~0 clearance.
  constexpr real min_target_fraction = 0.01;

  if (!opt->constraints.external_collisions && !opt->constraints.self_collisions) {
    opt->pinned_self_clearance.fill(buffer);
    for (auto& clearances: opt->pinned_col_clearance)
      clearances.fill(buffer);
    return;
  }
  const real min_target = buffer * min_target_fraction;
  // Measured in the buffered geometry, so the endpoint's planning clearance is distance + buffer.
  const auto endpoint_target = [&](real buffered_distance) { return std::clamp(buffered_distance + buffer, min_target, buffer); };

  ManipulatorTempData manip_data;
  Array               endpoint_position(opt->manip.n_joints);
  for (int endpoint_index = 0; endpoint_index < 2; endpoint_index++) {
    for (int joint = 0; joint < opt->manip.n_joints; joint++)
      endpoint_position[joint] = opt->task(joint, endpoint_index == 0 ? 0 : 3);
    forward_kinematics(opt->manip, manip_data, endpoint_position);
    compute_collision_model(opt->manip, manip_data);

    opt->pinned_self_clearance[endpoint_index] = opt->constraints.self_collisions
                                                         ? endpoint_target(min(get_internal_collisions(opt->manip, manip_data)))
                                                         : buffer;
    for (int capsule_id = 0; capsule_id < opt->manip._n_caps; capsule_id++) {
      if (!opt->constraints.external_collisions) {
        opt->pinned_col_clearance[endpoint_index][capsule_id] = buffer;
        continue;
      }
      const auto& capsule          = manip_data.capsule_list[capsule_id];
      real        closest_distance = INF_REAL;
      for (const auto& box: opt->world.boxes)
        closest_distance = std::min(closest_distance, distance(capsule, box));
      for (const auto& world_capsule: opt->world.capsules)
        closest_distance = std::min(closest_distance, distance(capsule, world_capsule));
      for (const auto& sphere: opt->world.spheres)
        closest_distance = std::min(closest_distance, distance(capsule, sphere));
      opt->pinned_col_clearance[endpoint_index][capsule_id] = endpoint_target(closest_distance);
    }
  }
}

/**
 * @brief Tightens the limits and collision geometry by success_tolerance, so that a solve accepted
 *        at success_tolerance (the solver's con_tol) satisfies the true limits.
 *
 * Limits shrink by a factor 1 + success_tolerance about their center; collision geometry grows by
 * collision_buffer. Where the pinned start/goal is closer to a limit than that, the first/last
 * segment's constraints target the endpoint instead (see position_constraint_params,
 * set_pinned_col_clearances).
 * Must be paired with restore_from_tolerance() before the Optimization is returned or reused.
 *
 * @return The true limits and geometry, for restore_from_tolerance().
 */
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
    opt->collision_scale          = 1.0;
    opt->pinned_tightening_active = false;
    return snap;
  }

  const real tol               = opt->success_tolerance;
  const real tightening_factor = 1 + tol; // divide each limit's half-range by this to shrink it by tol

  // The unit change that makes the dimensionless success_tolerance land on the buffer.
  // A non-positive buffer means no unit change.
  const real buffer    = opt->collision_buffer > 0 ? opt->collision_buffer : tol;
  opt->collision_scale = tol / buffer;

  // An unbounded joint stays unbounded.
  for (int joint = 0; joint < opt->manip.n_joints; joint++) {
    if (opt->manip.position_max[joint] == INF_REAL || opt->manip.position_min[joint] == -INF_REAL)
      continue;
    const real center              = (opt->manip.position_max[joint] + opt->manip.position_min[joint]) / 2;
    const real half_range          = (opt->manip.position_max[joint] - opt->manip.position_min[joint]) / 2 / tightening_factor;
    opt->manip.position_max[joint] = center + half_range;
    opt->manip.position_min[joint] = center - half_range;
  }
  for (int j = 0; j < opt->manip.n_joints; j++) {
    opt->manip.velocity_max[j] /= tightening_factor;
    opt->manip.acceleration_max[j] /= tightening_factor;
    opt->manip.torque_max[j] /= tightening_factor;
  }
  opt->manip.tool_speed_max /= tightening_factor;

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

  set_pinned_col_clearances(opt, buffer);
  opt->pinned_tightening_active = true;

  return snap;
}

/** @brief Undoes tighten_for_success_tolerance(): puts back the true limits and geometry. */
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
  opt->world                    = snap.world;
  opt->collision_scale          = 1.0; // raw metres again outside the tightening window
  opt->pinned_tightening_active = false;
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

  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);

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
  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);
  nlopt_res = nlopt_set_lower_bounds(o, lower_bounds.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bounds(o, upper_bounds.data);
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
      x = init_guess(opt);
      clamp_to_bounds(x, lower_bounds, upper_bounds); // NLopt refuses a start outside the bounds
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
              lower_bounds.data,
              upper_bounds.data,
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
      restore_from_tolerance(&opt_val_more, tolerance_snapshot); // validate on true geometry, not tightened
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < 0;

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

  opt->guess = start_guess;                        // reset to original
  restore_from_tolerance(opt, tolerance_snapshot); // undo tighten_for_success_tolerance(): caller keeps the real limits

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

  opt->lower_bounds        = Array(opt->x_len(), -HUGE_VAL);
  opt->upper_bounds        = Array(opt->x_len(), HUGE_VAL);
  opt->lower_bounds.back() = opt->min_duration;
  opt->upper_bounds.back() = opt->max_duration;

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

  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);

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
  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);
  nlopt_res = nlopt_set_lower_bounds(o, lower_bounds.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bounds(o, upper_bounds.data);
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
      x = init_guess_segments(opt);
      clamp_to_bounds(x, lower_bounds, upper_bounds); // NLopt refuses a start outside the bounds
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
              lower_bounds.data,
              upper_bounds.data,
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
      restore_from_tolerance(&opt_val_more, tolerance_snapshot); // validate on true geometry, not tightened
      opt_val_more.set_bspline(bspline_val_more);
      n_con_with_segments(&opt_val_more);
      Array  constraints_more_points(opt_val_more.constraints.n_constraints);
      Matrix gradient;
      constraints_and_gradients_with_segments(x, opt_val_more, constraints_more_points, gradient);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < 0;

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

  opt->guess = start_guess;                        // reset to original
  restore_from_tolerance(opt, tolerance_snapshot); // undo tighten_for_success_tolerance(): caller keeps the real limits

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

inline Result optimize_with_broadphase_impl(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
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

  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_with_broadphase;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_with_broadphase, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);
  nlopt_res = nlopt_set_lower_bounds(o, lower_bounds.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bounds(o, upper_bounds.data);
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
      x = init_guess_segments(opt);
      clamp_to_bounds(x, lower_bounds, upper_bounds); // NLopt refuses a start outside the bounds
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
              lower_bounds.data,
              upper_bounds.data,
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
      constraints_and_gradients_with_broadphase(x, *opt, constraints_points, gradient);
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance * 2;
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
      restore_from_tolerance(&opt_val_more, tolerance_snapshot); // validate on true geometry, not tightened
      opt_val_more.set_bspline(bspline_val_more);
      n_con_with_segments(&opt_val_more);
      Array  constraints_more_points(opt_val_more.constraints.n_constraints);
      Matrix gradient;
      constraints_and_gradients_with_broadphase(x, opt_val_more, constraints_more_points, gradient);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < 0;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess;                        // reset to original
  restore_from_tolerance(opt, tolerance_snapshot); // undo tighten_for_success_tolerance(): caller keeps the real limits

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

inline Result optimize_with_double_broadphase_impl(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
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

  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_with_double_broadphase;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_with_double_broadphase, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);
  nlopt_res = nlopt_set_lower_bounds(o, lower_bounds.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bounds(o, upper_bounds.data);
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
      x = init_guess_segments(opt);
      clamp_to_bounds(x, lower_bounds, upper_bounds); // NLopt refuses a start outside the bounds
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
              lower_bounds.data,
              upper_bounds.data,
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
      constraints_and_gradients_with_double_broadphase(x, *opt, constraints_points, gradient);
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance * 2;
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
      restore_from_tolerance(&opt_val_more, tolerance_snapshot); // validate on true geometry, not tightened
      opt_val_more.set_bspline(bspline_val_more);
      n_con_with_segments(&opt_val_more);
      Array  constraints_more_points(opt_val_more.constraints.n_constraints);
      Matrix gradient;
      constraints_and_gradients_with_double_broadphase(x, opt_val_more, constraints_more_points, gradient);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < 0;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess;                        // reset to original
  restore_from_tolerance(opt, tolerance_snapshot); // undo tighten_for_success_tolerance(): caller keeps the real limits

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

  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);

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
  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);
  nlopt_res = nlopt_set_lower_bounds(o, lower_bounds.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bounds(o, upper_bounds.data);
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
      x = init_guess(opt);
      clamp_to_bounds(x, lower_bounds, upper_bounds); // NLopt refuses a start outside the bounds
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
              lower_bounds.data,
              upper_bounds.data,
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
      restore_from_tolerance(&opt_val_more, tolerance_snapshot); // validate on true geometry, not tightened
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < 0;

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

  opt->guess = start_guess;                        // reset to original
  restore_from_tolerance(opt, tolerance_snapshot); // undo tighten_for_success_tolerance(): caller keeps the real limits

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

  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);

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
  Array lower_bounds(n), upper_bounds(n);
  decision_bounds(opt, lower_bounds, upper_bounds);
  nlopt_res = nlopt_set_lower_bounds(o, lower_bounds.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bounds(o, upper_bounds.data);
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
      x = init_guess(opt);
      clamp_to_bounds(x, lower_bounds, upper_bounds); // NLopt refuses a start outside the bounds
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
              lower_bounds.data,
              upper_bounds.data,
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
      restore_from_tolerance(&opt_val_more, tolerance_snapshot); // validate on true geometry, not tightened
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < 0;

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

  opt->guess = start_guess;                        // reset to original
  restore_from_tolerance(opt, tolerance_snapshot); // undo tighten_for_success_tolerance(): caller keeps the real limits

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
    case OptimizationMethod::broadphase:
      return optimize_with_broadphase_impl(opt, output_steps_ms);
    case OptimizationMethod::double_broadphase:
      return optimize_with_double_broadphase_impl(opt, output_steps_ms);
  }
  return optimize_with_segments_impl(opt, output_steps_ms); // unreachable
}

} // namespace blast
