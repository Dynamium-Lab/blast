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
  Optimization* opt;
  Array         x;
  Array         x0;
  nlopt_result  nlopt_exit_criteria = NLOPT_SUCCESS;
  int           num_eval            = 0;
  int           num_tries           = 0;
  Trajectory    trajectory;

  Result() = delete;

  explicit Result(Optimization* optim) {
    opt = optim;
  }
};

inline Optimization::Optimization(const Manipulator& new_manip, const Matrix& new_task) :
    manip(new_manip),
    bspline(new_manip.n_joints),
    task(new_task),
    custom_data(nullptr) {
  // Default values
  guess.type   = Guess::random;
  guess.n_shot = 100;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.tcp_speed    = true;

  objective.K_time = 1.0;
}

inline Optimization::Optimization(const Manipulator& new_manip, const Matrix& new_task, const Bspline& new_bspline) :
    manip(new_manip),
    bspline(new_bspline),
    task(new_task),
    custom_data(nullptr) {
  guess.type   = Guess::random;
  guess.n_shot = 100;

  constraints.position     = true;
  constraints.velocity     = true;
  constraints.acceleration = true;
  constraints.tcp_speed    = true;

  objective.K_time = 1.0;
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

inline void Optimization::set_task(Matrix new_task) {
  task = std::move(new_task);
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

  if (opt->constraints.tcp_speed)
    opt->constraints.n_constraints += n_points;

  if (opt->constraints.external_collisions)
    opt->constraints.n_constraints += opt->constraints.n_collision_constraints;
  if (opt->constraints.self_collisions) {
    // if constexpr (std::is_same_v<T_manip, GenericManipulator>) {
    //     // auto col_matrix = opt->manip._collision_matrix;
    //     for (int i = 0; i < opt->manip._collision_matrix.cols; i++) {
    //         for (int j = i+1; j < opt->manip._n_caps; j++) {
    //             opt->manip._n_internal_collisions += opt->manip._collision_matrix(j, i) == 0 ? 0 : 1;
    //         }
    //         opt->manip._n_internal_collisions += opt->manip._collision_base[i] == 0 ? 0 : 1;
    //     }
    // }
    // else {
    opt->constraints.n_constraints += opt->manip._n_internal_collisions * n_points;
    // }
  }

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
      if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
        task(i, 0) = tmp_value;
      } else {
        // Try updating end value
        tmp_value = task(i, 3);
        if (task(i, 3) < task(i, 0)) {
          tmp_value += 2 * PI;
        } else {
          tmp_value -= 2 * PI;
        }
        if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
          task(i, 0) = tmp_value;
        } else {
          break; // Nothing to be done
        }
      }
    }
  }
  opt->task = task;
}

inline Result optimize(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
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

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, 0.001);
  Array x_tol(n, 0.000001);

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
  stop.maxeval    = 100000;
  stop.maxtime    = 5000;
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
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
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
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      is_valid = max(constraints_points) < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      is_valid_more = max(constraints_more_points) < opt->success_tolerance;

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

  opt->guess = start_guess; // reset to original

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

inline void initialize_optimization_with_segments(Optimization* opt) {
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
      if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
        task(i, 0) = tmp_value;
      } else {
        // Try updating end value
        tmp_value = task(i, 3);
        if (task(i, 3) < task(i, 0)) {
          tmp_value += 2 * PI;
        } else {
          tmp_value -= 2 * PI;
        }
        if (tmp_value > opt->manip.pmin[i] && tmp_value < opt->manip.pmax[i]) {
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
  opt->constraints.n_constraints = (opt->manip.n_joints * 4) * ((int) opt->bspline.n_ctrl - (int) opt->bspline.p);
}

inline Result optimize_with_segments(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization_with_segments(opt);
  n_con_with_segments(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, 0.001);
  Array x_tol(n, 0.000001);

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
  stop.maxeval    = 100000;
  stop.maxtime    = 5000;
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
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
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
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints_with_segments(x, *opt, constraints_points);
      is_valid = max(constraints_points) < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con_with_segments(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints_with_segments(x, opt_val_more, constraints_more_points);
      is_valid_more = max(constraints_more_points) < opt->success_tolerance;

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

  opt->guess = start_guess; // reset to original

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

inline Result optimize_dev(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
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

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, 0.001);
  Array x_tol(n, 0.000001);

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
  stop.maxeval    = 100000;
  stop.maxtime    = 5000;
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
  fc.mf     = nlopt_constraints_dev;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_dev, opt, con_tol.data);
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
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
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
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      is_valid = max(constraints_points) < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      is_valid_more = max(constraints_more_points) < opt->success_tolerance;

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

  opt->guess = start_guess; // reset to original

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

inline Result optimize_dev_new(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
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

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, 0.001);
  Array x_tol(n, 0.000001);

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
  stop.maxeval    = 100000;
  stop.maxtime    = 5000;
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
  fc.mf     = nlopt_constraints_dev_new;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_dev_new, opt, con_tol.data);
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
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
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
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      is_valid = max(constraints_points) < opt->success_tolerance;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      is_valid_more = max(constraints_more_points) < opt->success_tolerance;

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

  opt->guess = start_guess; // reset to original

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


} // namespace blast
