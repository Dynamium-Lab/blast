// example_05_two_stage_optimization_selectable.cpp
//
// Same two-stage pipeline as example_04, with a selectable secondary objective:
//   1. Minimize trajectory duration.
//   2. Minimize either the total sampled joint velocity or the total sampled
//      joint acceleration while preserving the stage-1 duration. Two
//      strategies are compared:
//        2a. Duration stays a decision variable, bounded by a
//            max-duration constraint.
//        2b. Duration is removed from the decision vector and fixed
//            to the stage-1 optimum (SmoothnessOptimization).
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_05_two_stage_optimization_selectable
// Run:
//   ./build/examples/example_05_two_stage_optimization_selectable [vel|acc]   (default: vel)

#define BLAST_USE_NATIVE_SQP
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <blast>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

using namespace blast;

// -------------------------------------------------------------------------
// Secondary objective selection.
// -------------------------------------------------------------------------
enum class SecondaryObjective { velocity_sum,
                                acceleration_sum };

static const char* objective_label(SecondaryObjective o) {
  return o == SecondaryObjective::velocity_sum ? "Velocity sum" : "Acceleration sum";
}

static const char* objective_suffix(SecondaryObjective o) {
  return o == SecondaryObjective::velocity_sum ? "vel" : "acc";
}

template<typename M>
static real sum_abs(const M& m) {
  real sum = 0.0;
  for (int i = 0; i < m.cols; ++i)
    for (int j = 0; j < m.rows; ++j)
      sum += std::abs(m(j, i));
  return sum;
}

// -------------------------------------------------------------------------
// Sum of absolute sampled joint velocities.
// Lower values indicate lower aggregate sampled joint speed, but do not
// guarantee smoothness in every sense (e.g., jerk or acceleration continuity).
// -------------------------------------------------------------------------
real objective_velocity_sum(Optimization* opt, Array* x) {
  opt->bspline.compute_trajectory(*x, opt->task);
  return sum_abs(opt->bspline.traj.vel);
}

// -------------------------------------------------------------------------
// Sum of absolute sampled joint accelerations.
// Lower values indicate lower aggregate sampled joint effort, but do not
// bound jerk.
// -------------------------------------------------------------------------
real objective_acceleration_sum(Optimization* opt, Array* x) {
  opt->bspline.compute_trajectory(*x, opt->task);
  return sum_abs(opt->bspline.traj.acc);
}

using ObjectiveFn = real (*)(Optimization*, Array*);

static ObjectiveFn objective_fn(SecondaryObjective o) {
  return o == SecondaryObjective::velocity_sum ? objective_velocity_sum : objective_acceleration_sum;
}

// -------------------------------------------------------------------------
// Stage-2a constraint: trajectory duration <= maximum duration.
// Constraint convention: returned value must be <= 0.
// -------------------------------------------------------------------------
Array constraint_max_duration(Optimization* opt, Array* x) {
  const auto* max_duration = static_cast<const real*>(opt->custom_data);
  return Array(1, x->back() - *max_duration);
}

// -------------------------------------------------------------------------
// Stage-2b: optimization with the duration fixed. The solver sees a reduced
// decision vector (n - 1 entries, time removed); the wrappers below expand it
// back to the full vector expected by the library.
// -------------------------------------------------------------------------
struct SmoothnessOptimization : Optimization {
  real total_time;

  SmoothnessOptimization(const Optimization& base, real new_total_time) : Optimization(base), total_time(new_total_time) {}

  SmoothnessOptimization(const Manipulator& new_manip, const Task& new_task, real new_total_time) : Optimization(new_manip, new_task), total_time(new_total_time) {}

  SmoothnessOptimization(const Manipulator& new_manip, const Task& new_task, const Bspline& new_bspline, real new_total_time) : Optimization(new_manip, new_task, new_bspline), total_time(new_total_time) {}

  void set_total_time(real new_total_time) { total_time = new_total_time; }
};

inline void expand_x(const SmoothnessOptimization* opt, unsigned n_red, const real* x_red, Array& x_full) {
  std::memcpy(x_full.data, x_red, n_red * sizeof(real));
  x_full[n_red] = opt->total_time;
}

inline blast_fn real nlopt_objective_smoothness(unsigned n, const real* x, real* grad, void* f_data) {
  auto* opt = (SmoothnessOptimization*) f_data;

  Array xf(n + 1);
  expand_x(opt, n, x, xf);

  if (!grad)
    return objective_function(n + 1, xf.data, nullptr, opt);

  Array gf(n + 1);
  real  f = objective_function(n + 1, xf.data, gf.data, opt);
  std::memcpy(grad, gf.data, n * sizeof(real)); // drop d/d(time)
  return f;
}

inline blast_fn void nlopt_constraints_smoothness(unsigned m, real* result, unsigned n, const real* x, real* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  PROFILE_FUNCTION;
#endif
  auto* opt = (SmoothnessOptimization*) f_data;

  Array xf(n + 1);
  expand_x(opt, n, x, xf);

  if (!grad) {
    nlopt_constraints(m, result, n + 1, xf.data, nullptr, opt);
    return;
  }

  // Full Jacobian is m x (n + 1); keep the first n columns of each row.
  Array gf((int) (m * (n + 1)));
  nlopt_constraints(m, result, n + 1, xf.data, gf.data, opt);
  for (unsigned i = 0; i < m; i++)
    std::memcpy(grad + i * n, gf.data + i * (n + 1), n * sizeof(real));
}

inline Result optimize_for_smoothness(SmoothnessOptimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) {
    print(opt->task);
    return result;
  }

  auto tolerance_snapshot = tighten_for_success_tolerance(opt); // restored below, next to opt->guess restore

  const auto n     = opt->bspline.x_len(opt->task);             // full length, time is last
  const auto n_red = n - 1;                                     // what the solver sees

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n_red, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n_red;
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

  Array ub(n_red, INF_REAL);
  Array lb(n_red, -INF_REAL);

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_smoothness;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n_red);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_smoothness, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, nlopt_objective_smoothness, opt);
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
  for (; try_count < opt->max_tries; try_count++) {
#if BLAST_TRACE_LEVEL >= 1
    PROFILE_SCOPE("Optimization");
#endif

    // initial guess (full vector), then strip the time entry
    Array x;
    Array x_red(n_red);
    {
#if BLAST_TRACE_LEVEL >= 1
      PROFILE_SCOPE("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
      std::memcpy(x_red.data, x.data, n_red * sizeof(real));
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
              n_red,
              nlopt_objective_smoothness,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x_red.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;
#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x_red.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // rebuild the full decision vector for validation and output
    std::memcpy(x.data, x_red.data, n_red * sizeof(real));
    x.back() = opt->total_time;

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
        opt->guess.type = Guess::random; // random time entry is discarded; total_time stays fixed
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

static void report(const char* name, const Result& r) {
  std::cout << "  " << name
            << "\n    Success: " << (r.success ? "yes" : "no")
            << "\n    Compute time (ms): " << r.compute_time
            << "\n    Evaluations: " << r.num_eval << '\n';
}

TEST_CASE("Smoothness optimization test + benchmark", "[Optimization]") {
  // -----------------------------------------------------------------------
  // Step 0 — Select the secondary objective.
  // -----------------------------------------------------------------------
  SecondaryObjective secondary = SecondaryObjective::acceleration_sum;

  const ObjectiveFn secondary_fn    = objective_fn(secondary);
  const char*       secondary_label = objective_label(secondary);

  // -----------------------------------------------------------------------
  // Step 1 — Define the robot, task, and B-spline representation.
  // -----------------------------------------------------------------------
  Manipulator ur5e   = make_UR5e();
  ur5e.base_position = {-0.5, -0.3, 0.35};
  ur5e.base_rotation = {-1, 0, 0, 0, -1, 0, 0, 0, 1};

  Array start = {-0.704868, 2.17507, 0.834305, 0.132228, 2.27799, 0.00328674};
  Array end   = {3.12394, -0.553571, -0.509807, 1.06314, 1.52902, -3.10058};
  Task  task  = Task::stop_to_stop(start, end);

  auto world = get_kitchen_no_doors();

  Optimization opt(ur5e, task);
  opt.method  = OptimizationMethod::with_segments;
  opt.bspline = Bspline(12, 100, 5, ur5e.n_joints);

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.external_collisions = true;
  opt.set_world(world);
  opt.success_tolerance = 0.01;

  // Using baseline method to have sensible computation time comparisons
  opt.method = OptimizationMethod::baseline;

  // -----------------------------------------------------------------------
  // Step 2 — Primary optimization: minimize trajectory duration.
  // -----------------------------------------------------------------------
  opt.objective.time_weight = 1.0;

  std::cout << "Stage 1: minimizing trajectory duration...\n";

  Result time_result = optimize(&opt);
  report("Time-optimal", time_result);

  if (!time_result.success || time_result.x.is_empty()) {
    std::cerr << "Stage 1 failed; cannot proceed to stage 2.\n";
  }

  const real optimal_duration = time_result.x.back();
  std::cout << "    Optimal duration (s): " << optimal_duration << '\n';

  // -----------------------------------------------------------------------
  // Step 3 — Shared configuration for both refinement strategies.
  //
  // Start from the stage-1 solution, drop the duration objective, and
  // minimize the selected secondary objective. Each strategy gets its own
  // copy so the two runs are fully independent.
  // -----------------------------------------------------------------------
  Optimization smooth_base    = opt;
  smooth_base.guess.type      = Guess::custom;
  smooth_base.guess.initial_x = time_result.x;

  smooth_base.objective.time_weight = 0.0;
  smooth_base.objective.add_custom_objective(secondary_fn, 1.0);

  // 2a: duration free but bounded by the stage-1 optimum.
  // The referenced value must remain alive throughout optimize().
  Optimization opt_bounded  = smooth_base;
  const real   max_duration = optimal_duration;
  opt_bounded.custom_data   = const_cast<real*>(&max_duration);
  opt_bounded.constraints.add_constraint(constraint_max_duration, 1);

  // 2b: duration removed from the decision vector, fixed to the optimum.
  SmoothnessOptimization opt_fixed(smooth_base, optimal_duration);

  // -----------------------------------------------------------------------
  // Step 4 — Secondary optimizations.
  // -----------------------------------------------------------------------
  std::cout << "\nStage 2: minimizing " << secondary_label << " at the primary-stage duration...\n";

  Result bounded_result = optimize(&opt_bounded);
  report("Bounded duration (constraint)", bounded_result);

  Result fixed_result = optimize_for_smoothness(&opt_fixed);
  report("Fixed duration (reduced vector)", fixed_result);

  // -----------------------------------------------------------------------
  // Step 5 — Compare the three solutions on the selected objective.
  // Evaluate on a scratch copy so the optimizers' bsplines are untouched.
  // -----------------------------------------------------------------------
  Optimization evaluator = smooth_base;

  const real obj_time    = secondary_fn(&evaluator, &time_result.x);
  const auto improvement = [&](real v) { return obj_time > 0.0 ? 100.0 * (obj_time - v) / obj_time : 0.0; };

  std::cout << "\n--- Pipeline results (" << secondary_label << ") ---\n"
            << "Time-optimal:\n"
            << "  Duration (s):      " << time_result.x.back() << '\n'
            << "  " << secondary_label << ": " << obj_time << '\n';

  const auto print_refined = [&](const char* name, Result& r) {
    std::cout << name << ":\n";
    if (!r.success || r.x.is_empty()) {
      std::cout << "  failed\n";
      return;
    }
    const real v = secondary_fn(&evaluator, &r.x);
    std::cout << "  Duration (s):      " << r.x.back() << '\n'
              << "  " << secondary_label << ": " << v << '\n'
              << "  Improvement (%):   " << improvement(v) << '\n'
              << "  Compute time (ms): " << r.compute_time << '\n'
              << "  Evaluations:       " << r.num_eval << '\n';
    CHECK(improvement(v) >= 0);
  };

  print_refined("Bounded duration", bounded_result);
  print_refined("Fixed duration", fixed_result);
}
