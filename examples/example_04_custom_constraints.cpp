// example_04_two_stage_optimization.cpp
//
// This example shows how to:
//   1. Chain two optimizations: first minimize the trajectory duration, then
//      refine the result with a secondary objective at (about) that duration
//   2. Write a custom objective function with an analytical gradient and
//      register it with objective.add_custom_objective_with_gradient()
//   3. Keep the duration of stage 2 at the stage-1 optimum in two ways, and
//      compare them:
//        A. pin it: min_duration = max_duration, so the solver cannot move T
//        B. cap it with a custom global constraint, registered with
//           constraints.add_global_constraint(); this also introduces the
//           custom constraint framework (global and per-point constraints)
//   4. Pass user data to custom callbacks through Optimization::custom_data
//   5. Warm-start the second stage from the first-stage solution
//      (Guess::custom)
//
// The secondary objective is selectable on the command line:
//   vel : minimize the sum of absolute sampled joint velocities
//   acc : minimize the sum of absolute sampled joint accelerations
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_04_custom_constraints
// Run:
//   ./build/examples/example_04_custom_constraints [vel|acc]   (default: acc)

#ifndef BLAST_USE_NATIVE_SQP
#define BLAST_USE_NATIVE_SQP // use the built-in SQP solver
#endif

#include <blast>
#include <cmath>
#include <iostream>
#include <string>

using namespace blast;

// -------------------------------------------------------------------------
// Custom objectives.
//
// A custom objective has the signature real(Optimization*, Array*):
//   opt : the optimization problem being solved (robot, task, B-spline, ...)
//   x   : the current decision vector — the B-spline control points,
//         followed by the trajectory duration as the last entry (x->back())
// and is registered with objective.add_custom_objective(); the solver then
// takes its gradient by finite differences, recomputing the objective once
// per decision variable.
//
// The objectives below provide their gradient instead: the signature
// real(Optimization*, Array* x, Array* grad) returns the value and, when
// grad is not null, writes df/dx into it. They are registered with
// objective.add_custom_objective_with_gradient(), and are preferable to
// finite differences, when available.
//
// compute_trajectory() evaluates the B-spline defined by x and stores the
// sampled positions, velocities and accelerations in opt->bspline.traj,
// one column per sample and one row per joint. With the control points c
// and the duration T, at sample i:
//   vel(j, i) = sum_k c(k, j) * basis_v(k, i) / T
//   acc(j, i) = sum_k c(k, j) * basis_a(k, i) / T^2
// The decision vector holds the free control points k = 3 .. n_ctrl - 4 of
// each joint, joint by joint, then T. In a stop-to-stop task the pinned
// control points do not depend on T, so a sum of |vel| scales as 1/T and a
// sum of |acc| as 1/T^2.
// -------------------------------------------------------------------------
template<typename M>
static real sum_abs(const M& m) {
  real sum = 0;
  for (int i = 0; i < m.cols; ++i)
    for (int j = 0; j < m.rows; ++j)
      sum += std::abs(m(j, i));
  return sum;
}

// Gradient of sum_ij |d(j, i)| where d(j, i) = sum_k c(k, j) * basis(k, i) / T^power.
// d(T) for a stop-to-stop task: d/dT of the sum is -power * sum / T.
template<typename M>
static void sum_abs_gradient(const Optimization& opt, const M& d, const Matrix& basis, int power, real sum, real T, Array& grad) {
  const int  n_free = (int) opt.bspline.n_ctrl - 6; // free control points per joint
  const real scale  = (real) 1 / std::pow(T, (real) power);
  for (int i = 0; i < d.cols; ++i)
    for (int j = 0; j < d.rows; ++j) {
      const real s = sign(d(j, i));
      if (s == 0)
        continue;
      for (int k = 3; k < (int) opt.bspline.n_ctrl - 3; ++k)
        grad[j * n_free + k - 3] += s * basis(k, i) * scale;
    }
  grad.back() = -(real) power * sum / T;
}

// Sum of absolute sampled joint velocities. Lower values mean lower aggregate
// joint speed; this says nothing about acceleration or jerk.
real objective_velocity_sum(Optimization* opt, Array* x, Array* grad) {
  opt->bspline.compute_trajectory(*x, opt->task);
  const real sum = sum_abs(opt->bspline.traj.vel);
  if (grad)
    sum_abs_gradient(*opt, opt->bspline.traj.vel, opt->bspline.basis_v, 1, sum, x->back(), *grad);
  return sum;
}

// Sum of absolute sampled joint accelerations. Lower values mean lower
// aggregate joint acceleration; this does not bound jerk.
real objective_acceleration_sum(Optimization* opt, Array* x, Array* grad) {
  opt->bspline.compute_trajectory(*x, opt->task);
  const real sum = sum_abs(opt->bspline.traj.acc);
  if (grad)
    sum_abs_gradient(*opt, opt->bspline.traj.acc, opt->bspline.basis_a, 2, sum, x->back(), *grad);
  return sum;
}

using ObjectiveFn = real (*)(Optimization*, Array*, Array*);

// -------------------------------------------------------------------------
// Custom constraints.
//
// Besides the built-in constraints (position, velocity, collisions, ...),
// every OptimizationMethod accepts user constraints. A custom constraint
// returns an Array of values, and each value must be <= 0 for the
// constraint to be satisfied. A custom constraint is registered together
// with the number of values it returns. There are two kinds:
//
//   - Global: Array(Optimization* opt, Array* x), evaluated once per
//     trajectory. Use it for properties of the whole trajectory (duration,
//     an integral, ...). It is called at the end of the constraint function,
//     so the trajectory can be queried at opt->bspline.trajectory.
//     Registered with constraints.add_global_constraint(fn, n).
//
//   - Per point: Array(Optimization* opt, Array* x, ManipulatorTempData* state,
//     u32 point), evaluated at every B-spline sample. Use it for properties
//     that must hold all along the trajectory. `state` holds the sample's
//     forward kinematics (state->p_j: joint positions in space, the last one
//     being the tool), collision model and joint efforts; the sample's joint
//     positions, velocities and accelerations are
//     opt->bspline.traj.pos/vel/acc.col(point). It does not require recompute
//     of the trajectory.
//     Registered with constraints.add_point_constraint(fn, n).

//     For example, to keep the tool above 0.2 m:
//
//       Array tool_above(Optimization* opt, Array*, ManipulatorTempData* s, u32) {
//         return Array(1, (real) 0.2 - s->p_j[opt->manip.n_joints].z);
//       }
//       opt.constraints.add_point_constraint(tool_above, 1);
//
//     The segment-based methods (with_segments, broadphase,
//     double_broadphase) keep the worst sample of each B-spline segment,
//     value by value; the point-based ones keep every sample. This is done
//     on blast's side and does not need to be accounted for in your custom
//     constraints declarations.
//
// Scaling matters: like the built-in limits, custom constraints are
// tightened by success_tolerance during the solve (each value is shifted up
// by it), so that a successful solve has every value <= 0. Write the
// values on a scale where success_tolerance (1% by default) is a sensible
// margin, e.g. divided by the limit as below. If one wishes to scale a custom
// constraint in the same way blast's constraints are, choose to scale to
// obtain the constraint value cannot go below -1.
//
// Here: the trajectory duration T must not exceed a cap. The duration is a
// property of the whole trajectory, so the constraint is global. The cap is
// not part of the problem definition, so it is passed in through
// opt->custom_data (a user-owned pointer, see Step 4B).
// -------------------------------------------------------------------------
Array constraint_max_duration(Optimization* opt, Array* x) {
  const auto* max_duration = static_cast<const real*>(opt->custom_data);
  return Array(1, (x->back() - *max_duration) / *max_duration); // relative: <= 0 when T <= cap
}

struct StageResult {
  std::string label;
  Result      result;
};

static void print_stage(const std::string& label, const Result& result) {
  std::cout << label << "\n"
            << "  Success:           " << (result.success ? "yes" : "no") << "\n"
            << "  Compute time (ms): " << result.compute_time << "\n"
            << "  Evaluations:       " << result.num_eval << "\n";
}

int main(int argc, char** argv) {
  // -----------------------------------------------------------------------
  // Step 0 — Select the secondary objective from the command line.
  // -----------------------------------------------------------------------
  ObjectiveFn secondary_fn    = objective_acceleration_sum;
  std::string secondary_label = "Acceleration sum";
  std::string secondary_tag   = "acc";

  if (argc > 1) {
    const std::string arg = argv[1];
    if (arg == "vel") {
      secondary_fn    = objective_velocity_sum;
      secondary_label = "Velocity sum";
      secondary_tag   = "vel";
    } else if (arg != "acc") {
      std::cerr << "Usage: " << argv[0] << " [vel|acc]\n";
      return 2;
    }
  }

  // -----------------------------------------------------------------------
  // Step 1 — Robot, task, and collision world.
  // The UR5e is mounted with a custom base pose inside a kitchen scene.
  // -----------------------------------------------------------------------
  Manipulator ur5e   = make_UR5e();
  ur5e.base_position = {-(real) 0.5, -(real) 0.3, (real) 0.35};
  ur5e.base_rotation = {-1, 0, 0, 0, -1, 0, 0, 0, 1}; // 180° about z

  Array start = {-(real) 0.704868, (real) 2.17507, (real) 0.834305, (real) 0.132228, (real) 2.27799, (real) 0.00328674};
  Array end   = {(real) 3.12394, -(real) 0.553571, -(real) 0.509807, (real) 1.06314, (real) 1.52902, -(real) 3.10058};
  Task  task  = Task::stop_to_stop(start, end);

  World world = get_kitchen_no_doors();

  // -----------------------------------------------------------------------
  // Step 2 — Create the Optimization problem.
  // Bspline(n_ctrl, n_points, degree, n_joints): 12 control points, 105
  // sample points, degree 5. The sample points are also where the custom
  // objectives above are evaluated.
  // -----------------------------------------------------------------------
  Optimization opt(ur5e, task);
  opt.bspline = Bspline(12, 105, 5, ur5e.n_joints);
  opt.set_world(world);

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.external_collisions = true;

  opt.success_tolerance = (real) 0.01;

  // -----------------------------------------------------------------------
  // Step 3 — Stage 1: minimize trajectory duration.
  // time_weight is the weight of the built-in duration objective.
  // The initial guess is the straight line from start to end in joint
  // space: deterministic, and it converges for this task. (The default,
  // Guess::random, gives a different start, and result, on every run.)
  // -----------------------------------------------------------------------
  opt.objective.time_weight = 1;
  opt.guess.type            = Guess::straight_line;

  std::cout << "Stage 1: minimizing trajectory duration...\n";
  Result time_result = optimize(&opt);
  print_stage("Stage 1 (time-optimal)", time_result);

  if (!time_result.success || time_result.x.is_empty()) {
    std::cerr << "Stage 1 failed; cannot proceed to stage 2.\n";
    return 1;
  }

  const real optimal_duration = time_result.x.back();
  std::cout << "  Duration (s):      " << optimal_duration << "\n";

  // -----------------------------------------------------------------------
  // Step 4 — Stage 2 setup, shared by both variants: refine the trajectory.
  //
  //   - Warm start: Guess::custom uses initial_x as the starting point, so
  //     the optimizer begins from the time-optimal solution, which is
  //     already feasible.
  //   - Objective: disable the duration term and register the selected
  //     custom objective with weight 1.0.
  //
  // Without the duration term, nothing stops the solver from slowing the
  // trajectory down to reduce the secondary objective, so stage 2 must
  // keep the duration in check. Two ways follow.
  // -----------------------------------------------------------------------
  Optimization refined = opt; // same method as stage 1 (with_segments by default)

  refined.guess.type      = Guess::custom;
  refined.guess.initial_x = time_result.x;

  refined.objective.time_weight = 0;
  refined.objective.add_custom_objective_with_gradient(secondary_fn, 1);

  // -----------------------------------------------------------------------
  // Step 4A — Pin the duration.
  //
  // T is the last decision variable, bounded by [min_duration,
  // max_duration] (seconds). Setting both bounds to the stage-1 optimum
  // pins T exactly: the solver only moves the control points. This is the
  // simplest way to refine a trajectory at a known duration, and needs no
  // custom constraint.
  // -----------------------------------------------------------------------
  Optimization pinned = refined;
  pinned.min_duration = optimal_duration;
  pinned.max_duration = optimal_duration;

  // -----------------------------------------------------------------------
  // Step 4B — Cap the duration with a custom global constraint.
  //
  // T stays free, but constraint_max_duration keeps it below a cap. The
  // cap is the stage-1 optimum plus 2% of slack: the stage-1 duration is
  // already the fastest feasible one, and the tightening (see the
  // constraint above) asks for T to be 1% below the cap, so a cap exactly
  // at the optimum would leave no feasible trajectory. The slack lets the
  // solver trade a little duration for a smoother trajectory.
  //
  // custom_data must point to the cap for as long as the solve runs. It is
  // a non-const pointer, hence the non-const variable.
  // -----------------------------------------------------------------------
  real duration_cap = optimal_duration * (1 + 2 * opt.success_tolerance);

  Optimization capped = refined;
  capped.custom_data  = &duration_cap;
  capped.constraints.add_global_constraint(constraint_max_duration, 1); // 1: the function returns one value

  // -----------------------------------------------------------------------
  // Step 5 — Stage 2: minimize the secondary objective, both variants.
  // -----------------------------------------------------------------------
  std::cout << "\nStage 2: minimizing " << secondary_label << "...\n";
  const StageResult stages[] = {
          {"Stage 2A (pinned duration)", optimize(&pinned)},
          {"Stage 2B (duration cap, custom constraint)", optimize(&capped)},
  };
  for (const auto& stage: stages)
    print_stage(stage.label, stage.result);

  // -----------------------------------------------------------------------
  // Step 6 — Compare the solutions on the secondary objective.
  // The objective is evaluated on a scratch copy of the problem so that the
  // B-splines stored in the optimizations are left untouched.
  // -----------------------------------------------------------------------
  Optimization evaluator = refined;
  const real   obj_time  = secondary_fn(&evaluator, &time_result.x, nullptr);

  std::cout << "\n--- Results (" << secondary_label << ") ---\n"
            << "Time-optimal:\n"
            << "  Duration (s):      " << time_result.x.back() << "\n"
            << "  " << secondary_label << ": " << obj_time << "\n";

  bool all_succeeded = true;
  for (const auto& stage: stages) {
    std::cout << stage.label << ":\n";
    if (stage.result.success && !stage.result.x.is_empty()) {
      Array      x           = stage.result.x;
      const real obj_refined = secondary_fn(&evaluator, &x, nullptr);
      const real improvement = obj_time > 0 ? 100 * (obj_time - obj_refined) / obj_time : 0;

      std::cout << "  Duration (s):      " << stage.result.x.back() << "\n"
                << "  " << secondary_label << ": " << obj_refined << "\n"
                << "  Improvement (%):   " << improvement << "\n";
    } else {
      std::cout << "  failed\n";
      all_succeeded = false;
    }
  }

  // -----------------------------------------------------------------------
  // Step 7 — Export the trajectories to CSV for plotting (optional).
  // -----------------------------------------------------------------------
  // const std::string dir = "../../../../examples/";

  // print_to_csv(time_result.trajectory, (dir + "time_traj.csv").c_str());
  // for (int i = 0; i < 2; i++)
  //   if (stages[i].result.success)
  //     print_to_csv(stages[i].result.trajectory, (dir + "refined_" + secondary_tag + (i == 0 ? "_pinned" : "_capped") + "_traj.csv").c_str());

  return all_succeeded ? 0 : 1;
}
