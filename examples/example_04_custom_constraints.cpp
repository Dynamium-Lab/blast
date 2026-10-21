// example_04_two_stage_optimization.cpp
//
// Demonstrates a two-stage trajectory optimization pipeline:
//   1. Minimize trajectory duration.
//   2. Minimize total sampled joint velocity while preserving the
//      duration found in stage 1.
//
// The second stage starts from the first-stage solution and adds a
// maximum-duration constraint. This separates the primary objective
// (execution time) from the secondary objective (motion refinement).
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_04_two_stage_optimization
// Run:
//   ./build/examples/example_04_two_stage_optimization

#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <cmath>
#include <iostream>

using namespace blast;

// -------------------------------------------------------------------------
// Secondary objective: sum of absolute sampled joint velocities.
// Lower values indicate lower aggregate sampled joint speed, but do not
// guarantee smoothness in every sense (e.g., jerk or acceleration continuity).
// -------------------------------------------------------------------------
real objective_velocity_sum(Optimization* opt, Array* x) {
  opt->bspline.compute_trajectory(*x, opt->task);

  real sum = 0.0;

  for (int i = 0; i < opt->bspline.traj.vel.cols; ++i) {
    for (int j = 0; j < opt->bspline.traj.vel.rows; ++j) {
      sum += std::abs(opt->bspline.traj.vel(j, i));
    }
  }

  return sum;
}

// -------------------------------------------------------------------------
// Secondary-stage constraint: trajectory duration <= maximum duration.
// Constraint convention: returned value must be <= 0.
// -------------------------------------------------------------------------
Array constraint_max_duration(Optimization* opt, Array* x) {
  const auto* max_duration = static_cast<const real*>(opt->custom_data);
  return Array(1, x->back() - *max_duration);
}

int main() {
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

  std::cout << "  Success: " << (time_result.success ? "yes" : "no")
            << "\n  Compute time (ms): " << time_result.compute_time
            << "\n  Evaluations: " << time_result.num_eval << '\n';

  if (!time_result.success || time_result.x.is_empty()) {
    std::cerr << "Stage 1 failed; cannot proceed to stage 2.\n";
    return 1;
  }

  const real optimal_duration = time_result.x.back();

  std::cout << "  Optimal duration (s): " << optimal_duration << '\n';

  // -----------------------------------------------------------------------
  // Step 3 — Configure the secondary optimization.
  //
  // Use the primary solution as the initial guess. Remove the duration
  // objective and replace it with the velocity objective. Constrain the
  // duration to the primary-stage optimum.
  // -----------------------------------------------------------------------
  opt.guess.type      = Guess::custom;
  opt.guess.initial_x = time_result.x;

  opt.objective.time_weight = 0.0;
  opt.objective.add_custom_objective(objective_velocity_sum, 1.0);

  // Keep the duration at or below the stage-1 optimum.
  // The referenced value must remain alive throughout optimize().
  const real max_duration = optimal_duration + 0.01;
  opt.custom_data         = const_cast<real*>(&max_duration);

  opt.constraints.add_constraint(constraint_max_duration, 1);

  // Use the baseline method for the custom-objective stage.
  opt.method = OptimizationMethod::baseline;

  // -----------------------------------------------------------------------
  // Step 4 — Secondary optimization: refine the trajectory.
  // -----------------------------------------------------------------------
  std::cout << "\nStage 2: minimizing joint-velocity sum "
               "at the primary-stage duration...\n";

  Result smooth_result = optimize(&opt);

  std::cout << "  Success: " << (smooth_result.success ? "yes" : "no")
            << "\n  Compute time (ms): " << smooth_result.compute_time
            << "\n  Evaluations: " << smooth_result.num_eval << '\n';

  if (!smooth_result.success || smooth_result.x.is_empty()) {
    std::cerr << "Stage 2 failed; retaining the stage-1 solution.\n";
    return 1;
  }

  // -----------------------------------------------------------------------
  // Step 5 — Compare the two solutions.
  // -----------------------------------------------------------------------
  const real velocity_before =
          objective_velocity_sum(&opt, &time_result.x);

  const real velocity_after =
          objective_velocity_sum(&opt, &smooth_result.x);

  std::cout << "\n--- Pipeline results ---\n"
            << "Stage 1 duration (s):       " << time_result.x.back() << '\n'
            << "Stage 2 duration (s):       " << smooth_result.x.back() << '\n'
            << "Velocity sum before:        " << velocity_before << '\n'
            << "Velocity sum after:         " << velocity_after << '\n'
            << "Velocity-sum improvement:   "
            << (velocity_before > 0.0
                        ? 100.0 * (velocity_before - velocity_after) / velocity_before
                        : 0.0)
            << "%\n";

  print_to_csv(time_result.trajectory, "../../../../examples/time_traj.csv");
  print_to_csv(smooth_result.trajectory, "../../../../examples/smooth_traj.csv");

  return 0;
}
