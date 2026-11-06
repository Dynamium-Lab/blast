// example_04_two_stage_optimization.cpp
//
// This example shows how to:
//   1. Chain two optimizations: first minimize the trajectory duration, then
//      refine the result with a secondary objective at that same duration
//   2. Write a custom objective function and register it with
//      objective.add_custom_objective()
//   3. Write a custom constraint function and register it with
//      constraints.add_constraint()
//   4. Pass user data to custom callbacks through Optimization::custom_data
//   5. Warm-start the second stage from the first-stage solution
//      (Guess::custom)
//
// The secondary objective is selectable on the command line:
//   vel : minimize the sum of absolute sampled joint velocities
//   acc : minimize the sum of absolute sampled joint accelerations
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_04_two_stage_optimization
// Run:
//   ./build/examples/example_04_two_stage_optimization [vel|acc]   (default: acc)

#define BLAST_USE_NATIVE_SQP // use the built-in SQP solver

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
//
// compute_trajectory() evaluates the B-spline defined by x and stores the
// sampled positions, velocities and accelerations in opt->bspline.traj,
// one column per sample and one row per joint.
// -------------------------------------------------------------------------
template<typename M>
static real sum_abs(const M& m) {
  real sum = 0.0;
  for (int i = 0; i < m.cols; ++i)
    for (int j = 0; j < m.rows; ++j)
      sum += std::abs(m(j, i));
  return sum;
}

// Sum of absolute sampled joint velocities. Lower values mean lower aggregate
// joint speed; this says nothing about acceleration or jerk.
real objective_velocity_sum(Optimization* opt, Array* x) {
  opt->bspline.compute_trajectory(*x, opt->task);
  return sum_abs(opt->bspline.traj.vel);
}

// Sum of absolute sampled joint accelerations. Lower values mean lower
// aggregate joint acceleration; this does not bound jerk.
real objective_acceleration_sum(Optimization* opt, Array* x) {
  opt->bspline.compute_trajectory(*x, opt->task);
  return sum_abs(opt->bspline.traj.acc);
}

using ObjectiveFn = real (*)(Optimization*, Array*);

// -------------------------------------------------------------------------
// Custom constraint: trajectory duration <= maximum duration.
//
// A custom constraint has the signature Array(Optimization*, Array*) and
// returns one value per constraint. Convention: each value must be <= 0 for
// the constraint to be satisfied.
//
// The maximum duration is not part of the problem definition, so it is
// passed in through opt->custom_data (a user-owned pointer, see Step 4).
// -------------------------------------------------------------------------
Array constraint_max_duration(Optimization* opt, Array* x) {
  const auto* max_duration = static_cast<const real*>(opt->custom_data);
  return Array(1, x->back() - *max_duration);
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
  ur5e.base_position = {-0.5, -0.3, 0.35};
  ur5e.base_rotation = {-1, 0, 0, 0, -1, 0, 0, 0, 1}; // 180° about z

  Array start = {-0.704868, 2.17507, 0.834305, 0.132228, 2.27799, 0.00328674};
  Array end   = {3.12394, -0.553571, -0.509807, 1.06314, 1.52902, -3.10058};
  Task  task  = Task::stop_to_stop(start, end);

  World world = get_kitchen_no_doors();

  // -----------------------------------------------------------------------
  // Step 2 — Create the Optimization problem.
  // Bspline(n_ctrl, n_points, degree, n_joints): 12 control points, 100
  // sample points, degree 5. The sample points are also where the custom
  // objectives above are evaluated.
  // -----------------------------------------------------------------------
  Optimization opt(ur5e, task);
  opt.bspline = Bspline(12, 100, 5, ur5e.n_joints);
  opt.set_world(world);

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.external_collisions = true;

  opt.success_tolerance = 0.01;

  // -----------------------------------------------------------------------
  // Step 3 — Stage 1: minimize trajectory duration.
  // time_weight is the weight of the built-in duration objective.
  // -----------------------------------------------------------------------
  opt.objective.time_weight = 1.0;

  std::cout << "Stage 1: minimizing trajectory duration...\n";
  Result time_result = optimize(&opt);

  std::cout << "  Success:           " << (time_result.success ? "yes" : "no") << "\n"
            << "  Compute time (ms): " << time_result.compute_time << "\n"
            << "  Evaluations:       " << time_result.num_eval << "\n";

  if (!time_result.success || time_result.x.is_empty()) {
    std::cerr << "Stage 1 failed; cannot proceed to stage 2.\n";
    return 1;
  }

  // Not const: custom_data is a non-const pointer. This variable must stay
  // alive for as long as stage 2 runs, since the constraint reads it.
  real optimal_duration = time_result.x.back();
  std::cout << "  Duration (s):      " << optimal_duration << "\n";

  // -----------------------------------------------------------------------
  // Step 4 — Stage 2 setup: refine at the stage-1 duration.
  //
  //   - Warm start: Guess::custom uses initial_x as the starting point, so
  //     the optimizer begins from the time-optimal solution, which is
  //     already feasible.
  //   - Objective: disable the duration term and register the selected
  //     custom objective with weight 1.0.
  //   - Constraint: the duration stays a decision variable, but the custom
  //     constraint keeps it from exceeding the stage-1 optimum. The second
  //     argument of add_constraint() is the number of values the function
  //     returns.
  // -----------------------------------------------------------------------
  Optimization refined = opt;
  refined.method       = OptimizationMethod::baseline; // custom only supported on baseline

  refined.guess.type      = Guess::custom;
  refined.guess.initial_x = time_result.x;

  refined.objective.time_weight = 0.0;
  refined.objective.add_custom_objective(secondary_fn, 1.0);

  refined.custom_data = &optimal_duration;
  refined.constraints.add_constraint(constraint_max_duration, 1);

  // -----------------------------------------------------------------------
  // Step 5 — Stage 2: minimize the secondary objective.
  // -----------------------------------------------------------------------
  std::cout << "\nStage 2: minimizing " << secondary_label
            << " at the stage-1 duration...\n";
  Result refined_result = optimize(&refined);

  std::cout << "  Success:           " << (refined_result.success ? "yes" : "no") << "\n"
            << "  Compute time (ms): " << refined_result.compute_time << "\n"
            << "  Evaluations:       " << refined_result.num_eval << "\n";

  // -----------------------------------------------------------------------
  // Step 6 — Compare both solutions on the secondary objective.
  // The objective is evaluated on a scratch copy of the problem so that the
  // B-splines stored in opt and refined are left untouched.
  // -----------------------------------------------------------------------
  Optimization evaluator = refined;
  const real   obj_time  = secondary_fn(&evaluator, &time_result.x);

  std::cout << "\n--- Results (" << secondary_label << ") ---\n"
            << "Time-optimal:\n"
            << "  Duration (s):      " << time_result.x.back() << "\n"
            << "  " << secondary_label << ": " << obj_time << "\n";

  std::cout << "Refined:\n";
  if (refined_result.success && !refined_result.x.is_empty()) {
    const real obj_refined = secondary_fn(&evaluator, &refined_result.x);
    const real improvement = obj_time > 0.0 ? 100.0 * (obj_time - obj_refined) / obj_time : 0.0;

    std::cout << "  Duration (s):      " << refined_result.x.back() << "\n"
              << "  " << secondary_label << ": " << obj_refined << "\n"
              << "  Improvement (%):   " << improvement << "\n";
  } else {
    std::cout << "  failed\n";
  }

  // -----------------------------------------------------------------------
  // Step 7 — Export the trajectories to CSV for plotting (optional).
  // -----------------------------------------------------------------------
  // const std::string dir = "../../../../examples/";

  // print_to_csv(time_result.trajectory, (dir + "time_traj.csv").c_str());
  // if (refined_result.success)
  //   print_to_csv(refined_result.trajectory, (dir + "refined_" + secondary_tag + "_traj.csv").c_str());

  return refined_result.success ? 0 : 1;
}
