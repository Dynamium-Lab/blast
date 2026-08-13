// example_ur5e_pso_comparison.cpp
//
// Runs the CPU SLSQP optimizer on the exact same UR5e task and B-spline
// resolution the GPU PSO implementation (blast-gpu repo) now uses, so the
// two can be compared fairly:
//   - Same robot: UR5e (6 joints) -- switched from the original Kinova Gen3
//     port since UR5e is blast's own primary/best-tested robot (used in
//     all 3 tutorial examples, the core optimization test suite, and the
//     published ICRA2026 benchmark -- Gen3 had no prior CPU trajectory-
//     optimization history in this codebase before this comparison work).
//   - Same task: start/end joint positions below, zero velocity and
//     acceleration at both ends (matches blastParallelPso.cuh exactly) --
//     reused from blast's own ICRA2026 UR5e benchmark task list
//     (examples/ICRA2026/ICRA2026_benchmark_speedups.cpp, get_UR5e_tasks()'s
//     first start/end pair) rather than inventing a new one.
//   - Same B-spline resolution: 24 control points, 76 evaluation points,
//     degree 5. 76 = 19*4 was chosen deliberately: with_segments divides
//     points evenly across (n_ctrl - degree) = 19 segments via integer
//     division (n_points / n_segments). A non-evenly-divisible resolution
//     silently leaves trailing points of every trajectory unchecked by the
//     optimizer during the search. PSO's own kernel checks every point
//     directly with no segment structure, so it was never affected by
//     this -- only the resolution constant needed to move for the CPU
//     side's with_segments method.
//
// Constraint scope note: this enables position, velocity, and
// ACCELERATION -- not torque. UR5e has a real, finite acceleration limit
// (13.96 rad/s^2), unlike Gen3 which has none, so acceleration is the
// natural third constraint here and is what blast's own canonical UR5e
// example (example_02_trajectory_optimization.cpp) already uses. This also
// matches PSO's GPU-side cuUR5eMultiTraj, which checks position/velocity/
// acceleration for the same reason (blast/gpu/manipulator.cuh).
//
// Optimization's constructor turns position/velocity/acceleration/
// tool_speed on by default (blast/optimization/optimization.hpp:47-53) --
// tool_speed must be explicitly forced off here (acceleration is already
// wanted, unlike the earlier Gen3 version where it had to be turned off).
// Collision is deliberately left off: the GPU/PSO code doesn't model
// collision at all for UR5e.
//
// Build:
//   cmake -B build && cmake --build build --config Release --target example_ur5e_pso_comparison
// Run:
//   ./build/examples/example_ur5e_pso_comparison

#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

int main() {
  using namespace blast;

  Manipulator ur5e = make_UR5e();

  // Same start/end joint positions as blastParallelPso.cuh's task, zero
  // velocity and acceleration at both ends (stop-to-stop). Taken from
  // examples/ICRA2026/ICRA2026_benchmark_speedups.cpp's get_UR5e_tasks().
  Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  Task  task  = Task::stop_to_stop(start, end);

  Optimization opt(ur5e, task);
  opt.bspline = Bspline(24, 76, 5, ur5e.n_joints); // matches PSO: ncontrol=24, points=76, degree=5

  opt.constraints.position     = true;
  opt.constraints.velocity     = true;
  opt.constraints.acceleration = true;
  opt.constraints.torque       = false;
  opt.constraints.tool_speed   = false;

  opt.success_tolerance = 0.01;
  opt.guess.type        = Guess::random;
  opt.method             = OptimizationMethod::with_segments; // default method

  const int runs         = 10;
  const int max_attempts = 5; // matches example_03_collision_avoidance's retry pattern for hard problems
  std::cout << "Running SLSQP on UR5e (matching PSO's task/resolution), " << runs
            << " runs, up to " << max_attempts << " retries each...\n\n";

  double total_ms = 0.0, min_ms = 1e18, max_ms = 0.0;
  int    successes = 0;

  for (int i = 0; i < runs; i++) {
    double run_total_ms = 0.0;
    Result result(&opt);
    int    attempt = 0;
    for (; attempt < max_attempts; attempt++) {
      opt.guess.type = Guess::random;
      result          = optimize(&opt);
      run_total_ms += result.compute_time;
      if (result.success)
        break;
    }

    total_ms += run_total_ms;
    min_ms    = std::min(min_ms, run_total_ms);
    max_ms    = std::max(max_ms, run_total_ms);
    if (result.success)
      successes++;

    std::cout << "run " << i << ": " << (result.success ? "success" : "FAILED")
              << "  attempts=" << std::min(attempt + 1, max_attempts)
              << "  total_time=" << run_total_ms << "ms"
              << "  evals(last)=" << result.num_eval
              << "  max_constraint=" << result.max_constraint_value
              << "  T=" << (result.x.is_empty() ? -1.0 : (double) result.x.back())
              << "\n";
  }

  std::cout << "\n--- Summary ---\n";
  std::cout << "Successes: " << successes << "/" << runs << "\n";
  std::cout << "avg: " << (total_ms / runs) << "ms  min: " << min_ms << "ms  max: " << max_ms << "ms\n";

  return 0;
}
