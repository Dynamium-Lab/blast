// SLSQP baseline across ALL 64 UR5e tasks from the lab's ICRA benchmark.
//
// Why this exists: example_ur5e_pso_comparison.cpp benchmarks SLSQP on task 0
// only, which is the task the GPU optimizer was tuned against -- and task 0
// later turned out to be the EASIEST of the 64 (its trajectory duration ranks
// 1st of 64; the median task is nearly 3x longer). The GPU side has since been
// validated on all 64, so the missing half of the comparison was SLSQP on the
// other 63.
//
// Configuration is deliberately identical to example_ur5e_pso_comparison.cpp
// so the numbers are directly comparable:
//   - position / velocity / acceleration only (matching what the GPU kernel
//     models; torque, tool speed and collisions are off on BOTH sides)
//   - Bspline(24, 76, 5) matching the GPU exactly
//   - success_tolerance 0.01, random guess, OptimizationMethod::with_segments
//     (the lab's own fastest of four methods)
//   - up to 5 random-restart attempts, times summed, matching the task-0 run
//
// Tasks come from ICRA2026_benchmark_speedups.cpp's get_UR5e_tasks(): the full
// 8x8 cross product of start and end poses, so task index = start*8 + end.
// blast's native SQP, matching ICRA2026_benchmark_speedups.cpp and
// example_ur5e_pso_comparison.cpp. Without this define the build falls back to
// NLopt's SLSQP, which is a different solver and not what the ICRA numbers
// were produced with.
#define BLAST_USE_NATIVE_SQP

#include <blast>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace blast;

static std::vector<Array> ur5e_start_positions() {
  return {
      {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883},
      {2.4316, 0.0965, 1.1824, 1.8626, -0.8567, 0.0250},
      {-0.704868, 2.17507, 0.834305, 0.132228, 2.27799, 0.00328674},
      {-0.819005, 2.50686, 0.344732, -2.85083, -2.39379, -2.93992},
      {2.04549, 0.675615, 0.137326, -0.814844, 0.449157, -3.08683},
      {2.3499, 1.44305, -1.11122, -0.333013, 0.751039, -3.09482},
      {2.36232, 1.06656, -0.132829, -0.933749, 0.788944, -3.13519},
      {2.22452, 1.15771, -0.467326, -0.690436, 0.641695, -3.12205},
  };
}

static std::vector<Array> ur5e_end_positions() {
  return {
      {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328},
      {0.1706, -2.5340, -1.2649, -2.4857, -1.4118, 2.9011},
      {0.3572, 2.9674, 0.4629, 2.8527, -1.2170, 3.1408},
      {-0.0436368, -2.65335, 0.593745, 2.05947, -1.61882, -3.12196},
      {-2.6558, -0.6703, 0.0599, -2.5311, -2.0486, -0.0127},
      {3.12394, -0.553571, -0.509807, 1.06314, 1.52902, -3.10058},
      {2.4552, -0.1792, -1.0268, -1.9354, -0.8723, -0.0198},
      {2.14686, -0.283198, -0.264971, -2.59327, -0.566506, -0.0169222},
  };
}

int main(int argc, char** argv) {
  // Control-point count. The ICRA benchmark uses 16; sample points are kept at
  // 4 per segment (segments = n_ctrl - degree) to match its convention.
  const int n_ctrl   = (argc > 1) ? atoi(argv[1]) : 24;
  const int pps      = (argc > 2) ? atoi(argv[2]) : 4;   // points per segment
  const int n_points = pps * (n_ctrl - 5);
  Manipulator ur5e = make_UR5e();

  const auto starts = ur5e_start_positions();
  const auto ends   = ur5e_end_positions();

  const int max_attempts = 5; // matches example_ur5e_pso_comparison.cpp

  std::cout << "n_ctrl=" << n_ctrl << " n_points=" << n_points << "\n";
  std::cout << "SLSQP on all " << (starts.size() * ends.size())
            << " ICRA UR5e tasks (pos/vel/acc only, with_segments, up to "
            << max_attempts << " attempts each)\n\n";
  std::cout << "task  result    attempts       time_ms          T\n";

  std::vector<double> all_T, all_ms;
  std::vector<int>    failed;
  double total_ms = 0.0;
  int    successes = 0;

  for (size_t si = 0; si < starts.size(); si++) {
    for (size_t ei = 0; ei < ends.size(); ei++) {
      const int task_idx = (int) (si * ends.size() + ei);

      Task         task = Task::stop_to_stop(starts[si], ends[ei]);
      Optimization opt(ur5e, task);
      opt.bspline = Bspline(n_ctrl, n_points, 5, ur5e.n_joints);

      opt.constraints.position     = true;
      opt.constraints.velocity     = true;
      opt.constraints.acceleration = true;
      opt.constraints.torque       = false;
      opt.constraints.tool_speed   = false;

      opt.success_tolerance = 0.01;
      opt.guess.type        = Guess::random;
      opt.method            = OptimizationMethod::with_segments;

      double run_ms  = 0.0;
      Result result(&opt);
      int    attempt = 0;
      for (; attempt < max_attempts; attempt++) {
        opt.guess.type = Guess::random;
        result         = optimize(&opt);
        run_ms += result.compute_time;
        if (result.success) break;
      }

      total_ms += run_ms;
      const double T = result.x.is_empty() ? -1.0 : (double) result.x.back();
      if (result.success) {
        successes++;
        all_T.push_back(T);
        all_ms.push_back(run_ms);
      } else {
        failed.push_back(task_idx);
      }

      printf("%4d  %-8s  %8d  %12.1f  %9.4f\n", task_idx,
             result.success ? "success" : "FAILED",
             std::min(attempt + 1, max_attempts), run_ms, T);
      fflush(stdout);
    }
  }

  std::cout << "\n--- Summary ---\n";
  std::cout << "successes: " << successes << " / " << (starts.size() * ends.size()) << "\n";
  if (!failed.empty()) {
    std::cout << "failed tasks:";
    for (int f : failed) std::cout << " " << f;
    std::cout << "\n";
  }
  std::cout << "total wall time for all tasks: " << (total_ms / 1000.0) << " s\n";

  if (!all_ms.empty()) {
    auto sorted_ms = all_ms;
    auto sorted_T  = all_T;
    std::sort(sorted_ms.begin(), sorted_ms.end());
    std::sort(sorted_T.begin(), sorted_T.end());
    double sum_ms = 0.0, sum_T = 0.0;
    for (double v : all_ms) sum_ms += v;
    for (double v : all_T) sum_T += v;
    printf("time_ms  mean=%.1f  median=%.1f  min=%.1f  max=%.1f\n",
           sum_ms / all_ms.size(), sorted_ms[sorted_ms.size() / 2], sorted_ms.front(), sorted_ms.back());
    printf("T        mean=%.4f  median=%.4f  min=%.4f  max=%.4f\n",
           sum_T / all_T.size(), sorted_T[sorted_T.size() / 2], sorted_T.front(), sorted_T.back());
  }

  return 0;
}
