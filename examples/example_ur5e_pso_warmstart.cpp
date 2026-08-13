// example_ur5e_pso_warmstart.cpp
//
// The hybrid approach: read the GPU PSO optimizer's best result and use it
// as SLSQP's initial guess instead of a random one, on the exact same UR5e
// task/resolution/constraints as example_ur5e_pso_comparison.cpp.
//
// Rationale: PSO's earlier failures to match SLSQP's precision came from
// its own update rule being far less informed than a gradient step. SLSQP's
// failures (without this warm start) came from landing in a bad random
// starting basin on a hard, high-dimensional problem. Handing SLSQP a
// starting point that's already close to a good, feasible answer (found
// cheaply and safely by PSO) should let it converge fast and reliably,
// combining PSO's speed at finding the right neighborhood with SLSQP's
// precision at refining within it.
//
// The current baseline blastParallelPso.cuh doesn't write its result to a
// file (that was a later addition that got reverted along with other
// tuning changes) -- to generate PSO_RESULT_PATH, capture stdout and
// extract the "Global Best Position:" line manually:
//   cd ~/Downloads/blast-gpu && ./build/GPU/gpu 2>&1 | tee pso_run.log
//   grep "Global Best Position:" pso_run.log | sed 's/Global Best Position: //' | tr ' ' '\n' | grep -v '^$' > pso_result.txt
// Then build and run this:
//   cmake --build build --target example_ur5e_pso_warmstart
//   ./build/examples/example_ur5e_pso_warmstart

#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <fstream>
#include <iostream>
#include <sstream>

// Hardcoded to this machine's layout (blast and blast-gpu as sibling
// directories under ~/Downloads) -- change if your setup differs.
static const char* PSO_RESULT_PATH = "/home/sannan/Downloads/blast-gpu/pso_result.txt";

int main() {
  using namespace blast;

  Manipulator ur5e = make_UR5e();

  // Same start/end joint positions as blastParallelPso.cuh's task -- must
  // match exactly, since the warm-start vector's meaning depends on it.
  Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  Task  task  = Task::stop_to_stop(start, end);

  Optimization opt(ur5e, task);
  opt.bspline = Bspline(24, 76, 5, ur5e.n_joints); // matches PSO: ncontrol=24, points=76, degree=5

  // Position/velocity/acceleration -- matches PSO's GPU-side
  // cuUR5eMultiTraj (blast/gpu/manipulator.cuh). tool_speed must be
  // explicitly forced off (Optimization's constructor turns it on by
  // default, blast/optimization/optimization.hpp:47-53); acceleration is
  // already wanted so it's left on.
  opt.constraints.position     = true;
  opt.constraints.velocity     = true;
  opt.constraints.acceleration = true;
  opt.constraints.torque       = false;
  opt.constraints.tool_speed   = false;

  opt.success_tolerance = 0.01;
  opt.method             = OptimizationMethod::with_segments;

  const u32 x_len = (u32) opt.x_len();

  // ---- Read PSO's result -----------------------------------------------
  std::ifstream in(PSO_RESULT_PATH);
  if (!in) {
    std::cerr << "Could not open " << PSO_RESULT_PATH
              << " -- see the comment at the top of this file for how to generate it.\n";
    return 1;
  }
  Array pso_x(x_len);
  for (u32 i = 0; i < x_len; i++) {
    if (!(in >> pso_x[i])) {
      std::cerr << "pso_result.txt has fewer than " << x_len << " values (expected exactly "
                << x_len << " for this problem's dimensionality) -- was it generated with the "
                << "same ncontrol/points/joints as this program?\n";
      return 1;
    }
  }
  std::cout << "Read PSO warm-start vector (" << x_len << " values), T=" << pso_x.back() << "\n\n";

  // ---- Attempt 1: SLSQP starting from PSO's answer ----------------------
  opt.guess = Guess(pso_x); // Guess::custom, initial_x = pso_x

  auto t0     = std::chrono::high_resolution_clock::now();
  Result result = optimize(&opt);
  auto t1     = std::chrono::high_resolution_clock::now();

  std::cout << "Warm-started from PSO: " << (result.success ? "success" : "FAILED")
            << "  time=" << result.compute_time << "ms"
            << "  evals=" << result.num_eval
            << "  max_constraint(coarse," << result.max_constraint_idx << ")=" << result.max_constraint_value
            << "  max_constraint(dense," << result.max_constraint_more_points_idx << ")=" << result.max_constraint_more_points_value
            << "  T=" << (result.x.is_empty() ? -1.0 : (double) result.x.back())
            << "\n";

  // ---- Fallback: if the warm start didn't converge, try random restarts,
  // same as example_03_collision_avoidance's pattern -- report both so we
  // can see whether the warm start actually helped or was a wash.
  if (!result.success) {
    std::cout << "\nWarm start did not succeed on the first attempt -- falling back to "
                 "random restarts (up to 5), for comparison:\n";
    for (int attempt = 0; attempt < 5 && !result.success; attempt++) {
      opt.guess.type = Guess::random;
      result          = optimize(&opt);
      std::cout << "  fallback attempt " << attempt << ": " << (result.success ? "success" : "FAILED")
                << "  time=" << result.compute_time << "ms"
                << "  max_constraint=" << result.max_constraint_value
                << "  T=" << (result.x.is_empty() ? -1.0 : (double) result.x.back())
                << "\n";
    }
  }

  std::cout << "\n--- Summary ---\n";
  std::cout << "Final: " << (result.success ? "success" : "FAILED")
            << "  T=" << (result.x.is_empty() ? -1.0 : (double) result.x.back()) << "\n";

  return 0;
}
