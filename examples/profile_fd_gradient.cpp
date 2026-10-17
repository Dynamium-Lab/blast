// profile_fd_gradient.cpp
//
// TEMPORARY diagnostic: runs the same collision-aware optimization as
// example_03_collision_avoidance.cpp, but resets and prints the
// fd_profile:: timing accumulators (see blast/optimization/constraints.hpp)
// so we can see how much of a solve's wall-clock time is actually spent in
// the finite-difference gradient sweeps (the GPU-parallelizable part) vs.
// everything else (trajectory generation, SLSQP's own solver math, etc.).

#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <iostream>

inline blast::Task make_UR5e_task() {
  blast::Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  blast::Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  return blast::Task::stop_to_stop(start, end);
}

int main() {
  using namespace blast;

  Manipulator ur5e = make_UR5e();
  Task        task = make_UR5e_task();

  World world;
  world.add_box(
          Vec3{0.4, 0.0, 0.6},
          Vec3{0.05, 0.3, 0.3},
          Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});

  Optimization opt(ur5e, task);
  opt.bspline = Bspline(16, 110, 5, ur5e.n_joints);
  opt.world   = world;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.self_collisions     = true;
  opt.constraints.external_collisions = true;

  opt.success_tolerance = 0.01;

  const int n_runs             = 10;
  double    total_compute_ms   = 0, total_fd_ms = 0;
  double    sum_torque_ms = 0, sum_tool_ms = 0, sum_self_col_ms = 0, sum_ext_col_ms = 0;
  int       successes = 0;

  for (int run = 0; run < n_runs; run++) {
    fd_profile::reset();
    opt.guess.type = Guess::random;
    Result result  = optimize(&opt);

    double fd_ms = (fd_profile::torque_ns + fd_profile::tool_speed_ns +
                     fd_profile::self_collisions_ns + fd_profile::external_collisions_ns) /
                   1.0e6;

    printf("run %2d: success=%s  compute_time=%8.3f ms  num_eval=%4d  FD-gradient time=%7.3f ms (%.1f%%)\n",
           run, result.success ? "yes" : "no ", (double) result.compute_time, result.num_eval,
           fd_ms, result.compute_time > 0 ? 100.0 * fd_ms / result.compute_time : 0.0);

    total_compute_ms += result.compute_time;
    total_fd_ms += fd_ms;
    sum_torque_ms += fd_profile::torque_ns / 1.0e6;
    sum_tool_ms += fd_profile::tool_speed_ns / 1.0e6;
    sum_self_col_ms += fd_profile::self_collisions_ns / 1.0e6;
    sum_ext_col_ms += fd_profile::external_collisions_ns / 1.0e6;
    successes += result.success;
  }

  printf("\n--- summary over %d runs ---\n", n_runs);
  printf("successes: %d/%d\n", successes, n_runs);
  printf("total compute time:     %8.3f ms\n", total_compute_ms);
  printf("total FD-gradient time: %8.3f ms\n", total_fd_ms);
  printf("FD-gradient share:      %8.1f%%\n", 100.0 * total_fd_ms / total_compute_ms);
  printf("  of which torque:              %8.3f ms\n", sum_torque_ms);
  printf("  of which tool_speed:          %8.3f ms\n", sum_tool_ms);
  printf("  of which self_collisions:     %8.3f ms\n", sum_self_col_ms);
  printf("  of which external_collisions: %8.3f ms\n", sum_ext_col_ms);

  return 0;
}
