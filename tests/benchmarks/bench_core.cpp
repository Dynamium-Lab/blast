// Microbenchmarks for Blast. Built only when BLAST_BUILD_BENCHMARKS=ON. Run the
// same binary built with BLAST_USE_DOUBLES=ON and with
// -DBLAST_USE_DOUBLES=OFF -DBLAST_USE_NATIVE_SQP=ON to compare double vs float.
//
// The vectorized kernels (dot/norm/sincos) show the raw SIMD win: float doubles
// the AVX2 lane width (8-wide Vec8f vs 4-wide Vec4d). The full trajectory
// optimization shows the end-to-end effect: there float is *not* faster, because
// the solve is dominated by the SQP solver's scalar dense linear algebra (and a
// coarser float epsilon costs extra solver iterations), not by these kernels.
#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#include <blast>

#include "manipulator/UR5e.hpp" // make_UR5e

#include <string>

using namespace blast;

int main() {
  const std::string prec = sizeof(real) == 8 ? "double" : "float";
  ankerl::nanobench::Bench bench;
  bench.title("Blast core (real = " + prec + ")");

  // --- Vectorized math kernels --------------------------------------------
  constexpr u32 N = 1024;
  Array         a(N), b(N), sines(N), cosines(N);
  for (u32 i = 0; i < N; i++) {
    a[i] = real(i) * real(0.001);
    b[i] = real(i) * real(0.002);
  }

  bench.run("dot (N=1024)", [&] {
    ankerl::nanobench::doNotOptimizeAway(dot(a, b));
  });

  bench.run("norm (N=1024)", [&] {
    ankerl::nanobench::doNotOptimizeAway(norm(a));
  });

  bench.run("sincos (N=1024)", [&] {
    sincos(a, sines, cosines);
    ankerl::nanobench::doNotOptimizeAway(sines[0]);
  });

  // --- Full trajectory optimization ---------------------------------------
  // UR5e stop-to-stop with the default pva + tool_speed constraints. Exercises
  // trajectory evaluation, constraints, gradients and the solver end-to-end.
  {
    Manipulator robot = make_UR5e();
    Array       q0    = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
    Array       q1    = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
    Task        task  = Task::stop_to_stop(q0, q1);

    // A solve is ~20 ms, so force several iterations per epoch for a stable mean
    // (the random initial guess also makes per-solve time vary run to run).
    bench.minEpochIterations(10).run("optimize UR5e stop-to-stop", [&] {
      Optimization opt(robot, task);
      opt.success_tolerance = real(0.01);
      Result result         = optimize(&opt);
      ankerl::nanobench::doNotOptimizeAway(result.success);
    });
  }

  return 0;
}
