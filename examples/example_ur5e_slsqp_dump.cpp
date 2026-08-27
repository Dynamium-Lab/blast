// Run SLSQP on one ICRA task and dump its solution vector, so the answer can
// be pushed through the SAME independent validator used on the GPU results.
//
// Written because the 64-task comparison produced a result that looks
// physically impossible: on task 26, joint 5 must travel 6.081 rad and the
// velocity limit is 3.1416 rad/s, so the duration cannot be below ~1.94 s even
// at full speed the whole way -- yet SLSQP reported T = 0.616 s. Either the
// answer is not feasible, or something about how the optimization vector is
// being read is wrong. This settles which.
//
//   example_ur5e_slsqp_dump <task_idx> [out_file]
#include <blast>

#include <cstdio>
#include <iostream>
#include <vector>

using namespace blast;

static const double S[8][6] = {
    {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883},
    {2.4316, 0.0965, 1.1824, 1.8626, -0.8567, 0.0250},
    {-0.704868, 2.17507, 0.834305, 0.132228, 2.27799, 0.00328674},
    {-0.819005, 2.50686, 0.344732, -2.85083, -2.39379, -2.93992},
    {2.04549, 0.675615, 0.137326, -0.814844, 0.449157, -3.08683},
    {2.3499, 1.44305, -1.11122, -0.333013, 0.751039, -3.09482},
    {2.36232, 1.06656, -0.132829, -0.933749, 0.788944, -3.13519},
    {2.22452, 1.15771, -0.467326, -0.690436, 0.641695, -3.12205},
};
static const double E[8][6] = {
    {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328},
    {0.1706, -2.5340, -1.2649, -2.4857, -1.4118, 2.9011},
    {0.3572, 2.9674, 0.4629, 2.8527, -1.2170, 3.1408},
    {-0.0436368, -2.65335, 0.593745, 2.05947, -1.61882, -3.12196},
    {-2.6558, -0.6703, 0.0599, -2.5311, -2.0486, -0.0127},
    {3.12394, -0.553571, -0.509807, 1.06314, 1.52902, -3.10058},
    {2.4552, -0.1792, -1.0268, -1.9354, -0.8723, -0.0198},
    {2.14686, -0.283198, -0.264971, -2.59327, -0.566506, -0.0169222},
};

int main(int argc, char** argv) {
  const int   task_idx = (argc > 1) ? atoi(argv[1]) : 26;
  const char* out_path = (argc > 2) ? argv[2] : "/home/sannan/Downloads/blast-gpu/pso_result.txt";

  Manipulator ur5e = make_UR5e();
  Array start(6), end(6);
  for (u32 j = 0; j < 6; j++) {
    start[j] = (real) S[task_idx / 8][j];
    end[j]   = (real) E[task_idx % 8][j];
  }

  printf("task %d\n", task_idx);
  printf("  joint displacements:");
  double dmax = 0.0;
  for (u32 j = 0; j < 6; j++) {
    const double d = std::fabs((double) end[j] - (double) start[j]);
    dmax           = std::max(dmax, d);
    printf(" %.3f", d);
  }
  printf("\n  largest = %.3f rad; at the 3.1416 rad/s velocity limit that needs >= %.3f s\n",
         dmax, dmax / 3.1416);

  Task         task = Task::stop_to_stop(start, end);
  Optimization opt(ur5e, task);
  opt.bspline = Bspline(24, 76, 5, ur5e.n_joints);

  opt.constraints.position     = true;
  opt.constraints.velocity     = true;
  opt.constraints.acceleration = true;
  opt.constraints.torque       = false;
  opt.constraints.tool_speed   = false;

  opt.success_tolerance = 0.01;
  opt.guess.type        = Guess::random;
  opt.method            = OptimizationMethod::with_segments;

  Result result(&opt);
  for (int attempt = 0; attempt < 5; attempt++) {
    opt.guess.type = Guess::random;
    result         = optimize(&opt);
    if (result.success) break;
  }

  printf("  SLSQP: success=%d  x.size=%u  T(x.back())=%.6f\n",
         (int) result.success, result.x.size,
         result.x.is_empty() ? -1.0 : (double) result.x.back());
  printf("  reported max_constraint=%.6g  (dense: %.6g)\n",
         (double) result.max_constraint_value, (double) result.max_constraint_more_points_value);

  if (result.x.is_empty()) { printf("  no solution to dump\n"); return 1; }

  FILE* f = fopen(out_path, "w");
  for (u32 i = 0; i < result.x.size; i++) fprintf(f, "%.17g\n", (double) result.x[i]);
  fclose(f);
  printf("  wrote %u values to %s\n", result.x.size, out_path);
  return 0;
}
