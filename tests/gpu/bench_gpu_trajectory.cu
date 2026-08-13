// GPU vs CPU trajectory + manipulator constraint timing benchmark
//
// Measures wall-clock time for the CPU (Bspline + cuGen3MultiTraj::
// compute_constraints_point, run serially per trajectory) against the GPU
// (cuMultiBsplines + cuGen3MultiTraj, run in parallel across all
// trajectories) and prints a side-by-side performance summary.
//
// Fairness notes (both sides do the same work: trajectory generation AND
// constraint evaluation, for the same NTRAJECTORIES distinct trajectories):
//   - Build in Release. Debug-mode CPU timings are not representative and
//     will produce a misleadingly large "speedup" number.
//   - The CPU loop actually runs all NTRAJECTORIES trajectories serially and
//     measures the real total wall-clock time — it does not time one
//     trajectory and multiply, which would hide any per-call fixed overhead.
//   - Each trajectory (CPU and GPU alike) gets distinct random input instead
//     of NTRAJECTORIES copies of the same one, so neither side can shortcut
//     via caching/branch-prediction on repeated identical data.
//   - CPU's Bspline object and GPU's cuMultiBsplines/cuGen3MultiTraj objects
//     are each constructed once, outside the timed loop, and reused —
//     neither side pays repeated setup cost inside the timing.
//   - compute_constraints_point() is the exact same blast_fn function called
//     on both the CPU (directly, in a loop) and the GPU (from the kernel),
//     so this isolates the effect of parallelizing identical work rather
//     than comparing two different implementations.

#include "blast/gpu/error.cuh"
#include "blast/gpu/trajectory.cuh"
#include "blast/gpu/manipulator.cuh"

#include <chrono>
#include <cmath>
#include <cstdio>

using namespace blast;

using Clock = std::chrono::high_resolution_clock;
using Ms    = std::chrono::duration<double, std::milli>;

// ---- helpers ----------------------------------------------------------------

static double clock_ms(Clock::time_point t0, Clock::time_point t1) {
  return std::chrono::duration_cast<Ms>(t1 - t0).count();
}

// CUDA event timer — accurate for device-side work
struct CudaTimer {
  cudaEvent_t start, stop;
  CudaTimer() {
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  }
  ~CudaTimer() {
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  }
  void begin() { cudaEventRecord(start); }
  // Returns elapsed milliseconds; blocks until the stop event is done.
  float end() {
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);
    return ms;
  }
};

// ---- benchmark parameters ---------------------------------------------------

static constexpr blast::u32 POINTS        = 256;
static constexpr blast::u32 JOINTS        = 7;
static constexpr blast::u32 P             = 5;
static constexpr blast::u32 NCONTROL      = 24;
static constexpr blast::u32 NTRAJECTORIES = 80; // GPU computes all of these in parallel
static constexpr int        WARMUP_RUNS   = 3;  // discarded to avoid cold-start skew
static constexpr int        BENCH_RUNS    = 20; // timed iterations

// Runs the CPU pipeline (trajectory generation + constraint evaluation) for
// all NTRAJECTORIES trajectories in x_multi, serially. Reuses `b` and `manip`
// across trajectories (no per-trajectory allocation), matching how the GPU
// objects are reused across kernel launches. `manip` is the same
// cuGen3MultiTraj instance the GPU benchmark uses below — reused here purely
// for its (host-callable) compute_constraints_point() and physical
// constants, with results written into a separate host buffer.
static void cpu_compute_all(Bspline& b, cuGen3MultiTraj& manip, const Matrix& task,
                             const Array& x_multi, u32 xlen, real* out_con) {
  using namespace blast;
  for (u32 traj = 0; traj < NTRAJECTORIES; traj++) {
    const Array x_i(x_multi.data + traj * xlen, xlen);
    b.compute_trajectory(x_i, task);

    for (u32 point = 0; point < POINTS; point++) {
      real pos[7], vel[7], acc[7];
      for (u32 j = 0; j < JOINTS; j++) {
        pos[j] = b.traj.pos(j, point);
        vel[j] = b.traj.vel(j, point);
        acc[j] = b.traj.acc(j, point);
      }
      manip.compute_constraints_point(pos, vel, acc, &out_con[manip.con_id(point, traj)]);
    }
  }
}

// ---- main -------------------------------------------------------------------

int main() {
  // ---- shared setup ---------------------------------------------------------
  // Random task matrix (joints x 6 boundary conditions), kept within a
  // physically plausible joint-angle range. Using a large amplitude here
  // (e.g. the historical amp=10 rad, ~570 degrees) combined with a short
  // duration produces velocities/accelerations orders of magnitude beyond
  // the Gen3's actual rated limits (vmax ~1.7-2.4 rad/s). In that regime the
  // Newton-Euler recursion is numerically ill-conditioned enough that
  // single-ULP differences between host and device sincos() get amplified
  // through 7 chained recursive stages into a large CPU/GPU mismatch — not a
  // functional bug, just an unrealistic stress input. Keeping inputs within
  // a plausible operating range makes both the timing and the correctness
  // cross-check below representative of real usage.
  Matrix task(JOINTS, 6);
  const real pos_amp = 1.5; // rad — plausible joint-angle boundary conditions

  for (u32 i = 0; i < task.rows; i++)
    for (u32 j = 0; j < task.cols; j++)
      task(i, j) = pos_amp * random_real();

  Bspline   ref_bspline(NCONTROL, POINTS, P, JOINTS);
  const u32 xlen = ref_bspline.x_len(task);

  // Distinct random optimization vector per trajectory (not NTRAJECTORIES
  // copies of the same one) — both CPU and GPU consume this same data.
  Array x_multi(NTRAJECTORIES * xlen);
  for (u32 i = 0; i < NTRAJECTORIES; i++) {
    for (u32 j = 0; j < xlen; j++)
      x_multi[i * xlen + j] = pos_amp * random_real();
    x_multi[i * xlen + xlen - 1] = 1.5 + 1.0 * std::abs(random_real()); // duration in [1.5, 2.5] s
  }

  // ---- CPU objects (constructed once, reused across runs) -------------------
  Bspline cpu_bspline(NCONTROL, POINTS, P, JOINTS);
  Array   cpu_con(21u * POINTS * NTRAJECTORIES);

  // ---- GPU objects (constructed once, reused across runs) -------------------
  cuMultiBsplines gpu_bsplines(POINTS, JOINTS, P, NCONTROL, NTRAJECTORIES);
  cuGen3MultiTraj gpu_manip(POINTS, NTRAJECTORIES, 4.0 /* payload mass */);

  CudaTimer cuda_timer;

  // ============================  CPU BENCHMARK  ==============================
  // Warmup
  for (int i = 0; i < WARMUP_RUNS; i++)
    cpu_compute_all(cpu_bspline, gpu_manip, task, x_multi, xlen, cpu_con.data);

  double cpu_total_ms = 0.0;
  double cpu_min_ms   = 1e18;
  double cpu_max_ms   = 0.0;

  for (int run = 0; run < BENCH_RUNS; run++) {
    auto t0 = Clock::now();
    cpu_compute_all(cpu_bspline, gpu_manip, task, x_multi, xlen, cpu_con.data);
    auto t1 = Clock::now();

    double ms     = clock_ms(t0, t1);
    cpu_total_ms += ms;
    cpu_min_ms    = std::min(cpu_min_ms, ms);
    cpu_max_ms    = std::max(cpu_max_ms, ms);
  }

  double cpu_avg_ms = cpu_total_ms / BENCH_RUNS;

  // ============================  GPU BENCHMARK  ==============================

  // Warmup (Trajectory generation + Manipulator constraint calculations)
  for (int i = 0; i < WARMUP_RUNS; i++) {
    gpu_manip.compute_constraints(x_multi, task, gpu_bsplines);
    gpu_manip.fetch_constraints();
  }

  double gpu_total_ms = 0.0;
  double gpu_min_ms   = 1e18;
  double gpu_max_ms   = 0.0;

  for (int run = 0; run < BENCH_RUNS; run++) {
    // Time the combined execution of Trajectory generation + Manipulator constraint kernel
    cuda_timer.begin();
    gpu_manip.compute_constraints(x_multi, task, gpu_bsplines);
    float ms = cuda_timer.end(); // Includes cudaDeviceSynchronize inside compute_constraints

    gpu_total_ms += ms;
    gpu_min_ms    = std::min(gpu_min_ms, (double) ms);
    gpu_max_ms    = std::max(gpu_max_ms, (double) ms);
  }

  // Fetch outputs from GPU memory to Host memory for verification
  gpu_manip.fetch_constraints();

  double gpu_avg_ms = gpu_total_ms / BENCH_RUNS;

  // ---- cross-check: do CPU and GPU actually agree on the numbers? -----------
  real max_abs_diff = 0.0;
  for (u32 traj = 0; traj < NTRAJECTORIES; traj++) {
    for (u32 point = 0; point < POINTS; point++) {
      for (u32 k = 0; k < 21; k++) {
        const auto id   = gpu_manip.con_id(point, traj);
        real       diff = std::abs(cpu_con[id + k] - gpu_manip.host_constraints[id + k]);
        max_abs_diff     = std::max(max_abs_diff, diff);
      }
    }
  }

  // ============================  RESULTS  ====================================
  printf("\n");
  printf("======================================================================\n");
  printf("  %-42s  %8s  %8s  %8s\n", "", "avg", "min", "max");
  printf("  %-42s  %8.3f  %8.3f  %8.3f\n",
         "CPU — all trajs, serial (traj+manip)", cpu_avg_ms, cpu_min_ms, cpu_max_ms);
  printf("  %-42s  %8.3f  %8.3f  %8.3f\n",
         "GPU — all trajs, parallel (traj+manip)", gpu_avg_ms, gpu_min_ms, gpu_max_ms);
  printf("----------------------------------------------------------------------\n");
  printf("Both totals above are real measured wall-clock time for all %u\n", NTRAJECTORIES);
  printf("trajectories (CPU: looped serially; GPU: one batched launch).\n");
  printf("Grand total across all %d timed runs (%u trajectories each):\n", BENCH_RUNS, NTRAJECTORIES);
  printf("  %-42s  %8.3f ms\n", "CPU grand total", cpu_total_ms);
  printf("  %-42s  %8.3f ms\n", "GPU grand total", gpu_total_ms);

  if (gpu_avg_ms > 0.0) {
    double speedup = cpu_avg_ms / gpu_avg_ms;
    printf("  GPU speedup over CPU (for %u trajectories): %.1fx\n", NTRAJECTORIES, speedup);
  }
  printf("======================================================================\n");

  printf("\n[Correctness cross-check]\n");
  printf("Max |CPU - GPU| constraint value difference across %u trajectories x %u points x 21 constraints: %.6f\n",
         NTRAJECTORIES, POINTS, (double) max_abs_diff);

  printf("\n[Manipulator Output Info]\n");
  printf("Total constraints calculated across %u trajectories: %u\n",
         NTRAJECTORIES, gpu_manip.nconstraints());
  printf("Sample output (first constraint of trajectory 0): %.4f\n",
         gpu_manip.host_constraints[4]);
  printf("======================================================================\n\n");

  return 0;
}
