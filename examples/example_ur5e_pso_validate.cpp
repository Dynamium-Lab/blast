// example_ur5e_pso_validate.cpp
//
// Takes PSO's raw answer (as PSO produced it -- no SLSQP refinement at all)
// and judges it against the exact same two-tier validation blast's own
// SLSQP path (optimize_with_segments_impl) holds itself to before it will
// call a result a "success":
//   - coarse check: constraints evaluated at the same 76-point/19-segment
//     resolution the optimizer itself uses.
//   - dense check: constraints re-evaluated at a much finer resolution
//     (~1 sample per millisecond of trajectory), to catch violations that
//     occur between the coarse check's sample points.
// SLSQP is required to pass both. This asks whether PSO's answer would
// clear the same bar -- PSO's own GPU kernel only ever does the coarse-
// equivalent check (once, in-kernel, non-segment-grouped), so this has
// never been checked before.
//
// Reads PSO's raw output vector (one number per line, x_size values, last
// value = duration T) from PSO_X_PATH.
//
// Build:
//   cmake --build build --target example_ur5e_pso_validate
// Run:
//   ./build/examples/example_ur5e_pso_validate [path/to/pso_x.txt]

#define BLAST_USE_NATIVE_SQP

#include <blast>
#include <fstream>
#include <iostream>

static const char* DEFAULT_PSO_X_PATH = "/home/sannan/Downloads/blast-gpu/pso_result.txt";

int main(int argc, char** argv) {
  using namespace blast;

  const char* PSO_X_PATH = argc > 1 ? argv[1] : DEFAULT_PSO_X_PATH;

  Manipulator ur5e = make_UR5e();

  // Defaults to ICRA task 0, the task this project was tuned against.
  // Optionally overridden by 12 extra arguments (6 start joint positions then
  // 6 end joint positions) so results for any of the benchmark's 64 tasks can
  // be validated:
  //     example_ur5e_pso_validate <result_file> s0..s5 e0..e5
  Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  if (argc >= 14) {
    for (u32 j = 0; j < 6; j++) {
      start[j] = (real) atof(argv[2 + j]);
      end[j]   = (real) atof(argv[8 + j]);
    }
    std::cout << "(task overridden from arguments)\n";
  }
  Task task = Task::stop_to_stop(start, end);

  Optimization opt(ur5e, task);
  opt.bspline = Bspline(24, 76, 5, ur5e.n_joints); // matches PSO exactly: ncontrol=24, points=76, degree=5

  // Same scope as PSO's GPU-side cuUR5eMultiTraj: position/velocity/
  // acceleration, not torque.
  opt.constraints.position     = true;
  opt.constraints.velocity     = true;
  opt.constraints.acceleration = true;
  opt.constraints.torque       = false;
  opt.constraints.tool_speed   = false;

  opt.success_tolerance = 0.01;
  opt.method             = OptimizationMethod::with_segments;

  const u32 n = (u32) opt.bspline.x_len(opt.task);

  Array x(n);
  {
    std::ifstream in(PSO_X_PATH);
    if (!in) {
      std::cerr << "Could not open " << PSO_X_PATH << "\n";
      return 1;
    }
    for (u32 i = 0; i < n; i++) {
      if (!(in >> x[i])) {
        std::cerr << "input file has fewer than " << n << " values (expected " << n << ")\n";
        return 1;
      }
    }
  }
  std::cout << "Read PSO answer (" << n << " values), T=" << x.back() << "\n\n";

  n_con_with_segments(&opt);

  // ---- Coarse check: same resolution PSO itself checks against ----------
  real coarse_max;
  int  coarse_idx;
  {
    Array  constraints_points(opt.constraints.n_constraints);
    Matrix gradient; // empty -- skips gradient computation, values only
    constraints_and_gradients_with_segments(x, opt, constraints_points, gradient);
    coarse_max = max(constraints_points);
    coarse_idx = (int) argmax(constraints_points);
  }

  // ---- Dense check: same "more points" validation SLSQP is held to ------
  real dense_max;
  int  dense_idx;
  {
    const u32 output_steps_ms = 1;
    u64       steps_ms        = (u64) std::ceil(x.back() * 1e3 / output_steps_ms);
    int       points_more     = (int) (steps_ms + 1);

    Bspline bspline_val_more(opt.bspline.n_ctrl, points_more, opt.bspline.degree, opt.manip.n_joints);
    bspline_val_more.compute_trajectory(x, opt.task);
    Optimization opt_val_more(opt);
    opt_val_more.set_bspline(bspline_val_more);
    n_con_with_segments(&opt_val_more);

    Array  constraints_more_points(opt_val_more.constraints.n_constraints);
    Matrix gradient;
    constraints_and_gradients_with_segments(x, opt_val_more, constraints_more_points, gradient);
    dense_max = max(constraints_more_points);
    dense_idx = (int) argmax(constraints_more_points);

    std::cout << "(dense check used " << points_more << " points, vs. " << opt.bspline.n_points << " coarse)\n\n";
  }

  const bool is_valid      = coarse_max < opt.success_tolerance * 2;
  const bool is_valid_more = dense_max < opt.success_tolerance * 2;

  std::cout << "Coarse (76pt):  max_constraint=" << coarse_max << "  idx=" << coarse_idx
            << "  " << (is_valid ? "PASS" : "FAIL") << "\n";
  std::cout << "Dense  (fine):  max_constraint=" << dense_max << "  idx=" << dense_idx
            << "  " << (is_valid_more ? "PASS" : "FAIL") << "\n\n";

  std::cout << "Would count as SLSQP-style success: " << ((is_valid && is_valid_more) ? "YES" : "NO") << "\n";

  return 0;
}
