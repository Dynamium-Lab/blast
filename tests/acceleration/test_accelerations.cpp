#include <blast>
#include "test_helper/test_helper.hpp"

#define BLAST_TRACE_LEVEL 3

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"

#include <fcntl.h> // _O_WRONLY
#include <io.h>    // _open, _dup, _dup2, _close

struct IOSilencer {
  int             saved_stdout_fd = -1;
  int             null_fd         = -1;
  std::streambuf* old_cout        = nullptr;
  std::streambuf* old_cerr        = nullptr;
  std::ofstream   null_stream;

  IOSilencer() {
    // Flush buffers
    std::cout.flush();
    std::cerr.flush();
    fflush(stdout);

    // Redirect std::cout and std::cerr to /dev/null
    null_stream.open("NUL"); // Windows null device
    old_cout = std::cout.rdbuf();
    old_cerr = std::cerr.rdbuf();
    std::cout.rdbuf(null_stream.rdbuf());
    std::cerr.rdbuf(null_stream.rdbuf());

    // Redirect printf (C stdout)
    saved_stdout_fd = _dup(_fileno(stdout));
    null_fd         = _open("NUL", _O_WRONLY);
    _dup2(null_fd, _fileno(stdout));
  }

  ~IOSilencer() {
    std::cout.rdbuf(old_cout);
    std::cerr.rdbuf(old_cerr);
    fflush(stdout);

    _dup2(saved_stdout_fd, _fileno(stdout));
    _close(null_fd);
    _close(saved_stdout_fd);
  }
};

namespace blast {

TEST_CASE("Test gradient", "[accelerations]") {

  ManipulatorTempData manip_data;

  // Optimization opt(get_generic_gen3(), get_gen3_task());
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 100, 5, 6);
  opt.bspline = bspline;

  opt.guess.type = Guess::custom;
  opt.guess.initial_x   = Array{0.152355, -0.597958, -0.756348, 0.558375, -0.946557, -0.10268, 0.901938, -0.145201, -0.436917, -0.827993,
                       -0.210473, -0.0664338, 0.370776, -0.29955, -0.139142, 0.757138, -0.157328, 0.909626, -0.214084, 0.745555,
                       -0.0416009, 0.375275, 0.126153, 0.402337, 0.723458, -0.070362, -0.444889, 0.784517, 0.655188, 0.216073,
                       0.889643, -0.969162, -0.0279569, 0.496376, -0.891977, 0.0247388, -0.332795, -0.239096, -0.418835, 0.922935,
                       0.124817, 0.34831, -0.00905602, -0.157159, 0.492551, 0.0325292, -0.299891, 0.828649, -0.909505, -0.294613,
                       -0.768063, -0.489629, 0.0918699, -0.603901, -0.391773, -0.919239, 0.472415, 0.591886, -0.494069, 0.0488416,
                       4.81184};
  //   opt.guess.initial_x   = Array(opt.bspline.x_len(opt.task), 2.0);
  // opt.guess.n_random_shots = 100;

  opt.constraints.show_info = true;

  opt.constraints.position     = true;
  opt.constraints.velocity     = true;
  opt.constraints.acceleration = true;
  opt.constraints.torque       = true;

  opt.constraints.tcp_speed           = false;
  opt.constraints.self_collisions     = false;
  opt.constraints.external_collisions = false;

  opt.max_tries         = 1;
  opt.success_tolerance = 0.01;


  auto results = optimize(&opt, OptimizationMethod::baseline);

  Optimization opt_acc(get_generic_Link6(), get_link6_task());
  opt_acc.set_bspline(bspline);
  opt_acc.guess.type = Guess::custom;
  opt_acc.guess.initial_x   = Array{0.152355, -0.597958, -0.756348, 0.558375, -0.946557, -0.10268, 0.901938, -0.145201, -0.436917, -0.827993,
                           -0.210473, -0.0664338, 0.370776, -0.29955, -0.139142, 0.757138, -0.157328, 0.909626, -0.214084, 0.745555,
                           -0.0416009, 0.375275, 0.126153, 0.402337, 0.723458, -0.070362, -0.444889, 0.784517, 0.655188, 0.216073,
                           0.889643, -0.969162, -0.0279569, 0.496376, -0.891977, 0.0247388, -0.332795, -0.239096, -0.418835, 0.922935,
                           0.124817, 0.34831, -0.00905602, -0.157159, 0.492551, 0.0325292, -0.299891, 0.828649, -0.909505, -0.294613,
                           -0.768063, -0.489629, 0.0918699, -0.603901, -0.391773, -0.919239, 0.472415, 0.591886, -0.494069, 0.0488416,
                           4.81184};

  opt_acc.constraints.show_info = true;

  opt_acc.constraints.position     = true;
  opt_acc.constraints.velocity     = true;
  opt_acc.constraints.acceleration = true;
  opt_acc.constraints.torque       = true;

  opt_acc.constraints.tcp_speed           = false;
  opt_acc.constraints.self_collisions     = false;
  opt_acc.constraints.external_collisions = false;

  opt_acc.max_tries         = 1;
  opt_acc.success_tolerance = 0.01;

  auto result_acc = optimize(&opt_acc, OptimizationMethod::with_analytical_pva);

  // ------------- Constraint function -------------------------
  auto joints = opt_acc.manip.n_joints;

  Array pos_plus(joints);
  Array vel_plus(joints);
  Array acc_plus(joints);
  Array torque_constraint(joints);
  Array torque_constraint_plus(joints);

  Matrix gradient_torque_coeffs(opt_acc.manip.n_joints, 3 * opt_acc.manip.n_joints * opt_acc.bspline.n_points);

  constexpr real eps = 1e-5;

  initialize_optimization(&opt_acc);
  n_con(&opt_acc);

  auto x = opt.guess.initial_x;

  opt_acc.bspline.compute_trajectory(x, opt_acc.task);

  u32 i = 1;

  auto pos = opt_acc.bspline.traj.pos.col(i);

  if (opt_acc.constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Dynamics");
#endif
    auto vel = opt_acc.bspline.traj.vel.col(i);
    auto acc = opt_acc.bspline.traj.acc.col(i);
    dynamics(opt_acc.manip, manip_data, vel, acc);
    for (int j = 0; j < joints; j++) {
      torque_constraint[j] = abs_constraint(1.01 * manip_data.efforts[j], opt_acc.manip.torque_max[j]);
    }

    // grad_coeffs pos
    for (u32 j = 0; j < joints; j++) {
      pos[j] += eps;
      forward_kinematics(opt_acc.manip, manip_data, pos);
      dynamics(opt_acc.manip, manip_data, vel, acc);
      pos[j] -= eps;

      for (u32 k = 0; k < joints; k++) {
        torque_constraint_plus[k]                     = abs_constraint(1.01 * manip_data.efforts[k], opt_acc.manip.torque_max[k]);
        gradient_torque_coeffs(k, j + 3 * joints * i) = (torque_constraint_plus[k] - torque_constraint[k]) / eps;
      }
    }

    forward_kinematics(opt_acc.manip, manip_data, pos);
    // grad_coeffs vel
    for (u32 j = 0; j < joints; j++) {
      vel[j] += eps;
      dynamics(opt_acc.manip, manip_data, vel, acc);
      vel[j] -= eps;

      for (u32 k = 0; k < joints; k++) {
        torque_constraint_plus[k]                              = abs_constraint(1.01 * manip_data.efforts[k], opt_acc.manip.torque_max[k]);
        gradient_torque_coeffs(k, j + joints + 3 * joints * i) = (torque_constraint_plus[k] - torque_constraint[k]) / eps;
      }
    }

    // grad_coeffs acc
    for (u32 j = 0; j < joints; j++) {
      acc[j] += eps;
      dynamics(opt_acc.manip, manip_data, vel, acc);
      acc[j] -= eps;

      for (u32 k = 0; k < joints; k++) {
        torque_constraint_plus[k]                                  = abs_constraint(1.01 * manip_data.efforts[k], opt_acc.manip.torque_max[k]);
        gradient_torque_coeffs(k, j + 2 * joints + 3 * joints * i) = (torque_constraint_plus[k] - torque_constraint[k]) / eps;
      }
    }
  }

  const auto xlen = opt_acc.bspline.x_len(opt_acc.task);

  u32   n_con_lb = 0;
  Array x_plus(xlen);

  u32 x_per_joint = (xlen - 1) / opt_acc.manip.n_joints; // = nctrl - 6 (skip first and last 3)
  u32 joint       = 0;

  u32 grad_idx       = 0;
  u32 constraint_idx = 0;
  u32 x_idx          = 0;

  u32 n_torque         = 0;
  u32 n_tcp_speed      = 0;
  u32 n_self_collision = 0;

  u32   n_con = opt_acc.constraints.n_constraints;
  Array grad(xlen * n_con);

  memcpy(x_plus.data, x.data, xlen * sizeof(real)); // todo: is this necessary ? done once, and x_plus is modified to original at the end of each iteration

  for (u32 j = 0; j < xlen - 1; j++) {              // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
    // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
    joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
    x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

    x_plus[j] += eps;

    // Array xv_plus;
    // xv_plus.alias(x_plus.data, xlen);
    opt_acc.bspline.compute_trajectory(x_plus, opt_acc.task);

    n_con_lb = ncon_lb_acc(&opt_acc, x_idx); // find the amount of constraints before the current point

    Array external_collisions(opt_acc.bspline.upper_bounds[x_idx] - opt_acc.bspline.lower_bounds[x_idx]);
    // todo: create alias matrix that points to grad
    // todo: can we change the order in which we store the gradients ?
    grad_idx       = n_con_lb * xlen + j;                                          // gradients are stored column-wise xlen * npoints
    constraint_idx = n_con_lb;
    for (u32 i = opt_acc.bspline.lower_bounds[x_idx]; i <= opt_acc.bspline.upper_bounds[x_idx]; i++) { // lb & ub are inclusive
      grad_idx += joint * xlen;
      constraint_idx += joint;
      n_torque         = 0;
      n_tcp_speed      = 0;
      n_self_collision = 0;
      if (opt_acc.constraints.torque) { // todo: add analytical gradients for torque
        n_torque = opt_acc.manip.n_joints;

        for (u32 k = 0; k < opt_acc.manip.n_joints; k++) {
          grad[grad_idx] = opt_acc.bspline.basis_p(x_idx, i) * gradient_torque_coeffs(k, joint + 3 * opt_acc.manip.n_joints * i) +
                           opt_acc.bspline.basis_v(x_idx, i) / opt_acc.bspline.traj.t.back() * gradient_torque_coeffs(k, joint + opt_acc.manip.n_joints + 3 * opt_acc.manip.n_joints * i) +
                           opt_acc.bspline.basis_a(x_idx, i) / (opt_acc.bspline.traj.t.back() * opt_acc.bspline.traj.t.back()) * gradient_torque_coeffs(k, joint + 2 * opt_acc.manip.n_joints + 3 * opt_acc.manip.n_joints * i);
          grad_idx += xlen;
          constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }
      }
    }
  }

  CHECK(grad_idx == 100);
};

TEST_CASE("Test gradient accuracy", "[accelerations]") {

  ManipulatorTempData manip_data;

  // Optimization opt(get_generic_gen3(), get_gen3_task());
  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 100, 5, 6);
  opt.bspline = bspline;

  opt.guess.type = Guess::custom;
  opt.guess.initial_x   = Array{0.152355, -0.597958, -0.756348, 0.558375, -0.946557, -0.10268, 0.901938, -0.145201, -0.436917, -0.827993,
                       -0.210473, -0.0664338, 0.370776, -0.29955, -0.139142, 0.757138, -0.157328, 0.909626, -0.214084, 0.745555,
                       -0.0416009, 0.375275, 0.126153, 0.402337, 0.723458, -0.070362, -0.444889, 0.784517, 0.655188, 0.216073,
                       0.889643, -0.969162, -0.0279569, 0.496376, -0.891977, 0.0247388, -0.332795, -0.239096, -0.418835, 0.922935,
                       0.124817, 0.34831, -0.00905602, -0.157159, 0.492551, 0.0325292, -0.299891, 0.828649, -0.909505, -0.294613,
                       -0.768063, -0.489629, 0.0918699, -0.603901, -0.391773, -0.919239, 0.472415, 0.591886, -0.494069, 0.0488416,
                       4.81184};
  //   opt.guess.initial_x   = Array(opt.bspline.x_len(opt.task), 2.0);
  // opt.guess.n_random_shots = 100;

  opt.constraints.show_info = true;

  opt.constraints.position     = true;
  opt.constraints.velocity     = true;
  opt.constraints.acceleration = true;
  opt.constraints.torque       = true;

  opt.constraints.tcp_speed           = false;
  opt.constraints.self_collisions     = false;
  opt.constraints.external_collisions = false;

  opt.max_tries         = 1;
  opt.success_tolerance = 0.01;


  auto results = optimize(&opt, OptimizationMethod::baseline);

  Optimization opt_acc(get_generic_Link6(), get_link6_task());
  opt_acc.set_bspline(bspline);
  opt_acc.guess.type = Guess::custom;
  opt_acc.guess.initial_x   = Array{0.152355, -0.597958, -0.756348, 0.558375, -0.946557, -0.10268, 0.901938, -0.145201, -0.436917, -0.827993,
                           -0.210473, -0.0664338, 0.370776, -0.29955, -0.139142, 0.757138, -0.157328, 0.909626, -0.214084, 0.745555,
                           -0.0416009, 0.375275, 0.126153, 0.402337, 0.723458, -0.070362, -0.444889, 0.784517, 0.655188, 0.216073,
                           0.889643, -0.969162, -0.0279569, 0.496376, -0.891977, 0.0247388, -0.332795, -0.239096, -0.418835, 0.922935,
                           0.124817, 0.34831, -0.00905602, -0.157159, 0.492551, 0.0325292, -0.299891, 0.828649, -0.909505, -0.294613,
                           -0.768063, -0.489629, 0.0918699, -0.603901, -0.391773, -0.919239, 0.472415, 0.591886, -0.494069, 0.0488416,
                           4.81184};

  opt_acc.constraints.show_info = true;

  opt_acc.constraints.position     = true;
  opt_acc.constraints.velocity     = true;
  opt_acc.constraints.acceleration = true;
  opt_acc.constraints.torque       = true;

  opt_acc.constraints.tcp_speed           = false;
  opt_acc.constraints.self_collisions     = false;
  opt_acc.constraints.external_collisions = false;

  opt_acc.max_tries         = 1;
  opt_acc.success_tolerance = 0.01;

  auto result_acc = optimize(&opt_acc, OptimizationMethod::with_analytical_pva);

  //   Optimization opt_acc_2(get_generic_Link6(), get_link6_task());
  //   opt_acc_2.set_bspline(bspline);
  //   opt_acc_2.guess.type = Guess::custom;
  //   opt_acc_2.guess.initial_x   = Array{0.152355, -0.597958, -0.756348, 0.558375, -0.946557, -0.10268, 0.901938, -0.145201, -0.436917, -0.827993,
  //                              -0.210473, -0.0664338, 0.370776, -0.29955, -0.139142, 0.757138, -0.157328, 0.909626, -0.214084, 0.745555,
  //                              -0.0416009, 0.375275, 0.126153, 0.402337, 0.723458, -0.070362, -0.444889, 0.784517, 0.655188, 0.216073,
  //                              0.889643, -0.969162, -0.0279569, 0.496376, -0.891977, 0.0247388, -0.332795, -0.239096, -0.418835, 0.922935,
  //                              0.124817, 0.34831, -0.00905602, -0.157159, 0.492551, 0.0325292, -0.299891, 0.828649, -0.909505, -0.294613,
  //                              -0.768063, -0.489629, 0.0918699, -0.603901, -0.391773, -0.919239, 0.472415, 0.591886, -0.494069, 0.0488416,
  //                              4.81184};

  //   opt_acc_2.constraints.show_info = true;

  //   opt_acc_2.constraints.position     = true;
  //   opt_acc_2.constraints.velocity     = true;
  //   opt_acc_2.constraints.acceleration = true;
  //   opt_acc_2.constraints.torque       = true;

  //   opt_acc_2.constraints.tcp_speed           = false;
  //   opt_acc_2.constraints.self_collisions     = false;
  //   opt_acc_2.constraints.external_collisions = false;

  //   opt_acc_2.max_tries         = 1;
  //   opt_acc_2.success_tolerance = 0.01;

  //   auto result_acc_2 = optimize(&opt_acc_2, OptimizationMethod::with_analytical_dynamics);

  CHECK(std::abs(results.x.back() - result_acc.x.back()) < 1e-5);

  u32 i = 0;
  for (u32 j = 0; j < results.opt->constraints.grad_list[i].rows; j++) {
    for (u32 k = 0; k < results.opt->constraints.grad_list[i].cols; k++) {
      CHECK(std::abs(results.opt->constraints.grad_list[i](j, k) - result_acc.opt->constraints.grad_list[i](j, k)) < 0.01);
      CHECK(std::abs(results.opt->constraints.constr_list[i][k] - result_acc.opt->constraints.constr_list[i][k]) < 0.01);
      //   CHECK(std::abs(results.opt->constraints.grad_list[i](j, k) - result_acc_2.opt->constraints.grad_list[i](j, k)) < 0.01);
      //   CHECK(std::abs(results.opt->constraints.constr_list[i][k] - result_acc_2.opt->constraints.constr_list[i][k]) < 0.01);
    }
  }
}

} // namespace blast
