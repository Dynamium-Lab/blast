//
// Created by nikos on 2025-09-05.
//
#include <blast>
#include "../tests/test_helper/test_helper.hpp"

#define BLAST_TRACE_LEVEL 0
// #define BLAST_USE_NATIVE_SQP

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"

#include <fcntl.h> // _O_WRONLY
#include <io.h>    // _open, _dup, _dup2, _close

using namespace blast;

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

struct Config {
  int    n_ctrl   = 0;
  int    n_points = 0;
  int    task_idx = 0;
  int    n_optim  = 0;
  Matrix task;
};

// We create a config list that will be used for all benchmarks
inline void fill_config_list(std::array<Config, 13>& config_list) {

  std::array<std::tuple<int, int>, 1> bspline_config_list = {std::make_tuple(16, 110)}; // These give rounded number of n_points_per_segment
  std::vector<Matrix>                 tasks               = get_Link6_demo1_tasks();    // 13 tasks

  int config_idx = 0;
  int task_id    = 0;

  for (auto& task: tasks) {
    for (auto& bspline_config: bspline_config_list) {
      auto [n_ctrl, n_points]          = bspline_config;
      config_list[config_idx].n_ctrl   = n_ctrl;
      config_list[config_idx].n_points = n_points;
      config_list[config_idx].task_idx = task_id;
      config_list[config_idx].n_optim  = 1;
      config_list[config_idx++].task   = task;
    }
    task_id++;
  }
}

// ------------------------- Accelerated functions --------------------------------
//
// inline u32 ncon_acc(const Optimization* opt, const int x_idx) {
//   const u32 n_points            = opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]; // + 1 since ub is inclusive
//   const u32 n_joints            = opt->manip.n_joints;
//   const u32 n_constraints_basic = n_points * n_joints;
//   u32       n_constraints       = 0;
//
//   if (opt->constraints.position)
//     n_constraints += n_constraints_basic;
//   if (opt->constraints.velocity)
//     n_constraints += n_constraints_basic;
//   if (opt->constraints.acceleration)
//     n_constraints += n_constraints_basic;
//   if (opt->constraints.torque)
//     n_constraints += n_constraints_basic;
//
//   if (opt->constraints.tcp_speed)
//     n_constraints += n_points;
//
//   if (opt->constraints.external_collisions)
//     n_constraints += opt->constraints.n_collision_constraints;
//   if (opt->constraints.self_collisions) {
//     n_constraints += opt->manip._n_internal_collisions * n_points;
//   }
//
//   for (auto& n_con: opt->constraints.n_extra_constraints)
//     n_constraints += n_con;
//
//   return n_constraints;
// }
//
// inline void compute_constraints_grad1(double* result, const Array& x, const u32 x_idx, Optimization* opt) {
//   double* moving_result = result;
//
//   opt->bspline.compute_trajectory(x, opt->task);
//
//   ObjMatrix<Capsule> capsules(opt->manip._n_caps, (u32) (opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx] + 1) / opt->constraints.n_collision_skip);
//   for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) {
// #if BLAST_TRACE_LEVEL >= 3
//     ZoneScopedN("ConstraintSinglePoint");
// #endif
//
//     // todo: reset frame???
//     // rotations_computed -> bool
//     // forward_kinematics_computed -> bool
//
//     ManipulatorTempData manip_data;
//
//
//     // This is calculated every loop
//     auto pos = opt->bspline.traj.pos.col(i);
//     Assert(pos.is_alias);
//
//     {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("Kinematics");
// #endif
//       forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
//     }
//
//     {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("Capsules");
// #endif
//       // todo: this cleaner
//       if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
//         compute_capsules(opt->manip, manip_data);
//         if (opt->constraints.external_collisions) {
//           if (i % opt->constraints.n_collision_skip == 0 && i != opt->bspline.n_points - 1) {
//             for (int j = 0; j < opt->manip._n_caps; j++) {
//               capsules(j, i / (int) opt->constraints.n_collision_skip) = manip_data.capsule_list[j];
//             }
//           }
//         }
//       }
//     }
//
//     if (opt->constraints.position) {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("Position");
// #endif
//       for (int j = 0; j < opt->manip.n_joints; j++) {
//         moving_result[0] = bound_constraint(opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
//         moving_result++;
//       }
//     }
//
//     if (opt->constraints.velocity) {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("Velocity");
// #endif
//       for (int j = 0; j < opt->manip.n_joints; j++) {
//         // moving_result[0] = abs_constraint(opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
//         moving_result[0] = std::abs(opt->bspline.traj.vel(j, i)) / opt->manip.vmax[j] - 1.0;
//         moving_result++;
//       }
//     }
//
//     if (opt->constraints.acceleration) {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("Acceleration");
// #endif
//       for (int j = 0; j < opt->manip.n_joints; j++) {
//         moving_result[0] = std::abs(opt->bspline.traj.acc(j, i)) / opt->manip.amax[j] - 1.0;
//         moving_result++;
//       }
//     }
//
//     if (opt->constraints.torque) {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("Dynamics");
// #endif
//       auto vel = opt->bspline.traj.vel.col(i);
//       auto acc = opt->bspline.traj.acc.col(i);
//       dynamics(opt->manip, manip_data, vel, acc); // fills _efforts
//
//       for (int j = 0; j < opt->manip.n_joints; j++) {
//         moving_result[0] = std::abs(manip_data.efforts[j]) / opt->manip.tau_max[j] - 1.0;
//         moving_result++;
//       }
//     }
//
//     if (opt->constraints.tcp_speed) {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("TCPSpeed");
// #endif
//       auto J_tool      = get_J_tool(opt, manip_data);
//       real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
//       moving_result[0] = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
//       moving_result++;
//     }
//
//     if (opt->constraints.self_collisions) {
// #if BLAST_TRACE_LEVEL >= 3
//       ZoneScopedN("SelfCollisions");
// #endif
//       auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
//       for (u32 j = 0; j < tmp_coll.size; j++)
//         moving_result[j] = -tmp_coll[j] + 0.01; //*std::abs(tmp_coll[j]);
//       moving_result += tmp_coll.size;
//     }
//   }
//
//   if (opt->constraints.external_collisions) {
// #if BLAST_TRACE_LEVEL >= 3
//     ZoneScopedN("ExternalCollisions");
// #endif
//     auto collisions = test_collisions(capsules, &(opt->world), opt->constraints.n_collision_constraints, opt->trajectory_start_time, opt->trajectory_start_time + x.back());
//     for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
//       moving_result[j] = -collisions[j] + 0.01; //*std::abs(collisions[j]);
//     }
//     moving_result += opt->constraints.n_collision_constraints;
//   }
//
//   {
// #if BLAST_TRACE_LEVEL >= 3
//     ZoneScopedN("CustomConstraints");
// #endif
//     for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) { // todo: split in extra_constraint per point or once
//       auto extra_constraint = opt->constraints.extra_constraints[i];
//       extra_constraint(moving_result, opt);
//       moving_result += opt->constraints.n_extra_constraints[i];
//     }
//   }
// }
//
// inline void nlopt_constraints_acc1(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
// #if BLAST_TRACE_LEVEL >= 1
//   ZoneScoped;
// #endif
//   auto* opt = (Optimization*) f_data;
//
//   Array xv;
//   xv.alias(x, x_len);
//   compute_constraints(result, xv, opt);
//
//   // gradients calculation
//   if (grad) {
//     memset(grad, 0, m * x_len * sizeof(real)); // zeros grad, since grad originally starts with -6e+66 ...
//
//     constexpr real eps = 1e-5;
//     Array          x_plus(x_len);
//
//     // lb & ub automatically calculated by bsplines
//     u32 x_per_joint = (x_len - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
//     u32 mult        = 0;
//
//     for (u32 j = 0; j < x_len - 1; j++) {     // last one is T
//       mult      = j / x_per_joint > mult ? mult + 1 : mult;
//       u32 x_idx = j - mult * x_per_joint + 3; // todo: fix for NaN in PVA of task
//
//       memcpy(x_plus.data, x, x_len * sizeof(real));
//       x_plus[j] += eps;
//
//       Array xv_plus;
//       xv_plus.alias(x_plus.data, x_len);
//
//       // compute constraints
//       auto  n_con = ncon_acc(opt, x_idx);
//       Array r_plus(n_con);
//       compute_constraints_grad1(r_plus.data, x_plus, x_idx, opt);
//       auto n_con_lb = ncon_lb_acc(opt, x_idx);
//
//       // todo: fix this loop
//       // for (u32 i = 0; i < m; i++)
//       for (u32 i = 0; i < n_con; i++)
//         grad[(n_con_lb + i) * x_len + j] = (r_plus[i] - result[n_con_lb + i]) / eps;
//     }
//     {
//       // last point T
//       u32 j = x_len - 1;
//       memcpy(x_plus.data, x, x_len * sizeof(real));
//       x_plus[j] += eps;
//
//       Array xv_plus;
//       xv_plus.alias(x_plus.data, x_len);
//
//       // compute constraints
//       Array r_plus(m);
//       compute_constraints(r_plus.data, x_plus, opt);
//
//       for (u32 i = 0; i < m; i++)
//         grad[i * x_len + j] = (r_plus[i] - result[i]) / eps;
//     }
//   }
// }
// inline blast_fn void constraints_and_gradients_acc1(const Array& x, Optimization& opt, Array& constraints, Matrix& grad) {
//   // constraints (p,v,a,tor) for each joint, for each segment
//   // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
//
//   // basis: n_ctrl x n_points
//   Assert(constraints.is_alias);
//   if (grad.size) {
//     Assert(grad.is_alias);
//     Assert(grad.rows == x.size);
//     Assert(grad.cols == constraints.size);
//   }
//
//   const int n_segments              = (int) opt.bspline.n_ctrl - (int) opt.bspline.p;
//   const int n_points_per_segment    = (int) opt.bspline.n_points / n_segments; // todo: check if fine? -> (Nikos) constructor that sets n_points from n_points_per_segment
//   const int n_joints                = (int) opt.manip.n_joints;
//   const int n_ctrl                  = (int) opt.bspline.n_ctrl;
//   const int x_len                   = (int) x.size;
//   const int n_constraints_per_point = (n_joints * 4); //  todo: remove hard-coded 4
//   Assert(constraints.size == n_segments * n_constraints_per_point * n_points_per_segment);
//
//   // limits
//   auto pmax    = opt.manip.pmax;
//   auto pmin    = opt.manip.pmin;
//   auto vmax    = opt.manip.vmax;
//   auto amax    = opt.manip.amax;
//   auto tau_max = opt.manip.tau_max;
//
//   for (int j = 0; j < n_joints; j++) {
//     // todo (andre): do this when initializing the optimization, not here
//     // todo: document the current behaviour in the API
//     //        (doesn't currently work if one is inf and the other is not)
//     if (pmax[j] == INF_REAL)  // note: replace INF_REAL with huge value
//       pmax[j] = 1e300;
//     if (pmin[j] == -INF_REAL) // note: replace -INF_REAL with huge negative value
//       pmin[j] = -1e300;
//   }
//
//   ManipulatorTempData manip_data;
//
//   opt.bspline.compute_trajectory(x, opt.task);
//
//
//   for (int segment = 0; segment < n_segments; segment++) {
// #if BLAST_TRACE_LEVEL >= 2
//     ZoneScopedN("All Segment Constraints");
// #endif
//     const int first_affected_control_point = std::max(3, segment);
//     const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.p);
//     const int n_affected_control_points    = last_affected_control_point - first_affected_control_point + 1; // note: affected_control_points are inclusive, so when we have last = 5, first = 3, we want 3 (5 - 3) + 1
//     const int start_point_for_segment      = segment * n_points_per_segment;                                 // note:
//     Assert(n_affected_control_points >= 3);
//     Assert(n_affected_control_points <= 6);
//
//     Matrix bp(&opt.bspline.basis_p(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//     Matrix bv(&opt.bspline.basis_v(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//     Matrix ba(&opt.bspline.basis_a(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//
//     // Initialize constraints with -INF_REAL
//     Array pos_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array vel_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array acc_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array tor_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//
//     for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) { // note:
//       auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment);              // note:
//       auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment);              // note:
//       auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment);              // note:
//
//       forward_kinematics(opt.manip, manip_data, p);
//       dynamics(opt.manip, manip_data, v, a);
//
//       for (int j = 0; j < n_joints; j++) {
//         // position
//         pos_constraints[j + n_joints * point_in_segment] = bound_constraint(p[j], pmin[j], pmax[j]);
//         // velocity
//         vel_constraints[j + n_joints * point_in_segment] = std::abs(v[j]) / vmax[j] - 1.0;
//         // acceleration
//         acc_constraints[j + n_joints * point_in_segment] = std::abs(a[j]) / amax[j] - 1.0;
//         // torque
//         tor_constraints[j + n_joints * point_in_segment] = std::abs(manip_data.efforts[j]) / tau_max[j] - 1.0;
//       }
//     }
//
//     // at this point we have max constraints for pos, vel, acc, tor, for this segment
//
//     // fill in the constraints for the current segment
//     // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
//     auto fill_idx = segment * n_constraints_per_point * n_points_per_segment;
//     for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
//
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = pos_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = vel_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = acc_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = tor_constraints[point_in_segment * n_joints + j];
//     }
//
//     // The gradient should be a (x_len)x(n_constraints) matrix that looks like this:
//     // [dp1/dx1, dp2/dx1, ..., dv1/dx1, dv2/dx1, ..., da1/dx1, da2/dx1, ..., dt1/dx1, dt2/dx1]
//     // [dp1/dx2, dp2/dx2, ..., dv1/dx2, dv2/dx2, ..., da1/dx2, da2/dx2, ..., dt1/dx2, dt2/dx2]
//     // [dp1/dx3, dp2/dx3, ..., dv1/dx3, dv2/dx3, ..., da1/dx3, da2/dx3, ..., dt1/dx3, dt2/dx3]
//     // [dp1/dx4, dp2/dx4, ..., dv1/dx4, dv2/dx4, ..., da1/dx4, da2/dx4, ..., dt1/dx4, dt2/dx4]
//     // [dp1/dx5, dp2/dx5, ..., dv1/dx5, dv2/dx5, ..., da1/dx5, da2/dx5, ..., dt1/dx5, dt2/dx5]
//     // [dp1/dx6, dp2/dx6, ..., dv1/dx6, dv2/dx6, ..., da1/dx6, da2/dx6, ..., dt1/dx6, dt2/dx6]
//     // [.....................]
//     // [.....................]
//     // [.....................]
//     // [.....................]
//     // [dp1/dT=0,dp2/dT=0,..., dv1/dT , dv1/dT , ..., da1/dT , da2/dT , ..., dt1/dT , dt2/dT ]
//     // *** per segment ***
//     // where x is the optimization vector
//     if (grad.size) {
//       // Matrix (alias) in which we can insert the gradient for the current segment
//       Matrix grad_matrix(&grad(0, segment * n_constraints_per_point * n_points_per_segment), x_len, n_constraints_per_point * n_points_per_segment);
//       Assert(grad_matrix.is_alias);
//
//       int            con = 0;
//       constexpr real eps = 1e-5;
//
//       for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
//         Matrix dtau_dx(x_len, n_joints); // we can store and reuse some of the values
//         Matrix dp_dx(x_len, n_joints);
//         Matrix dv_dx(x_len, n_joints);
//         Matrix da_dx(x_len, n_joints);
//         {
//           auto x_idx      = first_affected_control_point - 3;
//           auto x_idx_skip = n_ctrl - 6;
//
//           for (int j = 0; j < n_joints; j++) {                    // These loops go from x_idx -> x_len - 2
//             for (int i = 0; i < n_affected_control_points; i++) { // ^^
//               x.data[x_idx + i] += eps;
//               opt.bspline.compute_trajectory(x, opt.task);
//               x.data[x_idx + i] -= eps;
//               // position
//               auto p                = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment); // note:
//               auto new_constraint_p = bound_constraint(p[j], pmin[j], pmax[j]);
//               dp_dx(x_idx + i, j)   = (new_constraint_p - pos_constraints[n_joints * point_in_segment + j]) / eps;
//               // velocity
//               auto v                = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment); // note:
//               auto new_constraint_v = std::abs(v[j]) / vmax[j] - 1.0;
//               dv_dx(x_idx + i, j)   = (new_constraint_v - vel_constraints[n_joints * point_in_segment + j]) / eps;
//               // acceleration
//               auto a                = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment); // note:
//               auto new_constraint_a = std::abs(a[j]) / amax[j] - 1.0;
//               da_dx(x_idx + i, j)   = (new_constraint_a - acc_constraints[n_joints * point_in_segment + j]) / eps;
//               // torque
//               forward_kinematics(opt.manip, manip_data, p);
//               dynamics(opt.manip, manip_data, v, a);
//               for (int joint = 0; joint < n_joints; joint++) { // This one inserts the gradient for each joint value (6 columns)
//                 const real new_constraint = std::abs(manip_data.efforts[joint]) / tau_max[joint] - 1;
//                 dtau_dx(x_idx + i, joint) = (new_constraint - tor_constraints[n_joints * point_in_segment + joint]) / eps;
//               }
//             }
//             x_idx += x_idx_skip;
//           }
//         }
//         { // d_dx for x_len -1
//           x.data[x_len - 1] += eps;
//           opt.bspline.compute_trajectory(x, opt.task);
//           x.data[x_len - 1] -= eps;
//
//           auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment); // note:
//           auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment); // note:
//           auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment); // note:
//           forward_kinematics(opt.manip, manip_data, p);
//           dynamics(opt.manip, manip_data, v, a);
//           for (int joint = 0; joint < n_joints; joint++) {
//             // position
//             auto new_constraint_p   = bound_constraint(p[joint], pmin[joint], pmax[joint]);
//             dp_dx(x_len - 1, joint) = (new_constraint_p - pos_constraints[n_joints * point_in_segment + joint]) / eps;
//             // velocity
//             auto new_constraint_v   = std::abs(v[joint]) / vmax[joint] - 1.0;
//             dv_dx(x_len - 1, joint) = (new_constraint_v - vel_constraints[n_joints * point_in_segment + joint]) / eps;
//             // acceleration
//             auto new_constraint_a   = std::abs(a[joint]) / amax[joint] - 1.0;
//             da_dx(x_len - 1, joint) = (new_constraint_a - acc_constraints[n_joints * point_in_segment + joint]) / eps;
//             // torque
//             const real new_constraint = std::abs(manip_data.efforts[joint]) / tau_max[joint] - 1;
//             dtau_dx(x_len - 1, joint) = (new_constraint - tor_constraints[n_joints * point_in_segment + joint]) / eps;
//           }
//         }
//         // position
//         for (int joint = 0; joint < n_joints; joint++) {
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           for (int i = 0; i < n_affected_control_points; i++)
//             fill_column[x_idx + i] = dp_dx(x_idx + i, joint);
//
//           // gradient w.r.t. T
//           fill_column.back() = dp_dx(x_len - 1, joint); // note:
//           con++;
//         }
//
//         // velocity
//         for (int joint = 0; joint < n_joints; joint++) {
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + n_joints + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           for (int i = 0; i < n_affected_control_points; i++)
//             fill_column[x_idx + i] = dv_dx(x_idx + i, joint);
//
//           // gradient w.r.t. T
//           fill_column.back() = dv_dx(x_len - 1, joint); // note:
//           con++;
//         }
//
//         // acceleration
//         for (int joint = 0; joint < n_joints; joint++) {
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + 2 * n_joints + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           for (int i = 0; i < n_affected_control_points; i++)
//             fill_column[x_idx + i] = da_dx(x_idx + i, joint);
//
//           // gradient w.r.t. T
//           fill_column.back() = da_dx(x_len - 1, joint); // note:
//           con++;
//         }
//
//         // torque
//         for (int joint = 0; joint < n_joints; joint++) {
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == n_constraints_per_point * point_in_segment + 3 * n_joints + joint); // note: Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           auto x_idx      = first_affected_control_point - 3;
//           auto x_idx_skip = n_ctrl - 6;
//
//           for (int j = 0; j < n_joints; j++) {
//             for (int i = 0; i < n_affected_control_points; i++) {
//               fill_column[x_idx + i] = dtau_dx(x_idx + i, joint);
//             }
//             x_idx += x_idx_skip;
//           }
//
//           // gradient w.r.t. T
//           fill_column.back() = dtau_dx(x_len - 1, joint); // note:
//           con++;
//         }
//       }
//     }
//   }
// }
//
// inline blast_fn void nlopt_constraints_acc1(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
// #if BLAST_TRACE_LEVEL >= 1
//   ZoneScoped;
// #endif
//   auto* opt = (Optimization*) f_data;
//
//   Array xv;
//   xv.alias(x, x_len);
//
//   Array constraints;
//   constraints.alias(result, m);
//
//   if (grad) {
//     memset(grad, 0, m * x_len * sizeof(real));
//     Matrix gradients;
//     gradients.alias(grad, x_len, m);
//     constraints_and_gradients_acc1(xv, *opt, constraints, gradients);
//   } else {
//     // no gradients requested
//     compute_constraints(constraints.data, xv, opt);
//   }
// }
//
// inline Result optimize_acc1(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
//   auto   T1 = get_tick_us();
//   Result result(opt); // todo: this is expensive
//
//   // Initialization
//   // configure_internal_data(opt); // todo: Ensure we can remove
//   initialize_optimization(opt);
//   n_con(opt);
//
//   // Initial validation
//   if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
//     print(opt->task);
//     return result;
//   }
//
//   const auto n = opt->bspline.x_len(opt->task);
//
//   Array con_tol(opt->constraints.n_constraints, 0.01);
//   Array x_tol(n, 0.001);
//
// #ifdef BLAST_USE_NATIVE_SQP
//   nlopt_stopping stop;
//   stop.n          = n;
//   stop.minf_max   = -HUGE_VAL;
//   stop.ftol_rel   = 0;
//   stop.ftol_abs   = 0.001;
//   stop.xtol_rel   = 0;
//   stop.xtol_abs   = x_tol.data;
//   stop.x_weights  = nullptr;
//   stop.nevals_p   = 0;
//   stop.maxeval    = 1000;
//   stop.maxtime    = 5;
//   stop.start      = nlopt_seconds();
//   stop.force_stop = false;
//   stop.stop_msg   = nullptr;
//
//   Array ub(n, INF_REAL);
//   Array lb(n, -INF_REAL);
//   ub.back() = 60.0;
//   lb.back() = 0.01;
//
//   nlopt_constraint fc{};
//   fc.m      = opt->constraints.n_constraints;
//   fc.f      = nullptr;
//   fc.mf     = nlopt_constraints_acc1;
//   fc.pre    = nullptr;
//   fc.f_data = opt;
//   fc.tol    = con_tol.data;
// #else
//   nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
//   nlopt_result nlopt_res;
//   nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc1, opt, con_tol.data);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_maxtime(o, 5);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_maxeval(o, 1000);
//   Assert(nlopt_res == NLOPT_SUCCESS);
// #endif
//
//
//   auto start_guess   = opt->guess; // save for restoration after restarts if necessary
//   int  try_count     = 0;
//   bool is_valid_more = false;
//   bool is_valid      = false;
//   for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
// #if BLAST_TRACE_LEVEL >= 1
//     ZoneScopedN("Optimization");
// #endif
//
//     // initial guess
//     Array x;
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Initial guess");
// #endif
//       x         = init_guess(opt);
//       result.x0 = x;
//     }
//
//     // launch optimization
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("NLopt optimization");
// #endif
//
//       double f = HUGE_VAL;
//       // note: can we initialize grad to 0 here
// #ifdef BLAST_USE_NATIVE_SQP
//       stop.nevals_p              = 0;
//       result.nlopt_exit_criteria = sqp(
//               opt->bspline.x_len(opt->task),
//               objective_function,
//               opt,
//               1,
//               &fc,
//               0,
//               nullptr,
//               lb.data,
//               ub.data,
//               x.data,
//               &f,
//               &stop);
//       result.num_eval = stop.nevals_p;
//
// #else
//       result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
//       result.num_eval            = nlopt_get_numevals(o);
// #endif
//     }
//
//     // validate solution
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Solution validation");
// #endif
//       Array constraints_points(opt->constraints.n_constraints);
//       compute_constraints(constraints_points.data, x, opt);
//       is_valid = max(constraints_points) < opt->success_tolerance;
//     }
//
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Solution validation (more points)");
// #endif
//       u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
//       x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
//       int points_more = (int) (steps_ms + 1);
//
//       Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
//       bspline_val_more.compute_trajectory(x, opt->task);
//       auto opt_val_more(*opt);
//       opt_val_more.set_bspline(bspline_val_more);
//       n_con(&opt_val_more);
//       Array constraints_more_points(opt_val_more.constraints.n_constraints);
//       compute_constraints(constraints_more_points.data, x, &opt_val_more);
//       is_valid_more = max(constraints_more_points) < opt_val_more.success_tolerance;
//
//       result.x = x;
//
//       if (is_valid && is_valid_more) {
//         result.trajectory = bspline_val_more.traj;
//         // break;
//       } else if (opt->guess.type != Guess::random && try_count == 0) {
//         opt->guess.type = Guess::random;
//       }
//     }
// #if BLAST_TRACE_LEVEL >= 1
//     FrameMark;
// #endif
//   }
//
//   opt->guess = start_guess; // reset to original
//
//   auto time = (real) (get_tick_us() - T1) / 1000.0;
//
//   // Output results
//   result.success       = is_valid && is_valid_more;
//   result.success_false = is_valid && !is_valid_more;
//   result.compute_time  = time;
//   result.opt           = opt;
//   result.num_tries     = try_count;
//
// #ifndef BLAST_USE_NATIVE_SQP
//   nlopt_destroy(o);
// #endif
//
//   return result;
// }

// acc 2
// inline blast_fn void constraints_and_gradients_acc2(const Array& x, Optimization& opt, Array& constraints, Matrix& grad) {
//   // constraints (p,v,a,tor) for each joint, for each segment
//   // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
//
//   // basis: n_ctrl x n_points
//   Assert(constraints.is_alias);
//   if (grad.size) {
//     Assert(grad.is_alias);
//     Assert(grad.rows == x.size);
//     Assert(grad.cols == constraints.size);
//   }
//
//   const int n_segments              = (int) opt.bspline.n_ctrl - (int) opt.bspline.p;
//   const int n_points_per_segment    = (int) opt.bspline.n_points / n_segments; // todo: check if fine? -> (Nikos) constructor that sets n_points from n_points_per_segment
//   const int n_joints                = (int) opt.manip.n_joints;
//   const int n_ctrl                  = (int) opt.bspline.n_ctrl;
//   const int x_len                   = (int) x.size;
//   const int n_constraints_per_point = (n_joints * 4); //  todo: remove hard-coded 4
//   Assert(constraints.size == n_segments * n_constraints_per_point * n_points_per_segment);
//
//   // limits
//   auto pmax    = opt.manip.pmax;
//   auto pmin    = opt.manip.pmin;
//   auto vmax    = opt.manip.vmax;
//   auto amax    = opt.manip.amax;
//   auto tau_max = opt.manip.tau_max;
//
//   for (int j = 0; j < n_joints; j++) {
//     // todo (andre): do this when initializing the optimization, not here
//     // todo: document the current behaviour in the API
//     //        (doesn't currently work if one is inf and the other is not)
//     if (pmax[j] == INF_REAL)  // note: replace INF_REAL with huge value
//       pmax[j] = 1e300;
//     if (pmin[j] == -INF_REAL) // note: replace -INF_REAL with huge negative value
//       pmin[j] = -1e300;
//   }
//
//   ManipulatorTempData manip_data;
//
//   opt.bspline.compute_trajectory(x, opt.task);
//
//
//   for (int segment = 0; segment < n_segments; segment++) {
// #if BLAST_TRACE_LEVEL >= 2
//     ZoneScopedN("All Segment Constraints");
// #endif
//     const int first_affected_control_point = std::max(3, segment);
//     const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.p);
//     const int n_affected_control_points    = last_affected_control_point - first_affected_control_point + 1; // note: affected_control_points are inclusive, so when we have last = 5, first = 3, we want 3 (5 - 3) + 1
//     const int start_point_for_segment      = segment * n_points_per_segment;                                 // note:
//     Assert(n_affected_control_points >= 3);
//     Assert(n_affected_control_points <= 6);
//
//     Matrix bp(&opt.bspline.basis_p(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//     Matrix bv(&opt.bspline.basis_v(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//     Matrix ba(&opt.bspline.basis_a(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//
//     // Initialize constraints with -INF_REAL
//     Array pos_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array vel_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array acc_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array tor_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//
//     for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) { // note:
//       auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment);              // note:
//       auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment);              // note:
//       auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment);              // note:
//
//       forward_kinematics(opt.manip, manip_data, p);
//       dynamics(opt.manip, manip_data, v, a);
//
//       for (int j = 0; j < n_joints; j++) {
//         // position
//         pos_constraints[j + n_joints * point_in_segment] = bound_constraint(p[j], pmin[j], pmax[j]);
//         // velocity
//         vel_constraints[j + n_joints * point_in_segment] = std::abs(v[j]) / vmax[j] - 1.0;
//         // acceleration
//         acc_constraints[j + n_joints * point_in_segment] = std::abs(a[j]) / amax[j] - 1.0;
//         // torque
//         tor_constraints[j + n_joints * point_in_segment] = std::abs(manip_data.efforts[j]) / tau_max[j] - 1.0;
//       }
//     }
//
//     // at this point we have max constraints for pos, vel, acc, tor, for this segment
//
//     // fill in the constraints for the current segment
//     // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
//     auto fill_idx = segment * n_constraints_per_point * n_points_per_segment;
//     for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
//
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = pos_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = vel_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = acc_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = tor_constraints[point_in_segment * n_joints + j];
//     }
//
//     // The gradient should be a (x_len)x(n_constraints) matrix that looks like this:
//     // [dp1/dx1, dp2/dx1, ..., dv1/dx1, dv2/dx1, ..., da1/dx1, da2/dx1, ..., dt1/dx1, dt2/dx1]
//     // [dp1/dx2, dp2/dx2, ..., dv1/dx2, dv2/dx2, ..., da1/dx2, da2/dx2, ..., dt1/dx2, dt2/dx2]
//     // [dp1/dx3, dp2/dx3, ..., dv1/dx3, dv2/dx3, ..., da1/dx3, da2/dx3, ..., dt1/dx3, dt2/dx3]
//     // [dp1/dx4, dp2/dx4, ..., dv1/dx4, dv2/dx4, ..., da1/dx4, da2/dx4, ..., dt1/dx4, dt2/dx4]
//     // [dp1/dx5, dp2/dx5, ..., dv1/dx5, dv2/dx5, ..., da1/dx5, da2/dx5, ..., dt1/dx5, dt2/dx5]
//     // [dp1/dx6, dp2/dx6, ..., dv1/dx6, dv2/dx6, ..., da1/dx6, da2/dx6, ..., dt1/dx6, dt2/dx6]
//     // [.....................]
//     // [.....................]
//     // [.....................]
//     // [.....................]
//     // [dp1/dT=0,dp2/dT=0,..., dv1/dT , dv1/dT , ..., da1/dT , da2/dT , ..., dt1/dT , dt2/dT ]
//     // *** per segment ***
//     // where x is the optimization vector
//     if (grad.size) {
//       // Matrix (alias) in which we can insert the gradient for the current segment
//       Matrix grad_matrix(&grad(0, segment * n_constraints_per_point * n_points_per_segment), x_len, n_constraints_per_point * n_points_per_segment);
//       Assert(grad_matrix.is_alias);
//
//       int con = 0;
//
//       for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
//         // note:      // positions
//         for (int joint = 0; joint < n_joints; joint++) {
//
//           // real pos = opt.bspline.traj.pos(joint, start_point_for_segment + point_in_segment);
//
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
//           Array bp_to_use(&bp(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Assert(bp_to_use.is_alias);
//
//           real coeff = 2.0 * std::abs(opt.bspline.traj.pos(joint, start_point_for_segment + point_in_segment)) / (pmax[joint] - pmin[joint]); // note:
//
//           for (int i = 0; i < n_affected_control_points; i++) {
//             fill_column[x_idx++] = bp_to_use[i] * coeff;
//           }
//
//           // note: dp/dT == 0
//           con++;
//         }
//
//         // velocities
//         auto one_over_T = 1 / opt.bspline.traj.t.back();
//         for (int joint = 0; joint < n_joints; joint++) {
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + n_joints + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           Array bv_to_use(&bv(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Assert(bv_to_use.is_alias);
//
//           real coeff = sign(opt.bspline.traj.vel(joint, start_point_for_segment + point_in_segment)) / vmax[joint] * one_over_T; // note:
//
//           for (int i = 0; i < n_affected_control_points; i++) {
//             fill_column[x_idx++] = bv_to_use[i] * coeff;
//           }
//
//           // dvj/dT = - (Cv + 1) / T
//           fill_column.back() = -(vel_constraints[n_joints * point_in_segment + joint] + 1) * one_over_T;
//
//           con++;
//         }
//
//         // accelerations
//         auto one_over_T2 = one_over_T * one_over_T;
//         for (int joint = 0; joint < n_joints; joint++) {
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + 2 * n_joints + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           Array ba_to_use(&ba(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Assert(ba_to_use.is_alias);
//
//           real coeff = sign(opt.bspline.traj.acc(joint, start_point_for_segment + point_in_segment)) / amax[joint] * one_over_T2; // note:
//
//           for (int i = 0; i < n_affected_control_points; i++) {
//             fill_column[x_idx++] = ba_to_use[i] * coeff;
//           }
//
//           // daj/dT = -2 * (Ca + 1) / T
//           fill_column.back() = -2 * (acc_constraints[n_joints * point_in_segment + joint] + 1) * one_over_T;
//
//           con++;
//         }
//
//         // torque
//         // [dt0/dp0, dt1/dp0, ..., dt4/dp0, dt5/dp0]
//         // [dt0/dp1, dt1/dp1, ..., dt4/dp1, dt5/dp1]
//         // [dt0/dp2, dt1/dp2, ..., dt4/dp2, dt5/dp2]
//         // [dt0/dp3, dt1/dp3, ..., dt4/dp3, dt5/dp3]
//         // [.
//         // [.
//         // [.
//         Matrix dtau_dx(x_len, n_joints); // we can store and reuse some of the values
//         auto   x_idx      = first_affected_control_point - 3;
//         auto   x_idx_skip = n_ctrl - 6;
//
//         constexpr real eps = 1e-5;
//
//         for (int j = 0; j < n_joints; j++) {                    // These loops go from x_idx -> x_len - 2
//           for (int i = 0; i < n_affected_control_points; i++) { // ^^
//             x.data[x_idx + i] += eps;
//             opt.bspline.compute_trajectory(x, opt.task);
//             x.data[x_idx + i] -= eps;
//
//             auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment); // note:
//             auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment); // note:
//             auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment); // note:
//             forward_kinematics(opt.manip, manip_data, p);
//             dynamics(opt.manip, manip_data, v, a);
//             for (int joint = 0; joint < n_joints; joint++) { // This one inserts the gradient for each joint value (6 columns)
//               const real new_constraint = std::abs(manip_data.efforts[joint]) / tau_max[joint] - 1;
//               dtau_dx(x_idx + i, joint) = (new_constraint - tor_constraints[n_joints * point_in_segment + joint]) / eps;
//             }
//           }
//           x_idx += x_idx_skip;
//         }
//         { // dtau_dx for x_len -1
//           x.data[x_len - 1] += eps;
//           opt.bspline.compute_trajectory(x, opt.task);
//           x.data[x_len - 1] -= eps;
//
//           auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment); // note:
//           auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment); // note:
//           auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment); // note:
//           forward_kinematics(opt.manip, manip_data, p);
//           dynamics(opt.manip, manip_data, v, a);
//           for (int joint = 0; joint < n_joints; joint++) {
//             const real new_constraint = std::abs(manip_data.efforts[joint]) / tau_max[joint] - 1;
//             dtau_dx(x_len - 1, joint) = (new_constraint - tor_constraints[n_joints * point_in_segment + joint]) / eps;
//           }
//         }
//
//         for (int joint = 0; joint < n_joints; joint++) {
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == n_constraints_per_point * point_in_segment + 3 * n_joints + joint); // note: Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           x_idx = first_affected_control_point - 3;
//
//           for (int j = 0; j < n_joints; j++) {
//             for (int i = 0; i < n_affected_control_points; i++) {
//               fill_column[x_idx + i] = dtau_dx(x_idx + i, joint);
//             }
//             x_idx += x_idx_skip;
//           }
//
//           // gradient w.r.t. T
//           fill_column.back() = dtau_dx(x_len - 1, joint); // note:
//           con++;
//         }
//       }
//     }
//   }
// }
//
// inline blast_fn void nlopt_constraints_acc2(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
// #if BLAST_TRACE_LEVEL >= 1
//   ZoneScoped;
// #endif
//   auto* opt = (Optimization*) f_data;
//
//   Array xv;
//   xv.alias(x, x_len);
//
//   Array constraints;
//   constraints.alias(result, m);
//
//   if (grad) {
//     memset(grad, 0, m * x_len * sizeof(real));
//     Matrix gradients;
//     gradients.alias(grad, x_len, m);
//     constraints_and_gradients_acc2(xv, *opt, constraints, gradients);
//   } else {
//     // no gradients requested
//     compute_constraints(constraints.data, xv, opt);
//   }
// }
//
// inline Result optimize_acc2(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
//   auto   T1 = get_tick_us();
//   Result result(opt); // todo: this is expensive
//
//   // Initialization
//   // configure_internal_data(opt); // todo: Ensure we can remove
//   initialize_optimization(opt);
//   n_con(opt);
//
//   // Initial validation
//   if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
//     print(opt->task);
//     return result;
//   }
//
//   const auto n = opt->bspline.x_len(opt->task);
//
//   Array con_tol(opt->constraints.n_constraints, 0.001);
//   Array x_tol(n, 0.000001);
//
// #ifdef BLAST_USE_NATIVE_SQP
//   nlopt_stopping stop;
//   stop.n          = n;
//   stop.minf_max   = -HUGE_VAL;
//   stop.ftol_rel   = 0;
//   stop.ftol_abs   = 0.0001;
//   stop.xtol_rel   = 0;
//   stop.xtol_abs   = x_tol.data;
//   stop.x_weights  = nullptr;
//   stop.nevals_p   = 0;
//   stop.maxeval    = 100000;
//   stop.maxtime    = 5000;
//   stop.start      = nlopt_seconds();
//   stop.force_stop = false;
//   stop.stop_msg   = nullptr;
//
//   Array ub(n, INF_REAL);
//   Array lb(n, -INF_REAL);
//   ub.back() = 60.0;
//   lb.back() = 0.01;
//
//   nlopt_constraint fc{};
//   fc.m      = opt->constraints.n_constraints;
//   fc.f      = nullptr;
//   fc.mf     = nlopt_constraints_acc2;
//   fc.pre    = nullptr;
//   fc.f_data = opt;
//   fc.tol    = con_tol.data;
// #else
//   nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
//   nlopt_result nlopt_res;
//   nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc2, opt, con_tol.data);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_maxtime(o, 5000);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_maxeval(o, 100000);
//   Assert(nlopt_res == NLOPT_SUCCESS);
// #endif
//
//
//   auto start_guess   = opt->guess; // save for restoration after restarts if necessary
//   int  try_count     = 0;
//   bool is_valid_more = false;
//   bool is_valid      = false;
//   for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
// #if BLAST_TRACE_LEVEL >= 1
//     ZoneScopedN("Optimization");
// #endif
//
//     // initial guess
//     Array x;
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Initial guess");
// #endif
//       x         = init_guess(opt);
//       result.x0 = x;
//     }
//
//     // launch optimization
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("NLopt optimization");
// #endif
//
//       double f = HUGE_VAL;
// #ifdef BLAST_USE_NATIVE_SQP
//       stop.nevals_p              = 0;
//       result.nlopt_exit_criteria = sqp(
//               opt->bspline.x_len(opt->task),
//               objective_function,
//               opt,
//               1,
//               &fc,
//               0,
//               nullptr,
//               lb.data,
//               ub.data,
//               x.data,
//               &f,
//               &stop);
//       result.num_eval = stop.nevals_p;
//
// #else
//       result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
//       result.num_eval            = nlopt_get_numevals(o);
// #endif
//     }
//
//     // validate solution
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Solution validation");
// #endif
//       Array constraints_points(opt->constraints.n_constraints);
//       compute_constraints(constraints_points.data, x, opt);
//       is_valid = max(constraints_points) < opt->success_tolerance;
//     }
//
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Solution validation (more points)");
// #endif
//       u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
//       x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
//       int points_more = (int) (steps_ms + 1);
//
//       Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
//       bspline_val_more.compute_trajectory(x, opt->task);
//       auto opt_val_more(*opt);
//       opt_val_more.set_bspline(bspline_val_more);
//       n_con(&opt_val_more);
//       Array constraints_more_points(opt_val_more.constraints.n_constraints);
//       compute_constraints(constraints_more_points.data, x, &opt_val_more);
//       is_valid_more = max(constraints_more_points) < opt->success_tolerance;
//
//       result.x = x;
//
//       if (is_valid && is_valid_more) {
//         result.trajectory = bspline_val_more.traj;
//         // break;
//       } else if (opt->guess.type != Guess::random && try_count == 0) {
//         opt->guess.type = Guess::random;
//       }
//     }
// #if BLAST_TRACE_LEVEL >= 1
//     FrameMark;
// #endif
//   }
//
//   opt->guess = start_guess; // reset to original
//
//   auto time = (real) (get_tick_us() - T1) / 1000.0;
//
//   // Output results
//   result.success       = is_valid && is_valid_more;
//   result.success_false = is_valid && !is_valid_more;
//   result.compute_time  = time;
//   result.opt           = opt;
//   result.num_tries     = try_count;
//
// #ifndef BLAST_USE_NATIVE_SQP
//   nlopt_destroy(o);
// #endif
//
//   return result;
// }

// acc 3
// inline blast_fn void constraints_and_gradients_acc3(const Array& x, Optimization& opt, Array& constraints, Matrix& grad) {
//   // constraints (p,v,a,tor) for each joint, for each segment
//   // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
//
//   // basis: n_ctrl x n_points
//   Assert(constraints.is_alias);
//   if (grad.size) {
//     Assert(grad.is_alias);
//     Assert(grad.rows == x.size);
//     Assert(grad.cols == constraints.size);
//   }
//
//   const int n_segments              = (int) opt.bspline.n_ctrl - (int) opt.bspline.p;
//   const int n_points_per_segment    = (int) opt.bspline.n_points / n_segments; // todo: check if fine? -> (Nikos) constructor that sets n_points from n_points_per_segment
//   const int n_joints                = (int) opt.manip.n_joints;
//   const int n_ctrl                  = (int) opt.bspline.n_ctrl;
//   const int x_len                   = (int) x.size;
//   const int n_constraints_per_point = (n_joints * 4); //  todo: remove hard-coded 4
//   Assert(constraints.size == n_segments * n_constraints_per_point * n_points_per_segment);
//
//   // limits
//   auto pmax    = opt.manip.pmax;
//   auto pmin    = opt.manip.pmin;
//   auto vmax    = opt.manip.vmax;
//   auto amax    = opt.manip.amax;
//   auto tau_max = opt.manip.tau_max;
//
//   for (int j = 0; j < n_joints; j++) {
//     // todo (andre): do this when initializing the optimization, not here
//     // todo: document the current behaviour in the API
//     //        (doesn't currently work if one is inf and the other is not)
//     if (pmax[j] == INF_REAL)  // note: replace INF_REAL with huge value
//       pmax[j] = 1e300;
//     if (pmin[j] == -INF_REAL) // note: replace -INF_REAL with huge negative value
//       pmin[j] = -1e300;
//   }
//
//   ManipulatorTempData manip_data;
//
//   opt.bspline.compute_trajectory(x, opt.task);
//
//
//   for (int segment = 0; segment < n_segments; segment++) {
// #if BLAST_TRACE_LEVEL >= 2
//     ZoneScopedN("All Segment Constraints");
// #endif
//     const int first_affected_control_point = std::max(3, segment);
//     const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.p);
//     const int n_affected_control_points    = last_affected_control_point - first_affected_control_point + 1; // note: affected_control_points are inclusive, so when we have last = 5, first = 3, we want 3 (5 - 3) + 1
//     const int start_point_for_segment      = segment * n_points_per_segment;                                 // note:
//     Assert(n_affected_control_points >= 3);
//     Assert(n_affected_control_points <= 6);
//
//     Matrix bp(&opt.bspline.basis_p(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//     Matrix bv(&opt.bspline.basis_v(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//     Matrix ba(&opt.bspline.basis_a(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
//
//     // Initialize constraints with -INF_REAL
//     Array pos_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array vel_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array acc_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//     Array tor_constraints(n_joints * n_points_per_segment, -INF_REAL);                            // note:
//
//     for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) { // note:
//       auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment);              // note:
//       auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment);              // note:
//       auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment);              // note:
//
//       forward_kinematics(opt.manip, manip_data, p);
//       dynamics(opt.manip, manip_data, v, a);
//
//       for (int j = 0; j < n_joints; j++) {
//         // position
//         pos_constraints[j + n_joints * point_in_segment] = bound_constraint(p[j], pmin[j], pmax[j]);
//         // velocity
//         vel_constraints[j + n_joints * point_in_segment] = std::abs(v[j]) / vmax[j] - 1.0;
//         // acceleration
//         acc_constraints[j + n_joints * point_in_segment] = std::abs(a[j]) / amax[j] - 1.0;
//         // torque
//         tor_constraints[j + n_joints * point_in_segment] = std::abs(manip_data.efforts[j]) / tau_max[j] - 1.0;
//       }
//     }
//
//     // at this point we have max constraints for pos, vel, acc, tor, for this segment
//
//     // fill in the constraints for the current segment
//     // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
//     auto fill_idx = segment * n_constraints_per_point * n_points_per_segment;
//     for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
//
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = pos_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = vel_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = acc_constraints[point_in_segment * n_joints + j];
//       for (int j = 0; j < n_joints; j++)
//         constraints[fill_idx++] = tor_constraints[point_in_segment * n_joints + j];
//     }
//
//     // The gradient should be a (x_len)x(n_constraints) matrix that looks like this:
//     // [dp1/dx1, dp2/dx1, ..., dv1/dx1, dv2/dx1, ..., da1/dx1, da2/dx1, ..., dt1/dx1, dt2/dx1]
//     // [dp1/dx2, dp2/dx2, ..., dv1/dx2, dv2/dx2, ..., da1/dx2, da2/dx2, ..., dt1/dx2, dt2/dx2]
//     // [dp1/dx3, dp2/dx3, ..., dv1/dx3, dv2/dx3, ..., da1/dx3, da2/dx3, ..., dt1/dx3, dt2/dx3]
//     // [dp1/dx4, dp2/dx4, ..., dv1/dx4, dv2/dx4, ..., da1/dx4, da2/dx4, ..., dt1/dx4, dt2/dx4]
//     // [dp1/dx5, dp2/dx5, ..., dv1/dx5, dv2/dx5, ..., da1/dx5, da2/dx5, ..., dt1/dx5, dt2/dx5]
//     // [dp1/dx6, dp2/dx6, ..., dv1/dx6, dv2/dx6, ..., da1/dx6, da2/dx6, ..., dt1/dx6, dt2/dx6]
//     // [.....................]
//     // [.....................]
//     // [.....................]
//     // [.....................]
//     // [dp1/dT=0,dp2/dT=0,..., dv1/dT , dv1/dT , ..., da1/dT , da2/dT , ..., dt1/dT , dt2/dT ]
//     // *** per segment ***
//     // where x is the optimization vector
//     if (grad.size) {
//       // Matrix (alias) in which we can insert the gradient for the current segment
//       Matrix grad_matrix(&grad(0, segment * n_constraints_per_point * n_points_per_segment), x_len, n_constraints_per_point * n_points_per_segment);
//       Assert(grad_matrix.is_alias);
//
//       int con = 0;
//
//       for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
//         // positions
//         for (int joint = 0; joint < n_joints; joint++) {
//
//           // real pos = opt.bspline.traj.pos(joint, start_point_for_segment + point_in_segment);
//
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
//           Array bp_to_use(&bp(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Assert(bp_to_use.is_alias);
//
//           real coeff = 2.0 * std::abs(opt.bspline.traj.pos(joint, start_point_for_segment + point_in_segment)) / (pmax[joint] - pmin[joint]); // note:
//
//           for (int i = 0; i < n_affected_control_points; i++) {
//             fill_column[x_idx++] = bp_to_use[i] * coeff;
//           }
//
//           // note: dp/dT == 0
//           con++;
//         }
//
//         // velocities
//         auto one_over_T = 1 / opt.bspline.traj.t.back();
//         for (int joint = 0; joint < n_joints; joint++) {
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + n_joints + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           Array bv_to_use(&bv(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Assert(bv_to_use.is_alias);
//
//           real coeff = sign(opt.bspline.traj.vel(joint, start_point_for_segment + point_in_segment)) / vmax[joint] * one_over_T; // note:
//
//           for (int i = 0; i < n_affected_control_points; i++) {
//             fill_column[x_idx++] = bv_to_use[i] * coeff;
//           }
//
//           // dvj/dT = - (Cv + 1) / T
//           fill_column.back() = -(vel_constraints[n_joints * point_in_segment + joint] + 1) * one_over_T;
//
//           con++;
//         }
//
//         // accelerations
//         auto one_over_T2 = one_over_T * one_over_T;
//         for (int joint = 0; joint < n_joints; joint++) {
//           // Array of the column where to put the gradient for the current constraint
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == (n_constraints_per_point * point_in_segment + 2 * n_joints + joint)); // note:
//           Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           // Which values in 'x' affect the current joint's position
//           auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
//           // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
//           x_idx += first_affected_control_point - 3;
//
//           Array ba_to_use(&ba(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Assert(ba_to_use.is_alias);
//
//           real coeff = sign(opt.bspline.traj.acc(joint, start_point_for_segment + point_in_segment)) / amax[joint] * one_over_T2; // note:
//
//           for (int i = 0; i < n_affected_control_points; i++) {
//             fill_column[x_idx++] = ba_to_use[i] * coeff;
//           }
//
//           // daj/dT = -2 * (Ca + 1) / T
//           fill_column.back() = -2 * (acc_constraints[n_joints * point_in_segment + joint] + 1) * one_over_T;
//
//           con++;
//         }
//
//         // torque
//         // [dt0/dp0, dt1/dp0, ..., dt4/dp0, dt5/dp0]
//         // [dt0/dp1, dt1/dp1, ..., dt4/dp1, dt5/dp1]
//         // [dt0/dp2, dt1/dp2, ..., dt4/dp2, dt5/dp2]
//         // [dt0/dp3, dt1/dp3, ..., dt4/dp3, dt5/dp3]
//         // [.
//         // [.
//         // [.
//         for (int joint = 0; joint < n_joints; joint++) {
//           Array fill_column = grad_matrix.col(con);
//           Assert(con == n_constraints_per_point * point_in_segment + 3 * n_joints + joint); // note: Assert(fill_column.size == x_len);
//           Assert(fill_column.is_alias);
//
//           constexpr real eps            = 1e-5;
//           auto           p              = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment); // note:
//           auto           v              = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment); // note:
//           auto           a              = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment); // note:
//           auto           old_constraint = tor_constraints[n_joints * point_in_segment + joint];                 // note:
//           auto           tau_max_now    = tau_max[joint];                                                       // note:
//
//           auto p_plus = p;
//           auto v_plus = v;
//           auto a_plus = a;
//
//           Array bp_to_use(&bp(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Array bv_to_use(&bv(first_affected_control_point, point_in_segment), n_affected_control_points);
//           Array ba_to_use(&ba(first_affected_control_point, point_in_segment), n_affected_control_points);
//
//           auto x_idx      = first_affected_control_point - 3;
//           auto x_idx_skip = n_ctrl - 6;
//
//           for (int j = 0; j < n_joints; j++) {
//
//             // partial derivative of torque constraints w.r.t. position
//             // finite difference on position
//             p_plus[j] += eps;
//             // compute the derivative of constraint(joint) w.r.t. theta(j). (remember, joint != j)
//             forward_kinematics(opt.manip, manip_data, p_plus);
//             dynamics(opt.manip, manip_data, v_plus, a_plus);
//             const real new_constraint_p = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // note: todo: remove 1.01 and use for validation only
//             const real dtau_dp          = (new_constraint_p - old_constraint) / eps;
//             // reset finite difference
//             p_plus[j] = p[j];
//
//             // note: reset forward kinematics because 'v' and 'a' don't change it
//             forward_kinematics(opt.manip, manip_data, p); // todo: precompute once
//
//             // partial derivative of torque constraints w.r.t. velocity
//             // finite difference on velocity
//             v_plus[j] += eps;
//             dynamics(opt.manip, manip_data, v_plus, a_plus);
//             const real new_constraint_v = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // note: todo: remove 1.01 and use for validation only
//             const real dtau_dv          = (new_constraint_v - old_constraint) / eps;
//             // reset finite difference
//             v_plus[j] = v[j];
//
//             // partial derivative of torque constraints w.r.t. acceleration
//             // finite difference on acceleration
//             a_plus[j] += eps;
//             dynamics(opt.manip, manip_data, v_plus, a_plus);
//             const real new_constraint_a = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // note: todo: remove 1.01 and use for validation only
//             const real dtau_da          = (new_constraint_a - old_constraint) / eps;
//             // reset finite difference
//             a_plus[j] = a[j];
//
//             // insert into the gradient
//             for (int i = 0; i < n_affected_control_points; i++) {
//               fill_column[x_idx + i] = bp_to_use[i] * dtau_dp +
//                                        bv_to_use[i] * dtau_dv * one_over_T +
//                                        ba_to_use[i] * dtau_da * one_over_T2;
//             }
//             x_idx += x_idx_skip;
//
//             // gradient w.r.t. T
//             fill_column.back() += dtau_dv * (-v[j] * one_over_T) + dtau_da * (-2 * a[j] * one_over_T); // note:
//           }
//           con++;
//         }
//       }
//     }
//   }
// }
//
// inline blast_fn void nlopt_constraints_acc3(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
// #if BLAST_TRACE_LEVEL >= 1
//   ZoneScoped;
// #endif
//   auto* opt = (Optimization*) f_data;
//
//   Array xv;
//   xv.alias(x, x_len);
//
//   Array constraints;
//   constraints.alias(result, m);
//
//   if (grad) {
//     memset(grad, 0, m * x_len * sizeof(real));
//     Matrix gradients;
//     gradients.alias(grad, x_len, m);
//     constraints_and_gradients_acc3(xv, *opt, constraints, gradients);
//   } else {
//     // no gradients requested
//     compute_constraints(constraints.data, xv, opt);
//   }
// }
//
// inline Result optimize_acc3(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
//   auto   T1 = get_tick_us();
//   Result result(opt); // todo: this is expensive
//
//   // Initialization
//   // configure_internal_data(opt); // todo: Ensure we can remove
//   initialize_optimization(opt);
//   n_con(opt);
//
//   // Initial validation
//   if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
//     print(opt->task);
//     return result;
//   }
//
//   const auto n = opt->bspline.x_len(opt->task);
//
//   Array con_tol(opt->constraints.n_constraints, 0.001);
//   Array x_tol(n, 0.000001);
//
// #ifdef BLAST_USE_NATIVE_SQP
//   nlopt_stopping stop;
//   stop.n          = n;
//   stop.minf_max   = -HUGE_VAL;
//   stop.ftol_rel   = 0;
//   stop.ftol_abs   = 0.0001;
//   stop.xtol_rel   = 0;
//   stop.xtol_abs   = x_tol.data;
//   stop.x_weights  = nullptr;
//   stop.nevals_p   = 0;
//   stop.maxeval    = 100000;
//   stop.maxtime    = 5000;
//   stop.start      = nlopt_seconds();
//   stop.force_stop = false;
//   stop.stop_msg   = nullptr;
//
//   Array ub(n, INF_REAL);
//   Array lb(n, -INF_REAL);
//   ub.back() = 60.0;
//   lb.back() = 0.01;
//
//   nlopt_constraint fc{};
//   fc.m      = opt->constraints.n_constraints;
//   fc.f      = nullptr;
//   fc.mf     = nlopt_constraints_acc3;
//   fc.pre    = nullptr;
//   fc.f_data = opt;
//   fc.tol    = con_tol.data;
// #else
//   nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
//   nlopt_result nlopt_res;
//   nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc3, opt, con_tol.data);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_maxtime(o, 5000);
//   Assert(nlopt_res == NLOPT_SUCCESS);
//   nlopt_res = nlopt_set_maxeval(o, 100000);
//   Assert(nlopt_res == NLOPT_SUCCESS);
// #endif
//
//
//   auto start_guess   = opt->guess; // save for restoration after restarts if necessary
//   int  try_count     = 0;
//   bool is_valid_more = false;
//   bool is_valid      = false;
//   for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
// #if BLAST_TRACE_LEVEL >= 1
//     ZoneScopedN("Optimization");
// #endif
//
//     // initial guess
//     Array x;
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Initial guess");
// #endif
//       x         = init_guess(opt);
//       result.x0 = x;
//     }
//
//     // launch optimization
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("NLopt optimization");
// #endif
//
//       double f = HUGE_VAL;
// #ifdef BLAST_USE_NATIVE_SQP
//       stop.nevals_p              = 0;
//       result.nlopt_exit_criteria = sqp(
//               opt->bspline.x_len(opt->task),
//               objective_function,
//               opt,
//               1,
//               &fc,
//               0,
//               nullptr,
//               lb.data,
//               ub.data,
//               x.data,
//               &f,
//               &stop);
//       result.num_eval = stop.nevals_p;
//
// #else
//       result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
//       result.num_eval            = nlopt_get_numevals(o);
// #endif
//     }
//
//     // validate solution
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Solution validation");
// #endif
//       Array constraints_points(opt->constraints.n_constraints);
//       compute_constraints(constraints_points.data, x, opt);
//       is_valid = max(constraints_points) < opt->success_tolerance;
//     }
//
//     {
// #if BLAST_TRACE_LEVEL >= 1
//       ZoneScopedN("Solution validation (more points)");
// #endif
//       u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
//       x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
//       int points_more = (int) (steps_ms + 1);
//
//       Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
//       bspline_val_more.compute_trajectory(x, opt->task);
//       auto opt_val_more(*opt);
//       opt_val_more.set_bspline(bspline_val_more);
//       n_con(&opt_val_more);
//       Array constraints_more_points(opt_val_more.constraints.n_constraints);
//       compute_constraints(constraints_more_points.data, x, &opt_val_more);
//       is_valid_more = max(constraints_more_points) < opt->success_tolerance;
//
//       result.x = x;
//
//       if (is_valid && is_valid_more) {
//         result.trajectory = bspline_val_more.traj;
//         // break;
//       } else if (opt->guess.type != Guess::random && try_count == 0) {
//         opt->guess.type = Guess::random;
//       }
//     }
// #if BLAST_TRACE_LEVEL >= 1
//     FrameMark;
// #endif
//   }
//
//   opt->guess = start_guess; // reset to original
//
//   auto time = (real) (get_tick_us() - T1) / 1000.0;
//
//   // Output results
//   result.success       = is_valid && is_valid_more;
//   result.success_false = is_valid && !is_valid_more;
//   result.compute_time  = time;
//   result.opt           = opt;
//   result.num_tries     = try_count;
//
// #ifndef BLAST_USE_NATIVE_SQP
//   nlopt_destroy(o);
// #endif
//
//   return result;
// }


// acc 1
inline u32 ncon_acc(const Optimization* opt, const int x_idx) {
  const u32 n_points            = opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]; // + 1 since ub is inclusive
  const u32 n_joints            = opt->manip.n_joints;
  const u32 n_constraints_basic = n_points * n_joints;
  u32       n_constraints       = 0;

  if (opt->constraints.position)
    n_constraints += n_constraints_basic;
  if (opt->constraints.velocity)
    n_constraints += n_constraints_basic;
  if (opt->constraints.acceleration)
    n_constraints += n_constraints_basic;
  if (opt->constraints.torque)
    n_constraints += n_constraints_basic;

  if (opt->constraints.tcp_speed)
    n_constraints += n_points;

  if (opt->constraints.self_collisions) {
    n_constraints += n_points;
  }
  if (opt->constraints.external_collisions)
    n_constraints += n_points * opt->manip._n_caps;

  for (auto& n_con: opt->constraints.n_extra_constraints)
    n_constraints += n_con;

  return n_constraints;
}

inline void compute_constraints_grad1(double* result, const Array& x, const u32 x_idx, Optimization* opt) {
  double* moving_result = result;

  opt->bspline.compute_trajectory(x, opt->task);

  // ObjMatrix<Capsule> capsules(opt->manip._n_caps, (u32) (opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx] + 1) / opt->constraints.n_collision_skip);
  for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    // todo: reset frame???
    // rotations_computed -> bool
    // forward_kinematics_computed -> bool

    ManipulatorTempData manip_data;


    // This is calculated every loop
    auto pos = opt->bspline.traj.pos.col(i);
    Assert(pos.is_alias);

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Capsules");
#endif
      // todo: this cleaner
      if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
        compute_capsules(opt->manip, manip_data);
        // if (opt->constraints.external_collisions) {
        //   if (i % opt->constraints.n_collision_skip == 0 && i / opt->constraints.n_collision_skip < capsules.cols) {
        //     for (int j = 0; j < opt->manip._n_caps; j++) {
        //       capsules(j, i / (int) opt->constraints.n_collision_skip) = manip_data.capsule_list[j];
        //     }
        //   }
        // }
      }
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = bound_constraint(opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        moving_result++;
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        // moving_result[0] = abs_constraint(opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
        moving_result[0] = std::abs(opt->bspline.traj.vel(j, i)) / opt->manip.vmax[j] - 1.0;
        moving_result++;
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = std::abs(opt->bspline.traj.acc(j, i)) / opt->manip.amax[j] - 1.0;
        moving_result++;
      }
    }

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      auto vel = opt->bspline.traj.vel.col(i);
      auto acc = opt->bspline.traj.acc.col(i);
      dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = std::abs(manip_data.efforts[j]) / opt->manip.tau_max[j] - 1.0;
        moving_result++;
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
      moving_result[0] = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
      moving_result++;
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = max(get_internal_collisions(opt->manip, manip_data));
      // for (u32 j = 0; j < tmp_coll.size; j++)
      moving_result[0] = -tmp_coll; //*std::abs(tmp_coll[j]);
      moving_result++;
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisions");
#endif
      auto collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
      // for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
      moving_result[0] = -collisions; //*std::abs(collisions[j]);
      // }
      moving_result++;
    }
  }


  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("CustomConstraints");
#endif
    for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) { // todo: split in extra_constraint per point or once
      auto extra_constraint = opt->constraints.extra_constraints[i];
      extra_constraint(moving_result, opt);
      moving_result += opt->constraints.n_extra_constraints[i];
    }
  }
}

inline void nlopt_constraints_acc1(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, x_len);
  compute_constraints(result, xv, opt);

  // gradients calculation
  if (grad) {
    memset(grad, 0, m * x_len * sizeof(real)); // zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps = 1e-5;
    Array          x_plus(x_len);

    // lb & ub automatically calculated by bspline
    u32 x_per_joint = (x_len - 1) / opt->manip.n_joints; // = n_ctrl - 6 (skip first and last 3)
    u32 mult        = 0;

    for (u32 j = 0; j < x_len - 1; j++) {     // last one is T
      mult      = j / x_per_joint > mult ? mult + 1 : mult;
      u32 x_idx = j - mult * x_per_joint + 3; // todo: fix for NaN in PVA of task

      memcpy(x_plus.data, x, x_len * sizeof(real));
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, x_len);

      // compute constraints
      auto  n_con = ncon_acc(opt, x_idx);
      Array r_plus(n_con);
      compute_constraints_grad1(r_plus.data, x_plus, x_idx, opt);
      auto n_con_lb = ncon_lb_acc(opt, x_idx);

      // todo: fix this loop
      // for (u32 i = 0; i < m; i++)
      for (u32 i = 0; i < n_con; i++)
        grad[(n_con_lb + i) * x_len + j] = (r_plus[i] - result[n_con_lb + i]) / eps;
    }
    {
      // last point T
      u32 j = x_len - 1;
      memcpy(x_plus.data, x, x_len * sizeof(real));
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, x_len);

      // compute constraints
      Array r_plus(m);
      compute_constraints(r_plus.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * x_len + j] = (r_plus[i] - result[i]) / eps;
    }
  }
}

// acc 2
inline u32 ncon_acc2(const Optimization* opt, const int x_idx) {
  const u32 n_points            = opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]; // + 1 since ub is inclusive
  const u32 n_joints            = opt->manip.n_joints;
  const u32 n_constraints_basic = n_points * n_joints;
  u32       n_constraints       = 0;

  if (opt->constraints.position)
    n_constraints += n_points;
  if (opt->constraints.velocity)
    n_constraints += n_points;
  if (opt->constraints.acceleration)
    n_constraints += n_points;
  if (opt->constraints.torque)
    n_constraints += n_constraints_basic;

  if (opt->constraints.tcp_speed)
    n_constraints += n_points;

  if (opt->constraints.self_collisions) {
    n_constraints += n_points;
  }
  if (opt->constraints.external_collisions)
    n_constraints += n_points * opt->manip._n_caps;

  for (auto& n_con: opt->constraints.n_extra_constraints)
    n_constraints += n_con;

  return n_constraints;
}

struct ConstraintPerPoint {
  Array  pos_constraint;
  Array  vel_constraint;
  Array  acc_constraint;
  Matrix tor_constraint;
  Array  tcp_constraint;
  Matrix collision_constraint;
  Array  self_collision_constraint;

  ConstraintPerPoint(int joints, int points, int n_capsules) {
    pos_constraint.resize(points);
    vel_constraint.resize(points);
    acc_constraint.resize(points);
    tor_constraint.resize(joints, points);
    tcp_constraint.resize(points);
    self_collision_constraint.resize(points);
    collision_constraint.resize(n_capsules, points);
  }
  ~ConstraintPerPoint() = default;
};

inline void compute_constraints_grad2(ConstraintPerPoint& constraints, const Array& x, const u32 x_idx, const u32 joint_idx, Optimization* opt) {

  opt->bspline.compute_trajectory(x, opt->task);

  auto n_joints = opt->manip.n_joints;

  // ObjMatrix<Capsule> capsules(opt->manip._n_caps, (u32) (opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx] + 1) / opt->constraints.n_collision_skip);
  for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    // todo: reset frame???
    // rotations_computed -> bool
    // forward_kinematics_computed -> bool

    ManipulatorTempData manip_data;


    // This is calculated every loop
    auto pos = opt->bspline.traj.pos.col(i);
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    Assert(pos.is_alias);

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Capsules");
#endif
      // todo: this cleaner
      // if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
      compute_capsules(opt->manip, manip_data);
      // }
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.pos_constraint[i - opt->bspline.lb[x_idx]] = bound_constraint(pos[joint_idx], opt->manip.pmin[joint_idx], opt->manip.pmax[joint_idx]);

      // }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.vel_constraint[i - opt->bspline.lb[x_idx]] = std::abs(vel[joint_idx]) / opt->manip.vmax[joint_idx] - 1.0;
      // }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.acc_constraint[i - opt->bspline.lb[x_idx]] = std::abs(acc[joint_idx]) / opt->manip.amax[joint_idx] - 1.0;

      // }
    }

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

      for (int j = 0; j < opt->manip.n_joints; j++) {
        constraints.tor_constraint(j, i - opt->bspline.lb[x_idx]) = std::abs(manip_data.efforts[j]) / opt->manip.tau_max[j] - 1.0;
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool                                            = get_J_tool(opt, manip_data);
      real tcp_speed                                         = norm(J_tool * opt->bspline.traj.vel.col(i));
      constraints.tcp_constraint[i - opt->bspline.lb[x_idx]] = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = max(-get_internal_collisions(opt->manip, manip_data));
      // for (u32 j = 0; j < tmp_coll.size; j++)
      constraints.self_collision_constraint[i - opt->bspline.lb[x_idx]] = tmp_coll; //*std::abs(tmp_coll[j]);
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisions");
#endif
      // auto collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
      Array max_col_constraints(opt->manip._n_caps, -INF_REAL);
      for (int capsule_id = 0; capsule_id < opt->manip._n_caps; capsule_id++) {
        real       dist_min = INF_REAL;
        const auto capsule  = manip_data.capsule_list[capsule_id];

        // check against boxes
        int count = 0;
        for (const auto& box: opt->world.boxes) {
          if (const auto dist = distance(capsule, box);
              dist < dist_min) {
            dist_min = dist;
          }
          count++;
        }

        // check against capsules
        count = 0;
        for (const auto caps: opt->world.capsules) {
          if (const auto dist = distance(capsule, caps);
              dist < dist_min) {
            dist_min = dist;
          }
          count++;
        }

        // check against spheres
        count = 0;
        for (const auto sphere: opt->world.spheres) {
          if (const auto dist = distance(capsule, sphere);
              dist < dist_min) {
            dist_min = dist;
          }
          count++;
        }

        dist_min = -dist_min; // negative distance is positive constraint

        // update worst position for the current capsule if necessary
        if (dist_min > max_col_constraints[capsule_id]) {
          max_col_constraints[capsule_id] = dist_min;
        }
      }
      for (int k = 0; k < opt->manip._n_caps; k++)
        constraints.collision_constraint(k, i - opt->bspline.lb[x_idx]) = max_col_constraints[k];
    }
  }
}

inline void nlopt_constraints_acc2(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);
  compute_constraints(result, xv, opt);

  // gradients calculation
  if (grad) {
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps = 1e-5;

    u32   n_con_lb = 0;
    Array x_plus(xlen);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx = 0;
    u32 x_idx    = 0;

    auto n_joints   = opt->manip.n_joints;
    auto n_capsules = opt->manip._n_caps;

    memcpy(x_plus.data, x, xlen * sizeof(real)); // todo: is this necessary ? done once, and x_plus is modified to original at the end of each iteration

    for (u32 j = 0; j < xlen - 1; j++) {
      // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
      // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
      joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
      x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

      x_plus[j] += eps;                                    // todo: add this is for extra constraints (tcp, collisions)
      opt->bspline.compute_trajectory(x_plus, opt->task);

      auto n_con_per_point = ncon_acc(opt, x_idx) / (opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]);

      auto               n_points = (opt->bspline.ub[x_idx] + 1 - opt->bspline.lb[x_idx]);
      ConstraintPerPoint constraint(n_joints, n_points, n_capsules);
      compute_constraints_grad2(constraint, x_plus, x_idx, joint, opt);

      n_con_lb = ncon_lb_acc(opt, x_idx); // find the amount of constraints before the current point

      grad_idx = n_con_lb * xlen + j;     // gradients are stored column-wise xlen * npoints
      for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) {
        // lb & ub are inclusive
        grad_idx += joint * xlen;

        // position
        grad[grad_idx] = (constraint.pos_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + joint]) / eps;
        grad_idx += n_joints * xlen;

        // velocity
        grad[grad_idx] = (constraint.vel_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + n_joints + joint]) / eps;
        grad_idx += n_joints * xlen;

        // acceleration
        grad[grad_idx] = (constraint.acc_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + 2 * n_joints + joint]) / eps;
        grad_idx += (n_joints - joint) * xlen;

        // torque
        for (int k = 0; k < n_joints; k++) {
          grad[grad_idx] = (constraint.tor_constraint(k, i - opt->bspline.lb[x_idx]) - result[i * n_con_per_point + 3 * n_joints + k]) / eps;
          grad_idx += xlen;
        }

        // tcp
        grad[grad_idx] = (constraint.tcp_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + 4 * n_joints]) / eps;
        grad_idx += xlen;

        // self col
        grad[grad_idx] = (constraint.self_collision_constraint[i - opt->bspline.lb[x_idx]] - result[i * n_con_per_point + 4 * n_joints + 1]) / eps;
        grad_idx += xlen;

        // col
        for (int k = 0; k < opt->manip._n_caps; k++) {
          grad[grad_idx] = (constraint.collision_constraint(k, i - opt->bspline.lb[x_idx]) - result[i * n_con_per_point + 4 * n_joints + 2 + k]) / eps;
          grad_idx += xlen;
        }
      }

      x_plus[j] = x[j];
    }

    {
      // last point T
      u32 j = xlen - 1;
      memcpy(x_plus.data, x, xlen * sizeof(real));
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, xlen);

      // compute constraints
      Array r_plus(m);
      compute_constraints(r_plus.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * xlen + j] = (r_plus[i] - result[i]) / eps;
    }
  }
}

inline Result optimize_acc2(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = 1000;
  stop.maxtime    = 30;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_acc2;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc2, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 30);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 1000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
      // note: can we initialize grad to 0 here
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      auto max_con = max(constraints_points);
      // cout << "max_con = " << max_con << endl;
      is_valid = max_con < opt->success_tolerance * 2;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more = max(constraints_more_points);
      // cout << "max_con_more = " << max_con_more << endl;
      is_valid_more = max_con_more < opt->success_tolerance * 2;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

// acc 4
template<bool is_grad> // note: n_collision_skip must be 1 for this to work !!!
blast_fn void compute_constraints_acc4(double* result, Array& gradient_coeffs, Matrix& dtau_dp, Matrix& dtau_dv, Matrix& dtau_da, Matrix& dtcp_dp, Matrix& dtcp_dv, Matrix& dselfcol_dp, std::vector<Matrix>& dcol_dp, const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  ZoneScoped;
#endif

  enum class CollisionObjectType {
    box,
    sphere,
    capsule,
  };
  struct CollisionEntities {
    // object in the world
    CollisionObjectType other_object_type = CollisionObjectType::box;
    union {
      Box     box{};
      Sphere  sphere;
      Capsule capsule;
    };

    // which point in time to look up the capsule
    int point_in_segment = 0;
  };

  double* moving_result = result;
  u32     grad_idx      = 0;

  constexpr real eps = 1e-5;

  auto joints     = opt->manip.n_joints;
  auto n_capsules = opt->manip._n_caps;

  std::array<CollisionEntities, MAX_CAPSULES> max_collision_entities{};

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // Lambda to process common bound constraint operations
  auto process_bound = [&](real value, real bound_max) {
    auto [constraint, gradient_coeff] = abs_constraint_dev<is_grad>(value, bound_max);
    *moving_result++                  = constraint;
    gradient_coeffs[grad_idx++]       = gradient_coeff;
  };

  for (u32 i = 0; i < opt->bspline.n_points; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    ManipulatorTempData manip_data;
    // Array               torque_constraint_plus(joints);

    auto pos = opt->bspline.traj.pos.col(i);
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    Assert(pos.is_alias);

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos);
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      for (int j = 0; j < joints; j++) {
        auto [constraint, gradient_coeff] = bound_constraint_dev<is_grad>(opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        *moving_result++                  = constraint;
        gradient_coeffs[grad_idx++]       = gradient_coeff;
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      for (int j = 0; j < joints; j++) {
        process_bound(opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      for (int j = 0; j < joints; j++) {
        process_bound(opt->bspline.traj.acc(j, i), opt->manip.amax[j]);
      }
    }

    Array torque_constraint(joints);
    real  tcp_constraint;
    real  self_collision_constraint;
    Array max_col_constraints(n_capsules, -INF_REAL);

    compute_capsules(opt->manip, manip_data);

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      dynamics(opt->manip, manip_data, vel, acc);
      for (int j = 0; j < joints; j++) {
        torque_constraint[j] = abs_constraint(manip_data.efforts[j], opt->manip.tau_max[j]);
        *moving_result++     = torque_constraint[j];
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * vel);
      tcp_constraint   = bound_constraint(tcp_speed, 0.0, opt->manip.tcp_max);
      *moving_result++ = tcp_constraint;
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      self_collision_constraint = max(-get_internal_collisions(opt->manip, manip_data));
      *moving_result++          = self_collision_constraint;
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisionsCalculate");
#endif
      // check every capsule with world
      for (int capsule_id = 0; capsule_id < n_capsules; capsule_id++) {
        real       dist_min = INF_REAL;
        const auto capsule  = manip_data.capsule_list[capsule_id];

        CollisionEntities collision_objects{};

        // check against boxes
        int count = 0;
        for (const auto& box: opt->world.boxes) {
          if (const auto dist = distance(capsule, box);
              dist < dist_min) {
            dist_min                            = dist;
            collision_objects.other_object_type = CollisionObjectType::box;
            collision_objects.box               = box;
            collision_objects.point_in_segment  = i;
          }
          count++;
        }

        // check against capsules
        count = 0;
        for (const auto caps: opt->world.capsules) {
          if (const auto dist = distance(capsule, caps);
              dist < dist_min) {
            dist_min                            = dist;
            collision_objects.other_object_type = CollisionObjectType::capsule;
            collision_objects.capsule           = capsule;
            collision_objects.point_in_segment  = i;
          }
          count++;
        }

        // check against spheres
        count = 0;
        for (const auto sphere: opt->world.spheres) {
          if (const auto dist = distance(capsule, sphere);
              dist < dist_min) {
            dist_min                            = dist;
            collision_objects.other_object_type = CollisionObjectType::sphere;
            collision_objects.sphere            = sphere;
            collision_objects.point_in_segment  = i;
          }
          count++;
        }

        dist_min = -dist_min; // negative distance is positive constraint

        // update worst position for the current capsule if necessary
        if (dist_min > max_col_constraints[capsule_id]) {
          max_col_constraints[capsule_id]    = dist_min;
          max_collision_entities[capsule_id] = collision_objects;
        }
      }
      for (int k = 0; k < n_capsules; k++)
        *moving_result++ = max_col_constraints[k];
    }

    if (is_grad) {
      Array pos_plus(joints);
      Array vel_plus(joints);
      Array acc_plus(joints);
      pos_plus = pos;
      vel_plus = vel;
      acc_plus = acc;

      // grad_coeffs pos
      for (u32 j = 0; j < joints; j++) {
        pos_plus[j] += eps;
        forward_kinematics(opt->manip, manip_data, pos_plus);

        if (opt->constraints.torque) {
          dynamics(opt->manip, manip_data, vel_plus, acc_plus);
          for (u32 k = 0; k < joints; k++) {
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
            dtau_dp(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        if (opt->constraints.tcp_speed) {
          auto J_tool_plus         = get_J_tool(opt, manip_data);
          real tcp_speed_plus      = norm(J_tool_plus * vel_plus);
          auto tcp_constraint_plus = bound_constraint(tcp_speed_plus, 0.0, opt->manip.tcp_max);
          dtcp_dp(j, i)            = (tcp_constraint_plus - tcp_constraint) / eps;
        }

        compute_capsules(opt->manip, manip_data);
        if (opt->constraints.self_collisions) {
          auto self_collision_constraint_plus = max(-get_internal_collisions(opt->manip, manip_data));
          dselfcol_dp(j, i)                   = (self_collision_constraint_plus - self_collision_constraint) / eps;
        }

        if (opt->constraints.external_collisions) {
          for (int capsule_id = 0; capsule_id < n_capsules; capsule_id++) {
            const auto  capsule = manip_data.capsule_list[capsule_id];
            real        distance_plus;
            const auto& objects = max_collision_entities[capsule_id];
            switch (objects.other_object_type) {
              case CollisionObjectType::box: {
                distance_plus = distance(capsule, objects.box);
                break;
              }
              case CollisionObjectType::capsule: {
                distance_plus = distance(capsule, objects.capsule);
                break;
              }
              case CollisionObjectType::sphere: {
                distance_plus = distance(capsule, objects.sphere);
                break;
              }
            }
            distance_plus = -distance_plus;
            // auto external_collisions_plus = -test_collisions_per_point(manip_data.capsule_list, &(opt->world));
            dcol_dp[i].resize(n_capsules, joints);
            dcol_dp[i](capsule_id, j) = (distance_plus - max_col_constraints[capsule_id]) / eps;
          }
        }

        pos_plus[j] = pos[j];
      }

      forward_kinematics(opt->manip, manip_data, pos_plus);
      // grad_coeffs vel
      for (u32 j = 0; j < joints; j++) {
        vel_plus[j] += eps;

        if (opt->constraints.torque) {
          dynamics(opt->manip, manip_data, vel_plus, acc_plus);
          for (u32 k = 0; k < joints; k++) {
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
            dtau_dv(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        if (opt->constraints.tcp_speed) {
          auto J_tool_plus         = get_J_tool(opt, manip_data);
          real tcp_speed_plus      = norm(J_tool_plus * vel_plus);
          auto tcp_constraint_plus = bound_constraint(tcp_speed_plus, 0.0, opt->manip.tcp_max);
          dtcp_dv(j, i)            = (tcp_constraint_plus - tcp_constraint) / eps;
        }

        vel_plus[j] = vel[j];
      }

      // grad_coeffs acc
      for (u32 j = 0; j < joints; j++) {
        acc_plus[j] += eps;

        if (opt->constraints.torque) {
          dynamics(opt->manip, manip_data, vel_plus, acc_plus);
          for (u32 k = 0; k < joints; k++) {
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.tau_max[k]);
            dtau_da(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        acc_plus[j] = acc[j];
      }
    }
  }

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("CustomConstraints");
#endif
    for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) {
      auto extra_constraint = opt->constraints.extra_constraints[i];
      extra_constraint(moving_result, opt);
      moving_result += opt->constraints.n_extra_constraints[i];
      grad_idx += opt->constraints.n_extra_constraints[i]; // todo: add analytical gradients
    }
  }
}

inline blast_fn void nlopt_constraints_acc4(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);

  int n_active_constraints = 0;
  if (opt->constraints.position)
    n_active_constraints++;
  if (opt->constraints.velocity)
    n_active_constraints++;
  if (opt->constraints.acceleration)
    n_active_constraints++;

  Array               gradient_coeffs(opt->manip.n_joints * n_active_constraints * opt->bspline.n_points); // todo: check performance
  Matrix              dtau_dp(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix              dtau_dv(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix              dtau_da(opt->manip.n_joints, opt->manip.n_joints * opt->bspline.n_points);           // todo: check performance
  Matrix              dtcp_dp(opt->manip.n_joints, opt->bspline.n_points);
  Matrix              dtcp_dv(opt->manip.n_joints, opt->bspline.n_points);
  Matrix              dselfcol_dp(opt->manip.n_joints, opt->bspline.n_points);
  std::vector<Matrix> dcol_dp; // (n_caps, n_joints)
  dcol_dp.resize(opt->bspline.n_points);
  // (opt->manip.n_joints, opt->bspline.n_points);
  if (!grad) {
    compute_constraints_acc4<false>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtcp_dp, dtcp_dv, dselfcol_dp, dcol_dp, xv, opt);
  }
  if (grad) {
    compute_constraints_acc4<true>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtcp_dp, dtcp_dv, dselfcol_dp, dcol_dp, xv, opt);
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps      = 1e-5;
    u32            n_con_lb = 0;
    Array          x_plus(xlen);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx       = 0;
    u32 constraint_idx = 0;
    u32 x_idx          = 0;

    for (u32 j = 0; j < xlen - 1; j++) { // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
      // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
      joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
      x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

      n_con_lb = ncon_lb_acc(opt, x_idx);                  // find the amount of constraints before the current point

      // todo: create alias matrix that points to grad
      // todo: can we change the order in which we store the gradients ?
      grad_idx       = n_con_lb * xlen + j;                                    // gradients are stored column-wise xlen * npoints
      constraint_idx = opt->manip.n_joints * n_active_constraints * opt->bspline.lb[x_idx];
      for (u32 i = opt->bspline.lb[x_idx]; i <= opt->bspline.ub[x_idx]; i++) { // lb & ub are inclusive
        grad_idx += joint * xlen;
        constraint_idx += joint;

        if (opt->constraints.position) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * gradient_coeffs[constraint_idx];
          grad_idx += opt->manip.n_joints * xlen; // increase index by the amount of joints * xlen
          constraint_idx += opt->manip.n_joints;
        }

        if (opt->constraints.velocity) {          // todo: basis_v / t once before this loop
          grad[grad_idx] = opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1] * gradient_coeffs[constraint_idx];
          grad_idx += opt->manip.n_joints * xlen; // increase index by the amount of joints * xlen
          constraint_idx += opt->manip.n_joints;
        }

        if (opt->constraints.acceleration) {                // todo: basis_a / t / t once before this loop
          grad[grad_idx] = opt->bspline.basis_a(x_idx, i) / (opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1] * opt->bspline.traj.t.data[opt->bspline.traj.t.size - 1]) * gradient_coeffs[constraint_idx];
          grad_idx += (opt->manip.n_joints - joint) * xlen; // increase index by the amount of (joints - current joint) * xlen
          constraint_idx += (opt->manip.n_joints - joint);
        }

        if (opt->constraints.torque) {

          for (u32 k = 0; k < opt->manip.n_joints; k++) {
            grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dtau_dp(k, joint + opt->manip.n_joints * i) +
                             opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * dtau_dv(k, joint + opt->manip.n_joints * i) +
                             opt->bspline.basis_a(x_idx, i) / (opt->bspline.traj.t.back() * opt->bspline.traj.t.back()) * dtau_da(k, joint + opt->manip.n_joints * i);
            grad_idx += xlen;
            // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
          }
        }

        if (opt->constraints.tcp_speed) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dtcp_dp(joint, i) +
                           opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * dtcp_dv(joint, i);
          grad_idx += xlen;
          // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }

        if (opt->constraints.self_collisions) {
          grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dselfcol_dp(joint, i);
          grad_idx += xlen;
          // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }

        if (opt->constraints.external_collisions) {
          for (int k = 0; k < opt->manip._n_caps; k++) {
            grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dcol_dp[i](k, joint);
            grad_idx += xlen;
          }
          // constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }
      }
    }

    {
      // last point T
      u32 j = xlen - 1;

      auto   n_joints             = opt->manip.n_joints;
      auto   constraint_per_point = m / opt->bspline.n_points;
      auto   one_over_T           = 1 / opt->bspline.traj.t.back();
      Matrix gradients;
      gradients.alias(grad, xlen, m);
      Array constraints;
      constraints.alias(result, m);

      for (int i = 1; i < opt->bspline.n_points; i++) {
        int    constraint_in_point_idx = 0;
        auto   vel                     = opt->bspline.traj.vel.col(i);
        auto   acc                     = opt->bspline.traj.acc.col(i);
        Matrix grad_point(&gradients(0, i * constraint_per_point), xlen, constraint_per_point);
        Array  constraint_point(&constraints[i * constraint_per_point], constraint_per_point);

        // dp/dT == 0
        constraint_in_point_idx += (int) n_joints;
        // dv/dT
        for (int k = 0; k < n_joints; k++) {
          grad_point.data[constraint_in_point_idx * xlen + j] = -(constraint_point[constraint_in_point_idx] + 1) * one_over_T;
          constraint_in_point_idx++;
        }
        // da_dT
        for (int k = 0; k < n_joints; k++) {
          grad_point.data[constraint_in_point_idx * xlen + j] = -2 * (constraint_point[constraint_in_point_idx] + 1) * one_over_T;
          constraint_in_point_idx++;
        }
        // dtau_dT
        for (int k = 0; k < n_joints; k++) {
          for (int l = 0; l < n_joints; l++) {
            grad_point.data[constraint_in_point_idx * xlen + j] += dtau_dv(k, l + n_joints * i) * (-vel[l] * one_over_T) + dtau_da(k, l + n_joints * i) * (-2 * acc[l] * one_over_T);
          }
          constraint_in_point_idx++;
        }
        // dtcp_dT
        {
          for (int k = 0; k < n_joints; k++) {
            grad_point.data[constraint_in_point_idx * xlen + j] += dtcp_dv(k, i) * (-vel[k] * one_over_T);
          }
          constraint_in_point_idx++; // unused, but added for uniformity
        }
      }
    }

    if (opt->constraints.show_info) { // when more info is needed per iteration
      Matrix gradients(xlen, m);
      Array  constraints(m);
      for (u32 j = 0; j < xlen; j++) {
        for (u32 i = 0; i < m; i++) {
          gradients(j, i) = grad[i * xlen + j];
          constraints[i]  = result[i];
        }
      }
      opt->constraints.grad_list.push_back(gradients);
      opt->constraints.constr_list.push_back(constraints);
      opt->constraints.x_list.push_back(xv);
    }
  }
}

inline Result optimize_acc4(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  initialize_optimization(opt);
  n_con(opt);

  // Initial validation
  if (!validate_task(opt)) { // todo: support validate_task when there are no capsules...
    print(opt->task);
    return result;
  }

  const auto n = opt->bspline.x_len(opt->task);

  Array con_tol(opt->constraints.n_constraints, opt->success_tolerance);
  Array x_tol(n, 0.001);

#ifdef BLAST_USE_NATIVE_SQP
  nlopt_stopping stop;
  stop.n          = n;
  stop.minf_max   = -HUGE_VAL;
  stop.ftol_rel   = 0;
  stop.ftol_abs   = 0.001;
  stop.xtol_rel   = 0;
  stop.xtol_abs   = x_tol.data;
  stop.x_weights  = nullptr;
  stop.nevals_p   = 0;
  stop.maxeval    = 100000;
  stop.maxtime    = 30;
  stop.start      = nlopt_seconds();
  stop.force_stop = false;
  stop.stop_msg   = nullptr;

  Array ub(n, INF_REAL);
  Array lb(n, -INF_REAL);
  ub.back() = 60.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_acc4;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc4, opt, con_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_min_objective(o, objective_function, opt);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_lower_bound(o, (int) n - 1, 0.01);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_upper_bound(o, (int) n - 1, 60.0);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_ftol_abs(o, 0.0001);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_xtol_abs(o, x_tol.data);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxtime(o, 5000);
  Assert(nlopt_res == NLOPT_SUCCESS);
  nlopt_res = nlopt_set_maxeval(o, 100000);
  Assert(nlopt_res == NLOPT_SUCCESS);
#endif


  auto start_guess   = opt->guess; // save for restoration after restarts if necessary
  int  try_count     = 0;
  bool is_valid_more = false;
  bool is_valid      = false;
  for (; try_count < opt->max_tries; try_count++) { // todo: add nlopt stop criteria to list, add max_time for full loop
#if BLAST_TRACE_LEVEL >= 1
    ZoneScopedN("Optimization");
#endif

    // initial guess
    Array x;
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Initial guess");
#endif
      x         = init_guess(opt);
      result.x0 = x;
    }

    // launch optimization
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("NLopt optimization");
#endif

      double f = HUGE_VAL;
#ifdef BLAST_USE_NATIVE_SQP
      stop.nevals_p              = 0;
      result.nlopt_exit_criteria = sqp(
              opt->bspline.x_len(opt->task),
              objective_function,
              opt,
              1,
              &fc,
              0,
              nullptr,
              lb.data,
              ub.data,
              x.data,
              &f,
              &stop);
      result.num_eval = stop.nevals_p;

#else
      result.nlopt_exit_criteria = nlopt_optimize(o, x.data, &f);
      result.num_eval            = nlopt_get_numevals(o);
#endif
    }

    // validate solution
    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation");
#endif
      Array constraints_points(opt->constraints.n_constraints);
      compute_constraints(constraints_points.data, x, opt);
      auto max_con = max(constraints_points);
      // cout << "max_con = " << max_con << endl;
      is_valid = max_con < opt->success_tolerance * 2;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.p, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more = max(constraints_more_points);
      // cout << "max_con_more = " << max_con_more << endl;
      is_valid_more = max_con_more < opt->success_tolerance * 2;

      result.x = x;

      if (is_valid && is_valid_more) {
        result.trajectory = bspline_val_more.traj;
        // break;
      } else if (opt->guess.type != Guess::random && try_count == 0) {
        opt->guess.type = Guess::random;
      }
    }
#if BLAST_TRACE_LEVEL >= 1
    FrameMark;
#endif
  }

  opt->guess = start_guess; // reset to original

  auto time = (real) (get_tick_us() - T1) / 1000.0;

  // Output results
  result.success       = is_valid && is_valid_more;
  result.success_false = is_valid && !is_valid_more;
  result.compute_time  = time;
  result.opt           = opt;
  result.num_tries     = try_count;

#ifndef BLAST_USE_NATIVE_SQP
  nlopt_destroy(o);
#endif

  return result;
}

TEST_CASE("Acceleration 1: Analytic Gradients for PVA", "[Paper2]") {
  std::array<Config, 13> config_list;
  fill_config_list(config_list);

  real eps = 0.01;

  for (int config_id = 0; config_id < config_list.size(); config_id++) {
    auto n_ctrl   = config_list[config_id].n_ctrl;
    auto n_points = config_list[config_id].n_points;
    auto task_idx = config_list[config_id].task_idx;
    auto n_optim  = config_list[config_id].n_optim;
    auto task     = config_list[config_id].task;
    // Optimization opt(get_generic_gen3(), get_gen3_task());
    Optimization opt(get_generic_Link6(), task);

    opt.world = get_lab_world();

    // Sleep(200);

    Bspline bspline(n_ctrl, n_points, 5, 6);
    opt.bspline = bspline;

    // opt.guess.type = Guess::custom;
    // opt.guess.x0   = Array(opt.bspline.x_len(opt.task), 1.5);
    // opt.guess.type = Guess::random;
    // opt.guess.n_shot = 100;

    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tcp_speed           = true;
    opt.constraints.self_collisions     = true;
    opt.constraints.external_collisions = true;

    opt.constraints.n_collision_skip = 1;

    opt.max_tries         = 1;
    opt.success_tolerance = 0.01;

    opt.guess.type = Guess::custom;
    opt.guess.x0   = guess_random((opt.bspline), opt.task);

    // opt.guess.x0 = x;

    cout << "Config ID:                  " << config_id << endl;
    cout << "Task id:                    " << task_idx << endl;
    cout << "n_ctrl:                     " << n_ctrl << endl;
    cout << "n_points:                   " << n_points << endl;

    Result result(&opt);
    Result result_acc1(&opt);
    Result result_acc2(&opt);
    Result result_acc3(&opt);
    for (int iter = 0; iter < n_optim; iter++) {
      {
        {
          IOSilencer _;
          result = optimize(&opt);
        }

        // {
        //   IOSilencer _;
        //   result_acc1 = optimize_acc1(&opt);
        // }

        {
          // IOSilencer _;
          // result_acc2 = optimize_acc2(&opt);
        }

        {
          IOSilencer _;
          result_acc3 = optimize_acc4(&opt);
        }
        // CHECK(result_acc1.success == result.success);
        // CHECK(result_acc2.success == result.success);
        CHECK(result_acc3.success == result.success);
        // CHECK(result_acc1.success_false == result.success_false);
        // CHECK(result_acc2.success_false == result.success_false);
        CHECK(result_acc3.success_false == result.success_false);
        // CHECK(std::abs(result.x.back() - result_acc1.x.back()) < 0.01);
        // CHECK(std::abs(result.x.back() - result_acc2.x.back()) < 0.01);
        CHECK(std::abs(result.x.back() - result_acc3.x.back()) < 0.01);
      }

      // cout << "result.x.bac(): " << result.x.back() << "result_acc.x.back(): " << result_acc.x.back() << endl;
      // CHECK(std::abs(result.x.back() - result_acc.x.back()) < eps);
    }
  }
};
