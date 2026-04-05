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
  int  n_ctrl   = 0;
  int  n_points = 0;
  int  task_idx = 0;
  int  n_optim  = 0;
  Task task;
};

// We create a config list that will be used for all benchmarks
inline void fill_config_list(std::array<Config, 13>& config_list) {

  std::array<std::tuple<int, int>, 1> bspline_config_list = {std::make_tuple(16, 110)}; // These give rounded number of n_points_per_segment
  std::vector<Task>                   tasks               = get_Link6_demo1_tasks();    // 13 tasks

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

// acc 1
inline u32 ncon_acc(const Optimization* opt, const int x_idx) {
  const u32 n_points            = opt->bspline.upper_bounds[x_idx] + 1 - opt->bspline.lower_bounds[x_idx]; // + 1 since ub is inclusive
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

  if (opt->constraints.tool_speed)
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
  Array  tool_constraint;
  Matrix collision_constraint;
  Array  self_collision_constraint;

  ConstraintPerPoint(int joints, int points, int n_capsules) {
    pos_constraint.resize(points);
    vel_constraint.resize(points);
    acc_constraint.resize(points);
    tor_constraint.resize(joints, points);
    tool_constraint.resize(points);
    self_collision_constraint.resize(points);
    collision_constraint.resize(n_capsules, points);
  }
  ~ConstraintPerPoint() = default;
};

inline void compute_constraints_grad1(ConstraintPerPoint& constraints, const Array& x, const u32 x_idx, const u32 joint_idx, Optimization* opt) {

  opt->bspline.compute_trajectory(x, opt->task);

  auto n_joints = opt->manip.n_joints;

  // ObjMatrix<Capsule> capsules(opt->manip._n_caps, (u32) (opt->bspline.upper_bounds[x_idx] - opt->bspline.lower_bounds[x_idx] + 1) / opt->constraints.n_collision_skip);
  for (u32 i = opt->bspline.lower_bounds[x_idx]; i <= opt->bspline.upper_bounds[x_idx]; i++) {
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
      constraints.pos_constraint[i - opt->bspline.lower_bounds[x_idx]] = bound_constraint(pos[joint_idx], opt->manip.position_min[joint_idx], opt->manip.position_max[joint_idx]);

      // }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.vel_constraint[i - opt->bspline.lower_bounds[x_idx]] = std::abs(vel[joint_idx]) / opt->manip.velocity_max[joint_idx] - 1.0;
      // }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      // for (int j = 0; j < opt->manip.n_joints; j++) {
      constraints.acc_constraint[i - opt->bspline.lower_bounds[x_idx]] = std::abs(acc[joint_idx]) / opt->manip.acceleration_max[joint_idx] - 1.0;

      // }
    }

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

      for (int j = 0; j < opt->manip.n_joints; j++) {
        constraints.tor_constraint(j, i - opt->bspline.lower_bounds[x_idx]) = std::abs(manip_data.efforts[j]) / opt->manip.torque_max[j] - 1.0;
      }
    }

    if (opt->constraints.tool_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ToolSpeed");
#endif
      auto J_tool     = get_J_tool(opt, manip_data);
      real tool_speed = norm(J_tool * opt->bspline.traj.vel.col(i));

      constraints.tool_constraint[i - opt->bspline.lower_bounds[x_idx]] = bound_constraint(tool_speed, 0.0, opt->manip.tool_speed_max);
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = max(-get_internal_collisions(opt->manip, manip_data));
      // for (u32 j = 0; j < tmp_coll.size; j++)
      constraints.self_collision_constraint[i - opt->bspline.lower_bounds[x_idx]] = tmp_coll; //*std::abs(tmp_coll[j]);
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
        constraints.collision_constraint(k, i - opt->bspline.lower_bounds[x_idx]) = max_col_constraints[k];
    }
  }
}

inline void nlopt_constraints_acc1(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);
  {

#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Constraints");
#endif
    compute_constraints(result, xv, opt);
  }

  // gradients calculation
  if (grad) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("Grad");
#endif
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

      x_plus[j] += eps;                                    // todo: add this is for extra constraints (tool, collisions)
      opt->bspline.compute_trajectory(x_plus, opt->task);

      auto n_con_per_point = ncon_acc(opt, x_idx) / (opt->bspline.upper_bounds[x_idx] + 1 - opt->bspline.lower_bounds[x_idx]);

      auto               n_points = (opt->bspline.upper_bounds[x_idx] + 1 - opt->bspline.lower_bounds[x_idx]);
      ConstraintPerPoint constraint(n_joints, n_points, n_capsules);
      compute_constraints_grad1(constraint, x_plus, x_idx, joint, opt);

      n_con_lb = ncon_lb_acc(opt, x_idx); // find the amount of constraints before the current point

      grad_idx = n_con_lb * xlen + j;     // gradients are stored column-wise xlen * npoints
      for (u32 i = opt->bspline.lower_bounds[x_idx]; i <= opt->bspline.upper_bounds[x_idx]; i++) {
        // lb & ub are inclusive
        grad_idx += joint * xlen;

        // position
        grad[grad_idx] = (constraint.pos_constraint[i - opt->bspline.lower_bounds[x_idx]] - result[i * n_con_per_point + joint]) / eps;
        grad_idx += n_joints * xlen;

        // velocity
        grad[grad_idx] = (constraint.vel_constraint[i - opt->bspline.lower_bounds[x_idx]] - result[i * n_con_per_point + n_joints + joint]) / eps;
        grad_idx += n_joints * xlen;

        // acceleration
        grad[grad_idx] = (constraint.acc_constraint[i - opt->bspline.lower_bounds[x_idx]] - result[i * n_con_per_point + 2 * n_joints + joint]) / eps;
        grad_idx += (n_joints - joint) * xlen;

        // torque
        for (int k = 0; k < n_joints; k++) {
          grad[grad_idx] = (constraint.tor_constraint(k, i - opt->bspline.lower_bounds[x_idx]) - result[i * n_con_per_point + 3 * n_joints + k]) / eps;
          grad_idx += xlen;
        }

        // tool
        grad[grad_idx] = (constraint.tool_constraint[i - opt->bspline.lower_bounds[x_idx]] - result[i * n_con_per_point + 4 * n_joints]) / eps;
        grad_idx += xlen;

        // self col
        grad[grad_idx] = (constraint.self_collision_constraint[i - opt->bspline.lower_bounds[x_idx]] - result[i * n_con_per_point + 4 * n_joints + 1]) / eps;
        grad_idx += xlen;

        // col
        for (int k = 0; k < opt->manip._n_caps; k++) {
          grad[grad_idx] = (constraint.collision_constraint(k, i - opt->bspline.lower_bounds[x_idx]) - result[i * n_con_per_point + 4 * n_joints + 2 + k]) / eps;
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

inline Result optimize_acc1(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  // initialize_optimization(opt);
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
  ub.back() = 30.0;
  lb.back() = 0.01;

  nlopt_constraint fc{};
  fc.m      = opt->constraints.n_constraints;
  fc.f      = nullptr;
  fc.mf     = nlopt_constraints_acc1;
  fc.pre    = nullptr;
  fc.f_data = opt;
  fc.tol    = con_tol.data;
#else
  nlopt_opt    o = nlopt_create(NLOPT_LD_SLSQP, n);
  nlopt_result nlopt_res;
  nlopt_res = nlopt_add_inequality_mconstraint(o, opt->constraints.n_constraints, nlopt_constraints_acc1, opt, con_tol.data);
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
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance * 2;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.degree, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance * 2;

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

// acc 2
template<bool is_grad> // note: n_collision_skip must be 1 for this to work !!!
blast_fn void compute_constraints_acc2(double* result, Array& gradient_coeffs, Matrix& dtau_dp, Matrix& dtau_dv, Matrix& dtau_da, Matrix& dtool_dp, Matrix& dtool_dv, Matrix& dselfcol_dp, std::vector<Matrix>& dcol_dp, const Array& x, Optimization* opt) {
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
#if BLAST_TRACE_LEVEL >= 2
    ZoneScopedN("ConstraintSinglePoint");
#endif

    ManipulatorTempData manip_data;
    // Array               torque_constraint_plus(joints);

    auto pos = opt->bspline.traj.pos.col(i);
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    Assert(pos.is_alias);

    Array torque_constraint(joints);
    real  tool_constraint;
    real  self_collision_constraint;
    Array max_col_constraints(n_capsules, -INF_REAL);

    {
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("Constraints");
#endif


      {
#if BLAST_TRACE_LEVEL >= 3
        ZoneScopedN("Constraints per point");
#endif


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
            auto [constraint, gradient_coeff] = bound_constraint_dev<is_grad>(opt->bspline.traj.pos(j, i), opt->manip.position_min[j], opt->manip.position_max[j]);
            *moving_result++                  = constraint;
            gradient_coeffs[grad_idx++]       = gradient_coeff;
          }
        }

        if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
          ZoneScopedN("Velocity");
#endif
          for (int j = 0; j < joints; j++) {
            process_bound(opt->bspline.traj.vel(j, i), opt->manip.velocity_max[j]);
          }
        }

        if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
          ZoneScopedN("Acceleration");
#endif
          for (int j = 0; j < joints; j++) {
            process_bound(opt->bspline.traj.acc(j, i), opt->manip.acceleration_max[j]);
          }
        }

        compute_capsules(opt->manip, manip_data);

        if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
          ZoneScopedN("Dynamics");
#endif
          dynamics(opt->manip, manip_data, vel, acc);
          for (int j = 0; j < joints; j++) {
            torque_constraint[j] = abs_constraint(manip_data.efforts[j], opt->manip.torque_max[j]);
            *moving_result++     = torque_constraint[j];
          }
        }

        if (opt->constraints.tool_speed) {
#if BLAST_TRACE_LEVEL >= 3
          ZoneScopedN("ToolSpeed");
#endif
          auto J_tool      = get_J_tool(opt, manip_data);
          real tool_speed  = norm(J_tool * vel);
          tool_constraint  = bound_constraint(tool_speed, 0.0, opt->manip.tool_speed_max);
          *moving_result++ = tool_constraint;
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
      }
    }
    if (is_grad) {
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("Grad");
#endif
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
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.torque_max[k]);
            dtau_dp(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        if (opt->constraints.tool_speed) {
          auto J_tool_plus          = get_J_tool(opt, manip_data);
          real tool_speed_plus      = norm(J_tool_plus * vel_plus);
          auto tool_constraint_plus = bound_constraint(tool_speed_plus, 0.0, opt->manip.tool_speed_max);
          dtool_dp(j, i)            = (tool_constraint_plus - tool_constraint) / eps;
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
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.torque_max[k]);
            dtau_dv(k, j + joints * i)  = (torque_constraint_plus - torque_constraint[k]) / eps;
          }
        }

        if (opt->constraints.tool_speed) {
          auto J_tool_plus          = get_J_tool(opt, manip_data);
          real tool_speed_plus      = norm(J_tool_plus * vel_plus);
          auto tool_constraint_plus = bound_constraint(tool_speed_plus, 0.0, opt->manip.tool_speed_max);
          dtool_dv(j, i)            = (tool_constraint_plus - tool_constraint) / eps;
        }

        vel_plus[j] = vel[j];
      }

      // grad_coeffs acc
      for (u32 j = 0; j < joints; j++) {
        acc_plus[j] += eps;

        if (opt->constraints.torque) {
          dynamics(opt->manip, manip_data, vel_plus, acc_plus);
          for (u32 k = 0; k < joints; k++) {
            auto torque_constraint_plus = abs_constraint(manip_data.efforts[k], opt->manip.torque_max[k]);
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

inline blast_fn void nlopt_constraints_acc2(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
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
  Matrix              dtool_dp(opt->manip.n_joints, opt->bspline.n_points);
  Matrix              dtool_dv(opt->manip.n_joints, opt->bspline.n_points);
  Matrix              dselfcol_dp(opt->manip.n_joints, opt->bspline.n_points);
  std::vector<Matrix> dcol_dp; // (n_caps, n_joints)
  dcol_dp.resize(opt->bspline.n_points);
  // (opt->manip.n_joints, opt->bspline.n_points);
  if (!grad) {

    compute_constraints_acc2<false>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtool_dp, dtool_dv, dselfcol_dp, dcol_dp, xv, opt);
  }
  if (grad) {
    compute_constraints_acc2<true>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtool_dp, dtool_dv, dselfcol_dp, dcol_dp, xv, opt);
    {
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("Grad Fill");
#endif
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
        grad_idx       = n_con_lb * xlen + j;                                                        // gradients are stored column-wise xlen * npoints
        constraint_idx = opt->manip.n_joints * n_active_constraints * opt->bspline.lower_bounds[x_idx];
        for (u32 i = opt->bspline.lower_bounds[x_idx]; i <= opt->bspline.upper_bounds[x_idx]; i++) { // lb & ub are inclusive
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

          if (opt->constraints.tool_speed) {
            grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dtool_dp(joint, i) +
                             opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * dtool_dv(joint, i);
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
          // dtool_dT
          {
            for (int k = 0; k < n_joints; k++) {
              grad_point.data[constraint_in_point_idx * xlen + j] += dtool_dv(k, i) * (-vel[k] * one_over_T);
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
}

inline Result optimize_acc2(Optimization* opt, u32 output_steps_ms = 1 /*ms*/) {
  auto   T1 = get_tick_us();
  Result result(opt); // todo: this is expensive

  // Initialization
  // configure_internal_data(opt); // todo: Ensure we can remove
  // initialize_optimization(opt);
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
  ub.back() = 30.0;
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
      auto max_con                = max(constraints_points);
      result.max_constraint_idx   = argmax(constraints_points);
      result.max_constraint_value = max_con;
      is_valid                    = max_con < opt->success_tolerance * 2;
    }

    {
#if BLAST_TRACE_LEVEL >= 1
      ZoneScopedN("Solution validation (more points)");
#endif
      u64 steps_ms    = (u64) (std::ceil(x.back() * 1e3 / output_steps_ms));
      x.back()        = (real) (std::ceil(x.back() * 1000.0 / output_steps_ms) * output_steps_ms) * 1e-3;
      int points_more = (int) (steps_ms + 1);

      Bspline bspline_val_more(opt->bspline.n_ctrl, points_more, opt->bspline.degree, opt->manip.n_joints); // todo: this is expensive
      bspline_val_more.compute_trajectory(x, opt->task);
      auto opt_val_more(*opt);
      opt_val_more.set_bspline(bspline_val_more);
      n_con(&opt_val_more);
      Array constraints_more_points(opt_val_more.constraints.n_constraints);
      compute_constraints(constraints_more_points.data, x, &opt_val_more);
      auto max_con_more                       = max(constraints_more_points);
      result.max_constraint_more_points_idx   = argmax(constraints_more_points);
      result.max_constraint_more_points_value = max_con_more;
      is_valid_more                           = max_con_more < opt->success_tolerance * 2;

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

TEST_CASE("Acceleration tests", "[Paper2]") {
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

    Bspline bspline(n_ctrl, n_points, 5, 6);
    opt.bspline = bspline;

    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tool_speed          = true;
    opt.constraints.self_collisions     = true;
    opt.constraints.external_collisions = true;

    opt.constraints.n_collision_skip = 1;

    opt.max_tries         = 1;
    opt.success_tolerance = 0.01;

    opt.guess.type      = Guess::custom;
    opt.guess.initial_x = guess_random((opt.bspline), opt.task);

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
          result = optimize(&opt, OptimizationMethod::baseline);
        }

        {
          IOSilencer _;
          result_acc1 = optimize_acc1(&opt);
        }

        {
          IOSilencer _;
          result_acc2 = optimize_acc2(&opt);
        }

        {
          IOSilencer _;
          result_acc3 = optimize(&opt);
        }
        CHECK(result_acc1.success == result.success);
        CHECK(result_acc2.success == result.success);
        CHECK(result_acc3.success == result.success);
        CHECK(result_acc1.success_false == result.success_false);
        CHECK(result_acc2.success_false == result.success_false);
        CHECK(result_acc3.success_false == result.success_false);
        CHECK(std::abs(result.x.back() - result_acc1.x.back()) < 0.01);
        CHECK(std::abs(result.x.back() - result_acc2.x.back()) < 0.01);
        CHECK(std::abs(result.x.back() - result_acc3.x.back()) < 0.01);
      }
    }
  }
};
