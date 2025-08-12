#pragma once

#include <blast>
#include "optimization/optimization.hpp"

namespace blast {


inline void ConstraintSelection::add_constraint(ConstraintFunction constraint, int n_con) {
  extra_constraints.push_back(constraint);
  n_extra_constraints.push_back(n_con);
}

inline real bound_constraint(const real& value, const real& value_min, const real& value_max) {
  real g1 = value_max == INF_REAL ? -1.0 : value_max == 0.0 ? value
                                                            : sign(value_max) * (value / value_max - 1);
  real g2 = value_min == -INF_REAL ? -1.0 : value_min == 0.0 ? -value // todo: wait, what!?
                                                             : sign(value_min) * (1 - value / value_min);
  return g1 > g2 ? g1
                 : g2;
}

// todo: replace by jacobian?
inline Matrix get_J_tool(const Optimization* opt) {
  std::vector<Vec3> r_tool(opt->manip.joints);
  r_tool[opt->manip.joints - 1] = opt->manip.dv[opt->manip.joints - 1];
  for (int i = (int) opt->manip.joints - 2; i >= 0; i--) {
    r_tool[i] = opt->manip.dv[i] + opt->manip._rotations[i + 1] * r_tool[i + 1];
  }

  for (int i = 0; i < opt->manip.joints; i++) {
    r_tool[i] = opt->manip._rotations_mult[i] * r_tool[i];
  }

  Matrix J_tool(3, opt->manip.joints);
  for (int i = 0; i < opt->manip.joints; i++) {
    // e = opt->manip._Q_mult[i] * opt->manip.ev[i]; // replaced e directly in function to skip copy
    Vec3 cr_tool = cross(opt->manip._rotations_mult[i] * opt->manip.ev[i], r_tool[i]);
    J_tool(0, i) = cr_tool.x;
    J_tool(1, i) = cr_tool.y;
    J_tool(2, i) = cr_tool.z;
  }

  return J_tool;
}

inline void compute_constraints(double* result, const Array& x, Optimization* opt) {
  double* moving_result = result;

  opt->bspline.compute_trajectory(x, opt->task);

  ObjMatrix<Capsule> capsules(opt->manip._n_caps, (u32) opt->bspline.n_points / opt->constraints.n_collision_skip);
  int                test_idx = 0; // todo: Comment this out

  // {
  //     blast_time_block("Points Constraints");
  for (int i = 0; i < opt->bspline.n_points; i++) {
    // This is calculated every loop
    auto p = opt->bspline.traj.pos.col(i);
    Assert(p.is_alias);

    forward_kinematics(opt->manip, p);                                                  // fills _Q, _Q_mult, and _p_j

    opt->manip.compute_capsules();                                                      // fills _capsule_list
    if (opt->constraints.external_collisions) {
      if (i % opt->constraints.n_collision_skip == 0 && i != opt->bspline.n_points - 1) { // todo: this cleaner
        for (int j = 0; j < opt->manip._n_caps; j++) {
          capsules(j, i / (int) opt->constraints.n_collision_skip) = opt->manip._capsule_list[j];
        }
      }
    }

    if (opt->constraints.position) {
      for (int j = 0; j < opt->manip.joints; j++) {
        moving_result[0] = bound_constraint(1.01 * opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        moving_result++;
        test_idx++;
      }
    }

    if (opt->constraints.velocity) {
      for (int j = 0; j < opt->manip.joints; j++) {
        moving_result[0] = bound_constraint(1.01 * opt->bspline.traj.vel(j, i), opt->manip.vmin[j], opt->manip.vmax[j]);
        moving_result++;
        test_idx++;
      }
    }

    if (opt->constraints.acceleration) {
      for (int j = 0; j < opt->manip.joints; j++) {
        moving_result[0] = bound_constraint(1.01 * opt->bspline.traj.acc(j, i), opt->manip.amin[j], opt->manip.amax[j]);
        moving_result++;
        test_idx++;
      }
    }

    if (opt->constraints.torque) {
      auto vel = opt->bspline.traj.vel.col(i);
      auto acc = opt->bspline.traj.acc.col(i);
      dynamics(opt->manip, vel, acc); // fills _efforts

      for (int j = 0; j < opt->manip.joints; j++) {
        moving_result[0] = bound_constraint(1.01 * opt->manip._efforts[j], opt->manip.tau_min[j], opt->manip.tau_max[j]);
        moving_result++;
        test_idx++;
      }
    }

    if (opt->constraints.tcp_speed) {
      auto J_tool      = get_J_tool(opt);
      real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
      moving_result[0] = bound_constraint(1.01 * tcp_speed, 0.0, opt->manip.tcp_max);
      moving_result++;
      test_idx++;
    }

    if (opt->constraints.self_collisions) {
      auto tmp_coll = opt->manip.get_internal_collisions();
      for (u32 j = 0; j < tmp_coll.size; j++)
        moving_result[j] = -tmp_coll[j] + 0.01; //*std::abs(tmp_coll[j]);
      moving_result += tmp_coll.size;
      test_idx += (int) tmp_coll.size;
    }
  }
  // }
  // {
  //     blast_time_block("External Collision Constraints");
  if (opt->constraints.external_collisions) {
    auto collisions = test_collisions(capsules, &(opt->world), opt->constraints.n_collision_constraints, opt->trajectory_start_time, opt->trajectory_start_time + x.back());
    for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
      moving_result[j] = -collisions[j] + 0.01; //*std::abs(collisions[j]);
    }
    moving_result += opt->constraints.n_collision_constraints;
  }
  // }

  for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) { // todo: split in extra_constraint per point or once
    auto extra_constraint = opt->constraints.extra_constraints[i];
    extra_constraint(moving_result, opt);
    moving_result += opt->constraints.n_extra_constraints[i];
  }
}

inline void nlopt_constraints(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, x_len);
  compute_constraints(result, xv, opt);

  // gradients calculation
  if (grad) {
    Array x_plus(x_len);
    Array r_plus(m);
    for (u32 j = 0; j < x_len; j++) {
      constexpr real eps = 1e-5;
      memcpy(x_plus.data, x, x_len * sizeof(real));
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, x_len);

      // compute constraints
      compute_constraints(r_plus.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * x_len + j] = (r_plus[i] - result[i]) / eps;
    }
    if (opt->constraints.show_info) { // when more info is needed per iteration
      Matrix gradients(x_len, m);
      Array  constraints(m);
      for (u32 j = 0; j < x_len; j++) {
        for (u32 i = 0; i < m; i++) {
          gradients(j, i) = grad[i * x_len + j];
          constraints[i]  = result[i];
        }
      }
      opt->constraints.grad_list.push_back(gradients);
      opt->constraints.constr_list.push_back(constraints);
      opt->constraints.x_list.push_back(xv);
    }
  }
}

inline bool validate_task(Optimization* opt) {
  // int n_constraints = opt->constraints.n_constraints;

  // todo: why copy??
  auto manip = opt->manip;

  Matrix pos(opt->manip.joints, 2);
  Matrix vel(opt->manip.joints, 2);
  Matrix acc(opt->manip.joints, 2);
  for (int i = 0; i < opt->task.rows; i++) {
    pos(i, 0) = opt->task(i, 0);
    pos(i, 1) = opt->task(i, 3);
    vel(i, 0) = opt->task(i, 1);
    vel(i, 1) = opt->task(i, 4);
    acc(i, 0) = opt->task(i, 2);
    acc(i, 1) = opt->task(i, 5);
  }
  real result;

  ObjMatrix<Capsule> capsules_begin(opt->manip._n_caps, 1);
  ObjMatrix<Capsule> capsules_end(opt->manip._n_caps, 1);
  for (int i = 0; i < 2; i++) {
    auto p = pos.col(i);
    Assert(p.is_alias);

    forward_kinematics(manip, p); // fills _Q, _Q_mult, and _p_j
    if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
      manip.compute_capsules();   // fills _capsule_list
      if (opt->constraints.external_collisions) {
        if (i == 0) {
          for (u32 j = 0; j < manip._n_caps; j++) {
            capsules_begin(j, 0) = manip._capsule_list[j];
          }
        }
        if (i == 1) {
          for (u32 j = 0; j < manip._n_caps; j++) {
            capsules_end(j, 0) = manip._capsule_list[j];
          }
        }
      }
      if (opt->constraints.self_collisions) {
        if (min(manip.get_internal_collisions()) < 0) { // min because collisions constraints are -d < 0
          std::cout << "Self-collision at start/end position." << std::endl;
          return false;
        }
      }
    }

    if (opt->constraints.position) {
      for (int j = 0; j < manip.joints; j++) {
        result = bound_constraint(pos(j, i), manip.pmin[j], manip.pmax[j]);
        if (result > 0) {
          std::cout << "Position outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (opt->constraints.velocity) {
      for (int j = 0; j < manip.joints; j++) {
        result = bound_constraint(vel(j, i), manip.vmin[j], manip.vmax[j]);
        if (result > 0) {
          std::cout << "Velocity outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (opt->constraints.acceleration) {
      for (int j = 0; j < manip.joints; j++) {
        result = bound_constraint(acc(j, i), manip.amin[j], manip.amax[j]);
        if (result > 0) {
          std::cout << "Acceleration outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (opt->constraints.torque) {
      auto vel_tmp = vel.col(i);
      auto acc_tmp = acc.col(i);
      dynamics(manip, vel_tmp, acc_tmp); // fills _efforts

      for (int j = 0; j < manip.joints; j++) {
        result = bound_constraint(manip._efforts[j], manip.tau_min[j], manip.tau_max[j]);
        if (result > 0) {
          std::cout << "Torque outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (opt->constraints.tcp_speed) {
      auto J_tool    = get_J_tool(opt);
      real tcp_speed = norm(J_tool * vel.col(i));
      result         = bound_constraint(tcp_speed, -INF_REAL, manip.tcp_max);
      if (result > 0) {
        std::cout << "TCP speed outside bounds." << std::endl;
        return false;
      }
    }
  }

  if (opt->constraints.external_collisions) {
    auto collisions_begin = test_collisions(capsules_begin, &(opt->world), opt->constraints.n_collision_constraints, opt->trajectory_start_time, opt->trajectory_start_time);
    if (min(collisions_begin) < 0) { // min because collisions constraints are -d < 0
      std::cout << "External collision at start position." << std::endl;
      return false;
    }
    // End of trajectory only takes static obstacles in account since we do not know the finish time
    World world_end;
    world_end.boxes     = opt->world.boxes;
    world_end.capsules  = opt->world.capsules;
    world_end.spheres   = opt->world.spheres;
    world_end.size      = world_end.boxes.size() + world_end.capsules.size() + world_end.spheres.size();
    auto collisions_end = test_collisions(capsules_end, &world_end, opt->constraints.n_collision_constraints, opt->trajectory_start_time, opt->trajectory_start_time);
    if (min(collisions_end) < 0) { // min because collisions constraints are -d < 0
      std::cout << "External collision at end position." << std::endl;
      return false;
    }
  }

  for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) { // todo: split in extra_constraint per point or once
    Array moving_result(opt->constraints.n_extra_constraints[i]);
    auto  extra_constraint = opt->constraints.extra_constraints[i];
    Array extra_result(opt->constraints.n_extra_constraints[i]);
    extra_constraint(extra_result.data, opt);
    if (max(extra_result) > 0) {
      std::cout << "Extra constraint " << i + 1 << " inadmissible at start/end position." << std::endl;
      return false;
    }
  }
  return true;
}

} // namespace blast
