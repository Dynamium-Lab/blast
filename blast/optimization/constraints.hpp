#pragma once

#include <blast>
#include "optimization/optimization.hpp"

namespace blast {


inline blast_fn void ConstraintSelection::add_constraint(ConstraintFunction constraint, int n_con) {
  extra_constraints.push_back(constraint);
  n_extra_constraints.push_back(n_con);
}

inline blast_fn real bound_constraint(const real& value, const real& value_min, const real& value_max) {
#if BLAST_TRACE_LEVEL >= 3
  ZoneScoped;
#endif
  // todo: remove INF_REAL from constraints at initialization
  if (value_max == INF_REAL || value_min == -INF_REAL) { // todo: fix for one is INF and not the other
    return -1.0;
  }
  real center     = (value_max + value_min) / 2;
  real half_range = (value_max - value_min) / 2;

  return std::abs(value - center) / half_range - 1.0;
}

inline blast_fn real abs_constraint(const real& value, const real& value_max) {
  return (std::abs(value) - value_max) / value_max;
}

inline blast_fn Matrix get_J_tool(const Optimization* opt, const ManipulatorTempData& temp) {
  std::vector<Vec3> r_tool(opt->manip.n_joints);
  r_tool[opt->manip.n_joints - 1] = opt->manip.dv[opt->manip.n_joints - 1];
  for (int i = (int) opt->manip.n_joints - 2; i >= 0; i--) {
    r_tool[i] = opt->manip.dv[i] + temp.rotations[i + 1] * r_tool[i + 1];
  }

  for (int i = 0; i < opt->manip.n_joints; i++) {
    r_tool[i] = temp.rotations_mult[i] * r_tool[i];
  }

  Matrix J_tool(3, opt->manip.n_joints);
  for (int i = 0; i < opt->manip.n_joints; i++) {
    // e = opt->manip._Q_mult[i] * opt->manip.ev[i]; // replaced e directly in function to skip copy
    Vec3 cr_tool = cross(temp.rotations_mult[i] * opt->manip.ev[i], r_tool[i]);
    J_tool(0, i) = cr_tool.x;
    J_tool(1, i) = cr_tool.y;
    J_tool(2, i) = cr_tool.z;
  }

  return J_tool;
}

inline blast_fn void constraints_with_segments(const Array& x, Optimization& opt, Array& constraints, Matrix& grad) {
  // constraints (p,v,a,tor) for each joint, for each segment
  // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]


  // basis: n_ctrl x n_points
  Assert(constraints.is_alias);
  if (grad.size) {
    Assert(grad.is_alias);
    Assert(grad.rows == x.size);
    Assert(grad.cols == constraints.size);
  }
  const int n_segments                = (int) opt.bspline.n_ctrl - (int) opt.bspline.p;
  const int n_points_per_segment      = (int) opt.bspline.n_points / n_segments; // todo: check if fine?
  const int n_joints                  = (int) opt.manip.n_joints;
  const int n_ctrl                    = (int) opt.bspline.n_ctrl;
  const int x_len                     = (int) x.size;
  const int n_constraints_per_segment = (n_joints * 4); //  todo: remove hard-coded 4
  Assert(constraints.size == n_segments * n_constraints_per_segment);

  // limits
  Array       pmax    = opt.manip.pmax;
  Array       pmin    = opt.manip.pmin;
  const Array vmax    = opt.manip.vmax;
  const Array amax    = opt.manip.amax;
  const Array tau_max = opt.manip.tau_max;

  ManipulatorTempData manip_data;
  std::vector<u8>     max_pos_indices(n_joints);
  std::vector<u8>     max_vel_indices(n_joints);
  std::vector<u8>     max_acc_indices(n_joints);
  std::vector<u8>     max_tor_indices(n_joints);
  Array               max_pos_constraints(n_joints);
  Array               max_vel_constraints(n_joints);
  Array               max_acc_constraints(n_joints);
  Array               max_tor_constraints(n_joints);


  opt.bspline.compute_trajectory(x, opt.task);


  for (int segment = 0; segment < n_segments; segment++) {
#if BLAST_TRACE_LEVEL >= 2
    ZoneScopedN("All Segment Constraints");
#endif
    const int first_affected_control_point = std::max(3, segment);
    const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.p);
    const int n_affected_control_points    = last_affected_control_point - first_affected_control_point;
    Assert(n_affected_control_points >= 3);
    Assert(n_affected_control_points <= 6);

    Matrix bp(&opt.bspline.basis_p(0, segment * n_points_per_segment), n_ctrl, n_points_per_segment);
    Matrix bv(&opt.bspline.basis_v(0, segment * n_points_per_segment), n_ctrl, n_points_per_segment);
    Matrix ba(&opt.bspline.basis_a(0, segment * n_points_per_segment), n_ctrl, n_points_per_segment);

    for (int point = segment * n_points_per_segment; point < (segment + 1) * n_points_per_segment; point++) {
      auto p = opt.bspline.traj.pos.col(point);
      auto v = opt.bspline.traj.vel.col(point);
      auto a = opt.bspline.traj.acc.col(point);

      forward_kinematics(opt.manip, manip_data, p);
      dynamics(opt.manip, manip_data, v, a);

      for (int j = 0; j < n_joints; j++) {
        // position
        if (auto c = bound_constraint(p[j], pmin[j], pmax[j]);
            c > max_pos_constraints[j]) {
          max_pos_constraints[j] = c;
          max_pos_indices[j]     = point;
        }
        // velocity
        if (auto c = std::abs(v[j]) / vmax[j] - 1.0;
            c > max_vel_constraints[j]) {
          max_vel_constraints[j] = c;
          max_vel_indices[j]     = point;
        }
        // acceleration
        if (auto c = std::abs(a[j]) / amax[j] - 1.0;
            c > max_acc_constraints[j]) {
          max_acc_constraints[j] = c;
          max_acc_indices[j]     = point;
        }
        // torque
        if (auto c = std::abs(manip_data.efforts[j]) / tau_max[j] - 1.0;
            c > max_tor_constraints[j]) {
          max_tor_constraints[j] = c;
          max_tor_indices[j]     = point;
        }
      }
    }

    // at this point we have max constraints for pos, vel, acc, tor, for this segment

    // fill in the constraints for the current segment
    // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
    auto fill_idx = segment * n_constraints_per_segment;
    std::copy_n(max_pos_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints; // note (andre): we can use the comma operator because we don't need the result of copy_n()
    std::copy_n(max_vel_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    std::copy_n(max_acc_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    std::copy_n(max_tor_constraints.data, n_joints, &constraints[fill_idx]);


    // The gradient should be a (x_len)x(n_constraints) matrix that looks like this:
    // [dp1/dx1, dp2/dx1, ..., dv1/dx1, dv2/dx1, ..., da1/dx1, da2/dx1, ..., dt1/dx1, dt2/dx1]
    // [dp1/dx2, dp2/dx2, ..., dv1/dx2, dv2/dx2, ..., da1/dx2, da2/dx2, ..., dt1/dx2, dt2/dx2]
    // [dp1/dx3, dp2/dx3, ..., dv1/dx3, dv2/dx3, ..., da1/dx3, da2/dx3, ..., dt1/dx3, dt2/dx3]
    // [dp1/dx4, dp2/dx4, ..., dv1/dx4, dv2/dx4, ..., da1/dx4, da2/dx4, ..., dt1/dx4, dt2/dx4]
    // [dp1/dx5, dp2/dx5, ..., dv1/dx5, dv2/dx5, ..., da1/dx5, da2/dx5, ..., dt1/dx5, dt2/dx5]
    // [dp1/dx6, dp2/dx6, ..., dv1/dx6, dv2/dx6, ..., da1/dx6, da2/dx6, ..., dt1/dx6, dt2/dx6]
    // [.....................]
    // [.....................]
    // [.....................]
    // [.....................]
    // [dp1/dT=0,dp2/dT=0,..., dv1/dT , dv1/dT , ..., da1/dT , da2/dT , ..., dt1/dT , dt2/dT ]
    // *** per segment ***
    // where x is the optimization vector
    if (grad.size) {
      // Matrix (alias) in which we can insert the gradient for the current segment
      Matrix grad_segment(&grad(0, segment * n_constraints_per_segment), x_len, n_constraints_per_segment);
      Assert(grad_segment.is_alias);

      int con = 0;

      // positions
      for (int joint = 0; joint < n_joints; joint++) {
        {
          // todo (andre): do this when initializing the optimization, not here
          // todo: document the current behaviour in the API
          //        (doesn't currently work if one is inf and the other is not)
          if (pmax[joint] == INF_REAL)  // note: replace INF_REAL with huge value
            pmax[joint] = 1e300;
          if (pmin[joint] == -INF_REAL) // note: replace -INF_REAL with huge negative value
            pmin[joint] = -1e300;
        }

        // Array of the column where to put the gradient for the current constraint
        Array fill_column = grad_segment.col(con);
        Assert(con == (segment * n_constraints_per_segment + joint));
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // Which values in 'x' affect the current joint's position
        auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
        // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
        x_idx += first_affected_control_point - 3;

        real center     = (pmax[joint] + pmin[joint]) / 2;
        real half_range = (pmax[joint] - pmin[joint]) / 2;

        // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
        Array bp_to_use(&bp(first_affected_control_point, max_pos_indices[joint]), n_affected_control_points);
        Assert(bp_to_use.is_alias);

        real coeff = 2.0 * std::abs(opt.bspline.traj.pos(joint, max_pos_indices[joint])) / (pmax[joint] - pmin[joint]);

        for (int i = 0; i < n_affected_control_points; i++) {
          fill_column[x_idx++] = bp_to_use[i] * coeff;
        }

        con++;
      }

      // velocities
      auto one_over_T = 1 / opt.bspline.traj.t.back();
      for (int joint = 0; joint < n_joints; joint++) {
        // Array of the column where to put the gradient for the current constraint
        Array fill_column = grad_segment.col(con);
        Assert(con == (segment * n_constraints_per_segment + n_joints + joint));
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // Which values in 'x' affect the current joint's position
        auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
        // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
        x_idx += first_affected_control_point - 3;

        Array bv_to_use(&bv(first_affected_control_point, max_vel_indices[joint]), n_affected_control_points);
        Assert(bv_to_use.is_alias);

        real coeff = sign(opt.bspline.traj.vel(joint, max_vel_indices[joint])) / vmax[joint] * one_over_T;

        for (int i = 0; i < n_affected_control_points; i++) {
          fill_column[x_idx++] = bv_to_use[i] * coeff;
        }
        con++;
      }

      // dvj/dT --> if T increases, abs(v) drops. If T == 1, T doubles, v halves. If T increases by 25%
      // todo: deal with T

      // accelerations
      auto one_over_T2 = one_over_T * one_over_T;
      for (int joint = 0; joint < n_joints; joint++) {
        // Array of the column where to put the gradient for the current constraint
        Array fill_column = grad_segment.col(con);
        Assert(con == (segment * n_constraints_per_segment + 2 * n_joints + joint));
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // Which values in 'x' affect the current joint's position
        auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
        // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
        x_idx += first_affected_control_point - 3;

        Array ba_to_use(&ba(first_affected_control_point, max_acc_indices[joint]), n_affected_control_points);
        Assert(ba_to_use.is_alias);

        real coeff = sign(opt.bspline.traj.acc(joint, max_acc_indices[joint])) / amax[joint] * one_over_T2;

        for (int i = 0; i < n_affected_control_points; i++) {
          fill_column[x_idx++] = ba_to_use[i] * coeff;
        }
        con++;
      }

      // todo: deal with T

      // torque
      real eps = 1e-5;

      // [dt0/dp0, dt1/dp0, ..., dt4/dp0, dt5/dp0]
      // [dt0/dp1, dt1/dp1, ..., dt4/dp1, dt5/dp1]
      // [dt0/dp2, dt1/dp2, ..., dt4/dp2, dt5/dp2]
      // [dt0/dp3, dt1/dp3, ..., dt4/dp3, dt5/dp3]
      Matrix grad_tau_p(n_joints, n_joints); // variation of torque[joint] in respect to delta p[joint]
      Matrix grad_tau_v(n_joints, n_joints); // variation of torque[joint] in respect to delta v[joint]
      Matrix grad_tau_a(n_joints, n_joints); // variation of torque[joint] in respect to delta a[joint]
      for (int joint = 0; joint < n_joints; joint++) {
        auto p_at_max_tor  = opt.bspline.traj.pos.col(max_tor_indices[joint]);
        auto v_at_max_tor  = opt.bspline.traj.vel.col(max_tor_indices[joint]);
        auto a_at_max_tore = opt.bspline.traj.acc.col(max_tor_indices[joint]);

        // variation in pos
        p_at_max_tor[joint] += eps;
        forward_kinematics(opt.manip, manip_data, p_at_max_tor);
        dynamics(opt.manip, manip_data, v_at_max_tor, a_at_max_tore);
        p_at_max_tor[joint] -= eps;

        for (int j = 0; j < n_joints; j++) {
          auto max_tor_plus    = (std::abs(manip_data.efforts[j])) / tau_max[j] - 1.0;
          grad_tau_p(j, joint) = (max_tor_plus - max_tor_constraints[j]) / eps;
        }

        // restore to original positions (only once, since vel and acc do not affect fk)
        forward_kinematics(opt.manip, manip_data, p_at_max_tor);

        // variation in vel
        v_at_max_tor[joint] += eps;
        dynamics(opt.manip, manip_data, v_at_max_tor, a_at_max_tore);
        v_at_max_tor[joint] -= eps;

        for (int j = 0; j < n_joints; j++) {
          auto max_tor_plus    = (std::abs(manip_data.efforts[j])) / tau_max[j] - 1.0;
          grad_tau_v(j, joint) = (max_tor_plus - max_tor_constraints[j]) / eps;
        }

        // variation in acc
        a_at_max_tore[joint] += eps;
        dynamics(opt.manip, manip_data, v_at_max_tor, a_at_max_tore);
        a_at_max_tore[joint] -= eps;

        for (int j = 0; j < n_joints; j++) {
          auto max_tor_plus    = (std::abs(manip_data.efforts[j])) / tau_max[j] - 1.0;
          grad_tau_a(j, joint) = (max_tor_plus - max_tor_constraints[j]) / eps;
        }
      }

      // Which values in 'x' affect the torque (these are the same for each joint)
      auto x_idx = first_affected_control_point - 3; // todo: does not work with tasks that don't impose p,v,a for every joint!!
      // skip the control points unaffected in the segment
      auto x_idx_skip = (n_ctrl - 6) - n_affected_control_points;
      for (int joint = 0; joint < n_joints; joint++) {
        // Array of the column where to put the gradient for the current constraint
        Array fill_column = grad_segment.col(con);
        Assert(con == (segment * n_constraints_per_segment + 2 * n_joints + joint));
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        Array bp_to_use(&bp(first_affected_control_point, max_tor_indices[joint]), n_affected_control_points);
        Assert(bp_to_use.is_alias);
        Array bv_to_use(&bv(first_affected_control_point, max_tor_indices[joint]), n_affected_control_points);
        Assert(bv_to_use.is_alias);
        Array ba_to_use(&ba(first_affected_control_point, max_tor_indices[joint]), n_affected_control_points);
        Assert(ba_to_use.is_alias);
        for (int j = 0; j < n_joints; j++) {
          // fill 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
          auto tmp_grad_p = grad_tau_p(j, joint);
          auto tmp_grad_v = grad_tau_v(j, joint);
          auto tmp_grad_a = grad_tau_a(j, joint);
          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx++] = bp_to_use[i] * tmp_grad_p +
                                   bv_to_use[i] * tmp_grad_v * one_over_T +
                                   ba_to_use[i] * tmp_grad_a * one_over_T2;
          }
          x_idx += x_idx_skip;
        }
        con++;
      }

      // todo: deal with T

      // --------------------- below are notes only -------------------------
      // dp1/dx (rowidx - rowidx+n_active_control_points), col_idx ---> fetch new basis function column (basis_p)
      // dp2/dx (rowidx - rowidx+n_active_control_points), col_idx+1 ---> fetch new basis function column (basis_p)

      // dv1/dx (rowidx - rowidx+n_active_control_points), col_idx+1 ---> fetch new basis function column (basis_v)

      // dt1/dx ??? 6xn_joints control points max


      // dp1/dc = basis_p.col(max_pos_indices[0])
    }
  }
}

inline blast_fn void compute_constraints(double* result, const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  ZoneScoped;
#endif

  double* moving_result = result;

  // todo: compute inside loop so that every worker can compute independently??
  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // todo: no copy every iteration
  ObjMatrix<Capsule> capsules(opt->manip._n_caps, (int) opt->bspline.n_points / (int) opt->constraints.n_collision_skip);

  // todo: compute n_con_per_point() so that every worker knows where to put their result

  for (int i = 0; i < opt->bspline.n_points; i++) {
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
        if (opt->constraints.external_collisions) {
          if (i % opt->constraints.n_collision_skip == 0 && i != opt->bspline.n_points - 1) {
            for (int j = 0; j < opt->manip._n_caps; j++) {
              capsules(j, i / (int) opt->constraints.n_collision_skip) = manip_data.capsule_list[j];
            }
          }
        }
      }
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Position");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = bound_constraint(1.01 * opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        moving_result++;
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        // moving_result[0] = abs_constraint(1.01 * opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
        moving_result[0] = (std::abs(1.01 * opt->bspline.traj.vel(j, i)) - opt->manip.vmax[j]) / opt->manip.vmax[j];
        moving_result++;
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        moving_result[0] = (std::abs(1.01 * opt->bspline.traj.acc(j, i)) - opt->manip.amax[j]) / opt->manip.amax[j];
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
        moving_result[0] = (std::abs(1.01 * manip_data.efforts[j]) - opt->manip.tau_max[j]) / opt->manip.tau_max[j];
        moving_result++;
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
      moving_result[0] = bound_constraint(1.01 * tcp_speed, 0.0, opt->manip.tcp_max);
      moving_result++;
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
      for (u32 j = 0; j < tmp_coll.size; j++)
        moving_result[j] = -tmp_coll[j] + 0.01; //*std::abs(tmp_coll[j]);
      moving_result += tmp_coll.size;
    }
  }

  if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ExternalCollisions");
#endif
    auto collisions = test_collisions(capsules, &(opt->world), opt->constraints.n_collision_constraints, opt->trajectory_start_time, opt->trajectory_start_time + x.back());
    for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
      moving_result[j] = -collisions[j] + 0.01; //*std::abs(collisions[j]);
    }
    moving_result += opt->constraints.n_collision_constraints;
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

inline blast_fn void nlopt_constraints(unsigned m, double* result, unsigned x_len, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, x_len);
  compute_constraints(result, xv, opt);

  // gradients calculation
  if (grad) {

    // todo: no construction??
    Array x_plus(x_len);
    Array r_plus(m);

    for (u32 j = 0; j < x_len; j++) {
      constexpr real eps = 1e-5;
      // todo: only copy value that changed last j
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

inline blast_fn bool validate_task(Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto& manip       = opt->manip;
  auto& constraints = opt->constraints;
  auto& task        = opt->task;

  ManipulatorTempData manip_data;

  Matrix pos(manip.n_joints, 2);
  Matrix vel(manip.n_joints, 2);
  Matrix acc(manip.n_joints, 2);
  for (int i = 0; i < task.rows; i++) {
    pos(i, 0) = task(i, 0);
    pos(i, 1) = task(i, 3);
    vel(i, 0) = task(i, 1);
    vel(i, 1) = task(i, 4);
    acc(i, 0) = task(i, 2);
    acc(i, 1) = task(i, 5);
  }
  real result;

  ObjMatrix<Capsule> capsules_begin(manip._n_caps, 1);
  ObjMatrix<Capsule> capsules_end(manip._n_caps, 1);
  for (int i = 0; i < 2; i++) {
    auto p = pos.col(i);
    Assert(p.is_alias);

    forward_kinematics(manip, manip_data, p); // fills _Q, _Q_mult, and _p_j
    if (constraints.external_collisions || constraints.self_collisions) {
      compute_capsules(manip, manip_data);
      if (constraints.external_collisions) {
        if (i == 0) {
          for (int j = 0; j < manip._n_caps; j++) {
            capsules_begin(j, 0) = manip_data.capsule_list[j];
          }
        }
        if (i == 1) {
          for (int j = 0; j < manip._n_caps; j++) {
            capsules_end(j, 0) = manip_data.capsule_list[j];
          }
        }
      }
      if (constraints.self_collisions) {
        if (auto tmp_coll = get_internal_collisions(manip, manip_data); min(tmp_coll) < 0) { // min because collisions constraints are -d < 0
          std::cout << "Self-collision at start/end position." << std::endl;
          return false;
        }
      }
    }

    if (constraints.position) {
      for (int j = 0; j < manip.n_joints; j++) {
        result = bound_constraint(pos(j, i), manip.pmin[j], manip.pmax[j]);
        if (result > 0) {
          std::cout << "Position outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (constraints.velocity) {
      for (int j = 0; j < manip.n_joints; j++) {
        result = abs_constraint(vel(j, i), manip.vmax[j]);
        if (result > 0) {
          std::cout << "Velocity outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (constraints.acceleration) {
      for (int j = 0; j < manip.n_joints; j++) {
        result = abs_constraint(acc(j, i), manip.amax[j]);
        if (result > 0) {
          std::cout << "Acceleration outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (constraints.torque) {
      auto vel_tmp = vel.col(i);
      auto acc_tmp = acc.col(i);
      dynamics(manip, manip_data, vel_tmp, acc_tmp); // fills _efforts

      for (int j = 0; j < manip.n_joints; j++) {
        result = abs_constraint(manip_data.efforts[j], manip.tau_max[j]);
        if (result > 0) {
          std::cout << "Torque outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (constraints.tcp_speed) {
      auto J_tool    = get_J_tool(opt, manip_data);
      real tcp_speed = norm(J_tool * vel.col(i));
      result         = bound_constraint(tcp_speed, -INF_REAL, manip.tcp_max);
      if (result > 0) {
        std::cout << "TCP speed outside bounds." << std::endl;
        return false;
      }
    }
  }

  if (constraints.external_collisions) {
    auto collisions_begin = test_collisions(capsules_begin, &(opt->world), constraints.n_collision_constraints, opt->trajectory_start_time, opt->trajectory_start_time);
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
    auto collisions_end = test_collisions(capsules_end, &world_end, constraints.n_collision_constraints, opt->trajectory_start_time, opt->trajectory_start_time);
    if (min(collisions_end) < 0) { // min because collisions constraints are -d < 0
      std::cout << "External collision at end position." << std::endl;
      return false;
    }
  }

  for (u32 i = 0; i < constraints.extra_constraints.size(); i++) { // todo: split in extra_constraint per point or once
    Array moving_result(constraints.n_extra_constraints[i]);
    auto  extra_constraint = constraints.extra_constraints[i];
    Array extra_result(constraints.n_extra_constraints[i]);
    extra_constraint(extra_result.data, opt);
    if (max(extra_result) > 0) {
      std::cout << "Extra constraint " << i + 1 << " inadmissible at start/end position." << std::endl;
      return false;
    }
  }
  return true;
}

// ------------------------- Accelerated functions --------------------------------
inline blast_fn std::tuple<real, real> abs_constraint_dev(const real& q, const real& q_max) {
  return std::make_tuple((std::abs(q) - q_max) / q_max, sign(q) * 1 / q_max);
}

template<bool is_grad>
blast_fn std::tuple<real, real> bound_constraint_dev(const real& q, const real& q_min, const real& q_max) {
  // todo: remove INF_REAL from constraints at initialization
  if (q_max == INF_REAL || q_min == -INF_REAL)
    return std::make_tuple(-1.0, 0.0);
  else {
    real c          = (q_max + q_min) / 2;
    real q_prime    = q - c;
    real b          = (q_max - q_min) / 2;
    real constraint = (std::abs(q_prime) - b) / b;

    // Conditionally compute gradient coefficient
    real gradient_coeff = 0.0;
    if constexpr (is_grad) {
      gradient_coeff = sign(q_prime) * 1 / b;
    }
    return std::make_tuple(constraint, gradient_coeff);
  }
}

template<bool is_grad>
blast_fn void compute_constraints_dev(double* result, Array& gradient_coeffs, const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  ZoneScoped;
#endif

  double* moving_result = result;
  u32     grad_idx      = 0;
  Array   external_collisions(opt->bspline.n_points);

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // Lambda to process common bound constraint operations
  auto process_bound = [&](real value, real bound_max) {
    auto [constraint, gradient_coeff] = abs_constraint_dev(1.01 * value, bound_max);
    *moving_result++                  = constraint;
    gradient_coeffs[grad_idx++]       = gradient_coeff;
  };

  for (u32 i = 0; i < opt->bspline.n_points; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    ManipulatorTempData manip_data;

    auto pos = opt->bspline.traj.pos.col(i);
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
      for (int j = 0; j < opt->manip.n_joints; j++) {
        auto [constraint, gradient_coeff] = bound_constraint_dev<is_grad>(1.01 * opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
        *moving_result++                  = constraint;
        gradient_coeffs[grad_idx++]       = gradient_coeff;
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Velocity");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        process_bound(opt->bspline.traj.vel(j, i), opt->manip.vmax[j]);
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Acceleration");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        process_bound(opt->bspline.traj.acc(j, i), opt->manip.amax[j]);
      }
    }


    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      auto vel = opt->bspline.traj.vel.col(i);
      auto acc = opt->bspline.traj.acc.col(i);
      dynamics(opt->manip, manip_data, vel, acc);
      for (int j = 0; j < opt->manip.n_joints; j++) {
        *moving_result++ = abs_constraint(1.01 * manip_data.efforts[j], opt->manip.tau_max[j]);
        grad_idx++; // todo: add analytical gradients
      }
    }

    if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
      *moving_result++ = bound_constraint(1.01 * tcp_speed, 0.0, opt->manip.tcp_max);
      grad_idx++; // todo: add analytical gradients
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Capsules");
#endif
      if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
        compute_capsules(opt->manip, manip_data);
      }
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
      for (u32 j = 0; j < tmp_coll.size; j++) {
        moving_result[j] = -tmp_coll[j] + 0.01;
      }
      moving_result += tmp_coll.size;
      grad_idx += tmp_coll.size; // todo: add analytical gradients
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisionsCalculate");
#endif
      external_collisions[i] = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
    }
  }

  if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ExternalCollisionsSort");
#endif
    Array dist_min(opt->constraints.n_collision_constraints, INF_REAL);
    for (u32 i = 0; i < opt->bspline.n_points; i++) {
      for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
        if (external_collisions[i] < dist_min[j]) {
          for (int k = (int) opt->constraints.n_collision_constraints - 1; k > j; k--) {
            dist_min[k] = dist_min[k - 1];
          }
          dist_min[j] = external_collisions[i];
          break;
        }
      }
    }
    for (u32 i = 0; i < opt->constraints.n_collision_constraints; i++) {
      moving_result[i] = -dist_min[i] + 0.01; //*std::abs(collisions[j]);
    }
    moving_result += opt->constraints.n_collision_constraints;
    grad_idx += opt->constraints.n_collision_constraints; // todo: add analytical gradients
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

template<bool is_grad>
blast_fn void compute_constraints_dev_new(double* result, Array& gradient_coeffs, Matrix& gradient_coeffs_torque, /*Array& gradient_coeffs_collisions,*/ const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  ZoneScoped;
#endif

  double* moving_result   = result;
  u32     grad_idx        = 0;
  u32     grad_torque_idx = 0;
  Array   external_collisions(opt->bspline.n_points);
  Array   external_collisions_plus(opt->bspline.n_points);

  auto joints = opt->manip.n_joints;

  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // Lambda to process common bound constraint operations
  auto process_bound = [&](real value, real bound_max) {
    auto [constraint, gradient_coeff] = abs_constraint_dev(1.01 * value, bound_max);
    *moving_result++                  = constraint;
    gradient_coeffs[grad_idx++]       = gradient_coeff;
  };

  for (u32 i = 0; i < opt->bspline.n_points; i++) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ConstraintSinglePoint");
#endif

    ManipulatorTempData manip_data;
    Array               pos_plus(joints);
    Array               vel_plus(joints);
    Array               acc_plus(joints);
    Array               torque_constraint(joints);
    Array               torque_constraint_plus(joints);

    constexpr real eps = 1e-5;

    auto pos = opt->bspline.traj.pos.col(i);
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
        auto [constraint, gradient_coeff] = bound_constraint_dev<is_grad>(1.01 * opt->bspline.traj.pos(j, i), opt->manip.pmin[j], opt->manip.pmax[j]);
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


    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Dynamics");
#endif
      auto vel = opt->bspline.traj.vel.col(i);
      auto acc = opt->bspline.traj.acc.col(i);
      dynamics(opt->manip, manip_data, vel, acc);
      for (int j = 0; j < joints; j++) {
        torque_constraint[j] = abs_constraint(1.01 * manip_data.efforts[j], opt->manip.tau_max[j]);
        *moving_result++     = torque_constraint[j];
      }

      // grad_coeffs pos
      for (u32 j = 0; j < joints; j++) {
        pos[j] += eps;
        forward_kinematics(opt->manip, manip_data, pos);
        dynamics(opt->manip, manip_data, vel, acc);
        pos[j] -= eps;

        for (u32 k = 0; k < joints; k++) {
          torque_constraint_plus[k]                     = abs_constraint(1.01 * manip_data.efforts[k], opt->manip.tau_max[k]);
          gradient_coeffs_torque(k, j + 3 * joints * i) = (torque_constraint_plus[k] - torque_constraint[k]) / eps;
        }
      }

      forward_kinematics(opt->manip, manip_data, pos);
      // grad_coeffs vel
      for (u32 j = 0; j < joints; j++) {
        vel[j] += eps;
        dynamics(opt->manip, manip_data, vel, acc);
        vel[j] -= eps;

        for (u32 k = 0; k < joints; k++) {
          torque_constraint_plus[k]                              = abs_constraint(1.01 * manip_data.efforts[k], opt->manip.tau_max[k]);
          gradient_coeffs_torque(k, j + joints + 3 * joints * i) = (torque_constraint_plus[k] - torque_constraint[k]) / eps;
        }
      }

      // grad_coeffs acc
      for (u32 j = 0; j < joints; j++) {
        acc[j] += eps;
        dynamics(opt->manip, manip_data, vel, acc);
        acc[j] -= eps;

        for (u32 k = 0; k < joints; k++) {
          torque_constraint_plus[k]                                  = abs_constraint(1.01 * manip_data.efforts[k], opt->manip.tau_max[k]);
          gradient_coeffs_torque(k, j + 2 * joints + 3 * joints * i) = (torque_constraint_plus[k] - torque_constraint[k]) / eps;
        }
      }
    }

    if (opt->constraints.tcp_speed) { // todo: no grad_coeffs yet
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("TCPSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
      *moving_result++ = bound_constraint(1.01 * tcp_speed, 0.0, opt->manip.tcp_max);
      grad_idx++; // todo: add analytical gradients
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("Capsules");
#endif
      if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
        compute_capsules(opt->manip, manip_data);
      }
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("SelfCollisions");
#endif
      auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
      for (u32 j = 0; j < tmp_coll.size; j++) {
        moving_result[j] = -tmp_coll[j] + 0.01;
      }
      moving_result += tmp_coll.size;
      grad_idx += tmp_coll.size; // todo: add analytical gradients
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      ZoneScopedN("ExternalCollisionsCalculate");
#endif
      external_collisions[i] = test_collisions_per_point(manip_data.capsule_list, &(opt->world));

      for (u32 j = 0; j < joints; j++) {
        pos[j] += eps;
        forward_kinematics(opt->manip, manip_data, pos);
        compute_capsules(opt->manip, manip_data);
        external_collisions_plus[i] = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
      }
    }
  }

  //   if (opt->constraints.external_collisions) {
  // #if BLAST_TRACE_LEVEL >= 3
  //     ZoneScopedN("ExternalCollisionsSort");
  // #endif
  //     Array dist_min(opt->constraints.n_collision_constraints, INF_REAL);
  //     Array dist_min_plus(opt->constraints.n_collision_constraints, INF_REAL);
  //     for (u32 i = 0; i < opt->bspline.n_points; i++) {
  //       for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
  //         if (external_collisions[i] < dist_min[j]) {
  //           for (int k = (int) opt->constraints.n_collision_constraints - 1; k > j; k--) {
  //             dist_min[k]      = dist_min[k - 1];
  //             dist_min_plus[k] = dist_min_plus[k - 1];
  //           }
  //           dist_min[j]      = external_collisions[i];
  //           dist_min_plus[j] = external_collisions_plus[i];
  //           break;
  //         }
  //       }
  //     }
  //     for (u32 i = 0; i < opt->constraints.n_collision_constraints; i++) {
  //       moving_result[i] = -dist_min[i] + 0.01; //*std::abs(collisions[j]);
  //     }
  //     moving_result += opt->constraints.n_collision_constraints;

  //     for (u32 i = 0; i < opt->constraints.n_collision_constraints; i++) {
  //       gradient_coeffs_collisions
  //     }
  //     grad_idx += opt->constraints.n_collision_constraints; // todo: add analytical gradients
  //   }
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

inline blast_fn void compute_constraints_per_point(double* result, real& external_collisions, const u32 i, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 3
  ZoneScopedN("ConstraintsPerPoint");
#endif
  double* moving_result = result;

  ManipulatorTempData manip_data;

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
    if (opt->constraints.self_collisions || opt->constraints.external_collisions) {
      compute_capsules(opt->manip, manip_data);
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
      moving_result[0] = (std::abs(1.01 * manip_data.efforts[j]) - opt->manip.tau_max[j]) / opt->manip.tau_max[j];
      moving_result++;
    }
  }

  if (opt->constraints.tcp_speed) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("TCPSpeed");
#endif
    auto J_tool      = get_J_tool(opt, manip_data);
    real tcp_speed   = norm(J_tool * opt->bspline.traj.vel.col(i));
    moving_result[0] = bound_constraint(1.01 * tcp_speed, 0.0, opt->manip.tcp_max);
    moving_result++;
  }

  if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("SelfCollisions");
#endif
    auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
    for (u32 j = 0; j < tmp_coll.size; j++)
      moving_result[j] = -tmp_coll[j] + 0.01; //*std::abs(tmp_coll[j]);
    moving_result += tmp_coll.size;
  }

  if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ExternalCollisionsCalculate");
#endif
    external_collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
  }
}

inline blast_fn u32 ncon_lb_acc(const Optimization* opt, const u32 x_idx) {
  const int n_points            = (int) opt->bspline.lb[x_idx];
  const int n_joints            = (int) opt->manip.n_joints;
  const int n_constraints_basic = n_points * n_joints;
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

  if (opt->constraints.external_collisions)
    n_constraints += opt->constraints.n_collision_constraints;
  if (opt->constraints.self_collisions) {
    n_constraints += opt->manip._n_internal_collisions * n_points;
  }

  for (auto& n_con: opt->constraints.n_extra_constraints)
    n_constraints += n_con;

  return n_constraints;
}

inline u32 ncon_extras(const Optimization* opt) {
  const int n_joints      = (int) opt->manip.n_joints;
  u32       n_constraints = 0;

  if (opt->constraints.torque)
    n_constraints += n_joints;

  if (opt->constraints.tcp_speed)
    n_constraints += 1;

  if (opt->constraints.external_collisions)
    n_constraints += opt->constraints.n_collision_constraints;
  if (opt->constraints.self_collisions) {
    n_constraints += opt->manip._n_internal_collisions;
  }

  for (auto& n_con: opt->constraints.n_extra_constraints)
    n_constraints += n_con;

  return n_constraints;
}

inline blast_fn void nlopt_constraints_dev(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);

  Array gradient_coeffs(m); // todo: check performance
  if (!grad) {
    compute_constraints_dev<false>(result, gradient_coeffs, xv, opt);
  }
  if (grad) {
    compute_constraints_dev<true>(result, gradient_coeffs, xv, opt);
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps      = 1e-5;
    auto           n_con    = ncon_extras(opt);
    u32            n_con_lb = 0;
    Array          x_plus(xlen);
    Array          r_plus(n_con);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx       = 0;
    u32 constraint_idx = 0;
    u32 x_idx          = 0;

    u32 n_torque         = 0;
    u32 n_tcp_speed      = 0;
    u32 n_self_collision = 0;

    memcpy(x_plus.data, x, xlen * sizeof(real)); // todo: is this necessary ? done once, and x_plus is modified to original at the end of each iteration

    for (u32 j = 0; j < xlen - 1; j++) {         // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
      // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
      joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
      x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

      x_plus[j] += eps;

      // Array xv_plus;
      // xv_plus.alias(x_plus.data, xlen);
      opt->bspline.compute_trajectory(x_plus, opt->task);

      n_con_lb = ncon_lb_acc(opt, x_idx); // find the amount of constraints before the current point

      Array external_collisions(opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx]);
      // todo: create alias matrix that points to grad
      // todo: can we change the order in which we store the gradients ?
      grad_idx       = n_con_lb * xlen + j;                                    // gradients are stored column-wise xlen * npoints
      constraint_idx = n_con_lb;
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
          constraint_idx += (opt->manip.n_joints - joint);  // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }

        if (opt->constraints.torque || opt->constraints.tcp_speed || opt->constraints.self_collisions || opt->constraints.external_collisions) {
          compute_constraints_per_point(r_plus.data, external_collisions[i], i, opt);
          n_torque         = 0;
          n_tcp_speed      = 0;
          n_self_collision = 0;
          if (opt->constraints.torque) { // todo: add analytical gradients for torque
            n_torque = opt->manip.n_joints;

            for (u32 k = 0; k < n_torque; k++) {
              grad[grad_idx] = (r_plus.data[k] - result[constraint_idx]) / eps;
              grad_idx += xlen;
              constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
            }
          }

          if (opt->constraints.tcp_speed) {
            n_tcp_speed = 1;

            grad[grad_idx] = (r_plus.data[n_torque] - result[constraint_idx]) / eps;
            grad_idx += xlen;
            constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
          }

          if (opt->constraints.self_collisions) {
            n_self_collision = opt->manip._n_internal_collisions;
            for (u32 k = n_torque + n_tcp_speed; k < n_torque + n_tcp_speed + n_self_collision; k++) {
              grad[grad_idx] = (r_plus.data[k] - result[constraint_idx]) / eps;
              grad_idx += xlen;
              constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
            }
          }
        }
      }
      // todo: add for non point dependant constraints (extras)
      if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
        ZoneScopedN("ExternalCollisionsSort");
#endif
        // todo: fix
        // Array dist_min(opt->constraints.n_collision_constraints, INF_REAL);
        // for (u32 i = 0; i < opt->bspline.n_points; i++) {
        //   for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
        //     if (external_collisions[i] < dist_min[j]) {
        //       for (int k = (int) opt->constraints.n_collision_constraints - 1; k > j; k--) {
        //         dist_min[k] = dist_min[k - 1];
        //       }
        //       dist_min[j] = external_collisions[i];
        //       break;
        //     }
        //   }
        // }
        // for (u32 i = 0; i < opt->constraints.n_collision_constraints; i++) {
        //   grad[grad_idx] = -dist_min[i] + 0.01;               //*std::abs(collisions[j]);
        // }
        // grad_idx += opt->constraints.n_collision_constraints; // todo: add analytical gradients
      }

      x_plus[j] -= eps;
    }

    {
      // last point T
      u32 j = xlen - 1;
      // memcpy(x_plus.data, x, xlen * sizeof(real)); todo: check if necessary
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, xlen);

      // compute constraints
      Array r_plus_T(m);
      compute_constraints(r_plus_T.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * xlen + j] = (r_plus_T[i] - result[i]) / eps;
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

inline blast_fn void nlopt_constraints_dev_new(unsigned m, double* result, unsigned xlen, const double* x, double* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  ZoneScoped;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, xlen);

  Array  gradient_coeffs(m);                                                                           // todo: check performance
  Matrix gradient_torque_coeffs(opt->manip.n_joints, 3 * opt->manip.n_joints * opt->bspline.n_points); // todo: check performance
  if (!grad) {
    compute_constraints_dev_new<false>(result, gradient_coeffs, gradient_torque_coeffs, xv, opt);
  }
  if (grad) {
    compute_constraints_dev_new<true>(result, gradient_coeffs, gradient_torque_coeffs, xv, opt);
    memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

    constexpr real eps      = 1e-5;
    auto           n_con    = ncon_extras(opt);
    u32            n_con_lb = 0;
    Array          x_plus(xlen);
    Array          r_plus(n_con);

    // lb & ub automatically calculated by bsplines
    u32 x_per_joint = (xlen - 1) / opt->manip.n_joints; // = nctrl - 6 (skip first and last 3)
    u32 joint       = 0;

    u32 grad_idx       = 0;
    u32 constraint_idx = 0;
    u32 x_idx          = 0;

    u32 n_torque         = 0;
    u32 n_tcp_speed      = 0;
    u32 n_self_collision = 0;

    memcpy(x_plus.data, x, xlen * sizeof(real)); // todo: is this necessary ? done once, and x_plus is modified to original at the end of each iteration

    for (u32 j = 0; j < xlen - 1; j++) {         // last one is T todo: maybe change to 2 for loops (joint & x_per_joint)
      // vector x is stored as (ctrl points for joint 0, ctrl points for joint 1, ...)
      joint = j / x_per_joint > joint ? joint + 1 : joint; // increase joint by 1 everytime we reach its ctrl points
      x_idx = j - joint * x_per_joint + 3;                 // todo: fix for NaN in PVA of task

      x_plus[j] += eps;

      // Array xv_plus;
      // xv_plus.alias(x_plus.data, xlen);
      opt->bspline.compute_trajectory(x_plus, opt->task);

      n_con_lb = ncon_lb_acc(opt, x_idx); // find the amount of constraints before the current point

      Array external_collisions(opt->bspline.ub[x_idx] - opt->bspline.lb[x_idx]);
      // todo: create alias matrix that points to grad
      // todo: can we change the order in which we store the gradients ?
      grad_idx       = n_con_lb * xlen + j;                                    // gradients are stored column-wise xlen * npoints
      constraint_idx = n_con_lb;
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
          constraint_idx += (opt->manip.n_joints - joint);  // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
        }

        if (opt->constraints.torque || opt->constraints.tcp_speed || opt->constraints.self_collisions || opt->constraints.external_collisions) {
          // compute_constraints_per_point(r_plus.data, external_collisions[i], i, opt);
          n_torque         = 0;
          n_tcp_speed      = 0;
          n_self_collision = 0;
          if (opt->constraints.torque) { // todo: add analytical gradients for torque
            n_torque = opt->manip.n_joints;

            for (u32 k = 0; k < opt->manip.n_joints; k++) {
              grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * gradient_torque_coeffs(k, joint + 3 * opt->manip.n_joints * i) +
                               opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * gradient_torque_coeffs(k, joint + opt->manip.n_joints + 3 * opt->manip.n_joints * i) +
                               opt->bspline.basis_a(x_idx, i) / (opt->bspline.traj.t.back() * opt->bspline.traj.t.back()) * gradient_torque_coeffs(k, joint + 2 * opt->manip.n_joints + 3 * opt->manip.n_joints * i);
              grad_idx += xlen;
              constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
            }
          }

          if (opt->constraints.tcp_speed) { // todo: not working
            n_tcp_speed = 1;

            grad[grad_idx] = (r_plus.data[n_torque] - result[constraint_idx]) / eps;
            grad_idx += xlen;
            constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
          }

          if (opt->constraints.self_collisions) { // todo: not working
            n_self_collision = opt->manip._n_internal_collisions;
            for (u32 k = n_torque + n_tcp_speed; k < n_torque + n_tcp_speed + n_self_collision; k++) {
              grad[grad_idx] = (r_plus.data[k] - result[constraint_idx]) / eps;
              grad_idx += xlen;
              constraint_idx++; // note: eventhough this gradient is not done analytically, we still need to increase the constraint index
            }
          }
        }
      }
      // todo: add for non point dependant constraints (extras)
      if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
        ZoneScopedN("ExternalCollisionsSort");
#endif
        // todo: fix
        // Array dist_min(opt->constraints.n_collision_constraints, INF_REAL);
        // for (u32 i = 0; i < opt->bspline.n_points; i++) {
        //   for (u32 j = 0; j < opt->constraints.n_collision_constraints; j++) {
        //     if (external_collisions[i] < dist_min[j]) {
        //       for (int k = (int) opt->constraints.n_collision_constraints - 1; k > j; k--) {
        //         dist_min[k] = dist_min[k - 1];
        //       }
        //       dist_min[j] = external_collisions[i];
        //       break;
        //     }
        //   }
        // }
        // for (u32 i = 0; i < opt->constraints.n_collision_constraints; i++) {
        //   grad[grad_idx] = -dist_min[i] + 0.01;               //*std::abs(collisions[j]);
        // }
        // grad_idx += opt->constraints.n_collision_constraints; // todo: add analytical gradients
      }

      x_plus[j] -= eps;
    }

    {
      // last point T
      u32 j = xlen - 1;
      // memcpy(x_plus.data, x, xlen * sizeof(real)); todo: check if necessary
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, xlen);

      // compute constraints
      Array r_plus_T(m);
      compute_constraints(r_plus_T.data, x_plus, opt);

      for (u32 i = 0; i < m; i++)
        grad[i * xlen + j] = (r_plus_T[i] - result[i]) / eps;
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

} // namespace blast
