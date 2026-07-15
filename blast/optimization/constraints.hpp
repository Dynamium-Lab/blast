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
  PROFILE_FUNCTION;
#endif
  // todo: remove INF_REAL from constraints at initialization
  // if (value_max == INF_REAL || value_min == -INF_REAL) { // todo: fix for one is INF and not the other
  //   return -1.0;
  // }
  const real center = (value_max + value_min) / 2;
  const real range  = value_max - value_min;

  auto result = 2 * std::abs(value - center) / range - 1.0;
  // auto result2 = std::abs(2 * value) / range - 1.0; // note: fixed, there was an error here todo: this is not true, unless min = -max
  // Assert(result == result2);
  // todo: if this always passes, only keep result2 because it is slightly faster

  return result;
}

inline blast_fn real abs_constraint(const real& value, const real& value_max) {
  return std::abs(value) / value_max - 1.0;
}

inline blast_fn Matrix get_J_tool(const Optimization* opt, const ManipulatorTempData& temp) {
  std::vector<Vec3> r_tool(opt->manip.n_joints);
  r_tool[opt->manip.n_joints - 1] = opt->manip.joint_offsets[opt->manip.n_joints - 1];
  for (int i = (int) opt->manip.n_joints - 2; i >= 0; i--) {
    r_tool[i] = opt->manip.joint_offsets[i] + temp.rotations[i + 1] * r_tool[i + 1];
  }

  for (int i = 0; i < opt->manip.n_joints; i++) {
    r_tool[i] = temp.rotations_mult[i] * r_tool[i];
  }

  Matrix J_tool(3, opt->manip.n_joints);
  for (int i = 0; i < opt->manip.n_joints; i++) {
    // e = opt->manip._Q_mult[i] * opt->manip.joint_axes[i]; // replaced e directly in function to skip copy
    Vec3 cr_tool = cross(temp.rotations_mult[i] * opt->manip.joint_axes[i], r_tool[i]);
    J_tool(0, i) = cr_tool.x;
    J_tool(1, i) = cr_tool.y;
    J_tool(2, i) = cr_tool.z;
  }

  return J_tool;
}

inline blast_fn void constraints_and_gradients_with_segments(const Array& x, Optimization& opt, Array& constraints, Matrix& grad) {
  // constraints (p,v,a,tor) for each joint, for each segment
  // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]

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

  // basis: n_ctrl x n_points
  // Assert(constraints.is_alias); todo: fix for validation
  if (grad.size) {
    Assert(grad.is_alias);
    Assert(grad.rows == x.size);
    Assert(grad.cols == constraints.size);
  }
  const int n_segments                = (int) opt.bspline.n_ctrl - (int) opt.bspline.degree;
  const int n_points_per_segment      = (int) opt.bspline.n_points / n_segments; // todo: check if fine?
  const int n_joints                  = (int) opt.manip.n_joints;
  const int n_ctrl                    = (int) opt.bspline.n_ctrl;
  const int x_len                     = (int) x.size;
  const int n_capsules                = opt.manip._n_caps;
  const int n_constraints_per_segment = opt.constraints.n_constraints_per_segment;
  Assert(constraints.size == n_segments * n_constraints_per_segment);

  const auto& world = opt.world;

  // limits
  auto pmax           = opt.manip.position_max;
  auto pmin           = opt.manip.position_min;
  auto vmax           = opt.manip.velocity_max;
  auto amax           = opt.manip.acceleration_max;
  auto tau_max        = opt.manip.torque_max;
  auto tool_speed_max = opt.manip.tool_speed_max;

  for (int j = 0; j < n_joints; j++) {
    // todo: document the current behaviour in the API
    //        (doesn't currently work if one is inf and the other is not)
    if (pmax[j] == INF_REAL)  // note: replace INF_REAL with huge value
      pmax[j] = 1e300;
    if (pmin[j] == -INF_REAL) // note: replace -INF_REAL with huge negative value
      pmin[j] = -1e300;
  }

  ManipulatorTempData                         manip_data;
  std::array<u8, MAX_JOINTS>                  max_pos_indices{};
  std::array<u8, MAX_JOINTS>                  max_vel_indices{};
  std::array<u8, MAX_JOINTS>                  max_acc_indices{};
  std::array<u8, MAX_JOINTS>                  max_tor_indices{};
  std::array<CollisionEntities, MAX_CAPSULES> max_collision_entities{};
  u8                                          max_internal_collision_index = 0;
  u8                                          max_tool_index               = 0;

  opt.bspline.compute_trajectory(x, opt.task);


  for (int segment = 0; segment < n_segments; segment++) {
#if BLAST_TRACE_LEVEL >= 2
    PROFILE_SCOPE("All Segment Constraints");
#endif
    const int first_affected_control_point = std::max(3, segment);
    const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.degree);
    const int n_affected_control_points    = last_affected_control_point - first_affected_control_point + 1; // note: affected_control_points are inclusive, so when we have last = 5, first = 3, we want 3 (5 - 3) + 1
    const int start_point_for_segment      = segment * n_points_per_segment;
    Assert(n_affected_control_points >= 3);
    Assert(n_affected_control_points <= 6);

    Matrix bp(&opt.bspline.basis_p(0, start_point_for_segment), n_ctrl, n_points_per_segment);
    Matrix bv(&opt.bspline.basis_v(0, start_point_for_segment), n_ctrl, n_points_per_segment);
    Matrix ba(&opt.bspline.basis_a(0, start_point_for_segment), n_ctrl, n_points_per_segment);

    Array max_pos_constraints(n_joints, -INF_REAL);
    Array max_vel_constraints(n_joints, -INF_REAL);
    Array max_acc_constraints(n_joints, -INF_REAL);
    Array max_tor_constraints(n_joints, -INF_REAL);
    real  max_tool_speed_constraints   = -INF_REAL;
    real  max_internal_col_constraints = -INF_REAL; // todo: worst or worst per capsule ?
    Array max_col_constraints(n_capsules, -INF_REAL);

    for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) {
      {
#if BLAST_TRACE_LEVEL >= 2
        PROFILE_SCOPE("All Point Constraints");
#endif
        auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment);
        auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment);
        auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment);

        forward_kinematics(opt.manip, manip_data, p);
        compute_capsules(opt.manip, manip_data);
        dynamics(opt.manip, manip_data, v, a);

        for (int j = 0; j < n_joints; j++) {
          // position
          if (opt.constraints.position) {

#if BLAST_TRACE_LEVEL >= 3
            PROFILE_SCOPE("Pos Constraints");
#endif

            if (const auto c = bound_constraint(p[j], pmin[j], pmax[j]);
                c > max_pos_constraints[j]) {
              max_pos_constraints[j] = c;
              max_pos_indices[j]     = point_in_segment;
            }
          }
          // velocity
          if (opt.constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
            PROFILE_SCOPE("Vel Constraints");
#endif
            if (const auto c = std::abs(v[j]) / vmax[j] - 1.0;
                c > max_vel_constraints[j]) {
              max_vel_constraints[j] = c;
              max_vel_indices[j]     = point_in_segment;
            }
          }
          // acceleration
          if (opt.constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
            PROFILE_SCOPE("Acc Constraints");
#endif
            if (const auto c = std::abs(a[j]) / amax[j] - 1.0;
                c > max_acc_constraints[j]) {
              max_acc_constraints[j] = c;
              max_acc_indices[j]     = point_in_segment;
            }
          }
          // torque
          if (opt.constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
            PROFILE_SCOPE("Tau Constraints");
#endif
            if (const auto c = std::abs(manip_data.efforts[j]) / tau_max[j] - 1.0;
                c > max_tor_constraints[j]) {
              max_tor_constraints[j] = c;
              max_tor_indices[j]     = point_in_segment;
            }
          }
        }

        // Tool speed
        if (opt.constraints.tool_speed) {
#if BLAST_TRACE_LEVEL >= 3
          PROFILE_SCOPE("Tool Constraints");
#endif
          const auto J_tool     = get_J_tool(&opt, manip_data); // todo: clean up get_J_tool to a get_tool_speed
          const auto tool_speed = norm(get_J_tool(&opt, manip_data) * v);
          if (const auto c = bound_constraint(tool_speed, 0.0, tool_speed_max);
              c > max_tool_speed_constraints) {
            max_tool_speed_constraints = c;
            max_tool_index             = point_in_segment;
          }
        }

        // self collision
        if (opt.constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
          PROFILE_SCOPE("Self Constraints");
#endif
          // check every internal collision
          if (const auto c = max(-get_internal_collisions(opt.manip, manip_data));
              c > max_internal_col_constraints) {
            max_internal_col_constraints = c;
            max_internal_collision_index = point_in_segment;
          }
        }

        // external collision
        if (opt.constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
          PROFILE_SCOPE("Ext Constraints");
#endif

          // check every capsule with world
          for (int capsule_id = 0; capsule_id < n_capsules; capsule_id++) {
            real       dist_min = INF_REAL;
            const auto capsule  = manip_data.capsule_list[capsule_id];

            CollisionEntities collision_objects{};

            // check against boxes
            int count = 0;
            for (const auto& box: world.boxes) {
              if (const auto dist = distance(capsule, box);
                  dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::box;
                collision_objects.box               = box;
                collision_objects.point_in_segment  = point_in_segment;
              }
              count++;
            }

            // check against capsules
            count = 0;
            for (const auto caps: world.capsules) {
              if (const auto dist = distance(capsule, caps);
                  dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::capsule;
                collision_objects.capsule           = capsule;
                collision_objects.point_in_segment  = point_in_segment;
              }
              count++;
            }

            // check against spheres
            count = 0;
            for (const auto sphere: world.spheres) {
              if (const auto dist = distance(capsule, sphere);
                  dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::sphere;
                collision_objects.sphere            = sphere;
                collision_objects.point_in_segment  = point_in_segment;
              }
              count++;
            }

            // --- Dynamic tests --- todo: check & fix gradients
            int  current_point = segment * n_points_per_segment + point_in_segment;
            int  max_point     = n_segments * n_points_per_segment - 1;                                            // todo: check -1 ?
            real current_time  = x.back() * ((real) current_point / (real) max_point) + opt.trajectory_start_time; // trajectory time * progression along trajectory

            count = 0;
            for (const auto& box: world.dynamic_boxes) {
              auto       current_box = box.lookup(current_time);
              const auto dist        = distance(capsule, current_box);
              if (dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::box; // todo: Add dynamic boxes? Not sure if necessary if we take the current box instead???
                collision_objects.box               = current_box;
                collision_objects.point_in_segment  = point_in_segment;
              }
              count++;
            }

            count = 0;
            for (const auto& caps: world.dynamic_capsules) {
              auto       current_caps = caps.lookup(current_time);
              const auto dist         = distance(capsule, current_caps);
              if (dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::capsule; // todo: Add dynamic capsules? Not sure if necessary if we take the current caps instead???
                collision_objects.capsule           = current_caps;
                collision_objects.point_in_segment  = point_in_segment;
              }
              count++;
            }

            count = 0;
            for (const auto& sphere: world.dynamic_spheres) {
              auto       current_sphere = sphere.lookup(current_time);
              const auto dist           = distance(capsule, current_sphere);
              if (dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::sphere; // todo: Add dynamic spheres? Not sure if necessary if we take the current sphere instead???
                collision_objects.sphere            = current_sphere;
                collision_objects.point_in_segment  = point_in_segment;
              }
              count++;
            }

            count = 0;
            for (const auto& door: world.dynamic_doors) {
              auto       current_door = door.lookup(current_time);
              const auto dist         = distance(capsule, current_door);
              if (dist < dist_min) {
                dist_min                            = dist;
                collision_objects.other_object_type = CollisionObjectType::box; // todo: Add dynamic doors? Not sure if necessary if we take the current door instead???
                collision_objects.box               = current_door;
                collision_objects.point_in_segment  = point_in_segment;
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
        }
      }
    }

    // at this point we have max constraints for pos, vel, acc, tor, and collisions for this segment

    // fill in the constraints for the current segment
    // [p1, p2,..., v1, v2,..., a1, a2,..., t1, t2,...]
    auto fill_idx = segment * n_constraints_per_segment;
    if (opt.constraints.position)
      std::copy_n(max_pos_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints; // note (andre): we can use the comma operator because we don't need the output of copy_n()
    if (opt.constraints.velocity)
      std::copy_n(max_vel_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    if (opt.constraints.acceleration)
      std::copy_n(max_acc_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    if (opt.constraints.torque)
      std::copy_n(max_tor_constraints.data, n_joints, &constraints[fill_idx]), fill_idx += n_joints;
    if (opt.constraints.tool_speed)
      constraints[fill_idx++] = max_tool_speed_constraints;
    if (opt.constraints.self_collisions)
      constraints[fill_idx++] = max_internal_col_constraints;
    if (opt.constraints.external_collisions)
      std::copy_n(max_col_constraints.data, n_capsules, &constraints[fill_idx]), fill_idx += n_capsules;


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
#if BLAST_TRACE_LEVEL >= 2
      PROFILE_SCOPE("Grad");
#endif
      // Matrix (alias) in which we can insert the gradient for the current segment
      Matrix grad_segment(&grad(0, segment * n_constraints_per_segment), x_len, n_constraints_per_segment);
      Assert(grad_segment.is_alias);

      int con     = 0;
      int con_idx = 0;

      // positions
      if (opt.constraints.position) {
        for (int joint = 0; joint < n_joints; joint++) {

          // real pos = opt.bspline.traj.pos(joint, max_pos_indices[joint]);

          // Array of the column where to put the gradient for the current constraint
          Array fill_column = grad_segment.col(con);
          Assert(con == joint);
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          // Which values in 'x' affect the current joint's position
          auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
          // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
          x_idx += first_affected_control_point - 3;

          // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
          Array bp_to_use(&bp(first_affected_control_point, max_pos_indices[joint]), n_affected_control_points);
          Assert(bp_to_use.is_alias);

          real coeff = 2.0 * sign(opt.bspline.traj.pos(joint, start_point_for_segment + max_pos_indices[joint]) - (pmax[joint] + pmin[joint]) / 2) / (pmax[joint] - pmin[joint]); // note:

          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx++] = bp_to_use[i] * coeff;
          }

          // note: dp/dT == 0
          con++;
        }
        con_idx += n_joints;
      }

      // velocities
      auto one_over_T = 1 / opt.bspline.traj.t.back();
      if (opt.constraints.velocity) {
        for (int joint = 0; joint < n_joints; joint++) {
          // Array of the column where to put the gradient for the current constraint
          Array fill_column = grad_segment.col(con);
          Assert(con == con_idx + joint);
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          // Which values in 'x' affect the current joint's position
          auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
          // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
          x_idx += first_affected_control_point - 3;

          Array bv_to_use(&bv(first_affected_control_point, max_vel_indices[joint]), n_affected_control_points);
          Assert(bv_to_use.is_alias);

          real coeff = sign(opt.bspline.traj.vel(joint, start_point_for_segment + max_vel_indices[joint])) / vmax[joint] * one_over_T;

          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx++] = bv_to_use[i] * coeff;
          }

          // dvj/dT = - (Cv + 1) / T
          fill_column.back() = -(max_vel_constraints[joint] + 1) * one_over_T;

          con++;
        }
        con_idx += n_joints;
      }

      // accelerations
      auto one_over_T2 = one_over_T * one_over_T;
      if (opt.constraints.acceleration) {
        for (int joint = 0; joint < n_joints; joint++) {
          // Array of the column where to put the gradient for the current constraint
          Array fill_column = grad_segment.col(con);
          Assert(con == con_idx + joint);
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          // Which values in 'x' affect the current joint's position
          auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
          // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
          x_idx += first_affected_control_point - 3;

          Array ba_to_use(&ba(first_affected_control_point, max_acc_indices[joint]), n_affected_control_points);
          Assert(ba_to_use.is_alias);

          real coeff = sign(opt.bspline.traj.acc(joint, start_point_for_segment + max_acc_indices[joint])) / amax[joint] * one_over_T2;

          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx++] = ba_to_use[i] * coeff;
          }

          // daj/dT = -2 * (Ca + 1) / T
          fill_column.back() = -2 * (constraints[con] + 1) * one_over_T;

          con++;
        }
        con_idx += n_joints;
      }

      // torque
      // [dt0/dp0, dt1/dp0, ..., dt4/dp0, dt5/dp0]
      // [dt0/dp1, dt1/dp1, ..., dt4/dp1, dt5/dp1]
      // [dt0/dp2, dt1/dp2, ..., dt4/dp2, dt5/dp2]
      // [dt0/dp3, dt1/dp3, ..., dt4/dp3, dt5/dp3]
      // [.
      // [.
      // [.
      if (opt.constraints.torque) {
        for (int joint = 0; joint < n_joints; joint++) {
          Array fill_column = grad_segment.col(con);
          Assert(con == con_idx + joint);
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          constexpr real eps            = BLAST_FD_STEP;
          const auto     point          = max_tor_indices[joint];
          const auto     p              = opt.bspline.traj.pos.col(start_point_for_segment + point);
          const auto     v              = opt.bspline.traj.vel.col(start_point_for_segment + point);
          const auto     a              = opt.bspline.traj.acc.col(start_point_for_segment + point);
          auto           old_constraint = max_tor_constraints[joint];
          auto           tau_max_now    = tau_max[joint];

          auto p_plus = p;
          auto v_plus = v;
          auto a_plus = a;

          Array bp_to_use(&bp(first_affected_control_point, point), n_affected_control_points);
          Array bv_to_use(&bv(first_affected_control_point, point), n_affected_control_points);
          Array ba_to_use(&ba(first_affected_control_point, point), n_affected_control_points);

          auto x_idx      = first_affected_control_point - 3;
          auto x_idx_skip = n_ctrl - 6;

          for (int j = 0; j < n_joints; j++) {

            // partial derivative of torque constraints w.r.t. position
            // finite difference on position
            p_plus[j] += eps;
            // compute the derivative of constraint(joint) w.r.t. theta(j). (remember, joint != j)
            forward_kinematics(opt.manip, manip_data, p_plus);
            dynamics(opt.manip, manip_data, v_plus, a_plus);
            const real new_constraint_p = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // todo: remove 1.01 and use for validation only
            const real dtau_dp          = (new_constraint_p - old_constraint) / eps;
            // reset finite difference
            p_plus[j] = p[j];

            // note: reset forward kinematics because 'v' and 'a' don't change it
            forward_kinematics(opt.manip, manip_data, p); // todo: precompute once

            // partial derivative of torque constraints w.r.t. velocity
            // finite difference on velocity
            v_plus[j] += eps;
            dynamics(opt.manip, manip_data, v_plus, a_plus);
            const real new_constraint_v = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // todo: remove 1.01 and use for validation only

            const real dtau_dv = (new_constraint_v - old_constraint) / eps;
            // reset finite difference
            v_plus[j] = v[j];

            // partial derivative of torque constraints w.r.t. acceleration
            // finite difference on acceleration
            a_plus[j] += eps;
            dynamics(opt.manip, manip_data, v_plus, a_plus);
            const real new_constraint_a = std::abs(manip_data.efforts[joint]) / tau_max_now - 1; // todo: remove 1.01 and use for validation only
            const real dtau_da          = (new_constraint_a - old_constraint) / eps;
            // reset finite difference
            a_plus[j] = a[j];

            // insert into the gradient
            for (int i = 0; i < n_affected_control_points; i++) {
              fill_column[x_idx + i] = bp_to_use[i] * dtau_dp +
                                       bv_to_use[i] * dtau_dv * one_over_T +
                                       ba_to_use[i] * dtau_da * one_over_T2;
            }
            x_idx += x_idx_skip;

            // gradient w.r.t. T
            fill_column.back() += dtau_dv * (-v[j] * one_over_T) + dtau_da * (-2 * a[j] * one_over_T);
          }
          con++; // note: moved out of the inner j loop (only changes at the end of all j torques per joint)
        }
        con_idx += n_joints;
      }

      // tool speed
      if (opt.constraints.tool_speed) {
        constexpr real eps    = BLAST_FD_STEP;
        const auto     point  = max_tool_index;
        const auto     p      = opt.bspline.traj.pos.col(start_point_for_segment + point);
        const auto     v      = opt.bspline.traj.vel.col(start_point_for_segment + point);
        auto           p_plus = p;
        auto           v_plus = v;

        auto x_idx      = first_affected_control_point - 3;
        auto x_idx_skip = n_ctrl - 6;

        Array fill_column = grad_segment.col(con);
        Assert(con == con_idx);
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
        Array bp_to_use(&bp(first_affected_control_point, point), n_affected_control_points);
        Array bv_to_use(&bv(first_affected_control_point, point), n_affected_control_points);
        Assert(bp_to_use.is_alias);
        Assert(bv_to_use.is_alias);
        for (int j = 0; j < n_joints; j++) {
          // todo: check if the joint can actually move the current capsule
          p_plus[j] += eps;

          // recompute tool_speed
          forward_kinematics(opt.manip, manip_data, p_plus);
          const auto J_tool_p         = get_J_tool(&opt, manip_data);
          const auto tool_speed_p     = norm(J_tool_p * v_plus);
          const auto new_constraint_p = bound_constraint(tool_speed_p, 0.0, tool_speed_max);
          // partial difference d(tool_speed)/dp
          const real dtool_speed_dp = (new_constraint_p - max_tool_speed_constraints) / eps;
          p_plus[j]                 = p[j]; // reset finite difference

          forward_kinematics(opt.manip, manip_data, p_plus);
          v_plus[j] += eps;
          const auto J_tool_v         = get_J_tool(&opt, manip_data);
          const auto tool_speed_v     = norm(J_tool_v * v_plus);
          const auto new_constraint_v = bound_constraint(tool_speed_v, 0.0, tool_speed_max);
          const real dtool_speed_dv   = (new_constraint_v - max_tool_speed_constraints) / eps;
          v_plus[j]                   = v[j];

          // fill the gradient for the control points that affect current joint 'j' d(collision)/dp * dp/d(control point)
          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx + i] = bp_to_use[i] * dtool_speed_dp +
                                     bv_to_use[i] * dtool_speed_dv * one_over_T;
          }
          x_idx += x_idx_skip;
          // gradient w.r.t. T
          fill_column.back() += dtool_speed_dv * (-v[j] * one_over_T); // todo: check this !!
        }
        con++;
        con_idx++;
      }

      // internal collisions
      if (opt.constraints.self_collisions) {
        constexpr real eps    = BLAST_FD_STEP;
        const auto     point  = max_internal_collision_index;
        const auto     p      = opt.bspline.traj.pos.col(start_point_for_segment + point);
        auto           p_plus = p;

        auto x_idx      = first_affected_control_point - 3;
        auto x_idx_skip = n_ctrl - 6;

        Array fill_column = grad_segment.col(con);
        Assert(con == con_idx);
        Assert(fill_column.size == x_len);
        Assert(fill_column.is_alias);

        // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
        Array bp_to_use(&bp(first_affected_control_point, point), n_affected_control_points);
        Assert(bp_to_use.is_alias);
        for (int j = 0; j < n_joints; j++) {
          // todo: check if the joint can actually move the current capsule
          p_plus[j] += eps;

          // recompute internal collisions at the worst point in segment
          forward_kinematics(opt.manip, manip_data, p_plus);
          compute_capsules(opt.manip, manip_data);
          const auto new_internal_collision_constraint = max(-get_internal_collisions(opt.manip, manip_data));
          // partial difference d(internal_collision)/dp
          const real dint_coll_dp = (new_internal_collision_constraint - max_internal_col_constraints) / eps;

          // fill the gradient for the control points that affect current joint 'j' d(collision)/dp * dp/d(control point)
          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx + i] = bp_to_use[i] * dint_coll_dp;
          }
          x_idx += x_idx_skip;

          p_plus[j] = p[j]; // reset finite difference
        }
        con++;
        con_idx++;
      }

      // collisions
      if (opt.constraints.external_collisions) {
        for (int capsule_id = 0; capsule_id < opt.manip._n_caps; capsule_id++) {
          constexpr real eps    = BLAST_FD_STEP;
          const auto     point  = max_collision_entities[capsule_id].point_in_segment;
          const auto     p      = opt.bspline.traj.pos.col(start_point_for_segment + point);
          auto           p_plus = p;

          auto x_idx      = first_affected_control_point - 3;
          auto x_idx_skip = n_ctrl - 6;

          Array fill_column = grad_segment.col(con);
          Assert(con == con_idx + capsule_id);
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
          Array bp_to_use(&bp(first_affected_control_point, point), n_affected_control_points);
          Assert(bp_to_use.is_alias);

          // finite difference w.r.t. joint positions then multiply by relevant basis functions.
          for (int j = 0; j < n_joints; j++) {
            // todo: check if the joint can actually move the current capsule
            p_plus[j] += eps;

            // recompute collision constraint, but only with the current capsule and the identified object.
            forward_kinematics(opt.manip, manip_data, p_plus);
            compute_capsules(opt.manip, manip_data);
            const auto capsule = manip_data.capsule_list[capsule_id];

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

            distance_plus = -distance_plus; // negative distance is positive constraint

            // partial difference d(collision)/dp
            const real dcoll_dp = (distance_plus - max_col_constraints[capsule_id]) / eps;

            // fill the gradient for the control points that affect current joint 'j' d(collision)/dp * dp/d(control point)
            for (int i = 0; i < n_affected_control_points; i++) {
              fill_column[x_idx + i] = bp_to_use[i] * dcoll_dp;
            }
            x_idx += x_idx_skip;

            p_plus[j] = p[j]; // reset finite difference
          }

          con++; // finished filling the column with the gradient of the collision of the current capsule w.r.t. each joint position
        }
      }
    }
  }
}

inline blast_fn void nlopt_constraints_with_segments(unsigned m, real* result, unsigned x_len, const real* x, real* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  PROFILE_FUNCTION;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, x_len);

  Array constraints;
  constraints.alias(result, m);

  Matrix gradients;
  if (grad) {
    memset(grad, 0, m * x_len * sizeof(real));
    gradients.alias(grad, x_len, m);
  }

  constraints_and_gradients_with_segments(xv, *opt, constraints, gradients);

  if (opt->constraints.collect_x_each_iteration) {
    opt->constraints.x_list.push_back(xv);
  }
}

inline blast_fn void compute_constraints(real* result, const Array& x, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 2
  PROFILE_FUNCTION;
#endif

  real* moving_result = result;

  const int n_capsules = opt->manip._n_caps;

  // todo: compute inside loop so that every worker can compute independently??
  {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // todo: no copy every iteration
  // ObjMatrix<Capsule> capsules(opt->manip._n_caps, (int) opt->bspline.n_points / (int) opt->constraints.n_collision_skip);

  // todo: compute n_con_per_point() so that every worker knows where to put their result

  for (int i = 0; i < opt->bspline.n_points; i++) {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("ConstraintSinglePoint");
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
      PROFILE_SCOPE("Kinematics");
#endif
      forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
    }

    {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("Capsules");
#endif
      // todo: this cleaner
      // if (opt->constraints.external_collisions || opt->constraints.self_collisions) {
      compute_capsules(opt->manip, manip_data);
      // if (opt->constraints.external_collisions) {
      //   if (i % opt->constraints.n_collision_skip == 0 && i / opt->constraints.n_collision_skip < capsules.cols) {
      //     for (int j = 0; j < opt->manip._n_caps; j++) {
      //       capsules(j, i / (int) opt->constraints.n_collision_skip) = manip_data.capsule_list[j];
      //     }
      //   }
      // }
      // }
    }

    if (opt->constraints.position) {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("Position");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        *moving_result++ = bound_constraint(opt->bspline.traj.pos(j, i), opt->manip.position_min[j], opt->manip.position_max[j]);
      }
    }

    if (opt->constraints.velocity) {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("Velocity");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        // moving_result[0] = abs_constraint(opt->bspline.traj.vel(j, i), opt->manip.velocity_max[j]);
        *moving_result++ = std::abs(opt->bspline.traj.vel(j, i)) / opt->manip.velocity_max[j] - 1.0;
      }
    }

    if (opt->constraints.acceleration) {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("Acceleration");
#endif
      for (int j = 0; j < opt->manip.n_joints; j++) {
        *moving_result++ = std::abs(opt->bspline.traj.acc(j, i)) / opt->manip.acceleration_max[j] - 1.0;
      }
    }

    if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("Dynamics");
#endif
      auto vel = opt->bspline.traj.vel.col(i);
      auto acc = opt->bspline.traj.acc.col(i);
      dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

      for (int j = 0; j < opt->manip.n_joints; j++) {
        *moving_result++ = std::abs(manip_data.efforts[j]) / opt->manip.torque_max[j] - 1.0;
      }
    }

    if (opt->constraints.tool_speed) {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("ToolSpeed");
#endif
      auto J_tool      = get_J_tool(opt, manip_data);
      real tool_speed  = norm(J_tool * opt->bspline.traj.vel.col(i));
      *moving_result++ = bound_constraint(tool_speed, 0.0, opt->manip.tool_speed_max);
    }

    if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("SelfCollisions");
#endif
      auto tmp_coll = max(-get_internal_collisions(opt->manip, manip_data));
      // for (u32 j = 0; j < tmp_coll.size; j++)
      *moving_result++ = tmp_coll; //*std::abs(tmp_coll[j]);
    }

    if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
      PROFILE_SCOPE("ExternalCollisions");
#endif
      Array max_col_constraints(n_capsules, -INF_REAL);
      for (int capsule_id = 0; capsule_id < n_capsules; capsule_id++) {
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
      for (u32 j = 0; j < n_capsules; j++)
        *moving_result++ = max_col_constraints[j];
    }
  }
  //   {
  // #if BLAST_TRACE_LEVEL >= 3
  //     PROFILE_SCOPE("CustomConstraints");
  // #endif
  //     for (u32 i = 0; i < opt->constraints.extra_constraints.size(); i++) { // todo: split in extra_constraint per point or once
  //       auto extra_constraint = opt->constraints.extra_constraints[i];
  //       extra_constraint(moving_result, opt);
  //       moving_result += opt->constraints.n_extra_constraints[i];
  //     }
  //   }
}

inline blast_fn void nlopt_constraints(unsigned m, real* result, unsigned x_len, const real* x, real* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  PROFILE_FUNCTION;
#endif
  auto* opt = (Optimization*) f_data;

  Array xv;
  xv.alias(x, x_len);
  {

#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("Constraints");
#endif
    compute_constraints(result, xv, opt);
  }

  // gradients calculation
  if (grad) {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("Grad");
#endif
    // todo: no construction??
    Array x_plus(x_len);
    Array r_plus(m);

    for (u32 j = 0; j < x_len; j++) {
      constexpr real eps = BLAST_FD_STEP;
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
    if (opt->constraints.collect_x_each_iteration) {
      opt->constraints.x_list.push_back(xv);
    }
  }
}

inline blast_fn bool validate_task(Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 1
  PROFILE_FUNCTION;
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
        result = bound_constraint(pos(j, i), manip.position_min[j], manip.position_max[j]);
        if (result > 0) {
          std::cout << "Position outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (constraints.velocity) {
      for (int j = 0; j < manip.n_joints; j++) {
        result = abs_constraint(vel(j, i), manip.velocity_max[j]);
        if (result > 0) {
          std::cout << "Velocity outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (constraints.acceleration) {
      for (int j = 0; j < manip.n_joints; j++) {
        result = abs_constraint(acc(j, i), manip.acceleration_max[j]);
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
        result = abs_constraint(manip_data.efforts[j], manip.torque_max[j]);
        if (result > 0) {
          std::cout << "Torque outside bounds." << std::endl;
          return false;
        }
      }
    }

    if (constraints.tool_speed) {
      auto J_tool     = get_J_tool(opt, manip_data);
      real tool_speed = norm(J_tool * vel.col(i));
      result          = bound_constraint(tool_speed, -INF_REAL, manip.tool_speed_max);
      if (result > 0) {
        std::cout << "Tool speed outside bounds." << std::endl;
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

template<bool is_grad>
inline blast_fn std::tuple<real, real> abs_constraint_analytical(const real& q, const real& q_max) {
  real constraint = (std::abs(q)) / q_max - 1.0;
  real gradient   = 0.0;
  if (is_grad)
    gradient = sign(q) * 1 / q_max;
  return std::make_tuple(constraint, gradient);
}

template<bool is_grad>
blast_fn std::tuple<real, real> bound_constraint_analytical(const real& q, const real& q_min, const real& q_max) {
  // todo: remove INF_REAL from constraints at initialization
  if (q_max == INF_REAL || q_min == -INF_REAL)
    return std::make_tuple(-1.0, 0.0);
  else {
    const real center = (q_max + q_min) / 2;
    const real range  = q_max - q_min;

    auto result = 2 * std::abs(q - center) / range - 1.0;

    // Conditionally compute gradient coefficient
    real gradient_coeff = 0.0;
    if constexpr (is_grad) {
      gradient_coeff = 2 * sign(q - center) / range;
    }
    return std::make_tuple(result, gradient_coeff);
  }
}

inline void compute_constraints_with_analytical_pva(ConstraintPerPoint& constraints, const Array& x, const u32 x_idx, const u32 joint_idx, Optimization* opt) {

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
      auto J_tool                                                       = get_J_tool(opt, manip_data);
      real tool_speed                                                   = norm(J_tool * opt->bspline.traj.vel.col(i));
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

template<bool is_grad> // note: n_collision_skip must be 1 for this to work !!!
blast_fn void compute_constraints_with_analytical_dynamics(real* result, Array& gradient_coeffs, Matrix& dtau_dp, Matrix& dtau_dv, Matrix& dtau_da, Matrix& dtool_dp, Matrix& dtool_dv, Matrix& dselfcol_dp, std::vector<Matrix>& dcol_dp, const Array& x, Optimization* opt) {
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

  real* moving_result = result;
  u32   grad_idx      = 0;

  constexpr real eps = BLAST_FD_STEP;

  auto joints     = opt->manip.n_joints;
  auto n_capsules = opt->manip._n_caps;

  std::array<CollisionEntities, MAX_CAPSULES> max_collision_entities{};

  opt->constraints.n_collision_skip = 1;
  {
#if BLAST_TRACE_LEVEL >= 3
    ZoneScopedN("ComputeTrajectory");
#endif
    opt->bspline.compute_trajectory(x, opt->task);
  }

  // Lambda to process common bound constraint operations
  auto process_bound = [&](real value, real bound_max) {
    auto [constraint, gradient_coeff] = abs_constraint_analytical<is_grad>(value, bound_max);
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
            auto [constraint, gradient_coeff] = bound_constraint_analytical<is_grad>(opt->bspline.traj.pos(j, i), opt->manip.position_min[j], opt->manip.position_max[j]);
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

inline blast_fn void compute_constraints_per_point(real* result, real& external_collisions, const u32 i, Optimization* opt) {
#if BLAST_TRACE_LEVEL >= 3
  PROFILE_SCOPE("ConstraintsPerPoint");
#endif
  real* moving_result = result;

  ManipulatorTempData manip_data;

  auto pos = opt->bspline.traj.pos.col(i);
  Assert(pos.is_alias);

  {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("Kinematics");
#endif
    forward_kinematics(opt->manip, manip_data, pos); // fills _Q, _Q_mult, and _p_j
  }

  {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("Capsules");
#endif
    if (opt->constraints.self_collisions || opt->constraints.external_collisions) {
      compute_capsules(opt->manip, manip_data);
    }
  }

  if (opt->constraints.torque) {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("Dynamics");
#endif
    auto vel = opt->bspline.traj.vel.col(i);
    auto acc = opt->bspline.traj.acc.col(i);
    dynamics(opt->manip, manip_data, vel, acc); // fills _efforts

    for (int j = 0; j < opt->manip.n_joints; j++) {
      moving_result[0] = (std::abs(manip_data.efforts[j]) - opt->manip.torque_max[j]) / opt->manip.torque_max[j];
      moving_result++;
    }
  }

  if (opt->constraints.tool_speed) {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("ToolSpeed");
#endif
    auto J_tool      = get_J_tool(opt, manip_data);
    real tool_speed  = norm(J_tool * opt->bspline.traj.vel.col(i));
    moving_result[0] = bound_constraint(tool_speed, 0.0, opt->manip.tool_speed_max);
    moving_result++;
  }

  if (opt->constraints.self_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("SelfCollisions");
#endif
    auto tmp_coll = get_internal_collisions(opt->manip, manip_data);
    for (u32 j = 0; j < tmp_coll.size; j++)
      moving_result[j] = -tmp_coll[j] + 0.01; //*std::abs(tmp_coll[j]);
    moving_result += tmp_coll.size;
  }

  if (opt->constraints.external_collisions) {
#if BLAST_TRACE_LEVEL >= 3
    PROFILE_SCOPE("ExternalCollisionsCalculate");
#endif
    external_collisions = test_collisions_per_point(manip_data.capsule_list, &(opt->world));
  }
}

inline blast_fn u32 ncon_lb_acc(const Optimization* opt, const u32 x_idx) {
  const int n_points            = (int) opt->bspline.lower_bounds[x_idx];
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

inline u32 ncon_extras(const Optimization* opt) {
  const int n_joints      = (int) opt->manip.n_joints;
  u32       n_constraints = 0;

  if (opt->constraints.torque)
    n_constraints += n_joints;

  if (opt->constraints.tool_speed)
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

inline u32 ncon_per_point_with_analytical_pva(const Optimization* opt, const int x_idx) {
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

inline blast_fn void nlopt_constraints_with_analytical_pva(unsigned m, real* result, unsigned xlen, const real* x, real* grad, void* f_data) {
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

    constexpr real eps = BLAST_FD_STEP;

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

      auto n_con_per_point = ncon_per_point_with_analytical_pva(opt, x_idx) / (opt->bspline.upper_bounds[x_idx] + 1 - opt->bspline.lower_bounds[x_idx]);

      auto               n_points = (opt->bspline.upper_bounds[x_idx] + 1 - opt->bspline.lower_bounds[x_idx]);
      ConstraintPerPoint constraint(n_joints, n_points, n_capsules);
      compute_constraints_with_analytical_pva(constraint, x_plus, x_idx, joint, opt);

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

inline blast_fn void nlopt_constraints_with_analytical_dynamics(unsigned m, real* result, unsigned xlen, const real* x, real* grad, void* f_data) {
#if BLAST_TRACE_LEVEL >= 1
  PROFILE_FUNCTION;
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
    compute_constraints_with_analytical_dynamics<false>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtool_dp, dtool_dv, dselfcol_dp, dcol_dp, xv, opt);
  }
  if (grad) {
    compute_constraints_with_analytical_dynamics<true>(result, gradient_coeffs, dtau_dp, dtau_dv, dtau_da, dtool_dp, dtool_dv, dselfcol_dp, dcol_dp, xv, opt);
    {
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("Grad Fill");
#endif
      memset(grad, 0, m * xlen * sizeof(real)); // note: zeros grad, since grad originally starts with -6e+66 ...

      constexpr real eps      = BLAST_FD_STEP;
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
            }
          }

          if (opt->constraints.tool_speed) {
            grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dtool_dp(joint, i) +
                             opt->bspline.basis_v(x_idx, i) / opt->bspline.traj.t.back() * dtool_dv(joint, i);
            grad_idx += xlen;
          }

          if (opt->constraints.self_collisions) {
            grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dselfcol_dp(joint, i);
            grad_idx += xlen;
          }

          if (opt->constraints.external_collisions) {
            for (int k = 0; k < opt->manip._n_caps; k++) {
              grad[grad_idx] = opt->bspline.basis_p(x_idx, i) * dcol_dp[i](k, joint);
              grad_idx += xlen;
            }
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

      if (opt->constraints.collect_x_each_iteration) { // when more info is needed per iteration
        opt->constraints.x_list.push_back(xv);
      }
    }
  }
}
} // namespace blast
