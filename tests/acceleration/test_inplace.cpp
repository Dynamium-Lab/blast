#include <blast>
#include "test_helper/test_helper.hpp"

#define BLAST_TRACE_LEVEL 3

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"

namespace blast {


TEST_CASE("Replace INF_REAL with huge value", "[acceleration]") {
  int n_iterations = 10000;

  // both pmax and pmin are +/-INF_REAL
  real pmax = INF_REAL;
  real pmin = -INF_REAL;

  pmax = 1e300;
  pmin = -1e300;

  real c          = (pmax + pmin) / 2;
  real b          = (pmax - pmin) / 2;
  real one_over_b = 1 / b;

  Array p(n_iterations);
  p = fill_random(p, 360);

  for (int i = 0; i < n_iterations; i++) {
    real q_prime    = p[i] - c;
    real constraint = (std::abs(q_prime) - b) / b;
    real gradient   = sign(q_prime) * 1 / b;

    CHECK(constraint == -1.0);
    CHECK(std::abs(gradient) <= 1e-300);
  }
};

TEST_CASE("Analytic dC/dX per segments tests", "[acceleration]") {
  int n_iterations = 1000;

  for (int iter = 0; iter < n_iterations; iter++) {

    Optimization opt(get_generic_Link6(), get_link6_task());

    opt.world = get_lab_world();

    Bspline bspline(16, 110, 5, 6);
    opt.bspline = bspline;

    auto x_len    = opt.bspline.x_len(opt.task);
    auto n_joints = opt.manip.n_joints;

    opt.guess.type = Guess::random;


    opt.constraints.position            = true;
    opt.constraints.velocity            = true;
    opt.constraints.acceleration        = true;
    opt.constraints.torque              = true;
    opt.constraints.tcp_speed           = false;
    opt.constraints.self_collisions     = false;
    opt.constraints.external_collisions = false;

    opt.max_tries         = 1;
    opt.success_tolerance = 0.01;

    real eps = 1e-5;

    initialize_optimization(&opt);
    n_con(&opt);

    auto x = init_guess(&opt);
    bspline.compute_trajectory(x, opt.task);

    Array result(opt.constraints.n_constraints);
    compute_constraints(result.data, x, &opt);

    Matrix grad(x_len, opt.constraints.n_constraints);

    Array r_plus(opt.constraints.n_constraints);

    Array x_plus(x_len);
    x_plus = x;

    int result_idx = 0;
    for (u32 j = 0; j < x_len; j++) {
      constexpr real eps = 1e-5;
      // todo: only copy value that changed last j
      x_plus[j] += eps;

      Array xv_plus;
      xv_plus.alias(x_plus.data, x_len);

      // compute constraints
      compute_constraints(r_plus.data, x_plus, &opt);

      for (u32 i = 0; i < opt.constraints.n_constraints; i++)
        grad(j, i) = (r_plus[i] - result[i]) / eps;

      x_plus[j] = x[j];
    }

    // New implementation

    const int n_segments                = (int) opt.bspline.n_ctrl - (int) opt.bspline.p;
    const int n_points_per_segment      = (int) opt.bspline.n_points / n_segments; // todo: check if fine?
    const int n_ctrl                    = (int) opt.bspline.n_ctrl;
    const int n_constraints_per_segment = (n_joints * 4);                          //  todo: remove hard-coded 4
    // Assert(constraints.size == n_segments * n_constraints_per_segment);

    Array  constraints(n_segments * n_constraints_per_segment);
    Matrix grad_analytic(x_len, n_segments * n_constraints_per_segment);

    // basis: n_ctrl x n_points
    // Assert(constraints.is_alias);
    if (grad_analytic.size) {
      // Assert(grad_analytic.is_alias);
      Assert(grad_analytic.rows == x.size);
      Assert(grad_analytic.cols == constraints.size);
    }

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


    opt.bspline.compute_trajectory(x, opt.task);


    for (int segment = 0; segment < n_segments; segment++) {
      Array max_pos_constraints(n_joints, -INF_REAL); // note:
      Array max_vel_constraints(n_joints, -INF_REAL); // note:
      Array max_acc_constraints(n_joints, -INF_REAL); // note:
      Array max_tor_constraints(n_joints, -INF_REAL); // note:
#if BLAST_TRACE_LEVEL >= 2
      ZoneScopedN("All Segment Constraints");
#endif
      const int first_affected_control_point = std::max(3, segment);
      const int last_affected_control_point  = std::min((n_ctrl - 1) - 3, segment + (int) opt.bspline.p);
      const int n_affected_control_points    = last_affected_control_point - first_affected_control_point + 1; // note: added +1 the bounds are inclusive are inclusive
      Assert(n_affected_control_points >= 3);
      Assert(n_affected_control_points <= 6);

      // Find the start_point_for_segment, so tha twe can access the right memory, and also iterate the points of the given segment
      const int start_point_for_segment = segment * n_points_per_segment; // note:

      // Matrices are always the size of n_ctrl * n_points_per_segment. They "slide" with the segment & do not deal with the rest
      Matrix bp(&opt.bspline.basis_p(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
      Matrix bv(&opt.bspline.basis_v(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:
      Matrix ba(&opt.bspline.basis_a(0, start_point_for_segment), n_ctrl, n_points_per_segment); // note:


      for (int point_in_segment = 0; point_in_segment < n_points_per_segment; point_in_segment++) { // note:
        auto p = opt.bspline.traj.pos.col(start_point_for_segment + point_in_segment);              // note:
        auto v = opt.bspline.traj.vel.col(start_point_for_segment + point_in_segment);              // note:
        auto a = opt.bspline.traj.acc.col(start_point_for_segment + point_in_segment);              // note:

        forward_kinematics(opt.manip, manip_data, p);
        dynamics(opt.manip, manip_data, v, a);

        for (int j = 0; j < n_joints; j++) {
          {
            // todo (andre): do this when initializing the optimization, not here
            // todo: document the current behaviour in the API
            //        (doesn't currently work if one is inf and the other is not)
            if (pmax[j] == INF_REAL)  // note: replace INF_REAL with huge value
              pmax[j] = 1e300;
            if (pmin[j] == -INF_REAL) // note: replace -INF_REAL with huge negative value
              pmin[j] = -1e300;
          }
          // position
          if (auto c = bound_constraint(1.01 * p[j], pmin[j], pmax[j]); // note:
              c > max_pos_constraints[j]) {
            max_pos_constraints[j] = c;
            max_pos_indices[j]     = point_in_segment; // note:
          }
          // velocity
          if (auto c = std::abs(1.01 * v[j]) / vmax[j] - 1.0; // note:
              c > max_vel_constraints[j]) {
            max_vel_constraints[j] = c;
            max_vel_indices[j]     = point_in_segment; // note:
          }
          // acceleration
          if (auto c = std::abs(1.01 * a[j]) / amax[j] - 1.0; // note:
              c > max_acc_constraints[j]) {
            max_acc_constraints[j] = c;
            max_acc_indices[j]     = point_in_segment; // note:
          }
          // torque
          if (auto c = std::abs(1.01 * manip_data.efforts[j]) / tau_max[j] - 1.0; // note:
              c > max_tor_constraints[j]) {
            max_tor_constraints[j] = c;
            max_tor_indices[j]     = point_in_segment; // note:
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
      if (grad_analytic.size) {
        // Matrix (alias) in which we can insert the gradient for the current segment
        Matrix grad_segment(&grad_analytic(0, segment * n_constraints_per_segment), x_len, n_constraints_per_segment);
        Assert(grad_segment.is_alias);

        int con = 0;

        // positions
        for (int joint = 0; joint < n_joints; joint++) {

          // Array of the column where to put the gradient for the current constraint
          Array fill_column = grad_segment.col(con);
          Assert(con == joint); // note:
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          // Which values in 'x' affect the current joint's position
          auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
          // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
          x_idx += first_affected_control_point - 3;

          // 3 to 6 basis functions depending on the segment (first and last 3 control points are not in x)
          Array bp_to_use(&bp(first_affected_control_point, max_pos_indices[joint]), n_affected_control_points);
          Assert(bp_to_use.is_alias);

          real coeff = 2.0 * std::abs(opt.bspline.traj.pos(joint, start_point_for_segment + max_pos_indices[joint])) / (pmax[joint] - pmin[joint]); // note:

          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx++] = bp_to_use[i] * coeff;
          }

          // note: dp/dT == 0
          con++;
        }

        // velocities
        auto one_over_T = 1 / opt.bspline.traj.t.back();
        for (int joint = 0; joint < n_joints; joint++) {
          // Array of the column where to put the gradient for the current constraint
          Array fill_column = grad_segment.col(con);
          Assert(con == n_joints + joint); // note:
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          // Which values in 'x' affect the current joint's position
          auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
          // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
          x_idx += first_affected_control_point - 3;

          Array bv_to_use(&bv(first_affected_control_point, max_vel_indices[joint]), n_affected_control_points);
          Assert(bv_to_use.is_alias);

          real coeff = sign(opt.bspline.traj.vel(joint, start_point_for_segment + max_vel_indices[joint])) / vmax[joint] * one_over_T; // note:

          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx++] = bv_to_use[i] * coeff;
          }

          // dvj/dT = - (Cv + 1) / T
          fill_column.back() = -(max_vel_constraints[joint] + 1) * one_over_T;

          con++;
        }

        // accelerations
        auto one_over_T2 = one_over_T * one_over_T;
        for (int joint = 0; joint < n_joints; joint++) {
          // Array of the column where to put the gradient for the current constraint
          Array fill_column = grad_segment.col(con);
          Assert(con == 2 * n_joints + joint); // note:
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          // Which values in 'x' affect the current joint's position
          auto x_idx = joint * (n_ctrl - 6); // todo: does not work with tasks that don't impose p,v,a for every joint!!
          // shift to the first affected control point keeping in mind that the first 3 are not in the optimization vector
          x_idx += first_affected_control_point - 3;

          Array ba_to_use(&ba(first_affected_control_point, max_acc_indices[joint]), n_affected_control_points);
          Assert(ba_to_use.is_alias);

          real coeff = sign(opt.bspline.traj.acc(joint, start_point_for_segment + max_acc_indices[joint])) / amax[joint] * one_over_T2; // note:

          for (int i = 0; i < n_affected_control_points; i++) {
            fill_column[x_idx++] = ba_to_use[i] * coeff;
          }

          // daj/dT = -2 * (Ca + 1) / T
          fill_column.back() = -2 * (max_acc_constraints[joint] + 1) * one_over_T;

          con++;
        }

        // torque
        // [dt0/dp0, dt1/dp0, ..., dt4/dp0, dt5/dp0]
        // [dt0/dp1, dt1/dp1, ..., dt4/dp1, dt5/dp1]
        // [dt0/dp2, dt1/dp2, ..., dt4/dp2, dt5/dp2]
        // [dt0/dp3, dt1/dp3, ..., dt4/dp3, dt5/dp3]
        // [.
        // [.
        // [.
        for (int joint = 0; joint < n_joints; joint++) {
          Array fill_column = grad_segment.col(con);
          Assert(con == 3 * n_joints + joint); // note:
          Assert(fill_column.size == x_len);
          Assert(fill_column.is_alias);

          constexpr real eps            = 1e-5;
          const auto     point          = max_tor_indices[joint];
          auto           p              = opt.bspline.traj.pos.col(start_point_for_segment + point); // note:
          auto           v              = opt.bspline.traj.vel.col(start_point_for_segment + point); // note:
          auto           a              = opt.bspline.traj.acc.col(start_point_for_segment + point); // note:
          auto           old_constraint = max_tor_constraints[joint];                                // note:
          auto           tau_max_now    = tau_max[joint];                                            // note:

          auto p_plus = p;
          auto v_plus = v;
          auto a_plus = a;

          Array bp_to_use(&bp(first_affected_control_point, point), n_affected_control_points);
          Array bv_to_use(&bv(first_affected_control_point, point), n_affected_control_points);
          Array ba_to_use(&ba(first_affected_control_point, point), n_affected_control_points);

          auto x_idx      = first_affected_control_point - 3;
          auto x_idx_skip = n_ctrl - 6;

          fill_column.back() = 0;
          for (int j = 0; j < n_joints; j++) {

            // partial derivative of torque constraints w.r.t. position
            // finite difference on position
            p_plus[j] += eps;
            // compute the derivative of constraint(joint) w.r.t. theta(j). (remember, joint != j)
            forward_kinematics(opt.manip, manip_data, p_plus);
            dynamics(opt.manip, manip_data, v_plus, a_plus);
            const real new_constraint_p = std::abs(1.01 * manip_data.efforts[joint]) / tau_max_now - 1; // note:
            const real dtau_dp          = (new_constraint_p - old_constraint) / eps;
            // reset finite difference
            p_plus[j] = p[j];

            // note: reset forward kinematics because 'v' and 'a' don't change it
            forward_kinematics(opt.manip, manip_data, p); // todo: precompute once

            // partial derivative of torque constraints w.r.t. velocity
            // finite difference on velocity
            v_plus[j] += eps;
            dynamics(opt.manip, manip_data, v_plus, a_plus);
            const real new_constraint_v = std::abs(1.01 * manip_data.efforts[joint]) / tau_max_now - 1; // note:
            const real dtau_dv          = (new_constraint_v - old_constraint) / eps;
            // reset finite difference
            v_plus[j] = v[j];

            // partial derivative of torque constraints w.r.t. acceleration
            // finite difference on acceleration
            a_plus[j] += eps;
            dynamics(opt.manip, manip_data, v_plus, a_plus);
            const real new_constraint_a = std::abs(1.01 * manip_data.efforts[joint]) / tau_max_now - 1; // note:
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
            // gradient w.r.t. T note: fix
            fill_column.back() += dtau_dv * (-v[j] * one_over_T) + dtau_da * (-2 * a[j] * one_over_T); // note:
          }
          con++;                                                                                       // note: moved out of the inner j loop (only changes at the end of all j torques per joint)
        }

        // Tests for grad of segment compared to the old compute_constraint function
        // position
        for (int joint = 0; joint < n_joints; joint++) {
          CHECK(is_close(grad_segment.col(joint), grad.col(n_constraints_per_segment * (start_point_for_segment + max_pos_indices[joint]) + joint), 0.02));
          CHECK(std::abs(max_pos_constraints[joint] - result[n_constraints_per_segment * (start_point_for_segment + max_pos_indices[joint]) + joint]) < eps);
        }
        // velocity
        for (int joint = 0; joint < n_joints; joint++) {
          bool condition_1 = is_close(grad_segment.col(n_joints + joint), grad.col(n_constraints_per_segment * (start_point_for_segment + max_vel_indices[joint]) + n_joints + joint), 0.02);
          bool condition_2 = std::abs(opt.bspline.traj.vel(joint, start_point_for_segment + max_vel_indices[joint])) < eps;
          CHECK((condition_1 || condition_2));

          CHECK(std::abs(max_vel_constraints[joint] - result[n_constraints_per_segment * (start_point_for_segment + max_vel_indices[joint]) + n_joints + joint]) < eps);
        }
        // acceleration
        for (int joint = 0; joint < n_joints; joint++) {
          bool condition_1 = is_close(grad_segment.col(2 * n_joints + joint), grad.col(n_constraints_per_segment * (start_point_for_segment + max_acc_indices[joint]) + 2 * n_joints + joint), 0.2);
          bool condition_2 = std::abs(opt.bspline.traj.acc(joint, start_point_for_segment + max_acc_indices[joint])) < eps;

          CHECK((condition_1 || condition_2));
          CHECK(std::abs(max_acc_constraints[joint] - result[n_constraints_per_segment * (start_point_for_segment + max_acc_indices[joint]) + 2 * n_joints + joint]) < eps);
        }
        // torque
        for (int joint = 0; joint < n_joints; joint++) {
          CHECK(is_close(grad_segment.col(3 * n_joints + joint), grad.col(n_constraints_per_segment * (start_point_for_segment + max_tor_indices[joint]) + 3 * n_joints + joint), 0.02));
          CHECK(std::abs(max_tor_constraints[joint] - result[n_constraints_per_segment * (start_point_for_segment + max_tor_indices[joint]) + 3 * n_joints + joint]) < eps);
        }
      }
    }
  }
};

TEST_CASE("Analytic dv/dT & da/dT tests", "[acceleration]") {
  int n_iterations = 1000;

  Optimization opt(get_generic_Link6(), get_link6_task());

  opt.world = get_lab_world();

  Bspline bspline(16, 100, 5, 6);
  opt.bspline = bspline;

  opt.guess.type = Guess::custom;
  opt.guess.x0   = Array(opt.bspline.x_len(opt.task), 2.0);
  // opt.guess.n_shot = 100;

  opt.constraints.position            = true;
  opt.constraints.velocity            = true;
  opt.constraints.acceleration        = true;
  opt.constraints.torque              = true;
  opt.constraints.tcp_speed           = false;
  opt.constraints.self_collisions     = false;
  opt.constraints.external_collisions = false;

  opt.max_tries         = 1;
  opt.success_tolerance = 0.01;

  real eps = 1e-5;

  auto xlen = opt.bspline.x_len(opt.task);

  for (int iter = 0; iter < n_iterations; iter++) {

    initialize_optimization(&opt);
    n_con(&opt);

    auto x = init_guess(&opt);
    bspline.compute_trajectory(x, opt.task);

    Array result(opt.constraints.n_constraints);
    compute_constraints(result.data, x, &opt);

    Array x_plus(xlen);
    x_plus = x;
    Array grad(opt.constraints.n_constraints);
    Array grad_analytic(opt.constraints.n_constraints);

    // last point T
    u32 j = xlen - 1;
    // memcpy(x_plus.data, x, xlen * sizeof(real)); todo: check if necessary
    x_plus[j] += eps;

    // compute constraints
    Array r_plus_T(opt.constraints.n_constraints);
    compute_constraints(r_plus_T.data, x_plus, &opt);

    for (u32 i = 0; i < opt.constraints.n_constraints; i++) {
      grad[i] = (r_plus_T[i] - result[i]) / eps;
    }

    int fill_idx = 24;
    for (int point = 1; point < opt.bspline.n_points; point++) {
      // Position is all zeros
      fill_idx += opt.manip.n_joints;
      // Velocity
      for (int joint = 0; joint < opt.manip.n_joints; joint++) {
        grad_analytic[fill_idx] = -(result[fill_idx] + 1) / x[xlen - 1];
        CHECK(std::abs(grad[fill_idx] - grad_analytic[fill_idx]) < 1e-4);
        fill_idx++;
      }
      // Acceleration
      for (int joint = 0; joint < opt.manip.n_joints; joint++) {
        grad_analytic[fill_idx] = -2.0 * (result[fill_idx] + 1) / x[xlen - 1];
        CHECK(std::abs(grad[fill_idx] - grad_analytic[fill_idx]) < 1e-4);
        fill_idx++;
      }

      for (int joint = 0; joint < opt.manip.n_joints; joint++) {
        fill_idx++;
      }
    }
  }
};

} // namespace blast
