#pragma once
#include <blast>

namespace blast {

inline blast_fn Bspline::Bspline(u32 n_control, u32 n_points, u32 P, u32 n_joints) :
    traj(n_points, n_joints),
    control(n_control, n_joints),
    basis_p(n_control, n_points),
    basis_v(n_control, n_points),
    basis_a(n_control, n_points),
    n_joints(n_joints),
    n_points(n_points),
    n_ctrl(n_control),
    degree(degree),
    lower_bounds(n_control),
    upper_bounds(n_control) {
  compute_basis();
}

inline blast_fn u32 Bspline::x_len(const Matrix& task) const {
  Assert(task.rows == n_joints);
  Assert(task.cols == 6);

  // minimum number of variables (number of free control points + total time)
  auto results = n_joints * (n_ctrl - 6) + 1;

  // add extra variables for less defined tasks
  for (u32 i = 0; i < task.size; i++)
    if (std::isnan(task.data[i]))
      results++;
  return results;
}

inline blast_fn void Bspline::compute_basis() {
  u32   m = n_ctrl + degree;
  Array knots(m + 1);
  {
    for (u32 i = m; i > m - degree - 1; i--)
      knots[i] = 1.0f;
    const real du = 1.0f / (real) (m + 1 - 2 * (degree + 1) + 1);
    for (u32 i = degree + 1; i < m - degree; i++)
      knots[i] = knots[i - 1] + du;
  }

  Array      N(m * (degree + 1)); // triangle basis function
  const real du = 1.0f / (n_points - 1);
  zero(basis_p);
  zero(basis_v);
  zero(basis_a);
  real* basis_p_col = basis_p.data;
  real* basis_v_col = basis_v.data;
  real* basis_a_col = basis_a.data;
  for (u32 point = 0; point < n_points; point++) {
    const real u = point * du;

    for (u32 i = 0; i < m; i++)
      N[i] = u >= knots[i] && u < knots[i + 1] ? 1.0f : 0.0f;
    if (point == n_points - 1)
      N[n_ctrl - 1] = 1.0f;
    for (u32 pi = 1; pi <= degree; pi++) {
      for (u32 i = 0; i < m - pi; i++) {
        if (knots[i + pi] != knots[i])
          N[m * pi + i] = N[m * (pi - 1) + i] * (u - knots[i]) / (knots[i + pi] - knots[i]);
        else
          N[m * pi + i] = 0.0f;
        if (knots[i + pi + 1] != knots[i + 1])
          N[m * pi + i] += N[m * (pi - 1) + i + 1] * (knots[i + pi + 1] - u) / (knots[i + pi + 1] - knots[i + 1]);
      }
    }

    // position basis functions
    for (u32 i = 0; i < n_ctrl; i++)
      basis_p_col[i] = N[m * degree + i];

    // velocity basis functions
    for (u32 i = 0; i < n_ctrl - 1; i++)
      basis_v_col[i] = -(real) degree * N[m * (degree - 1) + i + 1] / (knots[i + degree + 1] - knots[i + 1]);
    for (u32 i = 1; i < n_ctrl; i++)
      basis_v_col[i] += (real) degree * N[m * (degree - 1) + i] / (knots[i + degree] - knots[i]);

    // acceleration basis functions
    for (u32 i = 0; i < n_ctrl - 2; i++)
      basis_a_col[i] = (degree * (degree - 1)) * N[m * (degree - 2) + i + 2] / ((knots[i + degree + 1] - knots[i + 1]) * (knots[i + degree + 1] - knots[i + 2]));
    for (u32 i = 1; i < n_ctrl - 1; i++)
      basis_a_col[i] -= (degree * (degree - 1)) * N[m * (degree - 2) + i + 1] * (1.0f / (knots[i + degree] - knots[i]) + 1.0f / (knots[i + degree + 1] - knots[i + 1])) / (knots[i + degree] - knots[i + 1]);
    for (u32 i = 2; i < n_ctrl; i++)
      basis_a_col[i] += (degree * (degree - 1)) * N[m * (degree - 2) + i] / ((knots[i + degree - 1] - knots[i]) * (knots[i + degree] - knots[i]));

    // increment pointers
    basis_p_col += n_ctrl;
    basis_v_col += n_ctrl;
    basis_a_col += n_ctrl;
  }

  // find lower bound & upper bound for each nctrl [lower_bounds, upper_bounds] todo: accelerate
  for (u32 i = 0; i < n_ctrl; i++) {
    // we do not test the first and last points
    for (u32 point = 1; point < n_points - 1; point++) {
      if (basis_p(i, point) != 0.0) {
        lower_bounds[i] = point;
        break;
      }
    }
    for (u32 point = 1; point < n_points - 1; point++) {
      if (basis_p(i, point) != 0.0) {
        upper_bounds[i] = point;
      }
    }
  }
}

inline blast_fn void Bspline::compute_basis_open() {
  u32 m = n_ctrl + degree;

  Array knots(m + 1);
  knots[0]     = 0.0;
  knots.back() = 1.0;
  {
    const real du = 1.0f / (real) m;
    for (u32 i = 1; i < m; i++)
      knots[i] = i * du;
  }

#ifdef BLAST_DEBUG
  print(knots);
#endif

  Matrix     N(m, (degree + 1)); // triangle basis function
  const real du = 1.0f / (n_points - 1);
  zero(basis_p);
  zero(basis_v);
  zero(basis_a);
  real* basis_p_col = basis_p.data;
  real* basis_v_col = basis_v.data;
  real* basis_a_col = basis_a.data;
  for (u32 point = 0; point < n_points; point++) {
    const real u = point * du;

    for (u32 i = 0; i < m; i++)
      N(i, 0) = u >= knots[i] && u <= knots[i + 1] ? 1.0f : 0.0f;
    for (u32 pi = 1; pi <= degree; pi++) {
      for (u32 i = 0; i < m - pi; i++) {
        if (knots[i + pi] != knots[i])
          N(i, pi) = N(i, pi - 1) * (u - knots[i]) / (knots[i + pi] - knots[i]);
        else
          N(i, pi) = 0.0f;
        if (knots[i + pi + 1] != knots[i + 1])
          N(i, pi) += N(i + 1, pi - 1) * (knots[i + pi + 1] - u) / (knots[i + pi + 1] - knots[i + 1]);
      }
    }

    // position basis functions
    for (u32 i = 0; i < n_ctrl; i++)
      basis_p_col[i] = N(i, degree);

    // velocity basis functions
    for (u32 i = 0; i < n_ctrl - 1; i++)
      basis_v_col[i] = -(real) degree * N(i + 1, degree - 1) / (knots[i + degree + 1] - knots[i + 1]);
    for (u32 i = 1; i < n_ctrl; i++)
      basis_v_col[i] += (real) degree * N(i, degree - 1) / (knots[i + degree] - knots[i]);

    // acceleration basis functions
    for (u32 i = 0; i < n_ctrl - 2; i++)
      basis_a_col[i] = (real) (degree * (degree - 1)) * N(i + 2, degree - 2) / ((knots[i + degree + 1] - knots[i + 1]) * (knots[i + degree + 1] - knots[i + 2]));
    for (u32 i = 1; i < n_ctrl - 1; i++)
      basis_a_col[i] -= (real) (degree * (degree - 1)) * N(i + 1, degree - 2) * (1.0f / (knots[i + degree] - knots[i]) + 1.0f / (knots[i + degree + 1] - knots[i + 1])) / (knots[i + degree] - knots[i + 1]);
    for (u32 i = 2; i < n_ctrl; i++)
      basis_a_col[i] += (degree * (degree - 1)) * N(i, degree - 2) / ((knots[i + degree - 1] - knots[i]) * (knots[i + degree] - knots[i]));

    // increment pointers
    basis_p_col += n_ctrl;
    basis_v_col += n_ctrl;
    basis_a_col += n_ctrl;
  }
}

inline blast_fn void Bspline::compute_control(const Array& x, const Matrix& task, real* dst) const {
  using std::isnan;
  Assert(n_ctrl >= 6);
  const real T  = x[x.size - 1];
  const real du = 1.0f / (n_ctrl - degree);
  const real T2 = T * T;

  u32  ctr_i = 0;
  u32  x_i   = 0;
  auto ctr   = dst;

  const real kv = T * du / degree;
  const real ka = 2 * T2 * du * du / (degree * (degree - 1));
  for (u32 joint = 0; joint < n_joints; joint++) {
    // Initial PVA
    const auto pi = task(joint, 0);
    const auto vi = task(joint, 1);
    const auto ai = task(joint, 2);
    ctr[ctr_i++]  = isnan(pi) ? x[x_i++] : pi;
    ctr[ctr_i++]  = isnan(vi) ? x[x_i++] : kv * vi + ctr[ctr_i - 1];
    ctr[ctr_i++]  = isnan(ai) ? x[x_i++] : ka * ai + 3 * ctr[ctr_i - 1] - 2 * ctr[ctr_i - 2];

    // From optimization vector
    for (u32 i = 3; i < n_ctrl - 3; i++)
      ctr[ctr_i++] = x[x_i++];

    // Final PVA
    const auto pf         = task(joint, 3);
    const auto vf         = task(joint, 4);
    const auto af         = task(joint, 5);
    const real pn         = isnan(pf) ? x[x_i++] : pf;
    const real pn_minus_1 = isnan(vf) ? x[x_i++] : pn - kv * vf;
    const real pn_minus_2 = isnan(af) ? x[x_i++] : ka * af - 2 * pn + 3 * pn_minus_1;
    ctr[ctr_i++]          = pn_minus_2;
    ctr[ctr_i++]          = pn_minus_1;
    ctr[ctr_i++]          = pn;
  }
}

/**
 * @brief Computes the trajectory of a B-spline curve.
 *
 * @param x An array of real numbers representing the knots of the B-spline curve.
 * @param task A matrix representing the control points of the B-spline curve.
 *
 * This function computes the trajectory of a B-spline curve given its knots and control points.
 * The knots are passed as an array `x` and the control points are passed as a matrix `task`.
 * The function first checks that the size of `x` is equal to the length of `task` along the x-axis,
 * and that the number of rows and columns in `task` are equal to the number of joints and 6, respectively.
 * Then, it calls the `compute_control` function to compute the control points of the B-spline curve.
 * Finally, it computes the position, velocity, and acceleration of each joint at each point along the trajectory
 * using the computed control points and the basis functions for position, velocity, and acceleration.
 */
inline blast_fn void Bspline::compute_trajectory(const Array& x, const Matrix& task) {
#if BLAST_TRACE_LEVEL >= 2
  PROFILE_FUNCTION;
#endif


  Assert(x.size == x_len(task));
  Assert(task.rows == n_joints);
  Assert(task.cols == 6);

  compute_control(x, task, control.data);

  const real T           = x.back();
  const real dt          = T / (n_points - 1);
  const real one_over_T  = 1 / T;
  const real one_over_T2 = one_over_T * one_over_T;

  for (u32 point = 0; point < n_points; point++) {
    traj.t[point] = dt * point;
    auto bp       = basis_p.col(point);
    auto bv       = basis_v.col(point);
    auto ba       = basis_a.col(point);
    for (u32 joint = 0; joint < n_joints; joint++) {
      auto c                 = control.col(joint);
      traj.pos(joint, point) = dot(c, bp);
      traj.vel(joint, point) = dot(c, bv) * one_over_T;
      traj.acc(joint, point) = dot(c, ba) * one_over_T2;
    }
  }
}

} // namespace blast
