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
    p(P),
    lb(n_control),
    ub(n_control) {
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
  u32   m = n_ctrl + p;
  Array knots(m + 1);
  {
    for (u32 i = m; i > m - p - 1; i--)
      knots[i] = 1.0f;
    const real du = 1.0f / (real) (m + 1 - 2 * (p + 1) + 1);
    for (u32 i = p + 1; i < m - p; i++)
      knots[i] = knots[i - 1] + du;
  }

  Array      N(m * (p + 1)); // triangle basis function
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
    for (u32 pi = 1; pi <= p; pi++) {
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
      basis_p_col[i] = N[m * p + i];

    // velocity basis functions
    for (u32 i = 0; i < n_ctrl - 1; i++)
      basis_v_col[i] = -(real) p * N[m * (p - 1) + i + 1] / (knots[i + p + 1] - knots[i + 1]);
    for (u32 i = 1; i < n_ctrl; i++)
      basis_v_col[i] += (real) p * N[m * (p - 1) + i] / (knots[i + p] - knots[i]);

    // acceleration basis functions
    for (u32 i = 0; i < n_ctrl - 2; i++)
      basis_a_col[i] = (p * (p - 1)) * N[m * (p - 2) + i + 2] / ((knots[i + p + 1] - knots[i + 1]) * (knots[i + p + 1] - knots[i + 2]));
    for (u32 i = 1; i < n_ctrl - 1; i++)
      basis_a_col[i] -= (p * (p - 1)) * N[m * (p - 2) + i + 1] * (1.0f / (knots[i + p] - knots[i]) + 1.0f / (knots[i + p + 1] - knots[i + 1])) / (knots[i + p] - knots[i + 1]);
    for (u32 i = 2; i < n_ctrl; i++)
      basis_a_col[i] += (p * (p - 1)) * N[m * (p - 2) + i] / ((knots[i + p - 1] - knots[i]) * (knots[i + p] - knots[i]));

    // increment pointers
    basis_p_col += n_ctrl;
    basis_v_col += n_ctrl;
    basis_a_col += n_ctrl;
  }

  // find lower bound & upper bound for each nctrl [lb, ub] todo: accelerate
  for (u32 i = 0; i < n_ctrl; i++) {
    // we do not test the first and last points
    for (u32 point = 1; point < n_points - 1; point++) {
      if (basis_p(i, point) != 0.0) {
        lb[i] = point;
        break;
      }
    }
    for (u32 point = 1; point < n_points - 1; point++) {
      if (basis_p(i, point) != 0.0) {
        ub[i] = point;
      }
    }
  }
}

inline int uniformClampedSpan(real u, int n_ctrl, int p) {
  if (u <= 0.0) {
    return p;
  } else if (u >= 1.0) {
    return n_ctrl - 1;
  } else {
    const int n_spans = n_ctrl - p;
    int       span    = p + int(u * n_spans);
    return std::min(span, n_ctrl - 1);
  }
}

// Computes Basis functions up to d, with d = 0 -> position, d = 1 -> velocity, etc.
inline blast_fn void Bspline::compute_basis_derivative(int d) {
  basis.resize(d+1, Matrix(n_ctrl, n_points));

  u32   m = n_ctrl + p;
  Array knots(m + 1);
  {
    for (u32 i = m; i > m - p - 1; i--)
      knots[i] = 1.0f;
    const real du = 1.0f / (real) (m + 1 - 2 * (p + 1) + 1);
    for (u32 i = p + 1; i < m - p; i++)
      knots[i] = knots[i - 1] + du;
  }

  Array      N(m * (p + 1)); // triangle basis function
  const real du = 1.0f / (n_points - 1);
  
  for (u32 point = 0; point < n_points; point++) {
    const real u = point * du;

    // note: could save if statements if we hard-coded span calculation since we know u is always within [0, 1]
    const int span  = uniformClampedSpan(u, n_ctrl, p);
    const int first = span - p;

    // --- Algorithm A2.2 (ndu table) ---
    Matrix ndu(p + 1, p + 1);
    Array  left(p + 1), right(p + 1);

    ndu(0, 0) = 1.0;

    for (int j = 1; j <= p; ++j) {
      left[j]  = u - knots[span + 1 - j];
      right[j] = knots[span + j] - u;

      real saved = 0.0;
      for (int r = 0; r < j; ++r) {
        ndu(j, r) = right[r + 1] + left[j - r];
        real temp = ndu(r, j - 1) / ndu(j, r);

        ndu(r, j) = saved + right[r + 1] * temp;
        saved     = left[j - r] * temp;
      }
      ndu(j, j) = saved;
    }

    // --- Derivative computation (Algorithm A2.3) ---
    Matrix ders(d + 1, p + 1);
    for (int j = 0; j <= p; ++j)
      ders(0, j) = ndu(j, p);

    // Working array a[2][p+1]
    Matrix a(2, p + 1);

    for (int r = 0; r <= p; ++r) {
      int s1 = 0, s2 = 1;
      a(0, 0) = 1.0;

      for (int k = 1; k <= d; ++k) {
        real d  = 0.0;
        int  rk = r - k;
        int  pk = p - k;

        int j1;
        int j2;

        if (r >= k) {
          a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
          d        = a(s2, 0) * ndu(rk, pk);
        }
        if (rk >= -1) {
          j1 = 1;
        } else {
          j1 = -rk;
        }

        if (r - 1 <= pk) {
          j2 = k - 1;
        } else {
          j2 = p - r;
        }

        for (int j = j1; j <= j2; j++) {
          a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
          d += a(s2, j) * ndu(rk + j, pk);
        }

        if (r <= pk) {
          a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
          d += a(s2, k) * ndu(r, pk);
        }

        ders(k, r) = d;
        std::swap(s1, s2);
      }
    }

    // --- Multiply by factorial terms ---
    real factor = real(p);
    for (int k = 1; k <= d; ++k) {
      for (int j = 0; j <= p; ++j)
        ders(k, j) *= factor;
      factor *= real(p - k);
    }

    // --- Scatter into global result ---
    for (int k = 0; k <= d; ++k)
      for (int j = 0; j <= p; ++j)
        basis[k](first + j, point) = ders(k, j);
  }
}

inline blast_fn void Bspline::compute_basis_open() {
  u32 m = n_ctrl + p;

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

  Matrix     N(m, (p + 1)); // triangle basis function
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
    for (u32 pi = 1; pi <= p; pi++) {
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
      basis_p_col[i] = N(i, p);

    // velocity basis functions
    for (u32 i = 0; i < n_ctrl - 1; i++)
      basis_v_col[i] = -(real) p * N(i + 1, p - 1) / (knots[i + p + 1] - knots[i + 1]);
    for (u32 i = 1; i < n_ctrl; i++)
      basis_v_col[i] += (real) p * N(i, p - 1) / (knots[i + p] - knots[i]);

    // acceleration basis functions
    for (u32 i = 0; i < n_ctrl - 2; i++)
      basis_a_col[i] = (real) (p * (p - 1)) * N(i + 2, p - 2) / ((knots[i + p + 1] - knots[i + 1]) * (knots[i + p + 1] - knots[i + 2]));
    for (u32 i = 1; i < n_ctrl - 1; i++)
      basis_a_col[i] -= (real) (p * (p - 1)) * N(i + 1, p - 2) * (1.0f / (knots[i + p] - knots[i]) + 1.0f / (knots[i + p + 1] - knots[i + 1])) / (knots[i + p] - knots[i + 1]);
    for (u32 i = 2; i < n_ctrl; i++)
      basis_a_col[i] += (p * (p - 1)) * N(i, p - 2) / ((knots[i + p - 1] - knots[i]) * (knots[i + p] - knots[i]));

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
  const real du = 1.0f / (n_ctrl - p);
  const real T2 = T * T;

  u32  ctr_i = 0;
  u32  x_i   = 0;
  auto ctr   = dst;

  const real kv = T * du / p;
  const real ka = 2 * T2 * du * du / (p * (p - 1));
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

// Algorithm A2.3 – Basis Function Derivatives (Piegl & Tiller)
// Returns the deriv_order-th derivative of all non-zero B-spline basis functions
// aligned with the control point indices.
// Input
// u : Point of evaluation, within [0, 1] (u)
// i : u is between knot i and knot i+1 (span)
// p : Bspline degree (p)
// n : Derivative degree (deriv_order)
// U : Knot vector (knots)
// Output
// ders : Basis functions for derivative (result)
inline Array BsplineDerivative_book(real u, int n_ctrl, int p, int deriv_order) {
  ZoneScoped;

  Array result(n_ctrl);

  if (deriv_order > p)
    return result;

  // --- Uniform clamped knot vector ---
  Array knots(n_ctrl + p + 1);
  for (int i = 0; i <= p; ++i)
    knots[i] = 0.0;
  for (int i = p + 1; i < n_ctrl; ++i)
    knots[i] = real(i - p) / real(n_ctrl - p);
  for (int i = n_ctrl; i <= n_ctrl + p; ++i)
    knots[i] = 1.0;

  const int span  = uniformClampedSpan(u, n_ctrl, p);
  const int first = span - p;

  // --- Algorithm A2.2 (ndu table) ---
  Matrix ndu(p + 1, p + 1);
  Array  left(p + 1), right(p + 1);

  ndu(0, 0) = 1.0;

  for (int j = 1; j <= p; ++j) {
    left[j]  = u - knots[span + 1 - j];
    right[j] = knots[span + j] - u;

    real saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu(j, r) = right[r + 1] + left[j - r];
      real temp = ndu(r, j - 1) / ndu(j, r);

      ndu(r, j) = saved + right[r + 1] * temp;
      saved     = left[j - r] * temp;
    }
    ndu(j, j) = saved;
  }

  // --- Derivative computation (Algorithm A2.3) ---
  Matrix ders(deriv_order + 1, p + 1);
  for (int j = 0; j <= p; ++j)
    ders(0, j) = ndu(j, p);

  // Working array a[2][p+1]
  Matrix a(2, p + 1);

  for (int r = 0; r <= p; ++r) {
    int s1 = 0, s2 = 1;
    a(0, 0) = 1.0;

    for (int k = 1; k <= deriv_order; ++k) {
      real d  = 0.0;
      int  rk = r - k;
      int  pk = p - k;

      int j1;
      int j2;

      if (r >= k) {
        a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
        d        = a(s2, 0) * ndu(rk, pk);
      }
      if (rk >= -1) {
        j1 = 1;
      } else {
        j1 = -rk;
      }

      if (r - 1 <= pk) {
        j2 = k - 1;
      } else {
        j2 = p - r;
      }

      for (int j = j1; j <= j2; j++) {
        a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
        d += a(s2, j) * ndu(rk + j, pk);
      }

      if (r <= pk) {
        a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
        d += a(s2, k) * ndu(r, pk);
      }

      ders(k, r) = d;
      std::swap(s1, s2);
    }
  }

  // --- Multiply by factorial terms ---
  real factor = real(p);
  for (int k = 1; k <= deriv_order; ++k) {
    for (int j = 0; j <= p; ++j)
      ders(k, j) *= factor;
    factor *= real(p - k);
  }

  // --- Scatter into global result ---
  for (int j = 0; j <= p; ++j)
    result[first + j] = ders(deriv_order, j);

  return result;
}

} // namespace blast
