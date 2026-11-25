#pragma once

#include <blast>

namespace blast {

// Create n pseudo-random arrays of d dimensions form : d-by-n
// (https://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/#GeneralizingGoldenRatio)
inline blast::Matrix lds_golden_ratio(const u32 n, const u32 d) {
  blast::Matrix result(d, n);

  blast::real g = (real) 1.6180339887498948482; // Recommanded by wikipedia article

  // precompute factors a
  std::vector<blast::real> a;
  a.reserve(d);
  a.push_back(1 / g);
  for (u32 i = 1; i < d; i++)
    a.push_back(a.at(i - 1) * a.at(0));

  for (u32 i = 0; i < n; i++) {
    for (u32 j = 0; j < d; j++) {
      result(j, i) = std::fmod((real) 0.5 + a[j] * (real) (i + 1), (real) 1);
    }
  }

  return result;
}

Array guess_random(Bspline& bspline, Matrix& task) {
  Array x(bspline.x_len(task));
  fill_random(x, 1);
  x.back() = std::abs(x.back()) * 2 + (real) 0.1;
  return x;
}

inline Matrix get_random_guesses(const u32 n, const u32 d) {
  Array  x(d);
  Matrix result(d, n);
  for (int i = 0; i < n; i++) {
    fill_random(x, 1);
    x.back() = std::abs(x.back()) * 5 + (real) 0.1;
    for (int j = 0; j < d; j++) {
      result(j, i) = x[j];
    }
  }
  return result;
}

inline Array guess_shot_mean_segments(Optimization* opt) { // no collisions
  Array best_x(opt->bspline.x_len(opt->task));
  real  best_val = INF_REAL;
  for (int idx_nshot = 0; idx_nshot < opt->guess.n_random_shots; idx_nshot++) {
    auto x = guess_random((opt->bspline), opt->task);
    opt->bspline.compute_trajectory(x, opt->task);
    Array  c1(opt->constraints.n_constraints);
    Matrix gradient;
    constraints_and_gradients_with_segments(x, *opt, c1, gradient); // todo: change for segments
    real r = 0;
    for (u32 i = 0; i < c1.size; i++)
      r += std::max({c1[i], (real) 0});
    Assert(!isnan(r));
    r = r * x.back(); // todo: Evaluate time estimate impact on trajectory
    if (r < best_val) {
      best_x   = x;
      best_val = r;
    }
  }
  return best_x;
}

inline Array get_best_x_segments(Optimization* opt) {
  real  best_f = INF_REAL;
  Array best_x(opt->guess.candidates.rows);
  for (u32 i = 0; i < opt->guess.candidates.cols; i++) {
    auto current_x = opt->guess.candidates.col(i);
    auto current_f = compute_objective(current_x, opt);

    best_x = current_f < best_f ? current_x : best_x;
    best_f = current_f < best_f ? current_f : best_f;
  }

  return best_x;
}

// Deterministic guess: straight line start->goal, time from manip's vel/accel limits
// (per-joint trapezoidal/triangular, max across joints). Requires fully specified task.
inline Array guess_straight_line(Optimization* opt) {
  Array     x(opt->bspline.x_len(opt->task));
  const u32 n_free = opt->bspline.n_ctrl - 6;
  Assert(x.size == (u32) opt->manip.n_joints * n_free + 1);

  real total_time = 0;
  u32  k          = 0;
  for (int j = 0; j < opt->manip.n_joints; j++) {
    const real start = opt->task(j, 0);
    const real goal  = opt->task(j, 3);
    const real dist  = std::abs(goal - start);

    const real vmax = opt->manip.velocity_max[j];
    const real amax = opt->manip.acceleration_max[j];
    if (vmax > 0 && amax > 0) {
      const real d_switch = vmax * vmax / amax; // dist covered while accelerating to vmax and back down
      const real t_j      = (dist >= d_switch) ? dist / vmax + vmax / amax : 2 * std::sqrt(dist / amax);
      total_time          = std::max(total_time, t_j);
    }

    for (u32 i = 0; i < n_free; i++) {
      const real a = (real) (i + 1) / (real) (n_free + 1);
      x[k++]       = start + a * (goal - start);
    }
  }
  x.back() = total_time > 0 ? total_time : 2;
  return x;
}

/**
 * @brief Initial guess that follows a joint-space waypoint path, e.g. from a sampling-based planner.
 *
 * Fits the B-spline (this optimizer's own knots and basis) to the path by least squares:
 * - The path is a polyline in joint space. Its length is the sum of the Euclidean joint-space
 *   distances between consecutive waypoints.
 * - At `n_samples` evenly spaced spline parameters u_k = k / (n_samples - 1) in [0, 1] (evenly
 *   spaced in time too, since t = u T), the target is the point a fraction u_k of the way along that
 *   length. E.g. the curve at u = 0.3 should be 30 % of the way along the path. This is not waypoint
 *   0.3 * (n - 1), so unevenly spaced waypoints do not skew the fit.
 * - The free control points minimize the sum of squared distances between the curve and these
 *   targets. The first and last three control points stay pinned, as Bspline::compute_control sets
 *   them for a stop-to-stop task, so the curve starts and ends at rest.
 * - The fitted control points sit just outside the path's corners, so the curve follows the path
 *   instead of cutting them.
 *
 * The duration T is the shortest for which every control point of the velocity and acceleration
 * curves (the derivative control polygons) is within the manipulator's limits, and at least
 * opt->min_duration. By the convex-hull property this is a conservative choice for the whole curve,
 * not only at the samples.
 *
 * Requires a stop-to-stop task with finite start and goal (zero boundary velocity and acceleration).
 * The path is expected to start at the task's start and end at its goal: the curve is pinned to the
 * task's, so a mismatch shows up as fit error near the ends. A path of zero length gives
 * guess_straight_line().
 *
 * @param opt        Optimization whose bspline, task and manipulator limits define the guess.
 * @param path       Waypoints as an n_joints x n_waypoints Matrix, one column per waypoint (the layout of
 *                   Bspline trajectories); at least two columns.
 * @param n_samples  Points the curve is fitted at; 0 (the default) means max(100, 4 * n_segments), with
 *                   n_segments = n_ctrl - degree. An explicit value is used as given; if it is too few to
 *                   determine the free control points (about n_ctrl - 4), a message is printed and the
 *                   default is used instead.
 * @return The optimization vector x (free control points, joint-major, then T), laid out as
 *         guess_straight_line() returns it.
 */
inline Array guess_from_path(Optimization* opt, const Matrix& path, u32 n_samples = 0) {
  const u32 n_ctrl   = opt->bspline.n_ctrl;
  const u32 degree   = opt->bspline.degree;
  const u32 n_joints = opt->bspline.n_joints;
  const u32 n_free   = n_ctrl - 6;
  Assert(n_ctrl > 6);
  Assert(path.rows == n_joints && path.cols >= 2);
  const u32 default_samples = std::max((u32) 100, 4 * (n_ctrl - degree));
  if (n_samples == 0)
    n_samples = default_samples;
  for (u32 j = 0; j < n_joints; j++) {
    Assert(std::isfinite(opt->task(j, 0)) && std::isfinite(opt->task(j, 3)));
    Assert(opt->task(j, 1) == 0 && opt->task(j, 2) == 0 && opt->task(j, 4) == 0 && opt->task(j, 5) == 0);
  }

  // cumulative[i]: distance along the path from its start to waypoint i (sum of Euclidean
  // joint-space segment lengths). This only places the targets; nothing is minimized here. The
  // square root is needed: squared lengths do not add up along a path.
  std::vector<real> cumulative(path.cols, 0);
  for (u32 i = 1; i < path.cols; i++) {
    real squared = 0;
    for (u32 j = 0; j < n_joints; j++) {
      const real d = path(j, i) - path(j, i - 1);
      squared += d * d;
    }
    cumulative[i] = cumulative[i - 1] + std::sqrt(squared);
  }
  const real length = cumulative.back();
  if (!(length > 0))
    return guess_straight_line(opt);

  // Full control polygon (n_ctrl x n_joints). The first and last three rows are pinned to the task's
  // start and goal, as compute_control sets them for a stop-to-stop task; the rows between are fitted below.
  Matrix control(n_ctrl, n_joints);
  for (u32 j = 0; j < n_joints; j++) {
    for (u32 i = 0; i < 3; i++) {
      control(i, j)              = opt->task(j, 0);
      control(n_ctrl - 1 - i, j) = opt->task(j, 3);
    }
  }

  // Least squares: over the free control points, minimize sum_k |C(u_k) - target_k|^2, C being the
  // curve. Normal equations (A^T A) p_j = A^T r_j, with one n_free x n_free matrix shared by every joint j:
  //   A_ki = N_{3+i}(u_k), the basis function of free control point i at sample k;
  //   r_kj = target_kj minus the pinned control points' contribution at u_k.
  // Basis values come from this optimizer's own Bspline, sampled at u_k = k / (count - 1).
  //
  // build_and_factor(count) builds the normal equations from `count` samples and factors them in place
  // (Cholesky, normal = L L^T in the lower triangle). It returns false when the samples are too few to
  // determine the free control points: a pivot then drops to rounding level relative to the largest
  // diagonal entry, while a well-posed fit keeps every pivot within a small factor of it.
  std::vector<real> normal, rhs, point(n_joints);
  const auto        build_and_factor = [&](u32 count) -> bool {
    if (count < std::max(n_free, (u32) 2))
      return false;
    normal.assign(n_free * n_free, 0);
    rhs.assign(n_free * n_joints, 0);
    const Bspline samples(n_ctrl, count, degree, n_joints);
    u32           segment = 0;
    for (u32 k = 0; k < count; k++) {
      // target_k: the point at distance u_k * length along the path, interpolated linearly within its
      // segment. u_k only grows, so the segment index only advances.
      const real distance = length * (real) k / (real) (count - 1);
      while (segment + 2 < path.cols && cumulative[segment + 1] < distance)
        segment++;
      const real span     = cumulative[segment + 1] - cumulative[segment];
      const real fraction = span > 0 ? std::clamp((distance - cumulative[segment]) / span, (real) 0, (real) 1) : (real) 0;
      for (u32 j = 0; j < n_joints; j++)
        point[j] = path(j, segment) + fraction * (path(j, segment + 1) - path(j, segment));

      for (u32 j = 0; j < n_joints; j++) {
        real residual = point[j];
        for (u32 i = 0; i < 3; i++)
          residual -= samples.basis_p(i, k) * control(i, j) + samples.basis_p(n_ctrl - 1 - i, k) * control(n_ctrl - 1 - i, j);
        for (u32 a = 0; a < n_free; a++)
          rhs[a * n_joints + j] += samples.basis_p(3 + a, k) * residual;
      }
      for (u32 a = 0; a < n_free; a++)
        for (u32 b = 0; b < n_free; b++)
          normal[a * n_free + b] += samples.basis_p(3 + a, k) * samples.basis_p(3 + b, k);
    }

    real diagonal_max = 0;
    for (u32 a = 0; a < n_free; a++)
      diagonal_max = std::max(diagonal_max, normal[a * n_free + a]);
    const real pivot_min = std::sqrt(std::numeric_limits<real>::epsilon()) * diagonal_max;
    for (u32 a = 0; a < n_free; a++) {
      for (u32 b = 0; b <= a; b++) {
        real sum = normal[a * n_free + b];
        for (u32 c = 0; c < b; c++)
          sum -= normal[a * n_free + c] * normal[b * n_free + c];
        if (a == b) {
          if (!(sum > pivot_min))
            return false;
          normal[a * n_free + a] = std::sqrt(sum);
        } else {
          normal[a * n_free + b] = sum / normal[b * n_free + b];
        }
      }
    }
    return true;
  };
  if (!build_and_factor(n_samples)) {
    printf("guess_from_path: n_samples %u too few for %u control points; using %u.\n", n_samples, n_ctrl, default_samples);
    const bool factored = build_and_factor(default_samples);
    Assert(factored);
    (void) factored;
  }

  std::vector<real> solution(n_free);
  for (u32 j = 0; j < n_joints; j++) {
    for (u32 a = 0; a < n_free; a++) { // L y = rhs
      real sum = rhs[a * n_joints + j];
      for (u32 c = 0; c < a; c++)
        sum -= normal[a * n_free + c] * solution[c];
      solution[a] = sum / normal[a * n_free + a];
    }
    for (u32 a = n_free; a-- > 0;) { // L^T p = y
      real sum = solution[a];
      for (u32 c = a + 1; c < n_free; c++)
        sum -= normal[c * n_free + a] * solution[c];
      solution[a] = sum / normal[a * n_free + a];
    }
    for (u32 a = 0; a < n_free; a++)
      control(3 + a, j) = solution[a];
  }

  // Duration. Velocity is itself a B-spline whose control points are the degree-scaled first
  // differences of the control polygon, divided by T; acceleration likewise uses second differences
  // over T^2. T is the smallest value that brings every one of those control points within its
  // limit (conservative by the convex-hull property).
  // Knots are clamped uniform, as Bspline::compute_basis builds them.
  const auto knot = [&](u32 i) -> real {
    if (i <= degree)
      return 0;
    if (i >= n_ctrl)
      return 1;
    return (real) (i - degree) / (real) (n_ctrl - degree);
  };
  real              duration = opt->min_duration;
  std::vector<real> velocity(n_ctrl - 1);
  for (u32 j = 0; j < n_joints; j++) {
    const real vmax = opt->manip.velocity_max[j];
    const real amax = opt->manip.acceleration_max[j];
    for (u32 i = 0; i + 1 < n_ctrl; i++) {
      velocity[i] = (real) degree * (control(i + 1, j) - control(i, j)) / (knot(i + degree + 1) - knot(i + 1));
      if (vmax > 0)
        duration = std::max(duration, std::abs(velocity[i]) / vmax);
    }
    for (u32 i = 0; i + 2 < n_ctrl; i++) {
      const real acceleration = (real) (degree - 1) * (velocity[i + 1] - velocity[i]) / (knot(i + degree + 1) - knot(i + 2));
      if (amax > 0)
        duration = std::max(duration, std::sqrt(std::abs(acceleration) / amax));
    }
  }

  Array x(opt->bspline.x_len(opt->task));
  Assert(x.size == n_joints * n_free + 1);
  u32 k = 0;
  for (u32 j = 0; j < n_joints; j++)
    for (u32 i = 3; i < n_ctrl - 3; i++)
      x[k++] = control(i, j);
  x.back() = duration;
  return x;
}

inline Array init_guess_segments(Optimization* opt) {
  Array x(opt->bspline.x_len(opt->task));
  switch (opt->guess.type) {
    case Guess::random: {
      x = guess_random(opt->bspline, opt->task);
      break;
    }
    case Guess::shotgun: {
      x = guess_shot_mean_segments(opt);
      break;
    }
    case Guess::straight_line: {
      x = guess_straight_line(opt);
      break;
    }
    case Guess::custom: {
      x = opt->guess.initial_x;
      break;
    }
    case Guess::from_list: {
      x = get_best_x_segments(opt);
      break;
    }
    default:
      Assert(false);
  }
  return x;
}

// todo: Remove collisions
inline Array guess_shot_mean(Optimization* opt) { // no collisions
  Array best_x(opt->bspline.x_len(opt->task));
  real  best_val = INF_REAL;
  for (int idx_nshot = 0; idx_nshot < opt->guess.n_random_shots; idx_nshot++) {
    auto x = guess_random((opt->bspline), opt->task);
    opt->bspline.compute_trajectory(x, opt->task);
    Array c1(opt->constraints.n_constraints);
    compute_constraints(c1.data, x, opt); // todo: change for segments
    real r = 0;
    for (u32 i = 0; i < c1.size; i++)
      r += std::max({c1[i], (real) 0});
    Assert(!isnan(r));
    r = r * x.back(); // todo: Evaluate time estimate impact on trajectory
    if (r < best_val) {
      best_x   = x;
      best_val = r;
    }
  }
  return best_x;
}

inline Array get_best_x(Optimization* opt) {
  real  best_f = INF_REAL;
  Array best_x(opt->guess.candidates.rows);
  for (u32 i = 0; i < opt->guess.candidates.cols; i++) {
    auto current_x = opt->guess.candidates.col(i);
    auto current_f = compute_objective(current_x, opt);

    best_x = current_f < best_f ? current_x : best_x;
    best_f = current_f < best_f ? current_f : best_f;
  }

  return best_x;
}

inline Array init_guess(Optimization* opt) {
  Array x(opt->bspline.x_len(opt->task));
  switch (opt->guess.type) {
    case Guess::random: {
      x = guess_random(opt->bspline, opt->task);
      break;
    }
    case Guess::shotgun: {
      x = guess_shot_mean(opt);
      break;
    }
    case Guess::straight_line: {
      x = guess_straight_line(opt);
      break;
    }
    case Guess::custom: {
      x = opt->guess.initial_x;
      break;
    }
    case Guess::from_list: {
      x = get_best_x(opt);
      break;
    }
    default:
      Assert(false);
  }
  return x;
}

} // namespace blast
