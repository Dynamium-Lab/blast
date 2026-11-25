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
    Assert(!std::isnan(r));
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
 * @brief Distance along a joint-space path from its first waypoint to each waypoint.
 *
 * Entry i is the sum of the Euclidean joint-space lengths of the segments before waypoint i, so the
 * last entry is the path's total length.
 *
 * @param path  Waypoints, n_joints x n_waypoints, one column per waypoint.
 */
inline Array cumulative_length(const Matrix& path) {
  const u32 n_joints    = path.rows;
  const u32 n_waypoints = path.cols;
  Array     cumulative(n_waypoints);
  for (u32 waypoint = 1; waypoint < n_waypoints; waypoint++) {
    real squared_segment_length = 0;
    for (u32 joint = 0; joint < n_joints; joint++) {
      const real joint_step = path(joint, waypoint) - path(joint, waypoint - 1);
      squared_segment_length += joint_step * joint_step;
    }
    cumulative[waypoint] = cumulative[waypoint - 1] + std::sqrt(squared_segment_length);
  }
  return cumulative;
}

/**
 * @brief The point at a given distance along a path, interpolated linearly within the segment that
 *        contains it.
 *
 * @param path        Waypoints, n_joints x n_waypoints.
 * @param cumulative  cumulative_length(path).
 * @param distance    Distance from the first waypoint, in [0, cumulative.back()].
 * @param segment     Index of the segment to start searching from; updated to the segment found. Calls
 *                    with growing distances can pass the same variable, so the search only moves forward.
 * @param point       Output, n_joints values (allocated by the caller).
 */
inline void point_at_distance_along_path(const Matrix& path, const Array& cumulative, real distance, u32& segment, Array& point) {
  const u32 n_joints    = path.rows;
  const u32 n_waypoints = path.cols;
  while (segment + 2 < n_waypoints && cumulative[segment + 1] < distance)
    segment++;
  const real segment_length = cumulative[segment + 1] - cumulative[segment];
  real       fraction       = 0; // of the way from waypoint `segment` to waypoint `segment + 1`
  if (segment_length > 0)
    fraction = std::clamp((distance - cumulative[segment]) / segment_length, (real) 0, (real) 1);
  for (u32 joint = 0; joint < n_joints; joint++)
    point[joint] = path(joint, segment) + fraction * (path(joint, segment + 1) - path(joint, segment));
}

/**
 * @brief Cholesky factorization matrix = L L^T of a symmetric square matrix, in place (L in the lower triangle).
 *
 * @return false if the matrix is not clearly positive definite: a pivot falls to rounding level relative
 *         to the largest diagonal entry, as when the samples are too few to determine the unknowns.
 */
inline bool cholesky_factor(Matrix& matrix) {
  const u32 size         = matrix.rows;
  real      diagonal_max = 0;
  for (u32 row = 0; row < size; row++)
    diagonal_max = std::max(diagonal_max, matrix(row, row));
  const real pivot_min = std::sqrt(std::numeric_limits<real>::epsilon()) * diagonal_max;

  for (u32 row = 0; row < size; row++) {
    for (u32 col = 0; col <= row; col++) {
      real value = matrix(row, col);
      for (u32 inner = 0; inner < col; inner++)
        value -= matrix(row, inner) * matrix(col, inner);
      if (row == col) {
        if (!(value > pivot_min))
          return false;
        matrix(row, row) = std::sqrt(value);
      } else {
        matrix(row, col) = value / matrix(col, col);
      }
    }
  }
  return true;
}

/**
 * @brief Solves L L^T x = b for each column b of `columns`, in place, given the factor L from cholesky_factor().
 *
 * @param factor   The factored square matrix (L in the lower triangle).
 * @param columns  On entry the right-hand sides b, one per column; on return the solutions x.
 */
inline void cholesky_solve(const Matrix& factor, Matrix& columns) {
  const u32 size = factor.rows;
  for (u32 column = 0; column < columns.cols; column++) {
    // Forward substitution: L y = b.
    for (u32 row = 0; row < size; row++) {
      for (u32 col = 0; col < row; col++)
        columns(row, column) -= factor(row, col) * columns(col, column);
      columns(row, column) /= factor(row, row);
    }
    // Back substitution: L^T x = y.
    for (u32 row = size; row-- > 0;) {
      for (u32 col = row + 1; col < size; col++)
        columns(row, column) -= factor(col, row) * columns(col, column);
      columns(row, column) /= factor(row, row);
    }
  }
}

/**
 * @brief Least-squares fit of the free control points of a stop-to-stop B-spline to a path.
 *
 * At `n_samples` evenly spaced parameters u_k = k / (n_samples - 1), the target is the point a fraction
 * u_k of the way along the path. Over the free control points (rows 3 to n_ctrl - 4 of `control`),
 * minimizes sum_k |C(u_k) - target_k|^2, C being the curve, through the normal equations
 * (A^T A) p_j = A^T r_j, with one matrix shared by every joint j:
 *   A_ki = N_{3+i}(u_k), the basis function of free control point i at sample k;
 *   r_kj = target_kj minus the pinned control points' contribution at u_k.
 * The basis values come from a Bspline with the same control points and degree, so the fit uses the
 * optimizer's own knots and basis.
 *
 * @param bspline     The optimizer's B-spline (n_ctrl, degree, n_joints).
 * @param path        Waypoints, n_joints x n_waypoints.
 * @param cumulative  cumulative_length(path).
 * @param n_samples   Number of samples.
 * @param control     n_ctrl x n_joints control polygon; the pinned rows are read, the free rows written.
 * @return false, leaving `control` unchanged, if `n_samples` samples are too few to determine the free
 *         control points.
 */
inline bool fit_free_control_points(const Bspline& bspline, const Matrix& path, const Array& cumulative, u32 n_samples, Matrix& control) {
  const u32 n_ctrl   = bspline.n_ctrl;
  const u32 n_joints = bspline.n_joints;
  const u32 n_free   = n_ctrl - 6;
  if (n_samples < std::max(n_free, (u32) 2))
    return false;

  const Bspline samples(n_ctrl, n_samples, bspline.degree, n_joints);
  Matrix        normal(n_free, n_free); // A^T A
  Matrix        rhs(n_free, n_joints);  // A^T r_j, one column per joint
  Array         target(n_joints);
  u32           segment = 0;
  for (u32 sample = 0; sample < n_samples; sample++) {
    const real distance = cumulative.back() * (real) sample / (real) (n_samples - 1);
    point_at_distance_along_path(path, cumulative, distance, segment, target);
    for (u32 joint = 0; joint < n_joints; joint++) {
      real residual = target[joint];
      for (u32 pinned = 0; pinned < 3; pinned++) {
        const u32 last = n_ctrl - 1 - pinned;
        residual -= samples.basis_p(pinned, sample) * control(pinned, joint) + samples.basis_p(last, sample) * control(last, joint);
      }
      for (u32 free_row = 0; free_row < n_free; free_row++)
        rhs(free_row, joint) += samples.basis_p(3 + free_row, sample) * residual;
    }
    for (u32 free_row = 0; free_row < n_free; free_row++)
      for (u32 free_col = 0; free_col < n_free; free_col++)
        normal(free_row, free_col) += samples.basis_p(3 + free_row, sample) * samples.basis_p(3 + free_col, sample);
  }

  if (!cholesky_factor(normal))
    return false;
  cholesky_solve(normal, rhs);
  for (u32 joint = 0; joint < n_joints; joint++)
    for (u32 free_row = 0; free_row < n_free; free_row++)
      control(3 + free_row, joint) = rhs(free_row, joint);
  return true;
}

/**
 * @brief Shortest duration for which a B-spline's velocity and acceleration control points are within limits.
 *
 * Velocity is itself a B-spline whose control points are the degree-scaled first differences of the
 * control polygon, divided by the duration T; acceleration likewise uses second differences over T^2.
 * Returns the smallest T, at least `min_duration`, that brings every one of those control points within
 * its joint's limit. By the convex-hull property this is conservative for the whole curve. Knots are
 * clamped uniform, as Bspline::compute_basis builds them. A limit of 0 or less is ignored.
 *
 * @param control  n_ctrl x n_joints control polygon.
 */
inline real duration_within_limits(const Matrix& control, u32 degree, const Manipulator& manip, real min_duration) {
  const u32  n_ctrl   = control.rows;
  const u32  n_joints = control.cols;
  const auto knot     = [&](u32 index) -> real {
    if (index <= degree)
      return 0;
    if (index >= n_ctrl)
      return 1;
    return (real) (index - degree) / (real) (n_ctrl - degree);
  };
  real  duration = min_duration;
  Array velocity_ctrl(n_ctrl - 1); // velocity control points times T
  for (u32 joint = 0; joint < n_joints; joint++) {
    const real velocity_max     = manip.velocity_max[joint];
    const real acceleration_max = manip.acceleration_max[joint];
    for (u32 ctrl = 0; ctrl + 1 < n_ctrl; ctrl++) {
      velocity_ctrl[ctrl] = (real) degree * (control(ctrl + 1, joint) - control(ctrl, joint)) / (knot(ctrl + degree + 1) - knot(ctrl + 1));
      if (velocity_max > 0)
        duration = std::max(duration, std::abs(velocity_ctrl[ctrl]) / velocity_max);
    }
    for (u32 ctrl = 0; ctrl + 2 < n_ctrl; ctrl++) {
      // acceleration control point times T^2
      const real acceleration_ctrl = (real) (degree - 1) * (velocity_ctrl[ctrl + 1] - velocity_ctrl[ctrl]) / (knot(ctrl + degree + 1) - knot(ctrl + 2));
      if (acceleration_max > 0)
        duration = std::max(duration, std::sqrt(std::abs(acceleration_ctrl) / acceleration_max));
    }
  }
  return duration;
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
  const Bspline& bspline  = opt->bspline;
  const u32      n_ctrl   = bspline.n_ctrl;
  const u32      n_joints = bspline.n_joints;
  Assert(n_ctrl > 6);
  Assert(path.rows == n_joints && path.cols >= 2);
  for (u32 joint = 0; joint < n_joints; joint++) {
    Assert(std::isfinite(opt->task(joint, 0)) && std::isfinite(opt->task(joint, 3)));
    Assert(opt->task(joint, 1) == 0 && opt->task(joint, 2) == 0 && opt->task(joint, 4) == 0 && opt->task(joint, 5) == 0);
  }
  const u32 default_samples = std::max((u32) 100, 4 * (n_ctrl - bspline.degree));
  if (n_samples == 0)
    n_samples = default_samples;

  // 1. Measure the path.
  const Array cumulative = cumulative_length(path);
  if (!(cumulative.back() > 0))
    return guess_straight_line(opt);

  // 2. Pin the first and last three control points to the task's start and goal, as
  //    Bspline::compute_control sets them for a stop-to-stop task.
  Matrix control(n_ctrl, n_joints);
  for (u32 joint = 0; joint < n_joints; joint++) {
    for (u32 pinned = 0; pinned < 3; pinned++) {
      control(pinned, joint)              = opt->task(joint, 0);
      control(n_ctrl - 1 - pinned, joint) = opt->task(joint, 3);
    }
  }

  // 3. Fit the free control points between them, falling back to the default sample count.
  if (!fit_free_control_points(bspline, path, cumulative, n_samples, control)) {
    printf("guess_from_path: n_samples %u too few for %u control points; using %u.\n", n_samples, n_ctrl, default_samples);
    const bool fitted = fit_free_control_points(bspline, path, cumulative, default_samples, control);
    Assert(fitted);
    (void) fitted;
  }

  // 4. Pick the duration.
  const real duration = duration_within_limits(control, bspline.degree, opt->manip, opt->min_duration);

  // 5. Pack the free control points (joint-major) and the duration, as guess_straight_line() does.
  Array x(bspline.x_len(opt->task));
  Assert(x.size == n_joints * (n_ctrl - 6) + 1);
  u32 index = 0;
  for (u32 joint = 0; joint < n_joints; joint++)
    for (u32 ctrl = 3; ctrl < n_ctrl - 3; ctrl++)
      x[index++] = control(ctrl, joint);
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
    Assert(!std::isnan(r));
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
