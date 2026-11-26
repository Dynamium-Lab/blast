#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include "test_helper.hpp"

using namespace blast;

// guess_from_path: least-squares fit of the B-spline curve to a joint-space waypoint path. Every
// expectation below is computed from the test data, not from the function under test.

namespace {

#if BLAST_USE_DOUBLES
constexpr double tolerance = 1e-9;
#else
constexpr double tolerance = 1e-3;
#endif

const Array start = {0.0, 0.4, -0.2, -0.4, 0.3, 0.0};
const Array goal  = {1.2, -0.3, 0.5, 0.2, -0.4, 0.8};

// The point a fraction `fraction` of the way along a polyline, where distance is the sum of Euclidean
// joint-space segment lengths (the same target guess_from_path fits to).
Array point_along(const std::vector<Array>& path, real fraction) {
  std::vector<real> cumulative(path.size(), 0);
  for (size_t i = 1; i < path.size(); i++) {
    real squared = 0;
    for (u32 j = 0; j < path[i].size; j++)
      squared += (path[i][j] - path[i - 1][j]) * (path[i][j] - path[i - 1][j]);
    cumulative[i] = cumulative[i - 1] + std::sqrt(squared);
  }
  const real distance = fraction * cumulative.back();
  size_t     segment  = 0;
  while (segment + 2 < path.size() && cumulative[segment + 1] < distance)
    segment++;
  const real span = cumulative[segment + 1] - cumulative[segment];
  const real t    = span > 0 ? std::clamp((distance - cumulative[segment]) / span, (real) 0, (real) 1) : (real) 0;
  Array      p(path[0].size);
  for (u32 j = 0; j < p.size; j++)
    p[j] = path[segment][j] + t * (path[segment + 1][j] - path[segment][j]);
  return p;
}

// Distance from a point to a polyline (closest point on any segment).
real distance_to_path(const std::vector<Array>& path, const Array& q) {
  real best = INF_REAL;
  for (size_t i = 0; i + 1 < path.size(); i++) {
    real dot = 0, squared = 0;
    for (u32 j = 0; j < q.size; j++) {
      dot += (q[j] - path[i][j]) * (path[i + 1][j] - path[i][j]);
      squared += (path[i + 1][j] - path[i][j]) * (path[i + 1][j] - path[i][j]);
    }
    const real t = squared > 0 ? std::clamp(dot / squared, (real) 0, (real) 1) : (real) 0;
    real       d = 0;
    for (u32 j = 0; j < q.size; j++) {
      const real e = q[j] - (path[i][j] + t * (path[i + 1][j] - path[i][j]));
      d += e * e;
    }
    best = std::min(best, std::sqrt(d));
  }
  return best;
}

// A path with a sharp corner: first the first three joints move, then the other three.
std::vector<Array> corner_path() {
  Array corner = start;
  for (u32 j = 0; j < 3; j++)
    corner[j] = goal[j];
  return {start, corner, goal};
}

// The same turn placed near the goal (or, mirrored, near the start): the curve bends hardest next to
// the pinned control points, where the duration bound uses the end knots.
std::vector<Array> late_corner_path() {
  Array corner = goal;
  for (u32 j = 3; j < 6; j++)
    corner[j] = start[j] + (real) 0.9 * (goal[j] - start[j]);
  return {start, corner, goal};
}
std::vector<Array> early_corner_path() {
  Array corner = start;
  for (u32 j = 0; j < 3; j++)
    corner[j] = start[j] + (real) 0.1 * (goal[j] - start[j]);
  return {start, corner, goal};
}

// The n_joints x n_waypoints Matrix guess_from_path takes, one column per waypoint.
Matrix as_matrix(const std::vector<Array>& path) {
  Matrix m(path[0].size, (u32) path.size());
  for (u32 i = 0; i < m.cols; i++)
    for (u32 j = 0; j < m.rows; j++)
      m(j, i) = path[i][j];
  return m;
}

// Normal equations: at the least-squares optimum over `n_samples` samples, sum_k N_i(u_k) (C(u_k) - path(u_k))
// = 0 for every free control point i. C is evaluated by Bspline::compute_trajectory and path(u_k) by the
// test's own interpolation, so the fit and the optimizer must agree on knots, basis and pinned control points.
void check_least_squares_optimum(const Optimization& opt, const std::vector<Array>& path, const Array& x, u32 n_samples) {
  Bspline check(opt.bspline.n_ctrl, n_samples, opt.bspline.degree, opt.bspline.n_joints);
  check.compute_trajectory(x, opt.task);
  for (u32 i = 3; i < opt.bspline.n_ctrl - 3; i++) {
    for (u32 j = 0; j < opt.bspline.n_joints; j++) {
      real gradient = 0, scale = 0;
      for (u32 k = 0; k < n_samples; k++) {
        const Array target = point_along(path, (real) k / (real) (n_samples - 1));
        gradient += check.basis_p(i, k) * (check.traj.pos(j, k) - target[j]);
        scale += check.basis_p(i, k) * std::abs(target[j]);
      }
      INFO("control point " << i << ", joint " << j);
      CHECK(std::abs(gradient) <= tolerance * std::max(scale, (real) 1));
    }
  }
}

Array sample(const Bspline& bspline, u32 point) {
  Array q(bspline.n_joints);
  for (u32 j = 0; j < bspline.n_joints; j++)
    q[j] = bspline.traj.pos(j, point);
  return q;
}

} // namespace

TEST_CASE("a two-point path gives a straight curve from start to goal", "[Optimization]") {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  const Array  x = guess_from_path(&opt, as_matrix({start, goal}));
  REQUIRE(x.size == opt.bspline.x_len(opt.task));
  opt.bspline.compute_trajectory(x, opt.task);
  for (u32 point = 0; point < opt.bspline.n_points; point++) {
    INFO("point " << point);
    // Collinear with start -> goal: one fraction f explains every joint.
    const real f = (opt.bspline.traj.pos(0, point) - start[0]) / (goal[0] - start[0]);
    for (u32 j = 0; j < start.size; j++)
      CHECK(opt.bspline.traj.pos(j, point) == Approx(start[j] + f * (goal[j] - start[j])).margin(tolerance));
  }
  for (u32 j = 0; j < start.size; j++) {
    CHECK(opt.bspline.traj.pos(j, 0) == Approx(start[j]).margin(tolerance));
    CHECK(opt.bspline.traj.pos(j, opt.bspline.n_points - 1) == Approx(goal[j]).margin(tolerance));
  }
}

TEST_CASE("the fit is the least-squares optimum in this optimizer's B-spline", "[Optimization]") {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  check_least_squares_optimum(opt, corner_path(), guess_from_path(&opt, as_matrix(corner_path()), 100), 100);
}

TEST_CASE("an explicit low sample count is used as given", "[Optimization]") {
  // 12 control points, degree 5: 10 samples still determine the 6 free control points.
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  check_least_squares_optimum(opt, corner_path(), guess_from_path(&opt, as_matrix(corner_path()), 10), 10);
}

TEST_CASE("a sample count too low for the fit falls back to the default", "[Optimization]") {
  // 12 control points, degree 5: 7 samples are more than the 6 free control points but still leave
  // the normal equations singular, so this exercises the pivot check, not only the count check.
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  const Array  low      = guess_from_path(&opt, as_matrix(corner_path()), 7);
  const Array  fallback = guess_from_path(&opt, as_matrix(corner_path()));
  REQUIRE(low.size == fallback.size);
  for (u32 i = 0; i < low.size; i++)
    CHECK(low[i] == fallback[i]);
}

TEST_CASE("the fitted curve follows a path with a corner much closer than the straight line", "[Optimization]") {
  Optimization             opt(make_UR5e(), Task::stop_to_stop(start, goal));
  const std::vector<Array> path          = corner_path();
  const auto               mean_distance = [&](const Array& x) {
    opt.bspline.compute_trajectory(x, opt.task);
    real sum = 0;
    for (u32 point = 0; point < opt.bspline.n_points; point++)
      sum += distance_to_path(path, sample(opt.bspline, point));
    return sum / (real) opt.bspline.n_points;
  };
  const real fitted   = mean_distance(guess_from_path(&opt, as_matrix(path)));
  const real straight = mean_distance(guess_straight_line(&opt));
  INFO("mean distance to the path: fitted " << fitted << ", straight line " << straight);
  CHECK(fitted < (real) 0.05 * straight);
}

TEST_CASE("the duration keeps velocity and acceleration within limits over the whole curve", "[Optimization]") {
  const Manipulator robot = make_UR5e();
  Optimization      opt(robot, Task::stop_to_stop(start, goal));
  for (const auto& [name, path]: {std::pair{"corner", corner_path()}, std::pair{"early corner", early_corner_path()},
                                  std::pair{"late corner", late_corner_path()}}) {
    INFO(name << " path");
    const Array x = guess_from_path(&opt, as_matrix(path));
    CHECK(x.back() >= opt.min_duration);
    // Rendered far denser than the optimizer samples: the bound holds on the continuum, not only at samples.
    Bspline dense(opt.bspline.n_ctrl, 4001, opt.bspline.degree, opt.bspline.n_joints);
    dense.compute_trajectory(x, opt.task);
    real highest_ratio = 0; // how close the guess comes to a limit, so the bound is not arbitrarily loose
    for (u32 j = 0; j < start.size; j++) {
      real peak_velocity = 0, peak_acceleration = 0;
      for (u32 point = 0; point < dense.n_points; point++) {
        peak_velocity     = std::max(peak_velocity, (real) std::abs(dense.traj.vel(j, point)));
        peak_acceleration = std::max(peak_acceleration, (real) std::abs(dense.traj.acc(j, point)));
      }
      INFO("joint " << j);
      CHECK(peak_velocity <= robot.velocity_max[j] * (1 + tolerance));
      CHECK(peak_acceleration <= robot.acceleration_max[j] * (1 + tolerance));
      highest_ratio = std::max({highest_ratio, peak_velocity / robot.velocity_max[j], peak_acceleration / robot.acceleration_max[j]});
    }
    INFO("highest ratio to a limit " << highest_ratio);
    CHECK(highest_ratio >= (real) 0.3); // the convex-hull bound is conservative: about 0.45 on these paths
  }
}

TEST_CASE("the default sample count scales with the number of segments", "[Optimization]") {
  // 40 control points, degree 5: 35 segments, so the default is max(100, 4 * 35) = 140 samples.
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  opt.set_bspline(Bspline(40, 150, 5, 6));
  const Array x = guess_from_path(&opt, as_matrix(corner_path()));
  REQUIRE(x.size == opt.bspline.x_len(opt.task));
  check_least_squares_optimum(opt, corner_path(), x, 140);
  opt.bspline.compute_trajectory(x, opt.task);
  for (u32 j = 0; j < start.size; j++) {
    CHECK(opt.bspline.traj.pos(j, 0) == Approx(start[j]).margin(tolerance));
    CHECK(opt.bspline.traj.pos(j, opt.bspline.n_points - 1) == Approx(goal[j]).margin(tolerance));
  }
}

#if !(defined(BLAST_DEBUG) || defined(_DEBUG)) // a debug build stops on the Assert instead, by design
TEST_CASE("a path with fewer than two waypoints gives the straight-line guess", "[Optimization]") {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  const Array  straight = guess_straight_line(&opt);
  for (u32 n_waypoints: {0u, 1u}) {
    INFO(n_waypoints << " waypoints");
    const Array fitted = guess_from_path(&opt, Matrix(start.size, n_waypoints));
    REQUIRE(fitted.size == straight.size);
    for (u32 i = 0; i < fitted.size; i++)
      CHECK(fitted[i] == straight[i]);
  }
}
#endif

TEST_CASE("a path of zero length gives the straight-line guess", "[Optimization]") {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, start));
  const Array  fitted   = guess_from_path(&opt, as_matrix({start, start}));
  const Array  straight = guess_straight_line(&opt);
  REQUIRE(fitted.size == straight.size);
  for (u32 i = 0; i < fitted.size; i++)
    CHECK(fitted[i] == straight[i]);
}

// cholesky_factor and cholesky_solve on their own, against matrices with a known answer.

namespace {
// A symmetric positive definite matrix: B B^T + I, B a fixed 4 x 4 matrix.
Matrix spd_matrix() {
  const real b[4][4] = {{2, -1, 0, 3}, {1, 4, -2, 0}, {0, 1, 3, -1}, {-2, 0, 1, 1}};
  Matrix     matrix(4, 4);
  for (u32 row = 0; row < 4; row++)
    for (u32 col = 0; col < 4; col++) {
      real value = row == col ? (real) 1 : (real) 0;
      for (u32 inner = 0; inner < 4; inner++)
        value += b[row][inner] * b[col][inner];
      matrix(row, col) = value;
    }
  return matrix;
}
} // namespace

TEST_CASE("cholesky_factor gives L with L L^T equal to the matrix", "[Optimization]") {
  const Matrix original = spd_matrix();
  Matrix       factor   = original;
  REQUIRE(cholesky_factor(factor));
  for (u32 row = 0; row < 4; row++) {
    for (u32 col = 0; col < 4; col++) {
      // L is the lower triangle of `factor`; its upper triangle is left as it was and must be ignored.
      real product = 0;
      for (u32 inner = 0; inner <= std::min(row, col); inner++)
        product += factor(row, inner) * factor(col, inner);
      INFO("row " << row << ", col " << col);
      CHECK(product == Approx(original(row, col)).margin(tolerance * 100));
    }
  }
}

TEST_CASE("cholesky_solve recovers a known solution, one per column", "[Optimization]") {
  const Matrix matrix = spd_matrix();
  Matrix       expected(4, 2);
  const real   values[2][4] = {{1, -2, (real) 0.5, 3}, {0, 1, -1, 2}};
  for (u32 column = 0; column < 2; column++)
    for (u32 row = 0; row < 4; row++)
      expected(row, column) = values[column][row];
  Matrix columns(4, 2); // right-hand sides b = A x
  for (u32 column = 0; column < 2; column++)
    for (u32 row = 0; row < 4; row++)
      for (u32 inner = 0; inner < 4; inner++)
        columns(row, column) += matrix(row, inner) * expected(inner, column);
  Matrix factor = matrix;
  REQUIRE(cholesky_factor(factor));
  cholesky_solve(factor, columns);
  for (u32 column = 0; column < 2; column++)
    for (u32 row = 0; row < 4; row++) {
      INFO("row " << row << ", column " << column);
      CHECK(columns(row, column) == Approx(expected(row, column)).margin(tolerance * 100));
    }
}

TEST_CASE("cholesky_factor rejects a singular matrix", "[Optimization]") {
  // v v^T has rank 1: positive semidefinite but singular, as when samples are too few.
  const real vector[3] = {1, 2, -1};
  Matrix     matrix(3, 3);
  for (u32 row = 0; row < 3; row++)
    for (u32 col = 0; col < 3; col++)
      matrix(row, col) = vector[row] * vector[col];
  CHECK_FALSE(cholesky_factor(matrix));
}
