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

Array sample(const Bspline& bspline, u32 point) {
  Array q(bspline.n_joints);
  for (u32 j = 0; j < bspline.n_joints; j++)
    q[j] = bspline.traj.pos(j, point);
  return q;
}

} // namespace

TEST_CASE("a two-point path gives a straight curve from start to goal", "[Optimization]") {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  const Array  x = guess_from_path(&opt, {start, goal});
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
  // Normal equations: at the optimum, sum_k N_i(u_k) (C(u_k) - path(u_k)) = 0 for every free control
  // point i. C is evaluated by Bspline::compute_trajectory and path(u_k) by the test's own arc-length
  // interpolation, so the fit and the optimizer must agree on knots, basis and pinned control points.
  Optimization             opt(make_UR5e(), Task::stop_to_stop(start, goal));
  const std::vector<Array> path      = corner_path();
  const u32                n_samples = 100;
  const Array              x         = guess_from_path(&opt, path, n_samples);
  Bspline                  check(opt.bspline.n_ctrl, n_samples, opt.bspline.degree, opt.bspline.n_joints);
  check.compute_trajectory(x, opt.task);
  for (u32 i = 3; i < opt.bspline.n_ctrl - 3; i++) {
    for (u32 j = 0; j < start.size; j++) {
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
  const real fitted   = mean_distance(guess_from_path(&opt, path));
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
    const Array x = guess_from_path(&opt, path);
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

TEST_CASE("the default sample count works for any number of control points", "[Optimization]") {
  // 30 control points: the automatic count, max(100, 4 * n_ctrl) = 120, keeps the fit well posed.
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, goal));
  opt.set_bspline(Bspline(30, 150, 5, 6));
  const Array x = guess_from_path(&opt, corner_path());
  REQUIRE(x.size == opt.bspline.x_len(opt.task));
  for (u32 i = 0; i < x.size; i++)
    CHECK(std::isfinite(x[i]));
  opt.bspline.compute_trajectory(x, opt.task);
  for (u32 j = 0; j < start.size; j++) {
    CHECK(opt.bspline.traj.pos(j, 0) == Approx(start[j]).margin(tolerance));
    CHECK(opt.bspline.traj.pos(j, opt.bspline.n_points - 1) == Approx(goal[j]).margin(tolerance));
  }
}

TEST_CASE("a path of zero length gives the straight-line guess", "[Optimization]") {
  Optimization opt(make_UR5e(), Task::stop_to_stop(start, start));
  const Array  fitted   = guess_from_path(&opt, {start, start});
  const Array  straight = guess_straight_line(&opt);
  REQUIRE(fitted.size == straight.size);
  for (u32 i = 0; i < fitted.size; i++)
    CHECK(fitted[i] == straight[i]);
}
