#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;

TEST_CASE("compute_5order_trajectory - boundary conditions for single joint", "[Trajectory][Polynomial]") {
  Matrix task(1, 6);
  task(0, 0) = 0.0; // start pos
  task(0, 1) = 0.0; // start vel
  task(0, 2) = 0.0; // start acc
  task(0, 3) = 1.0; // end pos
  task(0, 4) = 0.0; // end vel
  task(0, 5) = 0.0; // end acc

  real T    = 2.0;
  auto traj = compute_5order_trajectory(T, task);

  constexpr real eps = 1e-6;
  const int      last = (int) traj.t.size - 1;
  CHECK(std::abs(traj.pos(0, 0) - 0.0) < eps);
  CHECK(std::abs(traj.vel(0, 0) - 0.0) < eps);
  CHECK(std::abs(traj.acc(0, 0) - 0.0) < eps);
  CHECK(std::abs(traj.pos(0, last) - 1.0) < eps);
  CHECK(std::abs(traj.vel(0, last) - 0.0) < eps);
  CHECK(std::abs(traj.acc(0, last) - 0.0) < eps);
}

TEST_CASE("compute_5order_trajectory - boundary conditions with nonzero start/end vel and acc", "[Trajectory][Polynomial]") {
  Matrix task(1, 6);
  task(0, 0) = -0.5; // start pos
  task(0, 1) = 0.2;  // start vel
  task(0, 2) = 0.1;  // start acc
  task(0, 3) = 0.5;  // end pos
  task(0, 4) = -0.1; // end vel
  task(0, 5) = 0.0;  // end acc

  real T    = 1.0;
  auto traj = compute_5order_trajectory(T, task);

  constexpr real eps = 1e-5;
  const int      last = (int) traj.t.size - 1;
  CHECK(std::abs(traj.pos(0, 0) - task(0, 0)) < eps);
  CHECK(std::abs(traj.vel(0, 0) - task(0, 1)) < eps);
  CHECK(std::abs(traj.acc(0, 0) - task(0, 2)) < eps);
  CHECK(std::abs(traj.pos(0, last) - task(0, 3)) < eps);
  CHECK(std::abs(traj.vel(0, last) - task(0, 4)) < eps);
  CHECK(std::abs(traj.acc(0, last) - task(0, 5)) < eps);
}

TEST_CASE("compute_5order_trajectory - output size matches expected point count", "[Trajectory][Polynomial]") {
  Matrix task(1, 6);
  task(0, 0) = 0.0; task(0, 1) = 0.0; task(0, 2) = 0.0;
  task(0, 3) = 1.0; task(0, 4) = 0.0; task(0, 5) = 0.0;

  real T       = 2.0;
  auto traj    = compute_5order_trajectory(T, task);
  u32  expected = (u32) std::ceil(T * 1000 + 1);
  CHECK(traj.t.size == expected);
  CHECK(traj.pos.rows == 1u);
  CHECK(traj.pos.cols == expected);
}

TEST_CASE("compute_5order_trajectory - time array starts at zero, is monotonic", "[Trajectory][Polynomial]") {
  Matrix task(1, 6);
  task(0, 0) = 0.0; task(0, 1) = 0.0; task(0, 2) = 0.0;
  task(0, 3) = 1.0; task(0, 4) = 0.0; task(0, 5) = 0.0;

  auto traj = compute_5order_trajectory(3.0, task);
  CHECK(traj.t[0] == Approx(0.0));
  for (u32 i = 1; i < traj.t.size; i++) {
    CHECK(traj.t[i] > traj.t[i - 1]);
  }
}

TEST_CASE("compute_5order_trajectory - trivial motion (start == end, all zero derivatives) stays constant", "[Trajectory][Polynomial]") {
  Matrix task(1, 6);
  task(0, 0) = 0.7; task(0, 1) = 0.0; task(0, 2) = 0.0;
  task(0, 3) = 0.7; task(0, 4) = 0.0; task(0, 5) = 0.0;

  auto traj = compute_5order_trajectory(1.0, task);
  for (u32 i = 0; i < traj.t.size; i++) {
    CHECK(std::abs(traj.pos(0, (int) i) - 0.7) < 1e-9);
    CHECK(std::abs(traj.vel(0, (int) i) - 0.0) < 1e-9);
    CHECK(std::abs(traj.acc(0, (int) i) - 0.0) < 1e-9);
  }
}

TEST_CASE("compute_5order_trajectory - multi-joint boundary conditions", "[Trajectory][Polynomial]") {
  const int n_joints = 6;
  Matrix    task(n_joints, 6);
  Array     p0 = {0.0, -1.57, 0.0, -1.57, 0.0, 0.0};
  Array     pf = {1.0, -1.0, 0.5, -1.0, 0.5, -0.5};

  for (int j = 0; j < n_joints; j++) {
    task(j, 0) = p0[j];
    task(j, 1) = 0.0;
    task(j, 2) = 0.0;
    task(j, 3) = pf[j];
    task(j, 4) = 0.0;
    task(j, 5) = 0.0;
  }

  auto      traj = compute_5order_trajectory(2.0, task);
  const int last = (int) traj.t.size - 1;
  for (int j = 0; j < n_joints; j++) {
    CHECK(std::abs(traj.pos(j, 0) - p0[j]) < 1e-6);
    CHECK(std::abs(traj.pos(j, last) - pf[j]) < 1e-6);
    CHECK(std::abs(traj.vel(j, 0) - 0.0) < 1e-6);
    CHECK(std::abs(traj.vel(j, last) - 0.0) < 1e-6);
  }
}
