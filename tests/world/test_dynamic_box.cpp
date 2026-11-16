#define CATCH_CONFIG_MAIN
#include <blast>
#include <catch2/catch.hpp>
#include <cmath>
#include "test_helper.hpp"

using namespace blast;

static DynamicBox two_keyframes(const Mat3& rotation_start, const Mat3& rotation_end) {
  DynamicBox box;
  box.n_points   = 2;
  box.start_time = 0;
  box.end_time   = 1;
  box.trajectory = {{{0, 0, 0}, {0.3, 0.2, 0.1}, rotation_start}, {{0, 0, 0}, {0.3, 0.2, 0.1}, rotation_end}};
  return box;
}

static void check_rotation(const Mat3& rotation) {
  CHECK(is_close(transpose(rotation) * rotation, eye(), test::approx_tol));
  const real determinant = rotation(0, 0) * (rotation(1, 1) * rotation(2, 2) - rotation(1, 2) * rotation(2, 1)) -
                           rotation(0, 1) * (rotation(1, 0) * rotation(2, 2) - rotation(1, 2) * rotation(2, 0)) +
                           rotation(0, 2) * (rotation(1, 0) * rotation(2, 1) - rotation(1, 1) * rotation(2, 0));
  CHECK(determinant == Approx(1).margin(test::approx_tol));
}

// Between keyframes the box must stay a rotated box: a linear blend of the two rotation matrices
// shrinks and shears it, so its bounding volume and its distance disagree (broadphase misses).
TEST_CASE("DynamicBox interpolates rotation along the shortest arc", "[World]") {
  const real pi = std::acos(-1.0);
  {
    const DynamicBox box = two_keyframes(eye(), rpy2rotation({0, 0, pi / 2}));
    check_rotation(box.lookup(0.5).rotation);
    CHECK(is_close(box.lookup(0.5).rotation, rpy2rotation({0, 0, pi / 4}), test::approx_tol)); // halfway: 45 deg about z
    CHECK(is_close(box.lookup(0).rotation, eye(), test::approx_tol));
    CHECK(is_close(box.lookup(1).rotation, rpy2rotation({0, 0, pi / 2}), test::approx_tol));
  }
  {
    const Mat3       start = rpy2rotation({0.3, -0.4, 1.1});
    const Mat3       end   = rpy2rotation({-0.2, 0.5, 1.9});
    const DynamicBox box   = two_keyframes(start, end);
    for (real fraction: {0.1, 0.25, 0.5, 0.75, 0.9})
      check_rotation(box.lookup(fraction).rotation);
    CHECK(is_close(box.lookup(1).rotation, end, test::approx_tol));
  }
}
