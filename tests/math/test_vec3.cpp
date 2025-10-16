//
// Created by nikos on 2025-10-15.
//
#define CATCH_CONFIG_MAIN

#include <blast>
#include "catch2/catch.hpp"

TEST_CASE("Constructor Test", "[Math]") {
  using namespace blast;
  real x = 1.0;
  real y = 2.0;
  real z = 3.0;
  Vec3 v(x, y, z);

  CHECK(v.x == x);
  CHECK(v.y == y);
  CHECK(v.z == z);

  Vec3 v2;
  CHECK(v2.x == 0);
  CHECK(v2.y == 0);
  CHECK(v2.z == 0);
}

TEST_CASE("Operator[] Test", "[Math]") { // NOTE: This should be removed...
  using namespace blast;
  real x = 1.0;
  real y = 2.0;
  real z = 3.0;
  Vec3 v(x, y, z);

  // todo: v[3], v[4], ... is unsafe and not checked. Should we remove this operator ?
  CHECK(v[0] == v.x);
  CHECK(v[1] == v.y);
  CHECK(v[2] == v.z);
}

TEST_CASE("Linear Algebra Test", "[Math]") {
  using namespace blast;
  real x = 1.0;
  real y = 2.0;
  real z = 3.0;
  Vec3 v(x, y, z);

  auto v_cr = cross(v, v);
  CHECK(v_cr == Vec3(0, 0, 0));

  auto v_dot = dot(v, v);
  CHECK(v_dot == 14.0);

  auto v_norm = norm(v);
  CHECK(std::abs(v_norm - 3.7416573868) < 1e-6);
}

TEST_CASE("Operators Test", "[Math]") {
  using namespace blast;
  real x_1 = 1.0;
  real y_1 = 2.0;
  real z_1 = 3.0;
  Vec3 v1(x_1, y_1, z_1);

  real x_2 = 2.0;
  real y_2 = 4.0;
  real z_2 = 6.0;
  Vec3 v2(x_2, y_2, z_2);

  real coeff = 2.0;

  CHECK(-v1 == Vec3(-1.0, -2.0, -3.0));
  CHECK(v2 - v1 == v1);
  CHECK(v1 + v2 == Vec3(3.0, 6.0, 9.0));
  CHECK(coeff * v1 == v2);
  CHECK(v1 * coeff == v2);
  CHECK(v2 / coeff == v1);

  // += and -=
  Vec3 v1_p(x_1, y_1, z_1);
  v1_p += v1;
  CHECK(v1_p == v2);
  Vec3 v1_m(x_2, y_2, z_2);
  v1_m -= v1;
  CHECK(v1_m == v1);

  // *= and /=
  Vec3 v1_mult(x_1, y_1, z_1);
  v1_mult *= coeff;
  CHECK(v1_mult == v2);
  Vec3 v1_div(x_2, y_2, z_2);
  v1_div /= coeff;
  CHECK(v1_div == v1);
}

TEST_CASE("Zero Test", "[Math]") {
  using namespace blast;
  real x = 1.0;
  real y = 2.0;
  real z = 3.0;
  Vec3 v(x, y, z);

  zero(v);
  CHECK(v.x == 0.0);
  CHECK(v.y == 0.0);
  CHECK(v.z == 0.0);
}

TEST_CASE("is_small & is_close Test", "[Math]") {
  using namespace blast;
  real x = 1e-06;
  real y = 1e-06;
  real z = 1e-06;
  Vec3 v(x, y, z);
  Vec3 v_zero(0.0, 0.0, 0.0);

  CHECK(is_small(v) == true);
  CHECK(is_close(v, v_zero) == true);
  v.x = 1e-05;
  CHECK(is_small(v) == false);
  CHECK(is_close(v, v_zero) == false);
  v.x = 1e-06;
  v.y = 1e-05;
  CHECK(is_small(v) == false);
  CHECK(is_close(v, v_zero) == false);
  v.y = 1e-06;
  v.z = 1e-05;
  CHECK(is_small(v) == false);
  CHECK(is_close(v, v_zero) == false);
  v.z = 1e-06;
}

TEST_CASE("Constant Test", "[Math]") {
  using namespace blast;
  real coeff = 2.0;
  Vec3 v;
  constant(v, coeff);

  CHECK(v.x == 2.0);
  CHECK(v.y == 2.0);
  CHECK(v.z == 2.0);
}
