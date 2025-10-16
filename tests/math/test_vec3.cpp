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

TEST_CASE("Operators", "[Math]") {
  using namespace blast;
  real x_1 = 1.0;
  real y_1 = 2.0;
  real z_1 = 3.0;
  Vec3 v1(x_1, y_1, z_1);

  real x_2 = 2.0;
  real y_2 = 4.0;
  real z_2 = 6.0;
  Vec3 v2(x_2, y_2, z_2);

  CHECK(v2 - v1 == v1);
}
