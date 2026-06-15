#define CATCH_CONFIG_MAIN

#include <blast>
#include "test_helper.hpp"
#include "catch2/catch.hpp"

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - default constructor zeroes all components", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v;
  CHECK(v.x == 0.0);
  CHECK(v.y == 0.0);
  CHECK(v.z == 0.0);
}

TEST_CASE("Vec3 - parameterized constructor stores x, y, z", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v(1.0, 2.0, 3.0);
  CHECK(v.x == 1.0);
  CHECK(v.y == 2.0);
  CHECK(v.z == 3.0);
}

TEST_CASE("Vec3 - copy constructor and assignment produce independent copies", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v(1.0, 2.0, 3.0);

  Vec3 v_copy(v);
  CHECK(v_copy.x == v.x);
  CHECK(v_copy.y == v.y);
  CHECK(v_copy.z == v.z);

  Vec3 v_assign = v;
  CHECK(v_assign.x == v.x);
  CHECK(v_assign.y == v.y);
  CHECK(v_assign.z == v.z);

  // Mutating the copy must not affect the original
  v_copy.x = 99.0;
  CHECK(v.x == 1.0);
}

// ---------------------------------------------------------------------------
// operator[]
// NOTE: unsafe — no bounds check in release builds.
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - operator[] maps indices 0,1,2 to x,y,z (mutable)", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v(1.0, 2.0, 3.0);
  CHECK(v[0] == v.x);
  CHECK(v[1] == v.y);
  CHECK(v[2] == v.z);
}

TEST_CASE("Vec3 - operator[] maps indices 0,1,2 to x,y,z (const)", "[Math][Vec3]") {
  using namespace blast;
  const Vec3 v(4.0, 5.0, 6.0);
  CHECK(v[0] == 4.0);
  CHECK(v[1] == 5.0);
  CHECK(v[2] == 6.0);
}

// ---------------------------------------------------------------------------
// Equality  (operator== delegates to is_close)
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - operator== treats vectors within epsilon as equal", "[Math][Vec3]") {
  using namespace blast;
  Vec3 a(1.0, 2.0, 3.0);
  Vec3 b(1.0 + 1e-7, 2.0, 3.0); // within default epsilon
  CHECK(a == b);
}

TEST_CASE("Vec3 - operator== treats vectors outside epsilon as not equal", "[Math][Vec3]") {
  using namespace blast;
  Vec3 a(1.0, 2.0, 3.0);
  Vec3 c(1.0 + 1e-4, 2.0, 3.0); // outside default epsilon
  CHECK_FALSE(a == c);
}

// ---------------------------------------------------------------------------
// Arithmetic operators
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - unary minus negates all components", "[Math][Vec3]") {
  using namespace blast;
  CHECK(-Vec3(1.0, 2.0, 3.0) == Vec3(-1.0, -2.0, -3.0));
  CHECK(-Vec3(0.0, 0.0, 0.0) == Vec3(0.0, 0.0, 0.0));
}

TEST_CASE("Vec3 - binary + and - are correct", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v1(1.0, 2.0, 3.0);
  Vec3 v2(2.0, 4.0, 6.0);

  CHECK(v1 + v2 == Vec3(3.0, 6.0, 9.0));
  CHECK(v2 - v1 == v1);
}

TEST_CASE("Vec3 - scalar multiplication and division are correct", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v1(1.0, 2.0, 3.0);
  Vec3 v2(2.0, 4.0, 6.0);
  real c = 2.0;

  CHECK(c * v1 == v2); // left scalar multiply
  CHECK(v1 * c == v2); // right scalar multiply
  CHECK(v2 / c == v1);
}

TEST_CASE("Vec3 - compound assignment operators are correct", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v1(1.0, 2.0, 3.0);
  Vec3 v2(2.0, 4.0, 6.0);
  real c = 2.0;

  Vec3 v_add(1.0, 2.0, 3.0);
  v_add += v1;
  CHECK(v_add == v2);

  Vec3 v_sub(2.0, 4.0, 6.0);
  v_sub -= v1;
  CHECK(v_sub == v1);

  Vec3 v_mul(1.0, 2.0, 3.0);
  v_mul *= c;
  CHECK(v_mul == v2);

  Vec3 v_div(2.0, 4.0, 6.0);
  v_div /= c;
  CHECK(v_div == v1);
}

// ---------------------------------------------------------------------------
// Cross product
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - cross product of a vector with itself is zero", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v(1.0, 2.0, 3.0);
  CHECK(cross(v, v) == Vec3(0.0, 0.0, 0.0));
}

TEST_CASE("Vec3 - cross product gives correct result for known inputs", "[Math][Vec3]") {
  using namespace blast;
  Vec3 a(1.0, 2.0, 3.0);
  Vec3 b(2.0, 3.0, 4.0);
  CHECK(cross(a, b) == Vec3(-1.0, 2.0, -1.0));
}

TEST_CASE("Vec3 - cross product is anticommutative", "[Math][Vec3]") {
  using namespace blast;
  Vec3 a(1.0, 2.0, 3.0);
  Vec3 b(2.0, 3.0, 4.0);
  CHECK(cross(a, b) == -cross(b, a));
}

TEST_CASE("Vec3 - cross product result is perpendicular to both inputs", "[Math][Vec3]") {
  using namespace blast;
  Vec3 a(1.0, 2.0, 3.0);
  Vec3 b(4.0, 5.0, 6.0);
  Vec3 c = cross(a, b);
  CHECK(std::abs(dot(c, a)) < blast::test::abs_tol);
  CHECK(std::abs(dot(c, b)) < blast::test::abs_tol);
}

// ---------------------------------------------------------------------------
// Dot product
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - dot product of a vector with itself equals squared norm", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v(1.0, 2.0, 3.0);
  CHECK(dot(v, v) == 14.0); // 1 + 4 + 9
}

TEST_CASE("Vec3 - dot product of orthogonal vectors is zero", "[Math][Vec3]") {
  using namespace blast;
  CHECK(dot(Vec3(1, 0, 0), Vec3(0, 1, 0)) == 0.0);
  CHECK(dot(Vec3(0, 1, 0), Vec3(0, 0, 1)) == 0.0);
  CHECK(dot(Vec3(1, 0, 0), Vec3(0, 0, 1)) == 0.0);
}

TEST_CASE("Vec3 - dot product gives correct result for known inputs", "[Math][Vec3]") {
  using namespace blast;
  // [1,2,3] · [4,5,6] = 4 + 10 + 18 = 32
  CHECK(dot(Vec3(1, 2, 3), Vec3(4, 5, 6)) == 32.0);
}

TEST_CASE("Vec3 - dot product is commutative", "[Math][Vec3]") {
  using namespace blast;
  Vec3 a(1.0, 2.0, 3.0);
  Vec3 b(4.0, 5.0, 6.0);
  CHECK(dot(a, b) == dot(b, a));
}

// ---------------------------------------------------------------------------
// Norm
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - norm gives correct result for known input", "[Math][Vec3]") {
  using namespace blast;
  CHECK(std::abs(norm(Vec3(1.0, 2.0, 3.0)) - 3.7416573868) < 1e-6);
}

TEST_CASE("Vec3 - norm of zero vector is zero", "[Math][Vec3]") {
  using namespace blast;
  CHECK(norm(Vec3(0.0, 0.0, 0.0)) == 0.0);
}

TEST_CASE("Vec3 - norm of unit axis vectors is one", "[Math][Vec3]") {
  using namespace blast;
  CHECK(std::abs(norm(Vec3(1, 0, 0)) - 1.0) < blast::test::abs_tol);
  CHECK(std::abs(norm(Vec3(0, 1, 0)) - 1.0) < blast::test::abs_tol);
  CHECK(std::abs(norm(Vec3(0, 0, 1)) - 1.0) < blast::test::abs_tol);
}

// ---------------------------------------------------------------------------
// zero()
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - zero() sets all components to 0.0", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v(1.0, 2.0, 3.0);
  zero(v);
  CHECK(v.x == 0.0);
  CHECK(v.y == 0.0);
  CHECK(v.z == 0.0);
}

TEST_CASE("Vec3 - zero() returns a reference to the same object", "[Math][Vec3]") {
  using namespace blast;
  Vec3  v(1.0, 2.0, 3.0);
  Vec3& ref = zero(v);
  CHECK(&ref == &v);
}

// ---------------------------------------------------------------------------
// constant()
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - constant() sets all components to the given value", "[Math][Vec3]") {
  using namespace blast;
  Vec3 v;

  constant(v, 2.0);
  CHECK(v.x == 2.0);
  CHECK(v.y == 2.0);
  CHECK(v.z == 2.0);

  constant(v, 0.0);
  CHECK(v.x == 0.0);
  CHECK(v.y == 0.0);
  CHECK(v.z == 0.0);

  constant(v, -5.0);
  CHECK(v.x == -5.0);
  CHECK(v.y == -5.0);
  CHECK(v.z == -5.0);
}

TEST_CASE("Vec3 - constant() returns a reference to the same object", "[Math][Vec3]") {
  using namespace blast;
  Vec3  v;
  Vec3& ref = constant(v, 7.0);
  CHECK(&ref == &v);
}

// ---------------------------------------------------------------------------
// is_small()
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - is_small() returns true when all components are within epsilon", "[Math][Vec3]") {
  using namespace blast;
  CHECK(is_small(Vec3(1e-6, 1e-6, 1e-6)) == true);
}

TEST_CASE("Vec3 - is_small() returns false when any component exceeds epsilon", "[Math][Vec3]") {
  using namespace blast;
  CHECK(is_small(Vec3(1e-5, 1e-6, 1e-6)) == false); // x exceeds
  CHECK(is_small(Vec3(1e-6, 1e-5, 1e-6)) == false); // y exceeds
  CHECK(is_small(Vec3(1e-6, 1e-6, 1e-5)) == false); // z exceeds
}

// ---------------------------------------------------------------------------
// is_close()
// ---------------------------------------------------------------------------

TEST_CASE("Vec3 - is_close() returns true when vectors are within epsilon", "[Math][Vec3]") {
  using namespace blast;
  CHECK(is_close(Vec3(1e-6, 1e-6, 1e-6), Vec3(0.0, 0.0, 0.0)) == true);

  Vec3 p(1.0, 2.0, 3.0);
  Vec3 q(1.0 + 1e-7, 2.0 + 1e-7, 3.0 + 1e-7);
  CHECK(is_close(p, q) == true);
}

TEST_CASE("Vec3 - is_close() returns false when any component exceeds epsilon", "[Math][Vec3]") {
  using namespace blast;
  Vec3 z(0.0, 0.0, 0.0);
  CHECK(is_close(Vec3(1e-5, 1e-6, 1e-6), z) == false); // x exceeds
  CHECK(is_close(Vec3(1e-6, 1e-5, 1e-6), z) == false); // y exceeds
  CHECK(is_close(Vec3(1e-6, 1e-6, 1e-5), z) == false); // z exceeds

  Vec3 p(1.0, 2.0, 3.0);
  Vec3 r(1.0 + 1e-4, 2.0, 3.0);
  CHECK(is_close(p, r) == false);
}
