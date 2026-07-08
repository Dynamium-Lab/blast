#define CATCH_CONFIG_MAIN

#include <blast>
#include "catch2/catch.hpp"

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - default constructor zeroes all elements", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m;
  for (int i = 0; i < 9; i++)
    CHECK(m[i] == 0.0);
}

TEST_CASE("Mat3 - 9-element constructor stores values in column-major order", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  // Column 0
  CHECK(m(0, 0) == 1.0);
  CHECK(m(1, 0) == 2.0);
  CHECK(m(2, 0) == 3.0);
  // Column 1
  CHECK(m(0, 1) == 4.0);
  CHECK(m(1, 1) == 5.0);
  CHECK(m(2, 1) == 6.0);
  // Column 2
  CHECK(m(0, 2) == 7.0);
  CHECK(m(1, 2) == 8.0);
  CHECK(m(2, 2) == 9.0);
}

TEST_CASE("Mat3 - copy constructor produces an independent copy", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 copy(m);
  CHECK(is_close(copy, m));
  copy[0] = 99.0;
  CHECK(m[0] == 1.0);
}

// ---------------------------------------------------------------------------
// Element access
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - operator()(row, col) provides mutable access", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m;
  m(1, 2) = 42.0;
  CHECK(m(1, 2) == 42.0);
}

TEST_CASE("Mat3 - operator()(row, col) provides const access", "[Math][Mat3]") {
  using namespace blast;
  const Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  CHECK(m(0, 0) == 1.0);
  CHECK(m(2, 2) == 9.0);
}

TEST_CASE("Mat3 - operator[] provides flat index access matching constructor order", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  for (int i = 0; i < 9; i++)
    CHECK(m[i] == static_cast<real>(i + 1));
}

// ---------------------------------------------------------------------------
// Equality
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - is_close() returns true when all elements are within epsilon", "[Math][Mat3]") {
  using namespace blast;
  Mat3 a(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 b(1 + 1e-6, 2 + 1e-6, 3 + 1e-6, 4 + 1e-6, 5 + 1e-6, 6 + 1e-6, 7 + 1e-6, 8 + 1e-6, 9 + 1e-6);
  CHECK(is_close(a, b));
}

TEST_CASE("Mat3 - is_close() returns false when any element difference exceeds epsilon", "[Math][Mat3]") {
  using namespace blast;
  Mat3 a(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 b(1 + 2e-5, 2, 3, 4, 5, 6, 7, 8, 9);
  CHECK_FALSE(is_close(a, b));
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - operator+ adds element-wise", "[Math][Mat3]") {
  using namespace blast;
  Mat3 a(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 b(9, 8, 7, 6, 5, 4, 3, 2, 1);
  Mat3 expected(10, 10, 10, 10, 10, 10, 10, 10, 10);
  CHECK(is_close(a + b, expected));
}

TEST_CASE("Mat3 - operator- subtracts element-wise", "[Math][Mat3]") {
  using namespace blast;
  Mat3 a(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 b(9, 8, 7, 6, 5, 4, 3, 2, 1);
  Mat3 expected(-8, -6, -4, -2, 0, 2, 4, 6, 8);
  CHECK(is_close(a - b, expected));
}

TEST_CASE("Mat3 - operator+= adds element-wise in place", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 n(1, 1, 1, 1, 1, 1, 1, 1, 1);
  m += n;
  Mat3 expected(2, 3, 4, 5, 6, 7, 8, 9, 10);
  CHECK(is_close(m, expected));
}

TEST_CASE("Mat3 - left scalar operator* scales all elements", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  real c = 2.0;
  Mat3 expected(2, 4, 6, 8, 10, 12, 14, 16, 18);
  CHECK(is_close(c * m, expected));
}

// ---------------------------------------------------------------------------
// Mat3 x Mat3 multiplication
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - operator* produces correct matrix product", "[Math][Mat3]") {
  using namespace blast;
  Mat3 a(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 b(9, 8, 7, 6, 5, 4, 3, 2, 1);
  // Computed by hand: data layout is column-major
  Mat3 expected(90, 114, 138, 54, 69, 84, 18, 24, 30);
  CHECK(is_close(a * b, expected));
}

TEST_CASE("Mat3 - multiplying by identity leaves matrix unchanged", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(2, 3, 1, 5, 4, 6, 8, 7, 9);
  CHECK(is_close(m * eye(), m));
}

TEST_CASE("Mat3 - operator*= with identity leaves matrix unchanged", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 orig(m);
  m *= eye();
  CHECK(is_close(m, orig));
}

// ---------------------------------------------------------------------------
// Mat3 x Vec3 multiplication
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - operator*(Mat3, Vec3) produces correct result", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Vec3 v(1, 2, 3);
  Vec3 result = m * v;
  CHECK(result == Vec3(30, 36, 42));
}

// ---------------------------------------------------------------------------
// Transpose
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - transpose() swaps rows and columns", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 t = transpose(m);
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < 3; c++)
      CHECK(t(r, c) == m(c, r));
}

TEST_CASE("Mat3 - transpose_inplace() produces the same result as transpose()", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3 expected = transpose(m);
  transpose_inplace(m);
  CHECK(is_close(m, expected));
}

TEST_CASE("Mat3 - double transpose returns the original matrix", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  CHECK(is_close(transpose(transpose(m)), m));
}

// ---------------------------------------------------------------------------
// Identity
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - eye() returns the 3x3 identity matrix", "[Math][Mat3]") {
  using namespace blast;
  Mat3 I = eye();
  CHECK(I(0, 0) == 1.0);
  CHECK(I(1, 1) == 1.0);
  CHECK(I(2, 2) == 1.0);
  CHECK(I(0, 1) == 0.0);
  CHECK(I(0, 2) == 0.0);
  CHECK(I(1, 0) == 0.0);
  CHECK(I(1, 2) == 0.0);
  CHECK(I(2, 0) == 0.0);
  CHECK(I(2, 1) == 0.0);
}

// ---------------------------------------------------------------------------
// Column access
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - col_copy() returns the correct column as a Vec3", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  CHECK(m.col_copy(0) == Vec3(1, 2, 3));
  CHECK(m.col_copy(1) == Vec3(4, 5, 6));
  CHECK(m.col_copy(2) == Vec3(7, 8, 9));
}

TEST_CASE("Mat3 - col() alias allows element write-through to matrix data", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  m.col(0)[1] = 99.0;
  CHECK(m(1, 0) == 99.0);
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

TEST_CASE("Mat3 - zero() sets all elements to 0.0", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  zero(m);
  CHECK(is_small(m));
}

TEST_CASE("Mat3 - zero() returns a reference to the same object", "[Math][Mat3]") {
  using namespace blast;
  Mat3  m(1, 2, 3, 4, 5, 6, 7, 8, 9);
  Mat3& ref = zero(m);
  CHECK(&ref == &m);
}

TEST_CASE("Mat3 - constant() sets all elements to the given value", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m;
  constant(m, 5.0);
  for (int i = 0; i < 9; i++)
    CHECK(m[i] == 5.0);
}

TEST_CASE("Mat3 - constant() returns a reference to the same object", "[Math][Mat3]") {
  using namespace blast;
  Mat3  m;
  Mat3& ref = constant(m, 3.0);
  CHECK(&ref == &m);
}

TEST_CASE("Mat3 - is_small() returns true when all elements are within epsilon", "[Math][Mat3]") {
  using namespace blast;
  Mat3 m(1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6);
  CHECK(is_small(m));
}

TEST_CASE("Mat3 - is_small() returns false when any element exceeds epsilon", "[Math][Mat3]") {
  using namespace blast;
  CHECK_FALSE(is_small(Mat3(2e-5, 0, 0, 0, 0, 0, 0, 0, 0)));
  CHECK_FALSE(is_small(Mat3(0, 0, 0, 0, 2e-5, 0, 0, 0, 0)));
  CHECK_FALSE(is_small(Mat3(0, 0, 0, 0, 0, 0, 0, 0, 2e-5)));
}
