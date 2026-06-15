#define CATCH_CONFIG_MAIN

#include <blast>
#include "test_helper.hpp"
#include "catch2/catch.hpp"

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - default constructor produces empty matrix", "[Math][Matrix]") {
  using namespace blast;
  Matrix m;
  CHECK(m.rows == 0u);
  CHECK(m.cols == 0u);
  CHECK(m.size == 0u);
  CHECK(m.data == nullptr);
  CHECK(m.is_alias == false);
}

TEST_CASE("Matrix - size constructor allocates and zeroes all elements", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(3u, 2u);
  CHECK(m.rows == 3u);
  CHECK(m.cols == 2u);
  CHECK(m.size == 6u);
  for (u32 r = 0; r < m.rows; r++)
    for (u32 c = 0; c < m.cols; c++)
      CHECK(m(r, c) == 0.0);
}

TEST_CASE("Matrix - size+value constructor fills all elements", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 3u, 5.0);
  CHECK(m.rows == 2u);
  CHECK(m.cols == 3u);
  for (u32 r = 0; r < m.rows; r++)
    for (u32 c = 0; c < m.cols; c++)
      CHECK(m(r, c) == 5.0);
}

TEST_CASE("Matrix - copy constructor produces a deep independent copy", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u, 1.0);
  Matrix copy(m);
  CHECK(copy.rows == 2u);
  CHECK(copy.cols == 2u);
  CHECK(copy.is_alias == false);
  CHECK(copy(0, 0) == 1.0);
  copy(0, 0) = 99.0;
  CHECK(m(0, 0) == 1.0);
}

TEST_CASE("Matrix - move constructor transfers ownership", "[Math][Matrix]") {
  using namespace blast;
  Matrix a(2u, 2u, 3.0);
  Matrix b(std::move(a));
  CHECK(b.rows == 2u);
  CHECK(b.cols == 2u);
  CHECK(b(0, 0) == 3.0);
  CHECK(a.rows == 0u);
  CHECK(a.cols == 0u);
}

// ---------------------------------------------------------------------------
// Element access
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - operator()(row, col) provides mutable access", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(3u, 3u);
  m(1, 2) = 42.0;
  CHECK(m(1, 2) == 42.0);
}

TEST_CASE("Matrix - operator()(row, col) provides const access", "[Math][Matrix]") {
  using namespace blast;
  Matrix        m(2u, 2u, 7.0);
  const Matrix& cm = m;
  CHECK(cm(0, 0) == 7.0);
  CHECK(cm(1, 1) == 7.0);
}

TEST_CASE("Matrix - address() returns pointer to the correct element", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u);
  m(1, 0) = 5.0;
  CHECK(m.address(1, 0) == &m(1, 0));
  CHECK(*m.address(1, 0) == 5.0);
}

TEST_CASE("Matrix - uses column-major storage layout", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u);
  m(0, 0) = 1.0;
  m(1, 0) = 2.0;
  m(0, 1) = 3.0;
  m(1, 1) = 4.0;
  // Column-major: data[row + rows*col]
  CHECK(m.data[0] == 1.0); // m(0,0)
  CHECK(m.data[1] == 2.0); // m(1,0)
  CHECK(m.data[2] == 3.0); // m(0,1)
  CHECK(m.data[3] == 4.0); // m(1,1)
}

// ---------------------------------------------------------------------------
// Dimensions
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - rows, cols, size fields reflect matrix dimensions", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(3u, 4u);
  CHECK(m.rows == 3u);
  CHECK(m.cols == 4u);
  CHECK(m.size == 12u);
}

// ---------------------------------------------------------------------------
// Equality
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - operator== returns true when all elements are within epsilon", "[Math][Matrix]") {
  using namespace blast;
  Matrix a(2u, 2u, 1.0);
  Matrix b(2u, 2u, 1.0 + 1e-6);
  CHECK(a == b);
}

TEST_CASE("Matrix - operator== returns false when any element difference exceeds epsilon", "[Math][Matrix]") {
  using namespace blast;
  Matrix a(2u, 2u, 1.0);
  Matrix b(2u, 2u, 1.0 + 2e-5);
  CHECK_FALSE(a == b);
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - unary minus negates all elements", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u, 3.0);
  Matrix neg = -m;
  CHECK(neg(0, 0) == -3.0);
  CHECK(neg(1, 1) == -3.0);
  CHECK(m(0, 0) == 3.0);
}

TEST_CASE("Matrix - operator+ adds element-wise", "[Math][Matrix]") {
  using namespace blast;
  Matrix a(2u, 2u, 1.0);
  Matrix b(2u, 2u, 2.0);
  Matrix result = a + b;
  CHECK(result(0, 0) == 3.0);
  CHECK(result(1, 1) == 3.0);
}

TEST_CASE("Matrix - operator- subtracts element-wise", "[Math][Matrix]") {
  using namespace blast;
  Matrix a(2u, 2u, 5.0);
  Matrix b(2u, 2u, 2.0);
  Matrix result = a - b;
  CHECK(result(0, 0) == 3.0);
  CHECK(result(1, 1) == 3.0);
}

TEST_CASE("Matrix - left scalar operator* scales all elements", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u, 3.0);
  real   c      = 2.0;
  Matrix result = c * m;
  CHECK(result(0, 0) == 6.0);
  CHECK(result(1, 1) == 6.0);
}

TEST_CASE("Matrix - operator/ divides all elements by scalar", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u, 6.0);
  real   c      = 2.0;
  Matrix result = m / c;
  CHECK(result(0, 0) == 3.0);
  CHECK(result(1, 1) == 3.0);
}

TEST_CASE("Matrix - operator*(Matrix, Matrix) produces correct product", "[Math][Matrix]") {
  using namespace blast;
  // A is 2x3, B is 3x2, C = A*B is 2x2
  Matrix A(2u, 3u);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(0, 2) = 3.0;
  A(1, 0) = 4.0;
  A(1, 1) = 5.0;
  A(1, 2) = 6.0;

  Matrix B(3u, 2u);
  B(0, 0) = 7.0;
  B(0, 1) = 8.0;
  B(1, 0) = 9.0;
  B(1, 1) = 10.0;
  B(2, 0) = 11.0;
  B(2, 1) = 12.0;

  Matrix C = A * B;
  CHECK(C.rows == 2u);
  CHECK(C.cols == 2u);
  CHECK(std::abs(C(0, 0) - 58.0) < blast::test::abs_tol);
  CHECK(std::abs(C(0, 1) - 64.0) < blast::test::abs_tol);
  CHECK(std::abs(C(1, 0) - 139.0) < blast::test::abs_tol);
  CHECK(std::abs(C(1, 1) - 154.0) < blast::test::abs_tol);
}

TEST_CASE("Matrix - operator*(Matrix, Array) produces correct product", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u);
  m(0, 0) = 1.0;
  m(0, 1) = 2.0;
  m(1, 0) = 3.0;
  m(1, 1) = 4.0;

  Array v{1.0, 1.0};
  Array result = m * v;
  CHECK(result.size == 2u);
  CHECK(std::abs(result[0] - 3.0) < blast::test::abs_tol);
  CHECK(std::abs(result[1] - 7.0) < blast::test::abs_tol);
}

// ---------------------------------------------------------------------------
// Transpose
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - transpose() swaps dimensions and elements", "[Math][Matrix]") {
  using namespace blast;
  Matrix A(2u, 3u);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(0, 2) = 3.0;
  A(1, 0) = 4.0;
  A(1, 1) = 5.0;
  A(1, 2) = 6.0;

  Matrix T = transpose(A);
  CHECK(T.rows == 3u);
  CHECK(T.cols == 2u);
  for (u32 r = 0; r < A.rows; r++)
    for (u32 c = 0; c < A.cols; c++)
      CHECK(std::abs(T(c, r) - A(r, c)) < blast::test::abs_tol);
}

// ---------------------------------------------------------------------------
// Identity
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - eye(n) returns the n×n identity matrix", "[Math][Matrix]") {
  using namespace blast;
  Matrix I = eye(3);
  CHECK(I.rows == 3u);
  CHECK(I.cols == 3u);
  for (u32 r = 0; r < 3u; r++)
    for (u32 c = 0; c < 3u; c++)
      CHECK(I(r, c) == (r == c ? 1.0 : 0.0));
}

TEST_CASE("Matrix - eye(n) * eye(n) equals eye(n)", "[Math][Matrix]") {
  using namespace blast;
  Matrix I = eye(3);
  CHECK(I * I == I);
}

TEST_CASE("Matrix - multiplying by identity leaves matrix unchanged", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u);
  m(0, 0) = 1.0;
  m(0, 1) = 2.0;
  m(1, 0) = 3.0;
  m(1, 1) = 4.0;
  CHECK(m * eye(2) == m);
}

// ---------------------------------------------------------------------------
// Column access
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - col() reads the correct column elements", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(3u, 2u);
  m(0, 0) = 1.0;
  m(1, 0) = 2.0;
  m(2, 0) = 3.0;
  m(0, 1) = 4.0;
  m(1, 1) = 5.0;
  m(2, 1) = 6.0;

  Array col0 = m.col(0);
  CHECK(col0.is_alias == true);
  CHECK(col0.size == 3u);
  CHECK(col0[0] == 1.0);
  CHECK(col0[1] == 2.0);
  CHECK(col0[2] == 3.0);
}

TEST_CASE("Matrix - col() element write-through modifies the matrix", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(3u, 2u);
  m.col(0)[1] = 77.0;
  CHECK(m(1, 0) == 77.0);
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - zero() sets all elements to 0.0", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u, 5.0);
  zero(m);
  CHECK(is_small(m));
}

TEST_CASE("Matrix - zero() returns a reference to the same object", "[Math][Matrix]") {
  using namespace blast;
  Matrix  m(2u, 2u, 5.0);
  Matrix& ref = zero(m);
  CHECK(&ref == &m);
}

TEST_CASE("Matrix - constant() sets all elements to the given value", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 3u);
  constant(m, 7.0);
  for (u32 r = 0; r < m.rows; r++)
    for (u32 c = 0; c < m.cols; c++)
      CHECK(m(r, c) == 7.0);
}

TEST_CASE("Matrix - constant() returns a reference to the same object", "[Math][Matrix]") {
  using namespace blast;
  Matrix  m(2u, 2u);
  Matrix& ref = constant(m, 1.0);
  CHECK(&ref == &m);
}

TEST_CASE("Matrix - is_small() returns true when all elements are within epsilon", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u, 1e-6);
  CHECK(is_small(m));
}

TEST_CASE("Matrix - is_small() returns false when any element exceeds epsilon", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u);
  m(0, 1) = 2e-5;
  CHECK_FALSE(is_small(m));
}

TEST_CASE("Matrix - is_close() returns true when matrices are within epsilon", "[Math][Matrix]") {
  using namespace blast;
  Matrix a(2u, 2u, 1.0);
  Matrix b(2u, 2u, 1.0 + 1e-6);
  CHECK(is_close(a, b));
}

TEST_CASE("Matrix - is_close() returns false when any element difference exceeds epsilon", "[Math][Matrix]") {
  using namespace blast;
  Matrix a(2u, 2u, 1.0);
  Matrix b(2u, 2u, 1.0 + 2e-5);
  CHECK_FALSE(is_close(a, b));
}

// ---------------------------------------------------------------------------
// Determinant
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - determinant() of a 2x2 matrix is correct", "[Math][Matrix]") {
  using namespace blast;
  Matrix m(2u, 2u);
  m(0, 0) = 1.0;
  m(0, 1) = 2.0;
  m(1, 0) = 3.0;
  m(1, 1) = 4.0;
  // det = 1*4 - 2*3 = -2
  CHECK(std::abs(determinant(m) - (-2.0)) < blast::test::abs_tol);
}

TEST_CASE("Matrix - determinant() of a 3x3 matrix is correct", "[Math][Matrix]") {
  using namespace blast;
  // [[1,2,3],[0,4,5],[1,0,6]] -> det = 22
  Matrix m(3u, 3u);
  m(0, 0) = 1.0;
  m(0, 1) = 2.0;
  m(0, 2) = 3.0;
  m(1, 0) = 0.0;
  m(1, 1) = 4.0;
  m(1, 2) = 5.0;
  m(2, 0) = 1.0;
  m(2, 1) = 0.0;
  m(2, 2) = 6.0;
  CHECK(std::abs(determinant(m) - 22.0) < blast::test::abs_tol);
}

// ---------------------------------------------------------------------------
// LU decomposition and solve
// ---------------------------------------------------------------------------

TEST_CASE("Matrix - solveLU() solves Ax=b via LU decomposition", "[Math][Matrix]") {
  using namespace blast;
  // A = [[2,1],[5,3]], b = {1,2}, solution x = {1,-1}
  Matrix A(2u, 2u);
  A(0, 0) = 2.0;
  A(0, 1) = 1.0;
  A(1, 0) = 5.0;
  A(1, 1) = 3.0;

  Array b{1.0, 2.0};
  Array x = solveLU(LU_decomp(A), b);

  CHECK(std::abs(x[0] - 1.0) < blast::test::abs_tol);
  CHECK(std::abs(x[1] - (-1.0)) < blast::test::abs_tol);
}

TEST_CASE("Matrix - solveLU() solution satisfies Ax = b", "[Math][Matrix]") {
  using namespace blast;
  Matrix A(2u, 2u);
  A(0, 0) = 2.0;
  A(0, 1) = 1.0;
  A(1, 0) = 5.0;
  A(1, 1) = 3.0;

  Array b{1.0, 2.0};
  Array x  = solveLU(LU_decomp(A), b);
  Array Ax = A * x;
  Array expected{1.0, 2.0};
  CHECK(Ax == expected);
}
