#define CATCH_CONFIG_MAIN

#include <blast>
#include "catch2/catch.hpp"

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

TEST_CASE("Array - default constructor produces empty array", "[Math][Array]") {
  using namespace blast;
  Array a;
  CHECK(a.size == 0u);
  CHECK(a.data == nullptr);
  CHECK(a.is_alias == false);
  CHECK(a.is_empty());
}

TEST_CASE("Array - size constructor allocates and zeroes elements", "[Math][Array]") {
  using namespace blast;
  Array a(4u);
  CHECK(a.size == 4u);
  CHECK(a.is_empty() == false);
  for (u32 i = 0; i < a.size; i++)
    CHECK(a[i] == 0.0);
}

TEST_CASE("Array - size+value constructor fills all elements", "[Math][Array]") {
  using namespace blast;
  Array a(3u, 7.0);
  CHECK(a.size == 3u);
  for (u32 i = 0; i < a.size; i++)
    CHECK(a[i] == 7.0);
}

TEST_CASE("Array - initializer list constructor stores elements in order", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  CHECK(a.size == 3u);
  CHECK(a[0] == 1.0);
  CHECK(a[1] == 2.0);
  CHECK(a[2] == 3.0);
}

TEST_CASE("Array - copy constructor produces a deep independent copy", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  Array b(a);
  CHECK(b.size == 3u);
  CHECK(b.is_alias == false);
  CHECK(b[0] == 1.0);
  b[0] = 99.0;
  CHECK(a[0] == 1.0);
}

TEST_CASE("Array - move constructor transfers ownership", "[Math][Array]") {
  using namespace blast;
  Array a(3u, 1.0);
  Array b(std::move(a));
  CHECK(b.size == 3u);
  CHECK(b[0] == 1.0);
  CHECK(a.size == 0u);
}

// ---------------------------------------------------------------------------
// Element access
// ---------------------------------------------------------------------------

TEST_CASE("Array - operator[] provides mutable element access", "[Math][Array]") {
  using namespace blast;
  Array a(3u);
  a[0] = 42.0;
  CHECK(a[0] == 42.0);
}

TEST_CASE("Array - operator[] provides const element access", "[Math][Array]") {
  using namespace blast;
  const Array a{1.0, 2.0, 3.0};
  CHECK(a[0] == 1.0);
  CHECK(a[2] == 3.0);
}

TEST_CASE("Array - back() returns and allows modifying the last element", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  CHECK(a.back() == 3.0);
  a.back() = 99.0;
  CHECK(a[2] == 99.0);
}

TEST_CASE("Array - back() const returns the last element", "[Math][Array]") {
  using namespace blast;
  const Array a{1.0, 2.0, 5.0};
  CHECK(a.back() == 5.0);
}

// ---------------------------------------------------------------------------
// Size and state
// ---------------------------------------------------------------------------

TEST_CASE("Array - is_empty() reflects whether the array has elements", "[Math][Array]") {
  using namespace blast;
  CHECK(Array().is_empty() == true);
  CHECK(Array(3u).is_empty() == false);
}

// ---------------------------------------------------------------------------
// Equality
// ---------------------------------------------------------------------------

TEST_CASE("Array - operator== returns true when arrays are within epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  Array b{1.0 + 1e-6, 2.0 + 1e-6, 3.0 + 1e-6};
  CHECK(a == b);
}

TEST_CASE("Array - operator== returns false when any element difference exceeds epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  Array b{1.0 + 2e-5, 2.0, 3.0};
  CHECK_FALSE(a == b);
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

TEST_CASE("Array - unary minus negates all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, -2.0, 3.0};
  Array result = -a;
  Array expected{-1.0, 2.0, -3.0};
  CHECK(result == expected);
}

TEST_CASE("Array - operator+ adds element-wise", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  Array b{4.0, 5.0, 6.0};
  Array expected{5.0, 7.0, 9.0};
  CHECK(a + b == expected);
}

TEST_CASE("Array - operator- subtracts element-wise", "[Math][Array]") {
  using namespace blast;
  Array a{4.0, 5.0, 6.0};
  Array b{1.0, 2.0, 3.0};
  Array expected{3.0, 3.0, 3.0};
  CHECK(a - b == expected);
}

TEST_CASE("Array - right scalar operator* scales all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  Array expected{2.0, 4.0, 6.0};
  CHECK(a * 2.0 == expected);
}

TEST_CASE("Array - left scalar operator* scales all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  Array expected{2.0, 4.0, 6.0};
  CHECK(2.0 * a == expected);
}

TEST_CASE("Array - operator/ divides all elements by scalar", "[Math][Array]") {
  using namespace blast;
  Array a{2.0, 4.0, 6.0};
  Array expected{1.0, 2.0, 3.0};
  CHECK(a / 2.0 == expected);
}

TEST_CASE("Array - operator*= scales in place", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  a *= 3.0;
  Array expected{3.0, 6.0, 9.0};
  CHECK(a == expected);
}

TEST_CASE("Array - operator/= divides in place", "[Math][Array]") {
  using namespace blast;
  Array a{4.0, 6.0, 8.0};
  a /= 2.0;
  Array expected{2.0, 3.0, 4.0};
  CHECK(a == expected);
}

// ---------------------------------------------------------------------------
// Dot product
// ---------------------------------------------------------------------------

TEST_CASE("Array - dot() computes correct inner product", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  Array b{4.0, 5.0, 6.0};
  CHECK(std::abs(dot(a, b) - 32.0) < 1e-9);
}

TEST_CASE("Array - dot(a, a) equals norm_sqr(a)", "[Math][Array]") {
  using namespace blast;
  Array a{3.0, 4.0};
  CHECK(std::abs(dot(a, a) - norm_sqr(a)) < 1e-9);
  CHECK(std::abs(dot(a, a) - 25.0) < 1e-9);
}

TEST_CASE("Array - dot() of orthogonal arrays is zero", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 0.0, 0.0};
  Array b{0.0, 1.0, 0.0};
  CHECK(std::abs(dot(a, b)) < 1e-9);
}

// ---------------------------------------------------------------------------
// Norms
// ---------------------------------------------------------------------------

TEST_CASE("Array - norm() computes the Euclidean norm", "[Math][Array]") {
  using namespace blast;
  Array a{3.0, 4.0};
  CHECK(std::abs(norm(a) - 5.0) < 1e-9);
}

TEST_CASE("Array - norm_sqr() computes the squared Euclidean norm", "[Math][Array]") {
  using namespace blast;
  Array a{3.0, 4.0};
  CHECK(std::abs(norm_sqr(a) - 25.0) < 1e-9);
}

TEST_CASE("Array - norm_1() computes the L1 norm", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, -2.0, 3.0};
  CHECK(std::abs(norm_1(a) - 6.0) < 1e-9);
}

TEST_CASE("Array - norm_inf() computes the L-infinity norm", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, -5.0, 3.0};
  CHECK(std::abs(norm_inf(a) - 5.0) < 1e-9);
}

// ---------------------------------------------------------------------------
// Stats
// ---------------------------------------------------------------------------

TEST_CASE("Array - sum() computes the sum of all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  CHECK(std::abs(sum(a) - 6.0) < 1e-9);
}

TEST_CASE("Array - mean() computes the arithmetic mean", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  CHECK(std::abs(mean(a) - 2.0) < 1e-9);
}

TEST_CASE("Array - min() returns the smallest element", "[Math][Array]") {
  using namespace blast;
  Array a{3.0, 1.0, 2.0};
  CHECK(std::abs(min(a) - 1.0) < 1e-9);
}

TEST_CASE("Array - max() returns the largest element", "[Math][Array]") {
  using namespace blast;
  Array a{3.0, 1.0, 2.0};
  CHECK(std::abs(max(a) - 3.0) < 1e-9);
}

TEST_CASE("Array - argmin() returns the index of the smallest element", "[Math][Array]") {
  using namespace blast;
  Array a{3.0, 1.0, 2.0};
  CHECK(argmin(a) == 1u);
}

TEST_CASE("Array - argmax() returns the index of the largest element", "[Math][Array]") {
  using namespace blast;
  Array a{3.0, 1.0, 2.0};
  CHECK(argmax(a) == 0u);
}

// ---------------------------------------------------------------------------
// Element-wise operations
// ---------------------------------------------------------------------------

TEST_CASE("Array - abs() computes absolute value of each element", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, -2.0, 3.0};
  Array expected{1.0, 2.0, 3.0};
  CHECK(abs(a) == expected);
}

TEST_CASE("Array - abs_inplace() modifies array in place", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, -2.0, -3.0};
  abs_inplace(a);
  Array expected{1.0, 2.0, 3.0};
  CHECK(a == expected);
}

TEST_CASE("Array - sqr() squares each element", "[Math][Array]") {
  using namespace blast;
  Array a{2.0, 3.0};
  Array expected{4.0, 9.0};
  CHECK(sqr(a) == expected);
}

TEST_CASE("Array - sqr_inplace() squares each element in place", "[Math][Array]") {
  using namespace blast;
  Array a{2.0, 3.0, 4.0};
  sqr_inplace(a);
  Array expected{4.0, 9.0, 16.0};
  CHECK(a == expected);
}

TEST_CASE("Array - clamp() clamps elements to [min, max]", "[Math][Array]") {
  using namespace blast;
  Array a{-1.0, 0.5, 2.0};
  Array result   = clamp(a, 0.0, 1.0);
  Array expected{0.0, 0.5, 1.0};
  CHECK(result == expected);
}

// ---------------------------------------------------------------------------
// Alias
// ---------------------------------------------------------------------------

TEST_CASE("Array - pointer constructor creates alias without ownership", "[Math][Array]") {
  using namespace blast;
  real  buf[3] = {1.0, 2.0, 3.0};
  Array a(buf, 3u);
  CHECK(a.is_alias == true);
  CHECK(a.size == 3u);
  CHECK(a[0] == 1.0);
}

TEST_CASE("Array - writing through alias modifies the original buffer", "[Math][Array]") {
  using namespace blast;
  real  buf[3] = {1.0, 2.0, 3.0};
  Array a(buf, 3u);
  a[0] = 99.0;
  CHECK(buf[0] == 99.0);
}

// ---------------------------------------------------------------------------
// Resize
// ---------------------------------------------------------------------------

TEST_CASE("Array - resize() increases size and preserves existing elements", "[Math][Array]") {
  using namespace blast;
  Array a(3u, 1.0);
  a.resize(5u);
  CHECK(a.size == 5u);
  CHECK(a[0] == 1.0);
  CHECK(a[2] == 1.0);
  CHECK(a[3] == 0.0);
  CHECK(a[4] == 0.0);
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

TEST_CASE("Array - zero() sets all elements to 0.0", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0, 3.0};
  zero(a);
  CHECK(is_small(a));
}

TEST_CASE("Array - zero() returns a reference to the same object", "[Math][Array]") {
  using namespace blast;
  Array  a{1.0, 2.0, 3.0};
  Array& ref = zero(a);
  CHECK(&ref == &a);
}

TEST_CASE("Array - constant() sets all elements to the given value", "[Math][Array]") {
  using namespace blast;
  Array a(3u);
  constant(a, 4.0);
  for (u32 i = 0; i < a.size; i++)
    CHECK(a[i] == 4.0);
}

TEST_CASE("Array - constant() returns a reference to the same object", "[Math][Array]") {
  using namespace blast;
  Array  a(3u);
  Array& ref = constant(a, 2.0);
  CHECK(&ref == &a);
}

TEST_CASE("Array - is_small() returns true when all elements are within epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1e-6, 1e-7, 1e-8};
  CHECK(is_small(a));
}

TEST_CASE("Array - is_small() returns false when any element exceeds epsilon", "[Math][Array]") {
  using namespace blast;
  CHECK_FALSE(is_small(Array{2e-5, 0.0, 0.0}));
  CHECK_FALSE(is_small(Array{0.0, 2e-5, 0.0}));
}

TEST_CASE("Array - is_close() returns true when arrays are within epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0};
  Array b{1.0 + 1e-6, 2.0 + 1e-6};
  CHECK(is_close(a, b));
}

TEST_CASE("Array - is_close() returns false when any element difference exceeds epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1.0, 2.0};
  Array b{1.0 + 2e-5, 2.0};
  CHECK_FALSE(is_close(a, b));
}
