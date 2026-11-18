#define CATCH_CONFIG_MAIN

#include <blast>
#include <cstdint>
#include "catch2/catch.hpp"
#include "test_helper.hpp"

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
    CHECK(a[i] == 0);
}

TEST_CASE("Array - size+value constructor fills all elements", "[Math][Array]") {
  using namespace blast;
  Array a(3u, 7.0);
  CHECK(a.size == 3u);
  for (u32 i = 0; i < a.size; i++)
    CHECK(a[i] == 7);
}

TEST_CASE("Array - initializer list constructor stores elements in order", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  CHECK(a.size == 3u);
  CHECK(a[0] == 1);
  CHECK(a[1] == 2);
  CHECK(a[2] == 3);
}

TEST_CASE("Array - copy constructor produces a deep independent copy", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  Array b(a);
  CHECK(b.size == 3u);
  CHECK(b.is_alias == false);
  CHECK(b[0] == 1);
  b[0] = 99;
  CHECK(a[0] == 1);
}

TEST_CASE("Array - move constructor transfers ownership", "[Math][Array]") {
  using namespace blast;
  Array a(3u, 1.0);
  Array b(std::move(a));
  CHECK(b.size == 3u);
  CHECK(b[0] == 1);
  CHECK(a.size == 0u);
}

// ---------------------------------------------------------------------------
// Element access
// ---------------------------------------------------------------------------

TEST_CASE("Array - operator[] provides mutable element access", "[Math][Array]") {
  using namespace blast;
  Array a(3u);
  a[0] = 42;
  CHECK(a[0] == 42);
}

TEST_CASE("Array - operator[] provides const element access", "[Math][Array]") {
  using namespace blast;
  const Array a{1, 2, 3};
  CHECK(a[0] == 1);
  CHECK(a[2] == 3);
}

TEST_CASE("Array - back() returns and allows modifying the last element", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  CHECK(a.back() == 3);
  a.back() = 99;
  CHECK(a[2] == 99);
}

TEST_CASE("Array - back() const returns the last element", "[Math][Array]") {
  using namespace blast;
  const Array a{1, 2, 5};
  CHECK(a.back() == 5);
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
  Array a{1, 2, 3};
  Array b{1.0 + 1e-6, 2.0 + 1e-6, 3.0 + 1e-6};
  CHECK(a == b);
}

TEST_CASE("Array - operator== returns false when any element difference exceeds epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  Array b{1.0 + 2e-5, 2, 3};
  CHECK_FALSE(a == b);
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

TEST_CASE("Array - unary minus negates all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1, -2, 3};
  Array result = -a;
  Array expected{-1, 2, -3};
  CHECK(result == expected);
}

TEST_CASE("Array - operator+ adds element-wise", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  Array b{4, 5, 6};
  Array expected{5, 7, 9};
  CHECK(a + b == expected);
}

TEST_CASE("Array - operator- subtracts element-wise", "[Math][Array]") {
  using namespace blast;
  Array a{4, 5, 6};
  Array b{1, 2, 3};
  Array expected{3, 3, 3};
  CHECK(a - b == expected);
}

TEST_CASE("Array - right scalar operator* scales all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  Array expected{2, 4, 6};
  CHECK(a * 2.0 == expected);
}

TEST_CASE("Array - left scalar operator* scales all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  Array expected{2, 4, 6};
  CHECK(2.0 * a == expected);
}

TEST_CASE("Array - operator/ divides all elements by scalar", "[Math][Array]") {
  using namespace blast;
  Array a{2, 4, 6};
  Array expected{1, 2, 3};
  CHECK(a / 2.0 == expected);
}

TEST_CASE("Array - operator*= scales in place", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  a *= 3;
  Array expected{3, 6, 9};
  CHECK(a == expected);
}

TEST_CASE("Array - operator/= divides in place", "[Math][Array]") {
  using namespace blast;
  Array a{4, 6, 8};
  a /= 2;
  Array expected{2, 3, 4};
  CHECK(a == expected);
}

// ---------------------------------------------------------------------------
// Dot product
// ---------------------------------------------------------------------------

TEST_CASE("Array - dot() computes correct inner product", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  Array b{4, 5, 6};
  CHECK(std::abs(dot(a, b) - 32.0) < blast::test::abs_tol);
}

TEST_CASE("Array - dot(a, a) equals norm_sqr(a)", "[Math][Array]") {
  using namespace blast;
  Array a{3, 4};
  CHECK(std::abs(dot(a, a) - norm_sqr(a)) < blast::test::abs_tol);
  CHECK(std::abs(dot(a, a) - 25.0) < blast::test::abs_tol);
}

TEST_CASE("Array - dot() of orthogonal arrays is zero", "[Math][Array]") {
  using namespace blast;
  Array a{1, 0, 0};
  Array b{0, 1, 0};
  CHECK(std::abs(dot(a, b)) < blast::test::abs_tol);
}

// ---------------------------------------------------------------------------
// Norms
// ---------------------------------------------------------------------------

TEST_CASE("Array - norm() computes the Euclidean norm", "[Math][Array]") {
  using namespace blast;
  Array a{3, 4};
  CHECK(std::abs(norm(a) - 5.0) < blast::test::abs_tol);
}

TEST_CASE("Array - norm_sqr() computes the squared Euclidean norm", "[Math][Array]") {
  using namespace blast;
  Array a{3, 4};
  CHECK(std::abs(norm_sqr(a) - 25.0) < blast::test::abs_tol);
}

TEST_CASE("Array - norm_1() computes the L1 norm", "[Math][Array]") {
  using namespace blast;
  Array a{1, -2, 3};
  CHECK(std::abs(norm_1(a) - 6.0) < blast::test::abs_tol);
}

TEST_CASE("Array - norm_inf() computes the L-infinity norm", "[Math][Array]") {
  using namespace blast;
  Array a{1, -5, 3};
  CHECK(std::abs(norm_inf(a) - 5.0) < blast::test::abs_tol);
}

// ---------------------------------------------------------------------------
// Stats
// ---------------------------------------------------------------------------

TEST_CASE("Array - sum() computes the sum of all elements", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  CHECK(std::abs(sum(a) - 6.0) < blast::test::abs_tol);
}

TEST_CASE("Array - mean() computes the arithmetic mean", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  CHECK(std::abs(mean(a) - 2.0) < blast::test::abs_tol);
}

TEST_CASE("Array - min() returns the smallest element", "[Math][Array]") {
  using namespace blast;
  Array a{3, 1, 2};
  CHECK(std::abs(min(a) - 1.0) < blast::test::abs_tol);
}

TEST_CASE("Array - max() returns the largest element", "[Math][Array]") {
  using namespace blast;
  Array a{3, 1, 2};
  CHECK(std::abs(max(a) - 3.0) < blast::test::abs_tol);
}

TEST_CASE("Array - argmin() returns the index of the smallest element", "[Math][Array]") {
  using namespace blast;
  Array a{3, 1, 2};
  CHECK(argmin(a) == 1u);
}

TEST_CASE("Array - argmax() returns the index of the largest element", "[Math][Array]") {
  using namespace blast;
  Array a{3, 1, 2};
  CHECK(argmax(a) == 0u);
}

// ---------------------------------------------------------------------------
// Element-wise operations
// ---------------------------------------------------------------------------

TEST_CASE("Array - abs() computes absolute value of each element", "[Math][Array]") {
  using namespace blast;
  Array a{1, -2, 3};
  Array expected{1, 2, 3};
  CHECK(abs(a) == expected);
}

TEST_CASE("Array - abs_inplace() modifies array in place", "[Math][Array]") {
  using namespace blast;
  Array a{1, -2, -3};
  abs_inplace(a);
  Array expected{1, 2, 3};
  CHECK(a == expected);
}

TEST_CASE("Array - sqr() squares each element", "[Math][Array]") {
  using namespace blast;
  Array a{2, 3};
  Array expected{4, 9};
  CHECK(sqr(a) == expected);
}

TEST_CASE("Array - sqr_inplace() squares each element in place", "[Math][Array]") {
  using namespace blast;
  Array a{2, 3, 4};
  sqr_inplace(a);
  Array expected{4, 9, 16};
  CHECK(a == expected);
}

TEST_CASE("Array - clamp() clamps elements to [min, max]", "[Math][Array]") {
  using namespace blast;
  Array a{-1, 0.5, 2};
  Array result = clamp(a, 0.0, 1.0);
  Array expected{0, 0.5, 1};
  CHECK(result == expected);
}

// ---------------------------------------------------------------------------
// Alias
// ---------------------------------------------------------------------------

TEST_CASE("Array - pointer constructor creates alias without ownership", "[Math][Array]") {
  using namespace blast;
  real  buf[3] = {1, 2, 3};
  Array a(buf, 3u);
  CHECK(a.is_alias == true);
  CHECK(a.size == 3u);
  CHECK(a[0] == 1);
}

TEST_CASE("Array - writing through alias modifies the original buffer", "[Math][Array]") {
  using namespace blast;
  real  buf[3] = {1, 2, 3};
  Array a(buf, 3u);
  a[0] = 99;
  CHECK(buf[0] == 99);
}

// ---------------------------------------------------------------------------
// Resize
// ---------------------------------------------------------------------------

TEST_CASE("Array - resize() increases size and preserves existing elements", "[Math][Array]") {
  using namespace blast;
  Array a(3u, 1.0);
  a.resize(5u);
  CHECK(a.size == 5u);
  CHECK(a[0] == 1);
  CHECK(a[2] == 1);
  CHECK(a[3] == 0);
  CHECK(a[4] == 0);
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

TEST_CASE("Array - zero() sets all elements to 0.0", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2, 3};
  zero(a);
  CHECK(is_small(a));
}

TEST_CASE("Array - zero() returns a reference to the same object", "[Math][Array]") {
  using namespace blast;
  Array  a{1, 2, 3};
  Array& ref = zero(a);
  CHECK(&ref == &a);
}

TEST_CASE("Array - constant() sets all elements to the given value", "[Math][Array]") {
  using namespace blast;
  Array a(3u);
  constant(a, 4.0);
  for (u32 i = 0; i < a.size; i++)
    CHECK(a[i] == 4);
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
  CHECK_FALSE(is_small(Array{2e-5, 0, 0}));
  CHECK_FALSE(is_small(Array{0, 2e-5, 0}));
}

TEST_CASE("Array - is_close() returns true when arrays are within epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2};
  Array b{1.0 + 1e-6, 2.0 + 1e-6};
  CHECK(is_close(a, b));
}

TEST_CASE("Array - is_close() returns false when any element difference exceeds epsilon", "[Math][Array]") {
  using namespace blast;
  Array a{1, 2};
  Array b{1.0 + 2e-5, 2};
  CHECK_FALSE(is_close(a, b));
}

// ---------------------------------------------------------------------------
// Memory
// ---------------------------------------------------------------------------

TEST_CASE("Array - an allocation whose byte size is not a multiple of the alignment", "[Math][Array]") {
  using namespace blast;
  // The rounding itself: glibc accepts an unrounded size, so only these checks catch losing it.
  CHECK(aligned_size(64, 24) == 64); // 3 doubles -> one aligned block
  CHECK(aligned_size(64, 64) == 64);
  CHECK(aligned_size(64, 65) == 128);

  Array a(3u); // 3 * sizeof(real) bytes, rounded up to ALIGN by the allocator
  REQUIRE(a.data != nullptr);
  CHECK(reinterpret_cast<std::uintptr_t>(a.data) % ALIGN == 0);
  a[0] = 1;
  a[1] = 2;
  a[2] = 3;
  CHECK(a[0] == 1);
  CHECK(a[2] == 3);
}

TEST_CASE("Array - resize grows an empty Array", "[Math][Array]") {
  using namespace blast;
  Array a;
  a.resize(3);
  REQUIRE(a.size == 3);
  CHECK(a[0] == 0);
  CHECK(a[2] == 0);
}

TEST_CASE("Array - copy-assigning an empty Array leaves no data behind", "[Math][Array]") {
  using namespace blast;
  Array a(3);
  Array empty;
  a = empty;
  CHECK(a.size == 0);
  CHECK(a.data == nullptr);
  a.resize(4); // used to free the old block a second time
  CHECK(a.size == 4);
  CHECK(a[3] == 0);
}

TEST_CASE("Array - assigning an empty list to an empty Array", "[Math][Array]") {
  using namespace blast;
  Array a;
  a = {};
  CHECK(a.size == 0);
}
