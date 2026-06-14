#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;

TEST_CASE("bound_constraint - center is fully satisfied (-1)", "[Optimization][Constraints]") {
  CHECK(bound_constraint(0.0, -1.0, 1.0) == Approx(-1.0));
  CHECK(bound_constraint(5.0, 0.0, 10.0) == Approx(-1.0));
}

TEST_CASE("bound_constraint - at boundary is zero", "[Optimization][Constraints]") {
  CHECK(bound_constraint(1.0, -1.0, 1.0) == Approx(0.0));
  CHECK(bound_constraint(-1.0, -1.0, 1.0) == Approx(0.0));
  CHECK(bound_constraint(10.0, 0.0, 10.0) == Approx(0.0));
}

TEST_CASE("bound_constraint - beyond boundary is positive (violated)", "[Optimization][Constraints]") {
  CHECK(bound_constraint(2.0, -1.0, 1.0) > 0.0);
  CHECK(bound_constraint(-2.0, -1.0, 1.0) > 0.0);
}

TEST_CASE("abs_constraint - zero is fully satisfied (-1)", "[Optimization][Constraints]") {
  CHECK(abs_constraint(0.0, 1.0) == Approx(-1.0));
  CHECK(abs_constraint(0.0, 5.0) == Approx(-1.0));
}

TEST_CASE("abs_constraint - at max is zero (at limit)", "[Optimization][Constraints]") {
  CHECK(abs_constraint(1.0, 1.0) == Approx(0.0));
  CHECK(abs_constraint(-1.0, 1.0) == Approx(0.0));
}

TEST_CASE("abs_constraint - beyond max is positive (violated)", "[Optimization][Constraints]") {
  CHECK(abs_constraint(2.0, 1.0) > 0.0);
  CHECK(abs_constraint(-2.0, 1.0) > 0.0);
}

TEST_CASE("abs_constraint_analytical<true> - value matches abs_constraint, gradient has correct sign", "[Optimization][Constraints]") {
  {
    auto [c, g] = abs_constraint_analytical<true>(2.0, 4.0);
    CHECK(c == Approx(abs_constraint(2.0, 4.0)));
    CHECK(g == Approx(1.0 / 4.0));
  }
  {
    auto [c, g] = abs_constraint_analytical<true>(-2.0, 4.0);
    CHECK(c == Approx(abs_constraint(-2.0, 4.0)));
    CHECK(g == Approx(-1.0 / 4.0));
  }
}

TEST_CASE("abs_constraint_analytical<false> - gradient is zero", "[Optimization][Constraints]") {
  auto [c, g] = abs_constraint_analytical<false>(2.0, 4.0);
  CHECK(c == Approx(abs_constraint(2.0, 4.0)));
  CHECK(g == 0.0);
}

TEST_CASE("bound_constraint_analytical<true> - value matches bound_constraint, gradient correct sign", "[Optimization][Constraints]") {
  {
    // value > center: gradient should be +2/range
    auto [c, g] = bound_constraint_analytical<true>(0.5, -1.0, 1.0);
    CHECK(c == Approx(bound_constraint(0.5, -1.0, 1.0)));
    CHECK(g == Approx(2.0 / 2.0));
  }
  {
    // value < center: gradient should be -2/range
    auto [c, g] = bound_constraint_analytical<true>(-0.5, -1.0, 1.0);
    CHECK(c == Approx(bound_constraint(-0.5, -1.0, 1.0)));
    CHECK(g == Approx(-2.0 / 2.0));
  }
}

TEST_CASE("bound_constraint_analytical<false> - gradient is zero", "[Optimization][Constraints]") {
  auto [c, g] = bound_constraint_analytical<false>(0.5, -1.0, 1.0);
  CHECK(c == Approx(bound_constraint(0.5, -1.0, 1.0)));
  CHECK(g == 0.0);
}

TEST_CASE("bound_constraint_analytical - INF limits return -1 (unconstrained)", "[Optimization][Constraints]") {
  {
    auto [c, g] = bound_constraint_analytical<true>(100.0, -INF_REAL, 1.0);
    CHECK(c == Approx(-1.0));
    CHECK(g == 0.0);
  }
  {
    auto [c, g] = bound_constraint_analytical<true>(100.0, -1.0, INF_REAL);
    CHECK(c == Approx(-1.0));
    CHECK(g == 0.0);
  }
}
