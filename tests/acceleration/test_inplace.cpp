#include <blast>
#include "test_helper/test_helper.hpp"

#define BLAST_TRACE_LEVEL 3

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

#include "tracy/Tracy.hpp"
#include "tracy/TracyC.h"

namespace blast {


TEST_CASE("replace INF_REAL with huge value", "[acceleratipon]") {
  int n_iterations = 10000;

  real pmax = INF_REAL;
  real pmin = -INF_REAL;

  pmax = 1e300;
  pmin = -1e300;

  real c          = (pmax + pmin) / 2;
  real b          = (pmax - pmin) / 2;
  real one_over_b = 1 / b;

  Array p(n_iterations);
  p = fill_random(p, 1e200);

  for (int i = 0; i < n_iterations; i++) {
    real q_prime    = p[i] - c;
    real constraint = (std::abs(q_prime) - b) / b;
    CHECK(constraint == -1.0);
  }
};

} // namespace blast
