#include <blast>
#include "test_helper/test_helper.hpp"

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

namespace blast {


// Testing by filling result instead of returning. No significant difference
inline void bound_constraint(real* result, const real& value, const real& lower_bound, const real& higher_bound) {
  real denom             = higher_bound - lower_bound;
  denom                  = denom == 0 ? 1 : 0; // todo: what??
  real constraint_lower  = (lower_bound == -INF_REAL) ? -1.0 : -(value - lower_bound) / denom;
  real constraint_higher = (higher_bound == INF_REAL) ? -1.0 : (value - higher_bound) / denom;
  result[0]              = constraint_higher >= constraint_lower ? constraint_higher : constraint_lower;
}

// No significant difference
// note : bounds constraint function changed since testing. It is now more robust and cannot be directly compared.
TEST_CASE("Bound constraint function benchmark", "[constraint]") {
  auto opt   = get_generic_link6_opt();
  auto manip = opt.manip;

  Array result(manip.joints);

  Array pos = get_Link6_singularity();

  Array expected_result = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

  BENCHMARK_ADVANCED("Written position function")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      for (int j = 0; j < manip.joints; j++) {
        real denom = manip.pmax[j] - manip.pmin[j];
        denom      = denom == 0 ? 1 : 0;
        // denom = manip.pmin[j] == 0 ? 1 : manip.pmin[j];
        real constraint_lower = (manip.pmin[j] == -INF_REAL) ? -1.0 : (pos[j] - manip.pmin[j]) / denom;
        // denom = manip.pmax[j] == 0 ? 1 : manip.pmax[j];
        real constraint_higher = (manip.pmax[j] == INF_REAL) ? -1.0 : (pos[j] - manip.pmax[j]) / denom;
        result[j]              = constraint_higher >= constraint_lower ? constraint_higher : constraint_lower;
      }
      return result;
    });
    CHECK(is_close(result, expected_result));
  };


  BENCHMARK_ADVANCED("Bounds constraint position function with pointer")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      for (int j = 0; j < manip.joints; j++) {
        bound_constraint(&(result[j]), pos[j], manip.pmin[j], manip.pmax[j]);
        // result[j] bound_constraint(pos[j], manip.pmin[j], manip.pmax[j]);
      }
      return result;
    });
    CHECK(is_close(result, expected_result));
  };

  BENCHMARK_ADVANCED("Bounds constraint position function")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] {
      for (int j = 0; j < manip.joints; j++) {
        // bound_constraint(result.data[j], pos[j], manip.pmin[j], manip.pmax[j]);
        result[j] = bound_constraint(pos[j], manip.pmin[j], manip.pmax[j]);
      }
      return result;
    });
    CHECK(is_close(result, expected_result));
  };
}

} // namespace blast
