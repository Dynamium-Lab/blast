#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include <blast>
#include "catch2/catch.hpp"
#include "test_helper/test_helper.hpp"

namespace blast {

// No significant difference
TEST_CASE("Benchmark of .at() vs [] for std::vector", "[general]") {
  std::vector<Capsule> capsules(250);
  Capsule              caps;
  caps.p1 = {1, 1, 1};
  caps.p2 = {2, 2, 2};
  caps.r  = 2.0;
  for (auto& capsule: capsules) {
    capsule = caps;
  }

  real result = 0;
  BENCHMARK(".at() performance") {
    result = 0;
    for (int i = 0; i < capsules.size(); i++) {
      result += capsules.at(i).r;
    }
    return result;
  };
  real at_result = result;

  BENCHMARK("[] performance") {
    result = 0;
    for (int i = 0; i < capsules.size(); i++) {
      result += capsules[i].r;
    }
    return result;
  };

  BENCHMARK("range based for loop performance") {
    result = 0;
    for (auto& capsule: capsules) {
      result += capsule.r;
    }
    return result;
  };
  CHECK(at_result == Approx(result));
}
} // namespace blast
