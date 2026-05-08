#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include <blast>
#include "catch2/catch.hpp"
#include "test_helper/test_functions.hpp"

using namespace blast;


TEST_CASE("BVH function tests", "[World]") {

}

TEST_CASE("SAP function tests", "[World]") {

}

TEST_CASE("Broadphase benchmarks", "[World]") {
  World world;

  BENCHMARK("BVH") {
    real dist = 0;
    return dist;
  };

  BENCHMARK("BVH") {
    real dist = 0;
    return dist;
  };
}