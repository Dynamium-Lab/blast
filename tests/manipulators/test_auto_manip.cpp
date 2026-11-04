#define CATCH_CONFIG_MAIN

#include <blast>
#include "catch2/catch.hpp"

TEST_CASE("make_Kinova_Gen3_auto() test", "[Generic]") {
  auto auto_manip    = blast::make_Kinova_Gen3_auto();
  auto generic_manip = make_Kinova_Gen3();
  CHECK(auto_manip == generic_manip);
}

TEST_CASE("make_Kinova_Link6_auto() test", "[Generic]") {
  auto auto_manip    = make_Kinova_Link6_auto();
  auto generic_manip = make_Kinova_Link6();
  CHECK(auto_manip == generic_manip);
}
