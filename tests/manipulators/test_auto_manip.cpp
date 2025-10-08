#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include <blast>
#include "test_helper/test_helper.hpp"
#include "test_helper/test_functions.hpp"



TEST_CASE("get_generic_gen3_auto() test", "[Generic]") {
    auto auto_manip = blast::get_generic_gen3_auto();
    auto generic_manip = get_generic_gen3();
    CHECK(is_close(auto_manip, generic_manip));
}

TEST_CASE("get_generic_link6_auto() test", "[Generic]") {
    auto auto_manip = get_generic_link6_auto();
    auto generic_manip = get_generic_Link6();
    CHECK(is_close(auto_manip, generic_manip));
}


