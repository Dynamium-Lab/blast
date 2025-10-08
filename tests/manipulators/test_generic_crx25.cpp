#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "blast_rush.h"
#include "test_helper/test_helper.hpp"
#include "test_helper/test_functions.hpp"

using namespace blast;


TEST_CASE("Fanuc Crx25ia forward_kinematics() test", "[Generic]") {
    int n_tests = 1;
    real epsilon = 1e-6;
    GenericManipulator manip = get_generic_fanuc_crx25ia();

    for (int i = 0; i < n_tests; i++) {
        Array test_position(6);
        fill_random(test_position, PI);
        print(test_position);
        std::cout << std::endl;
        forward_kinematics(manip, test_position);
        int stop = 0;
    }
}