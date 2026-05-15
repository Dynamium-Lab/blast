#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>
#include "manipulator/UR5e.hpp"
#include "test_helper/test_functions.hpp"
#include "test_helper/test_helper.hpp"


TEST_CASE("Conversions Q & rx, ry, rz", "[Misc]") {
  using namespace blast;
  int                 n_tests = 1;
  auto                manip   = make_UR5e();
  ManipulatorTempData data;

  for (u32 i = 0; i < n_tests; i++) {

    Array joint_pos(6);

    fill_random(joint_pos, PI);

    forward_kinematics(manip, data, joint_pos);
    auto rxyz_rot = rotation2rpy(data.rotations_mult[5]);
    auto rot_rxyz = rpy2rotation(rxyz_rot);

    CHECK(is_close(data.rotations_mult[5], rot_rxyz, 0.001));
  }
}
