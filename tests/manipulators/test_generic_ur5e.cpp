#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;

TEST_CASE("UR5e forward_kinematics() test", "[Generic]") {
  Manipulator generic_manip = make_UR5e();

  ManipulatorTempData data;

  //  @ home position
  Vec3  expected_tool_position = {0.0, -0.233, 1.08};
  Array test_position          = {0.0, -PI / 2, 0.0, -PI / 2, 0.0, 0.0};
  forward_kinematics(generic_manip, data, test_position);
  CHECK(is_close(data.p_j[generic_manip.n_joints], expected_tool_position, 1e-2));

  print(data.p_j[generic_manip.n_joints]);

  //  @ 2nd position
  Vec3  expected_tool_position_2 = {-0.1179, -0.436, 0.152};
  Array test_position_2          = {-1.54, -1.83, -2.28, -0.59, 1.60, 0.023};
  forward_kinematics(generic_manip, data, test_position_2);

  print(data.p_j[generic_manip.n_joints]);
  std::cout << std::endl;
  CHECK(is_close(data.p_j[generic_manip.n_joints], expected_tool_position_2, 1e-2));
}
