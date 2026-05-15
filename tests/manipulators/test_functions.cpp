#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>
#include "manipulator/UR5e.hpp"
#include "test_helper/test_functions.hpp"
#include "test_helper/test_helper.hpp"

// todo: fix inverse_kinematics
/*TEST_CASE("Inverse kinematics zero position test", "[Manipulator]") {
using namespace blast;
    auto manip = get_generic_Link6();
    Array goal_joint_pos = {0, -0.349, 1.92, 0, 0.698, 0}; // home
    // Array goal_joint_pos(6);
    // fill_random(goal_joint_pos, PI);
    Array joint_pos(6);

    // for (int i = 0; i < joint_pos.size; i++) {
    //     joint_pos[i] = goal_joint_pos[i] + 0.1;
    // }

    forward_kinematics(manip, goal_joint_pos);
    auto goal_pose = manip.get_tool_pose();

    auto joint_sol = inverse_kinematics(manip, goal_pose, joint_pos);
    auto joint_sol_nlopt = inverse_kinematics_nlopt(manip, goal_pose, joint_pos);
    forward_kinematics(manip, joint_sol);
    auto actual_pose = manip.get_tool_pose();
    Array position_diff(3);
    Vec3 goal_rpy;
    Vec3 actual_rpy;
    for (int i = 0; i < 3; i++) {
        position_diff[i] = (goal_pose[i] - actual_pose[i]);
        goal_rpy[i] = goal_pose[i+3];
        actual_rpy[i] = actual_pose[i+3];
    }
    auto goal_rotation = rpy2rotation(goal_rpy);
    auto actual_rotation = rpy2rotation(actual_rpy);
    auto rotation_diff = goal_rotation - actual_rotation;
    CHECK(is_close(goal_pose, actual_pose));
}

TEST_CASE("Inverse kinematics randomized test", "[Manipulator]") {
using namespace blast;
    auto manip = get_generic_Link6();
    Array goal_pose(6);
    Array joint_pos(6);

    goal_pose = {get_random(), get_random(), get_random(), PI*get_random(), PI*get_random(), PI*get_random()};
    fill_random(joint_pos, PI);

    auto joint_sol = inverse_kinematics(manip, goal_pose, joint_pos);
    forward_kinematics(manip, joint_sol);
    CHECK(is_close(goal_pose, manip.get_tool_pose()));
}*/
