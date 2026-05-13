
#pragma once
// Universal Robots UR5e definition for Blast examples.
// This header is part of the examples folder — it has no dependency on
// tests/test_helper or any other internal test infrastructure.

#include <blast>

namespace blast {

// Returns a fully configured UR5e Manipulator.
// The robot is placed at the world origin; call-sites may adjust
// manip.base_position and manip.base_rotation as needed.
inline Manipulator make_ur5e() {
  u32 joints = 6;

  // --- Joint limits ---
  ManipulatorLimits limits;
  limits.position_max     = {6.2832, 6.2832, 3.1416, 6.2832, 6.2832, 6.2832}; // rad
  limits.position_min     = {-6.2832, -6.2832, -3.1416, -6.2832, -6.2832, -6.2832};
  limits.velocity_max     = {3.1416, 3.1416, 3.1416, 3.1416, 3.1416, 3.1416}; // rad/s
  limits.acceleration_max = {13.96, 13.96, 13.96, 13.96, 13.96, 13.96};       // rad/s²
  limits.torque_max       = {150.0, 150.0, 150.0, 28.0, 28.0, 28.0};          // N·m
  limits.tool_speed_max   = 2.0;                                              // m/s

  // --- Kinematics (modified DH) ---
  ManipulatorKinematics kinematics;
  kinematics.joint_offsets = {
          Vec3{0, 0, 0},        // joint 0 → joint 1
          {-0.425, 0, 0},       // joint 1 → joint 2
          {-0.3922, 0, 0.1333}, // joint 2 → joint 3
          {0, -0.0997, 0},      // joint 3 → joint 4
          {0, 0.0996, 0},       // joint 4 → joint 5
          {0, 0, 0},            // joint 5 → To
  };
  kinematics.joint_axes = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
  };
  kinematics.first_joint_position = {0.0, 0.0, 0.1625}; // base joint height
  kinematics.static_rotations[0]  = {-1, 0, 0, 0, -1, 0, 0, 0, 1};
  kinematics.static_rotations[1]  = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.static_rotations[2]  = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.static_rotations[3]  = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.static_rotations[4]  = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.static_rotations[5]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};

  // --- Dynamics ---
  ManipulatorDynamics dynamics;
  dynamics.link_masses     = {3.7, 8.393, 2.275, 1.219, 1.219, 0.1879}; // kg
  dynamics.inertia_tensors = {
          Mat3{0.0103, 0, 0, 0, 0.0103, 0, 0, 0, 0.0067},
          {0.1339, 0, 0, 0, 0.1339, 0, 0, 0, 0.0151},
          {0.0312, 0, 0, 0, 0.0312, 0, 0, 0, 0.0041},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0001, 0, 0, 0, 0.0001, 0, 0, 0, 0.0001},
  };
  dynamics.cog_offsets = {
          Vec3{0, 0, 0}, // center of mass offsets
          {-0.2125, 0, 0.138},
          {-0.1961, 0, 0.007},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, -0.0229},
  };

  // --- Collision geometry (7 capsules + 1 base sphere) ---
  ManipulatorCapsules collisions;

  Sphere base;
  base.center               = {0, 0, 0.0325};
  base.radius               = 0.09;
  collisions.base_sphere    = base;
  collisions.collision_base = {0, 0, 0, 1, 1, 1, 1};

  collisions.collision_matrix.resize(7, 7);
  collisions.collision_matrix(3, 0) = 1;
  collisions.collision_matrix(4, 0) = 1;
  collisions.collision_matrix(5, 0) = 1;
  collisions.collision_matrix(6, 0) = 1;
  collisions.collision_matrix(4, 1) = 1;
  collisions.collision_matrix(5, 1) = 1;
  collisions.collision_matrix(6, 1) = 1;
  collisions.collision_matrix(6, 3) = 1;

  CollisionModelCapsule capsule;
  capsule.joint_frame = 0;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, -0.15, 0};
  capsule.radius      = 0.09;
  collisions.capsule_list.push_back(capsule);
  capsule.joint_frame = 1;
  capsule.p1          = {-0.42, 0, 0.1375};
  capsule.p2          = {0, 0, 0.1375};
  capsule.radius      = 0.06;
  collisions.capsule_list.push_back(capsule);
  capsule.joint_frame = 2;
  capsule.p1          = {0, 0, 0.02};
  capsule.p2          = {0, 0, 0.18};
  capsule.radius      = 0.065;
  collisions.capsule_list.push_back(capsule);
  capsule.joint_frame = 2;
  capsule.p1          = {-0.373156, 0, 0.00850418};
  capsule.p2          = {0.000440611, 0.000121443, 0.00850418};
  capsule.radius      = 0.05;
  collisions.capsule_list.push_back(capsule);
  capsule.joint_frame = 3;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, 0, -0.155};
  capsule.radius      = 0.0425;
  collisions.capsule_list.push_back(capsule);
  capsule.joint_frame = 3;
  capsule.p1          = {0, 0.03, 0};
  capsule.p2          = {0, -0.1, 0};
  capsule.radius      = 0.0425;
  collisions.capsule_list.push_back(capsule);
  capsule.joint_frame = 5;
  capsule.p1          = {0, 0, -0.14};
  capsule.p2          = {0, 0, -0.01};
  capsule.radius      = 0.038;
  collisions.capsule_list.push_back(capsule);

  return Manipulator(joints, limits, kinematics, &dynamics, &collisions);
}

// Returns a representative stop-to-stop task for the UR5e.
inline Task make_ur5e_task() {
  Array start = {1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883};
  Array end   = {2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328};
  return Task::stop_to_stop(start, end);
}

} // namespace blast
