#pragma once
#include <blast>

namespace blast {

Manipulator make_UR5e() {
  // Manipulator
  blast::u32 joints = 6;

  // limits
  blast::ManipulatorLimits limits;
  limits.position_max = {6.283200, 6.283200, 3.141600, 6.283200, 6.283200, 6.283200};
  limits.position_min = {-6.283200, -6.283200, -3.141600, -6.283200, -6.283200, -6.283200};
  limits.velocity_max = {3.141600, 3.141600, 3.141600, 3.141600, 3.141600, 3.141600};
  // limits.vmin = {-3.141600, -3.141600, -3.141600, -3.141600, -3.141600, -3.141600};
  limits.acceleration_max = {13.96, 13.96, 13.96, 13.96, 13.96, 13.96};
  // limits.amin = {-13.96, -13.96, -13.96, -13.96, -13.96, -13.96};
  limits.torque_max = {150.000000, 150.000000, 150.000000, 28.000000, 28.000000, 28.000000};
  // limits.tau_min = {-150.000000, -150.000000, -150.000000, -28.000000, -28.000000, -28.000000};
  limits.tool_speed_max = 2.0; // todo: verify this value

  blast::ManipulatorKinematics kinematics;
  kinematics.joint_offsets = {
          Vec3{0, 0, 0},
          {-0.425, 0, 0},
          {-0.3922, 0, 0.1333},
          {0, -0.0997, -0},
          {0, 0.0996, -0},
          {0, 0, 0} // to end effector
  };
  kinematics.joint_axes = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}};
  // kinematics.static_rotations.resize(6);
  kinematics.static_rotations[0] = {-1.000000, 0.000000, 0.000000, -0.000000, -1.000000, 0.000000, 0.000000, -0.000000, 1.000000}; // modified
  kinematics.static_rotations[1] = {1.000000, 0.000000, 0.000000, -0.000000, -0.000000, 1.000000, 0.000000, -1.000000, -0.000000};
  kinematics.static_rotations[2] = {1.000000, 0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000000, 1.000000};
  kinematics.static_rotations[3] = {1.000000, 0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000000, 1.000000};
  kinematics.static_rotations[4] = {1.000000, 0.000000, 0.000000, -0.000000, -0.000000, 1.000000, 0.000000, -1.000000, -0.000000};
  kinematics.static_rotations[5] = {1.000000, -0.000000, 0.000000, 0.000000, -0.000000, -1.000000, 0.000000, 1.000000, -0.000000};
  // kinematics.static_rotations[6] = {1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000};
  kinematics.first_joint_position = {0.000000, 0.000000, 0.162500};

  blast::ManipulatorDynamics dynamics;
  dynamics.link_masses     = {3.700000, 8.393000, 2.275000, 1.219000, 1.219000, 0.187900};
  dynamics.inertia_tensors = {
          Mat3{0.0103, 0, 0, 0, 0.0103, 0, 0, 0, 0.0067},
          {0.1339, 0, 0, 0, 0.1339, 0, 0, 0, 0.0151},
          {0.0312, 0, 0, 0, 0.0312, 0, 0, 0, 0.0041},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0001, 0, 0, 0, 0.0001, 0, 0, 0, 0.0001}};
  dynamics.cog_offsets = {
          Vec3{0, 0, 0},
          {-0.2125, 0, 0.138},
          {-0.1961, 0, 0.007},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, -0.0229}};

  blast::ManipulatorCapsules collisions;

  blast::Sphere base;
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

  // mirrored just for consistency, technically not used
  collisions.collision_matrix(0, 3) = 1;
  collisions.collision_matrix(0, 4) = 1;
  collisions.collision_matrix(0, 5) = 1;
  collisions.collision_matrix(0, 6) = 1;

  collisions.collision_matrix(1, 4) = 1;
  collisions.collision_matrix(1, 5) = 1;
  collisions.collision_matrix(1, 6) = 1;

  collisions.collision_matrix(3, 6) = 1;

  blast::CollisionModelCapsule capsule;
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

  blast::Manipulator generic_manip(joints, limits, kinematics, &dynamics, &collisions);
  return generic_manip;
}

} // namespace blast
