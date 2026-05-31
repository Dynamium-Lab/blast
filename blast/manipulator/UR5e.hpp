#pragma once
#include <blast>

namespace blast {

inline Manipulator make_UR5e() {
  // Manipulator
  u32 joints = 6;

  // limits
  ManipulatorLimits limits;
  limits.position_max     = {6.283200f, 6.283200f, 3.141600f, 6.283200f, 6.283200f, 6.283200f};       // rad
  limits.position_min     = {-6.283200f, -6.283200f, -3.141600f, -6.283200f, -6.283200f, -6.283200f}; // rad
  limits.velocity_max     = {3.141600f, 3.141600f, 3.141600f, 3.141600f, 3.141600f, 3.141600f};       // rad/s
  limits.acceleration_max = {13.96f, 13.96f, 13.96f, 13.96f, 13.96f, 13.96f};                         // rad/s^2
  limits.torque_max       = {150.0f, 150.0f, 150.0f, 28.0f, 28.0f, 28.0f};                            // Nm
  limits.tool_speed_max   = 1.0f;                                                                     // m/s

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.joint_offsets = {
          Vec3{0, 0, 0},
          {-0.425f, 0, 0},
          {-0.3922f, 0, 0.1333f},
          {0, -0.0997f, -0},
          {0, 0.0996f, -0},
          {0, 0, 0}}; // vector to next joint
  kinematics.joint_axes = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.static_rotations[0]  = {-1, 0, 0, -0, -1, 0, 0, -0, 1};
  kinematics.static_rotations[1]  = {1, 0, 0, -0, -0, 1, 0, -1, -0};
  kinematics.static_rotations[2]  = {1, 0, 0, -0, 1, 0, 0, -0, 1};
  kinematics.static_rotations[3]  = {1, 0, 0, -0, 1, 0, 0, -0, 1};
  kinematics.static_rotations[4]  = {1, 0, 0, -0, -0, 1, 0, -1, -0};
  kinematics.static_rotations[5]  = {1, -0, 0, 0, -0, -1, 0, 1, -0};
  kinematics.first_joint_position = {0, 0, 0.162500f};

  // dynamic properties
  ManipulatorDynamics dynamics;
  dynamics.link_masses = {
          3.700000f,
          8.393000f,
          2.275000f,
          1.219000f,
          1.219000f,
          0.187900f}; // link masses
  dynamics.inertia_tensors = {
          Mat3{0.0103f, 0, 0, 0, 0.0103f, 0, 0, 0, 0.0067f},
          {0.1339f, 0, 0, 0, 0.1339f, 0, 0, 0, 0.0151f},
          {0.0312f, 0, 0, 0, 0.0312f, 0, 0, 0, 0.0041f},
          {0.0026f, 0, 0, 0, 0.0026f, 0, 0, 0, 0.0022f},
          {0.0026f, 0, 0, 0, 0.0026f, 0, 0, 0, 0.0022f},
          {0.0001f, 0, 0, 0, 0.0001f, 0, 0, 0, 0.0001f}}; // Inertial tensors
  dynamics.cog_offsets = {
          Vec3{0, 0, 0},
          {-0.2125f, 0, 0.138f},
          {-0.1961f, 0, 0.007f},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, -0.0229f}}; // centers of mass

  ManipulatorCapsules collisions;

  Sphere base;
  base.center               = {0, 0, 0.0325f};
  base.radius               = 0.09f;
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

  CollisionModelCapsule capsule;
  capsule.joint_frame = 0;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, -0.15f, 0};
  capsule.radius      = 0.09f;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 1;
  capsule.p1          = {-0.42f, 0, 0.1375f};
  capsule.p2          = {0, 0, 0.1375f};
  capsule.radius      = 0.06f;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {0, 0, 0.02f};
  capsule.p2          = {0, 0, 0.18f};
  capsule.radius      = 0.065f;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {-0.373156f, 0, 0.00850418f};
  capsule.p2          = {0.000440611f, 0.000121443f, 0.00850418f};
  capsule.radius      = 0.05f;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, 0, -0.155f};
  capsule.radius      = 0.0425f;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, 0.03f, 0};
  capsule.p2          = {0, -0.1f, 0};
  capsule.radius      = 0.0425f;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 5;
  capsule.p1          = {0, 0, -0.14f};
  capsule.p2          = {0, 0, -0.01f};
  capsule.radius      = 0.038f;
  collisions.capsule_list.push_back(capsule);

  Manipulator ur5e(joints, limits, kinematics, &dynamics, &collisions);
  return ur5e;
}

} // namespace blast
