#pragma once
#include <blast>

namespace blast {

inline Manipulator make_Kinova_Link6() {
  // Manipulator
  u32 joints = 6;

  // limits
  ManipulatorLimits limits;
  limits.position_max = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad
  limits.position_min = -limits.position_max;                                         // rad
  limits.velocity_max = {3.4907f, 3.4907f, 3.4907f, 5.5851f, 5.5851f, 5.5851f};       // rad/s
  // 10.472 equals to the value of 600 deg/s^2 in rad/s^2
  limits.acceleration_max = {10.472f, 10.472f, 10.472f, 10.472f, 10.472f, 10.472f}; // rad/s^2
  limits.torque_max       = {210.0f, 210.0f, 210.0f, 100.0f, 100.0f, 100.0f};       // Nm
  limits.tool_speed_max   = 1.0f;                                                   // m/s

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.joint_offsets = {
          Vec3{0.11024f, -0.06926f, -0.1375f},
          {0, 0.4850f, 0},
          {0, -0.15216f, -0.0917f},
          {0, -0.06296f, -0.22275f},
          {0.08703f, 0.0860f, -0.07692f},
          {0, 0, -0.0920f}}; // vector to next joint
  kinematics.joint_axes = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.first_joint_position = {0.0, 0.0, 0.0530f};
  kinematics.base_position        = {0.0, 0.0, 0.0};
  kinematics.base_rotation        = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.static_rotations[0]  = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.static_rotations[1]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.static_rotations[2]  = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.static_rotations[3]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.static_rotations[4]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.static_rotations[5]  = {0, 0, 1, 0, 1, 0, -1, 0, 0};
  kinematics.static_rotations[6]  = {1, 0, 0, 0, -1, 0, 0, 0, -1}; // rotation matrix to tool

  // dynamic properties
  ManipulatorDynamics dynamics;
  dynamics.link_masses = {
          4.8257f,
          5.9860f,
          3.4159f,
          2.0849f,
          2.0076f,
          1.5193f}; // link masses
  dynamics.inertia_tensors = {
          Mat3{0.0192746f, -0.00239802f, -0.00896331f, -0.00239802f, 0.03087806f, 0.0016298f, -0.00896331f, 0.0016298f, 0.02134949f},
          {0.25899206f, -2.89E-05f, -1.23E-06f, -2.89E-05f, 0.01755445f, -0.02128064f, -1.23E-06f, -0.02128064f, 0.25291674f},
          {0.01742043f, -3.55E-06f, 8.4E-07f, -3.55E-06f, 0.01119175f, 0.00518163f, 8.4E-07f, 0.00518163f, 0.01212876f},
          {0.02454276f, 2.61E-06f, 1.799E-05f, 2.61E-06f, 0.02385702f, 0.00315758f, 1.799E-05f, 0.00315758f, 0.00294903f},
          {0.00734684f, 0.00124927f, -0.00090156f, 0.00124927f, 0.00464684f, -0.00236128f, -0.00090156f, -0.00236128f, 0.00589508f},
          {0.00390762f, -1.13E-06f, 1.16E-06f, -1.13E-06f, 0.00390722f, -2.21E-05f, 1.16E-06f, -2.21E-05f, 0.0013928f}}; // Inertial tensors
  dynamics.cog_offsets = {
          Vec3{0.03930119f, -0.00705889f, -0.08462154f},
          {2.53E-06f, 0.18829586f, -0.03988382f},
          {4.64E-06f, -0.02451414f, -0.02997969f},
          {-0.00010793f, -0.01056422f, -0.08091102f},
          {0.01243595f, 0.03284165f, -0.04091434f},
          {0.0f, 0.00050624f, -0.00388589f}}; // centers of mass

  // Collision model
  ManipulatorCapsules collisions;
  Sphere              base;
  base.center               = {0, 0, 0.0}; // because this is relative to p_base and p_base is {0, 0, 0.053}
  base.radius               = 0.2375f;
  collisions.base_sphere    = base;
  collisions.collision_base = {0, 0, 0, 1, 1, 1};

  collisions.collision_matrix.resize(6, 6);
  collisions.collision_matrix(5, 0) = 1;
  collisions.collision_matrix(4, 1) = 1;
  collisions.collision_matrix(5, 1) = 1;

  // Collision model
  CollisionModelCapsule model_caps;

  // Capsule 1
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, 0, -0.065f};
  model_caps.p2          = {0, 0, 0.045f};
  model_caps.radius      = 0.065f;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, 0, -0.08f};
  model_caps.p2          = {0, 0.485f, -0.08f};
  model_caps.radius      = 0.065f;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0, -0.065f};
  model_caps.p2          = {0, 0, 0.085f};
  model_caps.radius      = 0.065f;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 4
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0.00695f, -0.0917f};
  model_caps.p2          = {0, -0.36805f, -0.0917f};
  model_caps.radius      = 0.061f;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 5
  model_caps.joint_frame = 4;
  model_caps.p1          = {0, 0, 0};
  model_caps.p2          = {0, 0, -0.08f};
  model_caps.radius      = 0.060f;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 6
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, 0.08583f};
  model_caps.p2          = {0, 0, -0.06417f};
  model_caps.radius      = 0.060f;
  collisions.capsule_list.push_back(model_caps);

  Manipulator link6(joints, limits, kinematics, &dynamics, &collisions);

  return link6;
}

} // namespace blast
