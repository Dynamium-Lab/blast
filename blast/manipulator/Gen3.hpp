#pragma once
#include <blast>

namespace blast {

Manipulator make_kinova_gen3() {
  // Manipulator
  u32 joints = 7;

  // limits
  ManipulatorLimits limits;
  limits.position_max = {INF_REAL, 2.25, INF_REAL, 2.58f, INF_REAL, 2.1f, INF_REAL}; // rad
  limits.position_min = -limits.position_max;

  limits.velocity_max = {1.745f, 1.745f, 1.745f, 1.745f, 2.443f, 2.443f, 2.443f};                   // rad/s
                                                                                                    //   limits.vmin = -limits.velocity_max;

  limits.acceleration_max = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad/s^2
                                                                                                    //   limits.amin = -limits.acceleration_max;

  limits.torque_max = {52, 52, 52, 52, 17, 17, 17};                                                 // Nm
                                                                                                    //   limits.tau_min = -limits.torque_max;

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.joint_offsets = {
          Vec3{0.0, 0.0054, -0.1284},
          {0.0, -0.2104, -0.0064},
          {0.0, -0.0064, -0.2104},
          {0.0, -0.2084, -0.0064},
          {0.0, 0.0, -0.1059},
          {0.0, -0.1059, 0.0},
          {0.0, 0.0, -0.0615}}; // vector to next joint
  kinematics.joint_axes = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.static_rotations[0]  = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.static_rotations[1]  = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.static_rotations[2]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.static_rotations[3]  = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.static_rotations[4]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.static_rotations[5]  = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.static_rotations[6]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.first_joint_position = {0, 0, 0.1564f};
  kinematics.base_position        = {0, 0, 0};
  kinematics.base_rotation        = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  // dynamic properties
  ManipulatorDynamics dynamics;

  dynamics.link_masses = {
          1.377f,
          1.1636f,
          1.1636f,
          0.93f,
          0.678f,
          0.678f,
          0.364f}; // link masses
  dynamics.inertia_tensors = {
          Mat3{0.004570f, 0.000001f, 0.000002f, 0.000001f, 0.004831f, 0.000448f, 0.000002f, 0.000448f, 0.001409f},
          {0.011088f, 0.000005f, 0.000000f, 0.000005f, 0.001072f, -0.000691f, 0.000000f, -0.000691f, 0.011255f},
          {0.010932f, 0.000000f, -0.000007f, 0.000000f, 0.011127f, 0.000606f, -0.000007f, 0.000606f, 0.001043f},
          {0.008147f, -0.000001f, 0.000000f, -0.000001f, 0.000631f, -0.000500f, 0.000000f, -0.000500f, 0.008316f},
          {0.001596f, 0.000000f, 0.000000f, 0.000000f, 0.001607f, 0.000256f, 0.000000f, 0.000256f, 0.000399f},
          {0.001641f, 0.000000f, 0.000000f, 0.000000f, 0.000410f, -0.000278f, 0.000000f, -0.000278f, 0.001641f},
          {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, -0.000002f, 0.000001f, -0.000002f, 0.000240f}}; // Inertial tensors
  dynamics.cog_offsets = {
          Vec3{-0.000023f, -0.010364f, -0.073360f},
          {-0.000044f, -0.099580f, -0.013278f},
          {-0.000044f, -0.006641f, -0.117892f},
          {-0.000018f, -0.075478f, -0.015006f},
          {0.000001f, -0.009432f, -0.063883f},
          {0.000001f, -0.045483f, -0.009650f},
          {-0.000093f, 0.000132f, -0.022905f}}; // centers of mass

  // create manipulator gen3
  Manipulator generic_manip(joints, limits, kinematics, &dynamics);

  // Collision model
  ManipulatorCapsules collisions;
  Sphere              sphere;
  sphere.center             = {0, 0, 0.1214};
  sphere.radius             = 0.14;
  collisions.base_sphere    = sphere;
  collisions.collision_base = {0, 0, 1};

  collisions.collision_matrix.resize(3, 3);
  collisions.collision_matrix(2, 0) = 1;

  // Collision model
  CollisionModelCapsule model_caps;

  // Capsule 1
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, 0.035, 0};
  model_caps.p2          = {0, -0.425, 0};
  model_caps.radius      = 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 3;
  model_caps.p1          = {0, 0, -0.025};
  model_caps.p2          = {0, -0.3, -0.01};
  model_caps.radius      = 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, -0.015};
  model_caps.p2          = {0.0, -0.15, -0.015};
  model_caps.radius      = 0.055;
  collisions.capsule_list.push_back(model_caps);

  generic_manip.set_capsules(collisions);

  return generic_manip;
}
} // namespace blast
