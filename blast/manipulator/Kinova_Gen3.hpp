#pragma once
#include <blast>

namespace blast {

inline Manipulator make_Kinova_Gen3() {
  // Manipulator
  u32 joints = 7;

  // limits
  ManipulatorLimits limits;
  limits.position_max     = {INF_REAL, (real) 2.25, INF_REAL, (real) 2.58, INF_REAL, (real) 2.1, INF_REAL};           // rad
  limits.position_min     = -limits.position_max;                                                   // rad
  limits.velocity_max     = {(real) 1.745, (real) 1.745, (real) 1.745, (real) 1.745, (real) 2.443, (real) 2.443, (real) 2.443};               // rad/s
  limits.acceleration_max = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad/s^2
  limits.torque_max       = {(real) 52, (real) 52, (real) 52, (real) 52, (real) 17, (real) 17, (real) 17};                      // Nm
  limits.tool_speed_max   = (real) 1;

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.joint_offsets = {
          Vec3{0, (real) 0.0054, -(real) 0.1284},
          {0, -(real) 0.2104, -(real) 0.0064},
          {0, -(real) 0.0064, -(real) 0.2104},
          {0, -(real) 0.2084, -(real) 0.0064},
          {0, 0, -(real) 0.1059},
          {0, -(real) 0.1059, 0},
          {0, 0, -(real) 0.0615}}; // vector to next joint
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
  kinematics.first_joint_position = {0, 0, (real) 0.1564};
  kinematics.base_position        = {0, 0, 0};
  kinematics.base_rotation        = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  // dynamic properties
  ManipulatorDynamics dynamics;
  dynamics.link_masses = {
          (real) 1.377,
          (real) 1.1636,
          (real) 1.1636,
          (real) 0.93,
          (real) 0.678,
          (real) 0.678,
          (real) 0.364}; // link masses
  dynamics.inertia_tensors = {
          Mat3{(real) 0.004570, (real) 0.000001, (real) 0.000002, (real) 0.000001, (real) 0.004831, (real) 0.000448, (real) 0.000002, (real) 0.000448, (real) 0.001409},
          {(real) 0.011088, (real) 0.000005, (real) 0, (real) 0.000005, (real) 0.001072, -(real) 0.000691, (real) 0, -(real) 0.000691, (real) 0.011255},
          {(real) 0.010932, (real) 0, -(real) 0.000007, (real) 0, (real) 0.011127, (real) 0.000606, -(real) 0.000007, (real) 0.000606, (real) 0.001043},
          {(real) 0.008147, -(real) 0.000001, (real) 0, -(real) 0.000001, (real) 0.000631, -(real) 0.000500, (real) 0, -(real) 0.000500, (real) 0.008316},
          {(real) 0.001596, (real) 0, (real) 0, (real) 0, (real) 0.001607, (real) 0.000256, (real) 0, (real) 0.000256, (real) 0.000399},
          {(real) 0.001641, (real) 0, (real) 0, (real) 0, (real) 0.000410, -(real) 0.000278, (real) 0, -(real) 0.000278, (real) 0.001641},
          {(real) 0.000214, (real) 0, (real) 0.000001, (real) 0, (real) 0.000223, -(real) 0.000002, (real) 0.000001, -(real) 0.000002, (real) 0.000240}}; // Inertial tensors
  dynamics.cog_offsets = {
          Vec3{-(real) 0.000023, -(real) 0.010364, -(real) 0.073360},
          {-(real) 0.000044, -(real) 0.099580, -(real) 0.013278},
          {-(real) 0.000044, -(real) 0.006641, -(real) 0.117892},
          {-(real) 0.000018, -(real) 0.075478, -(real) 0.015006},
          {(real) 0.000001, -(real) 0.009432, -(real) 0.063883},
          {(real) 0.000001, -(real) 0.045483, -(real) 0.009650},
          {-(real) 0.000093, (real) 0.000132, -(real) 0.022905}}; // centers of mass

  // create manipulator gen3
  Manipulator gen3(joints, limits, kinematics, &dynamics);

  // Collision model
  ManipulatorCapsules collisions;
  Sphere              sphere;
  sphere.center             = {0, 0, (real) 0.1214};
  sphere.radius             = (real) 0.14;
  collisions.base_sphere    = sphere;
  collisions.collision_base = {0, 0, 1};

  collisions.collision_matrix.resize(3, 3);
  collisions.collision_matrix(2, 0) = 1;

  // mirrored just for consistency, technically not used
  collisions.collision_matrix(0, 2) = 1;

  // Collision model
  CollisionModelCapsule model_caps;

  // Capsule 1
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, (real) 0.035, 0};
  model_caps.p2          = {0, -(real) 0.425, 0};
  model_caps.radius      = (real) 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 3;
  model_caps.p1          = {0, 0, -(real) 0.025};
  model_caps.p2          = {0, -(real) 0.3, -(real) 0.01};
  model_caps.radius      = (real) 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, -(real) 0.015};
  model_caps.p2          = {0, -(real) 0.15, -(real) 0.015};
  model_caps.radius      = (real) 0.055;
  collisions.capsule_list.push_back(model_caps);

  gen3.set_capsules(collisions);

  return gen3;
}

inline host_fn Array get_gen3_home() {
  return wrap2pi(deg2rad(Array({0, 15, 180, -130, 0, 55, 90})));
}

} // namespace blast
