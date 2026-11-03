#pragma once
#include <blast>

namespace blast {

inline Manipulator make_UR5e() {
  // Manipulator
  u32 joints = 6;

  // limits
  ManipulatorLimits limits;
  limits.position_max     = {(real) 6.283200, (real) 6.283200, (real) 3.141600, (real) 6.283200, (real) 6.283200, (real) 6.283200};       // rad
  limits.position_min     = {-(real) 6.283200, -(real) 6.283200, -(real) 3.141600, -(real) 6.283200, -(real) 6.283200, -(real) 6.283200}; // rad
  limits.velocity_max     = {(real) 3.141600, (real) 3.141600, (real) 3.141600, (real) 3.141600, (real) 3.141600, (real) 3.141600};       // rad/s
  limits.acceleration_max = {(real) 13.96, (real) 13.96, (real) 13.96, (real) 13.96, (real) 13.96, (real) 13.96};                         // rad/s^2
  limits.torque_max       = {(real) 150, (real) 150, (real) 150, (real) 28, (real) 28, (real) 28};                                        // Nm
  limits.tool_speed_max   = (real) 1;                                                                                                     // m/s

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.joint_offsets = {
          Vec3{0, 0, 0},
          {-(real) 0.425, 0, 0},
          {-(real) 0.3922, 0, (real) 0.1333},
          {0, -(real) 0.0997, 0},
          {0, (real) 0.0996, 0},
          {0, 0, 0}}; // vector to next joint
  kinematics.joint_axes = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.base_position        = {0.0, 0.0, 0.0};
  kinematics.base_rotation        = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.static_rotations[0]  = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.static_rotations[1]  = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.static_rotations[2]  = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.static_rotations[3]  = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.static_rotations[4]  = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.static_rotations[5]  = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.first_joint_position = {0, 0, (real) 0.162500};

  // dynamic properties
  ManipulatorDynamics dynamics;
  dynamics.link_masses = {
          (real) 3.700000,
          (real) 8.393000,
          (real) 2.275000,
          (real) 1.219000,
          (real) 1.219000,
          (real) 0.187900}; // link masses
  dynamics.inertia_tensors = {
          Mat3{(real) 0.0103, 0, 0, 0, (real) 0.0103, 0, 0, 0, (real) 0.0067},
          {(real) 0.1339, 0, 0, 0, (real) 0.1339, 0, 0, 0, (real) 0.0151},
          {(real) 0.0312, 0, 0, 0, (real) 0.0312, 0, 0, 0, (real) 0.0041},
          {(real) 0.0026, 0, 0, 0, (real) 0.0026, 0, 0, 0, (real) 0.0022},
          {(real) 0.0026, 0, 0, 0, (real) 0.0026, 0, 0, 0, (real) 0.0022},
          {(real) 0.0001, 0, 0, 0, (real) 0.0001, 0, 0, 0, (real) 0.0001}}; // Inertial tensors
  dynamics.cog_offsets = {
          Vec3{0, 0, 0},
          {-(real) 0.2125, 0, (real) 0.138},
          {-(real) 0.1961, 0, (real) 0.007},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, -(real) 0.0229}}; // centers of mass

  ManipulatorCapsules collisions;

  Sphere base;
  base.center               = {0, 0, (real) 0.0325};
  base.radius               = (real) 0.09;
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
  capsule.p2          = {0, -(real) 0.15, 0};
  capsule.radius      = (real) 0.09;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 1;
  capsule.p1          = {-(real) 0.42, 0, (real) 0.1375};
  capsule.p2          = {0, 0, (real) 0.1375};
  capsule.radius      = (real) 0.06;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {0, 0, (real) 0.02};
  capsule.p2          = {0, 0, (real) 0.18};
  capsule.radius      = (real) 0.065;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {-(real) 0.373156, 0, (real) 0.00850418};
  capsule.p2          = {(real) 0.000440611, (real) 0.000121443, (real) 0.00850418};
  capsule.radius      = (real) 0.05;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, 0, -(real) 0.155};
  capsule.radius      = (real) 0.0425;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, (real) 0.03, 0};
  capsule.p2          = {0, -(real) 0.1, 0};
  capsule.radius      = (real) 0.0425;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 5;
  capsule.p1          = {0, 0, -(real) 0.14};
  capsule.p2          = {0, 0, -(real) 0.01};
  capsule.radius      = (real) 0.038;
  collisions.capsule_list.push_back(capsule);

  Manipulator ur5e(joints, limits, kinematics, &dynamics, &collisions);
  return ur5e;
}

} // namespace blast
