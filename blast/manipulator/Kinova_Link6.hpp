#pragma once
#include <blast>

namespace blast {

inline Manipulator make_Kinova_Link6() {
  // Manipulator
  u32 joints = 6;

  // limits
  ManipulatorLimits limits;
  limits.position_max = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL};                               // rad
  limits.position_min = -limits.position_max;                                                                       // rad
  limits.velocity_max = {(real) 3.4907, (real) 3.4907, (real) 3.4907, (real) 5.5851, (real) 5.5851, (real) 5.5851}; // rad/s
  // 10.472 equals to the value of 600 deg/s^2 in rad/s^2
  limits.acceleration_max = {(real) 10.472, (real) 10.472, (real) 10.472, (real) 10.472, (real) 10.472, (real) 10.472}; // rad/s^2
  limits.torque_max       = {(real) 210, (real) 210, (real) 210, (real) 100, (real) 100, (real) 100};                   // Nm
  limits.tool_speed_max   = (real) 1;                                                                                   // m/s

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.joint_offsets = {
          Vec3{(real) 0.11024, -(real) 0.06926, -(real) 0.1375},
          {0, (real) 0.4850, 0},
          {0, -(real) 0.15216, -(real) 0.0917},
          {0, -(real) 0.06296, -(real) 0.22275},
          {(real) 0.08703, (real) 0.0860, -(real) 0.07692},
          {0, 0, -(real) 0.0920}}; // vector to next joint
  kinematics.joint_axes = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.first_joint_position = {0.0, 0.0, (real) 0.0530};
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
          (real) 4.8257,
          (real) 5.9860,
          (real) 3.4159,
          (real) 2.0849,
          (real) 2.0076,
          (real) 1.5193}; // link masses
  dynamics.inertia_tensors = {
          Mat3{(real) 0.0192746, -(real) 0.00239802, -(real) 0.00896331, -(real) 0.00239802, (real) 0.03087806, (real) 0.0016298, -(real) 0.00896331, (real) 0.0016298, (real) 0.02134949},
          {(real) 0.25899206, -(real) 2.89E-05, -(real) 1.23E-06, -(real) 2.89E-05, (real) 0.01755445, -(real) 0.02128064, -(real) 1.23E-06, -(real) 0.02128064, (real) 0.25291674},
          {(real) 0.01742043, -(real) 3.55E-06, (real) 8.4E-07, -(real) 3.55E-06, (real) 0.01119175, (real) 0.00518163, (real) 8.4E-07, (real) 0.00518163, (real) 0.01212876},
          {(real) 0.02454276, (real) 2.61E-06, (real) 1.799E-05, (real) 2.61E-06, (real) 0.02385702, (real) 0.00315758, (real) 1.799E-05, (real) 0.00315758, (real) 0.00294903},
          {(real) 0.00734684, (real) 0.00124927, -(real) 0.00090156, (real) 0.00124927, (real) 0.00464684, -(real) 0.00236128, -(real) 0.00090156, -(real) 0.00236128, (real) 0.00589508},
          {(real) 0.00390762, -(real) 1.13E-06, (real) 1.16E-06, -(real) 1.13E-06, (real) 0.00390722, -(real) 2.21E-05, (real) 1.16E-06, -(real) 2.21E-05, (real) 0.0013928}}; // Inertial tensors
  dynamics.cog_offsets = {
          Vec3{(real) 0.03930119, -(real) 0.00705889, -(real) 0.08462154},
          {(real) 2.53E-06, (real) 0.18829586, -(real) 0.03988382},
          {(real) 4.64E-06, -(real) 0.02451414, -(real) 0.02997969},
          {-(real) 0.00010793, -(real) 0.01056422, -(real) 0.08091102},
          {(real) 0.01243595, (real) 0.03284165, -(real) 0.04091434},
          {(real) 0, (real) 0.00050624, -(real) 0.00388589}}; // centers of mass

  // Collision model
  ManipulatorCapsules collisions;
  Sphere              base;
  base.center               = {0, 0, 0.0}; // because this is relative to p_base and p_base is {0, 0, 0.053}
  base.radius               = (real) 0.2375;
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
  model_caps.p1          = {0, 0, -(real) 0.065};
  model_caps.p2          = {0, 0, (real) 0.045};
  model_caps.radius      = (real) 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, 0, -(real) 0.08};
  model_caps.p2          = {0, (real) 0.485, -(real) 0.08};
  model_caps.radius      = (real) 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0, -(real) 0.065};
  model_caps.p2          = {0, 0, (real) 0.085};
  model_caps.radius      = (real) 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 4
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, (real) 0.00695, -(real) 0.0917};
  model_caps.p2          = {0, -(real) 0.36805, -(real) 0.0917};
  model_caps.radius      = (real) 0.061;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 5
  model_caps.joint_frame = 4;
  model_caps.p1          = {0, 0, 0};
  model_caps.p2          = {0, 0, -(real) 0.08};
  model_caps.radius      = (real) 0.060;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 6
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, (real) 0.08583};
  model_caps.p2          = {0, 0, -(real) 0.06417};
  model_caps.radius      = (real) 0.060;
  collisions.capsule_list.push_back(model_caps);

  Manipulator link6(joints, limits, kinematics, &dynamics, &collisions);

  return link6;
}

} // namespace blast
