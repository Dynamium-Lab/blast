#pragma once
#include "blast"

namespace blast {


// todo: Add camera capsule
inline Manipulator get_generic_Link6_fixed() {
  // Manipulator
  u32 joints = 6;
  // limits
  ManipulatorLimits limits;
  limits.pmax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad
  limits.pmin = -limits.pmax;

  limits.vmax = {0.9 * 3.4907f, 0.9 * 3.4907f, 0.9 * 3.4907f, 0.9 * 5.5851f, 0.9 * 5.5851f, 0.9 * 5.5851f};                               // rad/s

  limits.amax = {0.9 * deg2rad(600), 0.9 * deg2rad(600), 0.9 * deg2rad(600), 0.9 * deg2rad(600), 0.9 * deg2rad(600), 0.9 * deg2rad(600)}; // rad/s^2

  limits.tau_max = {210, 210, 210, 100, 100, 100};                                                                                        // Nm

  limits.tcp_max = 2.0;

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.dv = {
          Vec3{0.11024, -0.06926, -0.1375},
          {0.0, 0.4850, 0.0},
          {0.0, -0.15216, -0.0917},
          {0.0, -0.06296, -0.22275},
          {0.08703, 0.0860, -0.07692},
          {0.0, 0.0, -0.0920}}; // vector to next joint
  kinematics.ev = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.p_j0   = {0.0, 0.0, 0.0530f};
  kinematics.p_base = {0.0, 0.0, 0.0};
  kinematics.Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  // kinematics.Q_static.resize(7);
  kinematics.Q_static[0] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.Q_static[1] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[2] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.Q_static[3] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[4] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[5] = {0, 0, 1, 0, 1, 0, -1, 0, 0};
  kinematics.Q_static[6] = {1, 0, 0, 0, -1, 0, 0, 0, -1};

  // dynamic properties
  ManipulatorDynamics dynamics;
  dynamics.m = {
          4.8257f,
          5.9860f,
          3.4159f,
          2.0849f,
          2.0076f,
          1.5193f}; // link masses
  dynamics.I = {
          Mat3{0.0192746f, -0.00239802f, -0.00896331f, -0.00239802f, 0.03087806f, 0.0016298f, -0.00896331f, 0.0016298f, 0.02134949f},
          {0.25899206f, -2.89E-05f, -1.23E-06f, -2.89E-05f, 0.01755445f, -0.02128064f, -1.23E-06f, -0.02128064f, 0.25291674f},
          {0.01742043f, -3.55E-06f, 8.4E-07f, -3.55E-06f, 0.01119175f, 0.00518163f, 8.4E-07f, 0.00518163f, 0.01212876f},
          {0.02454276f, 2.61E-06f, 1.799E-05f, 2.61E-06f, 0.02385702f, 0.00315758f, 1.799E-05f, 0.00315758f, 0.00294903f},
          {0.00734684f, 0.00124927f, -0.00090156f, 0.00124927f, 0.00464684f, -0.00236128f, -0.00090156f, -0.00236128f, 0.00589508f},
          {0.00390762f, -1.13E-06f, 1.16E-06f, -1.13E-06f, 0.00390722f, -2.21E-05f, 1.16E-06f, -2.21E-05f, 0.0013928f}}; // Inertial tensors
  dynamics.av = {
          Vec3{0.03930119f, -0.00705889f, -0.08462154f},
          {2.53E-06f, 0.18829586f, -0.03988382f},
          {4.64E-06f, -0.02451414f, -0.02997969f},
          {-0.00010793f, -0.01056422f, -0.08091102f},
          {0.01243595f, 0.03284165f, -0.04091434f},
          {0.0f, 0.00050624f, -0.00388589f}}; // centers of mass

  // capsules & internal collision data
  // Collision model
  ManipulatorCapsules collisions;
  Sphere              base;
  base.c                    = {0, 0, 0.0}; // because this is relative to p_base and p_base is {0, 0, 0.053}
  base.r                    = 0.2375;
  collisions.base_sphere    = base;
  collisions.collision_base = {0, 0, 0, 1, 1, 1};
  // collisions.collision_base = {0, 0, 0, 1, 1, 1, 1};  // includes camera capsule

  collisions.collision_matrix.resize(6, 6);
  // collisions.collision_matrix.resize(7, 7); // includes camera capsule
  collisions.collision_matrix(5, 0) = 1;
  // collisions.collision_matrix(6, 0) = 1; // includes camera capsule
  collisions.collision_matrix(4, 1) = 1;
  collisions.collision_matrix(5, 1) = 1;
  // collisions.collision_matrix(6, 1) = 1; // includes camera capsule

  // collisions.collision_matrix.resize(6, 6);
  // collisions.collision_matrix(5, 0) = 1;
  // collisions.collision_matrix(4, 1) = 1;
  // collisions.collision_matrix(5, 1) = 1;

  // Collision model
  CollisionModelCapsule model_caps;

  // Capsule 1
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, 0, -0.065};
  model_caps.p2          = {0, 0, 0.045};
  model_caps.r           = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 1;
  // model_caps.p1 = {0, 0, -0.065};
  model_caps.p1 = {0, 0, -0.08};
  // model_caps.p2 = {0, 0.485, -0.065};
  model_caps.p2 = {0, 0.485, -0.08};
  model_caps.r  = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0, -0.065};
  model_caps.p2          = {0, 0, 0.085};
  model_caps.r           = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 4
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0.00695, -0.0917};
  model_caps.p2          = {0, -0.36805, -0.0917};
  model_caps.r           = 0.061;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 5
  model_caps.joint_frame = 4;
  model_caps.p1          = {0, 0, 0};
  model_caps.p2          = {0, 0, -0.08};
  model_caps.r           = 0.060;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 6
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, 0.08583};
  model_caps.p2          = {0, 0, -0.06417};
  model_caps.r           = 0.060;
  collisions.capsule_list.push_back(model_caps);

  // // Capsule 7
  // model_caps.joint_frame = 5;
  // model_caps.p1 = {0, 0.02125, -0.007};
  // // model_caps.p2 = {0, 0.02125, -0.013};
  // model_caps.p2 = {0, 0.02125, 0.143};
  // model_caps.r = 0.085;
  // collisions.capsule_list.push_back(model_caps); // todo: change camera capsule

  // create manipulator link6
  Manipulator link6(joints, limits, kinematics, &dynamics, &collisions);

  return link6;
}

// note: p_base and Q_base are set in respect to lab's layout (in Link6 frame)
Manipulator get_generic_gen3_fixed() { // todo: fix capsules, not working
  // Manipulator
  u32 joints = 7;

  // limits
  ManipulatorLimits limits;
  limits.pmax = {INF_REAL, 2.25, INF_REAL, 2.58f, INF_REAL, 2.1f, INF_REAL}; // rad
  limits.pmin = -limits.pmax;

  // limits.vmax = {1.745f, 1.745f, 1.745f, 1.745f, 2.443f, 2.443f, 2.443f}; // rad/s (used for demo 2)
  limits.vmax = {1.39f, 1.39f, 1.39f, 1.39f, 1.22f, 1.22f, 1.22f}; // rad/s (updated for demo 3)

  // limits.amax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad/s^2
  limits.amax = {5.2, 5.2, 5.2, 5.2, 10.0, 10.0, 10.0}; // rad/s^2

  limits.tau_max = {52, 52, 52, 52, 17, 17, 17};        // Nm

  limits.tcp_max = 1.0;                                 // todo: validate with webapp (Bruno from Kinova said 1.0 could be adequate), webapp says 0.5
  // limits.tcp_max = 0.5; // todo: validate with webapp (Bruno from Kinova said 1.0 could be adequate), webapp says 0.5

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.dv = {
          Vec3{0.0, 0.0054, -0.1284},
          {0.0, -0.2104, -0.0064},
          {0.0, -0.0064, -0.2104},
          {0.0, -0.2084, -0.0064},
          {0.0, 0.0, -0.1059},
          {0.0, -0.1059, 0.0},
          {0.0, 0.0, -0.0615 /*- 0.164*/} // todo: why did we subtract this? Ans : We subtracted because we used the gripper
  }; // vector to next joint
  kinematics.ev = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  // kinematics.Q_static.resize(7);
  kinematics.Q_static[0] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.Q_static[1] = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.Q_static[2] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[3] = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.Q_static[4] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[5] = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.Q_static[6] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.p_j0        = {0, 0, 0.1564f};
  kinematics.p_base      = {0, 0, 0};
  kinematics.Q_base      = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  // dynamic properties
  ManipulatorDynamics dynamics;

  dynamics.m = {
          1.377f,
          1.1636f,
          1.1636f,
          0.93f,
          0.678f,
          0.678f,
          0.364f}; // link masses
  dynamics.I = {
          Mat3{0.004570f, 0.000001f, 0.000002f, 0.000001f, 0.004831f, 0.000448f, 0.000002f, 0.000448f, 0.001409f},
          {0.011088f, 0.000005f, 0.000000f, 0.000005f, 0.001072f, -0.000691f, 0.000000f, -0.000691f, 0.011255f},
          {0.010932f, 0.000000f, -0.000007f, 0.000000f, 0.011127f, 0.000606f, -0.000007f, 0.000606f, 0.001043f},
          {0.008147f, -0.000001f, 0.000000f, -0.000001f, 0.000631f, -0.000500f, 0.000000f, -0.000500f, 0.008316f},
          {0.001596f, 0.000000f, 0.000000f, 0.000000f, 0.001607f, 0.000256f, 0.000000f, 0.000256f, 0.000399f},
          {0.001641f, 0.000000f, 0.000000f, 0.000000f, 0.000410f, -0.000278f, 0.000000f, -0.000278f, 0.001641f},
          {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, -0.000002f, 0.000001f, -0.000002f, 0.000240f}}; // Inertial tensors
  dynamics.av = {
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
  sphere.c                  = {0, 0, 0.035};
  sphere.r                  = 0.14;
  collisions.base_sphere    = sphere;
  collisions.collision_base = {0, 0, 1};

  collisions.collision_matrix.resize(3, 3);
  collisions.collision_matrix(2, 0) = 1;
  collisions.collision_matrix(0, 2) = 1;

  // Collision model
  CollisionModelCapsule model_caps;

  // Capsule 1
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, 0.035, 0};
  model_caps.p2          = {0, -0.425, 0};
  model_caps.r           = 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 3;
  model_caps.p1          = {0, 0, -0.025};
  model_caps.p2          = {0, -0.3, -0.01};
  model_caps.r           = 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, -0.015};
  model_caps.p2          = {0.0, -0.15, -0.015};
  model_caps.r           = 0.055;
  collisions.capsule_list.push_back(model_caps);

  generic_manip.set_capsules(collisions);

  return generic_manip;
}

// todo: Change camera capsule
inline host_fn Manipulator get_generic_Link6() {
  // Manipulator
  constexpr u32 joints = 6;
  // limits
  ManipulatorLimits limits;
  limits.pmax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad
  limits.pmin = -limits.pmax;

  limits.vmax = {3.4907f, 3.4907f, 3.4907f, 5.5851f, 5.5851f, 5.5851f};                               // rad/s
                                                                                                      //   limits.vmin = -limits.vmax;

  limits.amax = {deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600)}; // rad/s^2
                                                                                                      //   limits.amin = -limits.amax;

  limits.tau_max = {210, 210, 210, 100, 100, 100};                                                    // Nm
                                                                                                      //   limits.tau_min = -limits.tau_max;

  limits.tcp_max = 2.0;

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.dv = {
          Vec3{0.11024, -0.06926, -0.1375},
          Vec3{0.0, 0.4850, 0.0},
          Vec3{0.0, -0.15216, -0.0917},
          Vec3{0.0, -0.06296, -0.22275},
          Vec3{0.08703, 0.0860, -0.07692},
          Vec3{0.0, 0.0, -0.0920}}; // vector to next joint
  kinematics.ev = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.p_j0        = {0.0, 0.0, 0.0530f};
  kinematics.p_base      = {0.0, 0.0, 0.0};
  kinematics.Q_base      = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  kinematics.Q_static[0] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.Q_static[1] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[2] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.Q_static[3] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[4] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[5] = {0, 0, 1, 0, 1, 0, -1, 0, 0};
  kinematics.Q_static[6] = {1, 0, 0, 0, -1, 0, 0, 0, -1};

  // dynamic properties
  ManipulatorDynamics dynamics;
  dynamics.m = {
          4.8257f,
          5.9860f,
          3.4159f,
          2.0849f,
          2.0076f,
          1.5193f}; // link masses
  dynamics.I = {
          Mat3{0.0192746f, -0.00239802f, -0.00896331f, -0.00239802f, 0.03087806f, 0.0016298f, -0.00896331f, 0.0016298f, 0.02134949f},
          {0.25899206f, -2.89E-05f, -1.23E-06f, -2.89E-05f, 0.01755445f, -0.02128064f, -1.23E-06f, -0.02128064f, 0.25291674f},
          {0.01742043f, -3.55E-06f, 8.4E-07f, -3.55E-06f, 0.01119175f, 0.00518163f, 8.4E-07f, 0.00518163f, 0.01212876f},
          {0.02454276f, 2.61E-06f, 1.799E-05f, 2.61E-06f, 0.02385702f, 0.00315758f, 1.799E-05f, 0.00315758f, 0.00294903f},
          {0.00734684f, 0.00124927f, -0.00090156f, 0.00124927f, 0.00464684f, -0.00236128f, -0.00090156f, -0.00236128f, 0.00589508f},
          {0.00390762f, -1.13E-06f, 1.16E-06f, -1.13E-06f, 0.00390722f, -2.21E-05f, 1.16E-06f, -2.21E-05f, 0.0013928f}}; // Inertial tensors
  dynamics.av = {
          Vec3{0.03930119f, -0.00705889f, -0.08462154f},
          {2.53E-06f, 0.18829586f, -0.03988382f},
          {4.64E-06f, -0.02451414f, -0.02997969f},
          {-0.00010793f, -0.01056422f, -0.08091102f},
          {0.01243595f, 0.03284165f, -0.04091434f},
          {0.0f, 0.00050624f, -0.00388589f}}; // centers of mass

  // capsules & internal collision data
  // Collision model
  ManipulatorCapsules collisions;
  Sphere              base;
  base.c                    = {0, 0, 0.0}; // because this is relative to p_base and p_base is {0, 0, 0.053}
  base.r                    = 0.2375;
  collisions.base_sphere    = base;
  collisions.collision_base = {0, 0, 0, 1, 1, 1, 1};

  collisions.collision_matrix.resize(7, 7);
  collisions.collision_matrix(5, 0) = 1;
  collisions.collision_matrix(6, 0) = 1;
  collisions.collision_matrix(4, 1) = 1;
  collisions.collision_matrix(5, 1) = 1;
  collisions.collision_matrix(6, 1) = 1;

  // collisions.collision_matrix.resize(6, 6);
  // collisions.collision_matrix(5, 0) = 1;
  // collisions.collision_matrix(4, 1) = 1;
  // collisions.collision_matrix(5, 1) = 1;

  // Collision model
  CollisionModelCapsule model_caps;

  // Capsule 1
  model_caps.joint_frame = 1;
  model_caps.p1          = {0, 0, -0.065};
  model_caps.p2          = {0, 0, 0.045};
  model_caps.r           = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 1;
  // model_caps.p1 = {0, 0, -0.065};
  model_caps.p1 = {0, 0, -0.08};
  // model_caps.p2 = {0, 0.485, -0.065};
  model_caps.p2 = {0, 0.485, -0.08};
  model_caps.r  = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0, -0.065};
  model_caps.p2          = {0, 0, 0.085};
  model_caps.r           = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 4
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0.00695, -0.0917};
  model_caps.p2          = {0, -0.36805, -0.0917};
  model_caps.r           = 0.061;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 5
  model_caps.joint_frame = 4;
  model_caps.p1          = {0, 0, 0};
  model_caps.p2          = {0, 0, -0.08};
  model_caps.r           = 0.060;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 6
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, 0.08583};
  model_caps.p2          = {0, 0, -0.06417};
  model_caps.r           = 0.060;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 7
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0.02125, -0.007};
  // model_caps.p2 = {0, 0.02125, -0.013};
  model_caps.p2 = {0, 0.02125, 0.143};
  model_caps.r  = 0.085;
  collisions.capsule_list.push_back(model_caps); // todo: change camera capsule

  // create manipulator link6
  Manipulator link6(joints, limits, kinematics, &dynamics, &collisions);

  return link6;
}

// note: p_base and Q_base are set in respect to lab's layout (in Link6 frame)
inline host_fn Manipulator get_generic_gen3() { // todo: fix capsules, not working
  // Manipulator
  u32 joints = 7;

  // limits
  ManipulatorLimits limits;
  limits.pmax = {INF_REAL, 2.25, INF_REAL, 2.58f, INF_REAL, 2.1f, INF_REAL}; // rad
  // limits.pmax = {INF_REAL, 4.5, INF_REAL, 5.16, INF_REAL, 4.2, INF_REAL}; // rad
  // limits.pmin = {-INF_REAL, 0, -INF_REAL, 0, -INF_REAL, 0, -INF_REAL}; // rad
  limits.pmin = -limits.pmax;

  limits.vmax = {1.745f, 1.745f, 1.745f, 1.745f, 2.443f, 2.443f, 2.443f};               // rad/s
                                                                                        //   limits.vmin = -limits.vmax;

  limits.amax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad/s^2
                                                                                        //   limits.amin = -limits.amax;

  limits.tau_max = {52, 52, 52, 52, 17, 17, 17};                                        // Nm
                                                                                        //   limits.tau_min = -limits.tau_max;

  limits.tcp_max = 0.5;                                                                 // todo: validate with webapp (Bruno from Kinova said 1.0 could be adequate), webapp says 0.5

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.dv = {
          Vec3{0.0, 0.0054, -0.1284},
          {0.0, -0.2104, -0.0064},
          {0.0, -0.0064, -0.2104},
          {0.0, -0.2084, -0.0064},
          {0.0, 0.0, -0.1059},
          {0.0, -0.1059, 0.0},
          {0.0, 0.0, -0.0615 /*- 0.164*/} // todo: why did we subtract this? Ans : We subtracted because we used the gripper
  }; // vector to next joint
  kinematics.ev = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint
  kinematics.Q_static[0] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  kinematics.Q_static[1] = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.Q_static[2] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[3] = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.Q_static[4] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.Q_static[5] = {1, 0, 0, 0, 0, 1, 0, -1, 0};
  kinematics.Q_static[6] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  kinematics.p_j0        = {0, 0, 0.1564f};
  kinematics.p_base      = {1.4, 0, 0};
  kinematics.Q_base      = {-1, 0, 0, 0, -1, 0, 0, 0, 1};

  // dynamic properties
  ManipulatorDynamics dynamics;

  dynamics.m = {
          1.377f,
          1.1636f,
          1.1636f,
          0.93f,
          0.678f,
          0.678f,
          0.364f}; // link masses
  dynamics.I = {
          Mat3{0.004570f, 0.000001f, 0.000002f, 0.000001f, 0.004831f, 0.000448f, 0.000002f, 0.000448f, 0.001409f},
          {0.011088f, 0.000005f, 0.000000f, 0.000005f, 0.001072f, -0.000691f, 0.000000f, -0.000691f, 0.011255f},
          {0.010932f, 0.000000f, -0.000007f, 0.000000f, 0.011127f, 0.000606f, -0.000007f, 0.000606f, 0.001043f},
          {0.008147f, -0.000001f, 0.000000f, -0.000001f, 0.000631f, -0.000500f, 0.000000f, -0.000500f, 0.008316f},
          {0.001596f, 0.000000f, 0.000000f, 0.000000f, 0.001607f, 0.000256f, 0.000000f, 0.000256f, 0.000399f},
          {0.001641f, 0.000000f, 0.000000f, 0.000000f, 0.000410f, -0.000278f, 0.000000f, -0.000278f, 0.001641f},
          {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, -0.000002f, 0.000001f, -0.000002f, 0.000240f}}; // Inertial tensors
  dynamics.av = {
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
  sphere.c = {0, 0, -0.075}; // todo: Verify this sphere
  // sphere.c = {0, 0, -0.1564f};
  sphere.r                  = 0.15;
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
  model_caps.r           = 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 3;
  model_caps.p1          = {0, 0, -0.025};
  model_caps.p2          = {0, -0.3, -0.01};
  model_caps.r           = 0.06;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, -0.015};
  model_caps.p2          = {0.0, -0.15, -0.015};
  model_caps.r           = 0.055;
  collisions.capsule_list.push_back(model_caps);

  generic_manip.set_capsules(collisions);

  return generic_manip;
}

// Warning: Will not do self-collisions
inline host_fn Manipulator get_generic_fanuc_crx25ia() {
  u32 joints = 6;

  ManipulatorLimits limits;
  limits.pmax = {3.1241, 3.1241, 4.6949, 3.2987, 3.1241, 3.9095};
  limits.pmin = {-3.1241, -3.1241, -4.6949, -3.2987, -3.1241, -3.9095};
  limits.vmax = {2.0944, 2.0944, 3.1416, 3.1416, 3.1416, 3.1416};
  //   limits.vmin    = {-2.0944, -2.0944, -3.1416, -3.1416, -3.1416, -3.1416};
  limits.tau_max = {0, 0, 0, 0, 0, 0};
  //   limits.tau_min = {-0, -0, -0, -0, -0, -0};
  limits.tcp_max = 0;

  ManipulatorKinematics kinematics;
  kinematics.dv = {
          Vec3{0, 0, 0},
          {0, -0.95, 0},
          {0, 0, 0},
          {-0.185, 0, -0.75},
          {-0.18, 0, 0},
          {0, 0, 0}};
  kinematics.ev = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}};
  kinematics.Q_static[0] = {1, 0, 0, 0, 1, 0, -0, 0, 1};
  kinematics.Q_static[1] = {1, 0, -0, 0, -3.67321e-06, -1, -0, 1, -3.67321e-06};
  kinematics.Q_static[2] = {-1, 0, 2.65359e-06, -0, 1, 0, -2.65359e-06, -0, -1};
  kinematics.Q_static[3] = {-3.67321e-06, 0, -1, -0, 1, 0, 1, -0, -3.67321e-06};
  kinematics.Q_static[4] = {-3.67321e-06, -0, 1, -0, 1, -0, -1, -0, -3.67321e-06};
  kinematics.Q_static[5] = {-3.67321e-06, 0, -1, -0, 1, 0, 1, -0, -3.67321e-06};
  kinematics.p_j0        = {0, 0, 0};
  kinematics.p_base      = {0, 0, 0};

  ManipulatorDynamics dynamics;
  dynamics.m = {0, 0, 0, 0, 0, 0};
  dynamics.I = {
          Mat3{0, 0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0, 0}};
  dynamics.av = {
          Vec3{0, 0, 0},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, 0}};

  Manipulator generic_manip(joints, limits, kinematics, &dynamics);

  ManipulatorCapsules collisions;

  CollisionModelCapsule capsule;
  capsule.joint_frame = 0;
  capsule.p1          = {-1.93343e-06, -0.0256679, -0.0423751};
  capsule.joint_frame = 0;
  capsule.p2          = {-1.93343e-06, -0.0256679, -0.0253869};
  capsule.r           = 0.15152;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 1;
  capsule.p1          = {-0.000221565, -0.897957, -0.264728};
  capsule.joint_frame = 1;
  capsule.p2          = {-0.000221565, -0.0264647, -0.264728};
  capsule.r           = 0.164111;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {-0.0821085, 2.98396e-06, 0.00941996};
  capsule.joint_frame = 2;
  capsule.p2          = {-0.00950472, 2.98396e-06, 0.00941996};
  capsule.r           = 0.106481;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {-0.00857209, 0, -0.717588};
  capsule.joint_frame = 3;
  capsule.p2          = {-0.00857209, 0, -0.280478};
  capsule.r           = 0.103252;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 4;
  capsule.p1          = {-0.0348868, 4.9632e-05, -0.00215328};
  capsule.joint_frame = 4;
  capsule.p2          = {0.00902324, 4.9632e-05, -0.00215328};
  capsule.r           = 0.0836233;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 5;
  capsule.p1          = {-0.0097164, 0.000276191, 0.0315251};
  capsule.joint_frame = 5;
  capsule.p2          = {0.0172386, 0.000276191, 0.0315251};
  capsule.r           = 0.0562761;
  collisions.capsule_list.push_back(capsule);

  ObjMatrix<u8> collision_matrix(collisions.capsule_list.size(), collisions.capsule_list.size());
  collision_matrix(2, 0)      = 1;
  collision_matrix(3, 0)      = 1;
  collision_matrix(4, 0)      = 1;
  collision_matrix(5, 0)      = 1;
  collision_matrix(3, 1)      = 1;
  collision_matrix(4, 1)      = 1;
  collision_matrix(5, 1)      = 1;
  collisions.collision_matrix = collision_matrix;

  Array collision_base(collisions.capsule_list.size());
  collision_base            = {0, 0, 1, 1, 1};
  collisions.collision_base = collision_base;

  generic_manip.set_capsules(collisions);

  // add end effector properties
  EndEffector ee;
  ee.dv_ee = {0, 0, 0.3};
  ee.Q_ee  = {0, 1, 0, 1, 0, 0, 0, 0, -1};

  generic_manip.add_end_effector(ee);

  return generic_manip;
}

inline host_fn World get_demo2_world() {
  World result;
  result.add_box({0.7, 0, 0.381}, {1.0, 0.75, 0.381}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});       // table 76 cm high
  result.add_box({0.67, -0.1475, 0.96}, {0.35, 0.025, 0.2}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // vertical plate (no coll)
  // result.add_box({0.7, 0, -0.55}, {1.0, 0.75, 0.4}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});                 // table 76 cm high
  // result.add_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // vertical plate (no coll)
  return result;
}

inline host_fn World get_lab_world() {
  World world;
  // Lab environment
  // add_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); //table
  world.add_box({0.7, 0, -0.55}, {1.0, 0.75, 0.4}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});                 // table
  world.add_box({-0.5654, -0.8145, 0.3248}, {0.381, 0.635, 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // computer and screens next to link6
  world.add_box({0.4506, -1.3479, 0.3248}, {0.635, 0.381, 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // computer and screens next to gen3lite
  world.add_box({0.0, 0.0, -0.4091}, {0.1459, 0.2018, 0.4091}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));     // link6 base
  world.add_box({0.6665, 1.1286, 0.0}, {0.508, 0.3175, 1.0457}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));    // UR5
  world.add_box({0, -0.75, 1.0}, {0.10, 0.15, 2.0}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));                // cables

  // DEMO 1 bin 1
  // world.add_box({0.4114, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // bottom
  // world.add_box({0.2114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // front
  // world.add_box({0.6114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // back
  // world.add_box({0.4114, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // left
  // world.add_box({0.4114, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // right

  // // DEMO 1 bin 2
  // world.add_box({0.8716, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // bottom
  // world.add_box({0.6716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // front
  // world.add_box({1.0716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1));  // back
  // world.add_box({0.8716, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // left
  // world.add_box({0.8716, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // right

  // DEMO 1 obstacle
  world.add_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)); // vertical plate (no coll)

  return world;
}

inline host_fn Matrix get_link6_task() {
  // Array pi = {-1.787, -0.370, -1.391, -1.766, 0.120, -0.569, 1.944};
  // Array pf = {-0.574, 0.421, 2.293, -2.161, 0.482, -0.740, -0.176};
  Array pi = deg2rad({-40.445762634277344,
                      -26.876392364501953,
                      83.60868835449219,
                      1.49383544921875,
                      19.951095581054688,
                      -42.22943115234375}); // w1 + 10
  Array pf = deg2rad({51.851436614990234,
                      -13.578636169433594,
                      107.87167358398438,
                      3.6194305419921875,
                      33.133209228515625,
                      51.21833801269531}); // wb1
  Array vi(6);
  Array vf(6);
  Array ai(6);
  Array af(6);

  Matrix task(6, 6);
  for (int i = 0; i < 6; i++) {
    task(i, 0) = pi[i];
    task(i, 1) = vi[i];
    task(i, 2) = ai[i];
    task(i, 3) = pf[i];
    task(i, 4) = vf[i];
    task(i, 5) = af[i];
  }

  return task;
}

inline host_fn Matrix get_gen3_task() {
  Array pi = wrap2pi(deg2rad({215.6, 282.3, 337.6, 325.2, 159.0, 21.7, 323.0}));
  Array pf = wrap2pi(deg2rad({114.3, 319.8, 29.1, 266.1, 243.0, 25.3, 164.3}));
  Array vi(7);
  Array vf(7);
  Array ai(7);
  Array af(7);

  Matrix task(7, 6);
  for (int i = 0; i < 7; i++) {
    task(i, 0) = pi[i];
    task(i, 1) = vi[i];
    task(i, 2) = ai[i];
    task(i, 3) = pf[i];
    task(i, 4) = vf[i];
    task(i, 5) = af[i];
  }

  return task;
}

inline host_fn std::vector<Matrix> get_Link6_demo1_tasks() {
  // Read tasks Demo 1 (13 tasks)
  // Array pos_zero = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Array pos_home = deg2rad({0.00177001953125,
                            -19.997268676757812,
                            110.00163269042969,
                            0.001434326171875,
                            39.99864196777344,
                            0.0013427734375});

  Array pos_wb1 = deg2rad({51.851436614990234,
                           -13.578636169433594,
                           107.87167358398438,
                           3.6194305419921875,
                           33.133209228515625,
                           51.21833801269531});
  Array pos_wb3 = deg2rad({53.879573822021484,
                           -6.246559143066406,
                           121.92112731933594,
                           5.29833984375,
                           34.81269836425781,
                           45.617919921875});
  Array pos_wb4 = deg2rad({62.90712356567383,
                           -2.7591018676757812,
                           119.23771667480469,
                           8.18115234375,
                           28.954452514648438,
                           47.103912353515625});

  Array pos_bb2 = deg2rad({27.443552017211914,
                           -54.93880081176758,
                           28.312477111816406,
                           1.807525634765625,
                           -4.7853546142578125,
                           29.731155395507812});
  Array pos_bb5 = deg2rad({32.459312438964844,
                           -41.956546783447266,
                           48.75336837768555,
                           4.1028289794921875,
                           1.5039520263671875,
                           15.067718505859375});
  Array pos_bb6 = deg2rad({29.169939041137695,
                           -44.1719856262207,
                           55.125301361083984,
                           0.5225830078125,
                           12.9627685546875,
                           25.127288818359375});

  Array pos_w1      = deg2rad({-40.19879913330078,
                               -34.87216567993164,
                               96.39874267578125,
                               1.4193572998046875,
                               41.67781066894531,
                               -41.464874267578125});
  Array pos_w1_10cm = deg2rad({-40.445762634277344,
                               -26.876392364501953,
                               83.60868835449219,
                               1.49383544921875,
                               19.951095581054688,
                               -42.22943115234375});

  Array pos_b2      = deg2rad({-23.84050750732422,
                               -23.07632064819336,
                               121.39236450195312,
                               -0.2214813232421875,
                               52.64436340332031,
                               -23.962554931640625});
  Array pos_b2_10cm = deg2rad({-24.826316833496094,
                               -9.07125473022461,
                               110.22100830078125,
                               1.2542877197265625,
                               30.40911865234375,
                               -27.223480224609375});

  Array pos_w3      = deg2rad({-34.8427734375,
                               -42.51811599731445,
                               81.55162811279297,
                               0.0153350830078125,
                               33.299774169921875,
                               -35.105316162109375});
  Array pos_w3_10cm = deg2rad({-33.93860626220703,
                               -35.56209945678711,
                               70.89104461669922,
                               -0.6509857177734375,
                               16.881134033203125,
                               -31.8890380859375});

  Array pos_w4      = deg2rad({-20.46697998046875,
                               -30.862834930419922,
                               105.23368835449219,
                               1.3421478271484375,
                               43.346038818359375,
                               -22.13421630859375});
  Array pos_w4_10cm = deg2rad({-19.595504760742188,
                               -19.332134246826172,
                               95.92303466796875,
                               -3.042205810546875,
                               23.505020141601562,
                               -17.61700439453125});

  Array pos_b5      = deg2rad({-30.808380126953125,
                               -53.44142532348633,
                               60.08732986450195,
                               0.0470733642578125,
                               20.517074584960938,
                               -33.375030517578125});
  Array pos_b5_10cm = deg2rad({-31.708511352539062,
                               -46.73225784301758,
                               46.94545364379883,
                               1.8437347412109375,
                               3.8740386962890625,
                               -30.906951904296875});

  Array pos_b6      = deg2rad({-16.236770629882812,
                               -40.21072769165039,
                               86.60179901123047,
                               -2.0247955322265625,
                               32.93904113769531,
                               -12.799804687});
  Array pos_b6_10cm = deg2rad({-17.104469299316406,
                               -30.512935638427734,
                               77.69818115234375,
                               0.307220458984375,
                               17.362060546875,
                               -12.498809814453125});

  // For every bloc, the sequence is  +10cm -> bloc -> +10cm -> dbox -> next bloc's +10cm

  std::vector<Matrix> tasks;

  // home to w1_10cm
  Matrix task_home_start(6, 6);

  // w1_10cm to w1 & w1 to w1_10cm & w1_10cm to wb1
  // Matrix task_w1_0(6, 6);
  // Matrix task_w1_1(6, 6);
  Matrix task_w1_wb1(6, 6);

  // wb1 to b2_10cm
  Matrix task_wb1_b2(6, 6);
  // b2_10cm to b2 & b2 to b2_10cm & b2_10cm to bb2
  // Matrix task_b2_0(6, 6);
  // Matrix task_b2_1(6, 6);
  Matrix task_b2_bb2(6, 6);

  // bb2 to w3_10cm
  Matrix task_bb2_w3(6, 6);
  // w3_10cm to w3 & w3 to w3_10cm & w3_10cm to wb3
  // Matrix task_w3_0(6, 6);
  // Matrix task_w3_1(6, 6);
  Matrix task_w3_wb3(6, 6);

  // wb3 to w4_10cm
  Matrix task_wb3_w4(6, 6);
  // w4_10cm to w4 & w4 to w4_10cm & w4_10cm to wb4
  // Matrix task_w4_0(6, 6);
  // Matrix task_w4_1(6, 6);
  Matrix task_w4_wb4(6, 6);

  // wb4 to b5_10cm
  Matrix task_wb4_b5(6, 6);
  // b5_10cm to b5 & b5 to b5_10cm & b5_10cm to bb5
  // Matrix task_b5_0(6, 6);
  // Matrix task_b5_1(6, 6);
  Matrix task_b5_bb5(6, 6);

  // bb5 to b6_10cm
  Matrix task_bb5_b6(6, 6);
  // b6_10cm to b6 & b6 to b6_10cm & b6_10cm to bb6
  // Matrix task_b6_0(6, 6);
  // Matrix task_b6_1(6, 6);
  Matrix task_b6_bb6(6, 6);

  // bb6 to home
  Matrix task_home_finish(6, 6);

  for (u32 i = 0; i < 6; i++) {
    // Home to start
    task_home_start(i, 0) = pos_home[i];
    task_home_start(i, 3) = pos_w1_10cm[i];

    // Bloc 1
    // task_w1_0(i, 0) = pos_w1_10cm[i];
    // task_w1_0(i, 3) = pos_w1[i];
    // task_w1_1(i, 0) = pos_w1[i];
    // task_w1_1(i, 3) = pos_w1_10cm[i];

    task_w1_wb1(i, 0) = pos_w1_10cm[i];
    task_w1_wb1(i, 3) = pos_wb1[i];

    // Bloc 2
    task_wb1_b2(i, 0) = pos_wb1[i];
    task_wb1_b2(i, 3) = pos_b2_10cm[i];

    // task_b2_0(i, 0) = pos_b2_10cm[i];
    // task_b2_0(i, 3) = pos_b2[i];
    // task_b2_1(i, 0) = pos_b2[i];
    // task_b2_1(i, 3) = pos_b2_10cm[i];

    task_b2_bb2(i, 0) = pos_b2_10cm[i];
    task_b2_bb2(i, 3) = pos_bb2[i];

    // Bloc 3
    task_bb2_w3(i, 0) = pos_bb2[i];
    task_bb2_w3(i, 3) = pos_w3_10cm[i];

    // task_w3_0(i, 0) = pos_w3_10cm[i];
    // task_w3_0(i, 3) = pos_w3[i];
    // task_w3_1(i, 0) = pos_w3[i];
    // task_w3_1(i, 3) = pos_w3_10cm[i];

    task_w3_wb3(i, 0) = pos_w3_10cm[i];
    task_w3_wb3(i, 3) = pos_wb3[i];

    // Bloc 4
    task_wb3_w4(i, 0) = pos_wb3[i];
    task_wb3_w4(i, 3) = pos_w4_10cm[i];

    // task_w4_0(i, 0) = pos_w4_10cm[i];
    // task_w4_0(i, 3) = pos_w4[i];
    // task_w4_1(i, 0) = pos_w4[i];
    // task_w4_1(i, 3) = pos_w4_10cm[i];

    task_w4_wb4(i, 0) = pos_w4_10cm[i];
    task_w4_wb4(i, 3) = pos_wb4[i];

    // Bloc 5
    task_wb4_b5(i, 0) = pos_wb4[i];
    task_wb4_b5(i, 3) = pos_b5_10cm[i];

    // task_b5_0(i, 0) = pos_b5_10cm[i];
    // task_b5_0(i, 3) = pos_b5[i];
    // task_b5_1(i, 0) = pos_b5[i];
    // task_b5_1(i, 3) = pos_b5_10cm[i];

    task_b5_bb5(i, 0) = pos_b5_10cm[i];
    task_b5_bb5(i, 3) = pos_bb5[i];

    // Bloc 6
    task_bb5_b6(i, 0) = pos_bb5[i];
    task_bb5_b6(i, 3) = pos_b6_10cm[i];

    // task_b6_0(i, 0) = pos_b6_10cm[i];
    // task_b6_0(i, 3) = pos_b6[i];
    // task_b6_1(i, 0) = pos_b6[i];
    // task_b6_1(i, 3) = pos_b6_10cm[i];

    task_b6_bb6(i, 0) = pos_b6_10cm[i];
    task_b6_bb6(i, 3) = pos_bb6[i];

    // Home to finish
    task_home_finish(i, 0) = pos_bb6[i];
    task_home_finish(i, 3) = pos_home[i];
  }
  // Sequence tasks (13 tasks)
  tasks.push_back(task_home_start); // 0
  // Bloc 1
  // tasks.push_back(task_w1_0);   // 1
  // tasks.push_back(task_w1_1);   // 2
  tasks.push_back(task_w1_wb1); // 3

  tasks.push_back(task_wb1_b2); // 4
  // Bloc 2
  // tasks.push_back(task_b2_0);   // 5
  // tasks.push_back(task_b2_1);   // 6
  tasks.push_back(task_b2_bb2); // 7

  tasks.push_back(task_bb2_w3); // 8
  // Bloc 3
  // tasks.push_back(task_w3_0);   // 9
  // tasks.push_back(task_w3_1);   // 10
  tasks.push_back(task_w3_wb3); // 11

  tasks.push_back(task_wb3_w4); // 12
  // Bloc 4
  // tasks.push_back(task_w4_0);   // 13
  // tasks.push_back(task_w4_1);   // 14
  tasks.push_back(task_w4_wb4); // 15

  tasks.push_back(task_wb4_b5); // 16
  // Bloc 5
  // tasks.push_back(task_b5_0);   // 17
  // tasks.push_back(task_b5_1);   // 18
  tasks.push_back(task_b5_bb5); // 19

  tasks.push_back(task_bb5_b6); // 20
  // Bloc 6
  // tasks.push_back(task_b6_0);        // 21
  // tasks.push_back(task_b6_1);        // 22
  tasks.push_back(task_b6_bb6);      // 23

  tasks.push_back(task_home_finish); // 24

  return tasks;
}

// HierarchicalTask get_demo2_task() {
//     // Initialization
//     MultiRobotTask task;
//
//     // --- Define manipulators ---
//     task.robots.reserve(2);
//     task.add_robot(get_generic_Link6());
//     task.add_robot(get_generic_gen3());
//
//     // --- Define world ---
//     World world;
//     add_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1}); // vertical plate
//     // add_box({0.7, 0.0, -0.8}, {1.0, 0.75, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1}); //table
//     // task.world = get_lab_world();
//     task.world = world;
//
//
//     // --- Define tasks ---
//     // Gripper toggles
//     Matrix gripper_toggle(13, 2);
//
//     // Link6 positions
//     Array pos_0_link6(6);
//
//     Array pos_w1_10cm = deg2rad({-40.4, -26.9, 83.6, 1.5, 20.0, -42.2});
//     Array pos_w1 = deg2rad({-40.2, -34.9, 96.4, 1.4, 41.7, -41.5});
//     Array pos_wb1 = deg2rad({51.9, -13.6, 107.9, 3.6, 33.1, 51.2});
//
//     Array pos_b2_10cm_link6 = deg2rad({-24.8, -9.1, 110.2, 1.254, 30.4, -27.2});
//     Array pos_b2_link6 = deg2rad({-23.8, -23.1, 121.4, -0.2, 52.6, -24.0});
//     Array pos_bb2_link6 = deg2rad({27.4, -54.9, 28.3, 1.8, -4.8, 29.7});
//
//     Array pos_w3_10cm = deg2rad({-33.9, -35.6, 70.9, -0.7, 16.9, -31.8});
//     Array pos_w3 = deg2rad({-34.8, -42.5, 81.6, 0.0, 33.3, -35.1});
//     Array pos_wb3 = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});
//
//     Array pos_w4_10cm = deg2rad({-20.0, -19.3, 95.9, -3.0, 23.5, -17.6});
//     Array pos_w4 = deg2rad({-20.5, -30.9, 105.2, 1.3, 43.3, -22.1});
//     Array pos_wb4 = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});
//
//     // Gen3 positions
//     Array pos_0_gen3(7);
//
//     Array pos_b2_10cm_gen3 = wrap2pi(deg2rad({36.1,  79.9,  114.8, 339.7, 41.9, 321.6, 116.2}));
//     Array pos_b2_gen3 = wrap2pi(deg2rad(     {37.6,  89.1,  114.6, 335.7, 42.7, 332.0, 119.4}));
//     Array pos_bb2_gen3 = wrap2pi(deg2rad(    {315.3, 33.9,  162.4, 264.0, 84.6, 306.9, 334.0}));
//
//     Array pos_b5_10cm = wrap2pi(deg2rad({72.4, 98.5,  98.5, 318.8, 42.0, 315.7, 135.7}));
//     Array pos_b5 = wrap2pi(deg2rad(     {73.3, 88.8,  98.4, 318.0, 43.4, 321.6, 142.4}));
//     Array pos_bb5 = wrap2pi(deg2rad(    {315.3, 33.9,  162.4, 264.0, 84.6, 306.9, 334.0}));
//
//     Array pos_b6_10cm = wrap2pi(deg2rad({69.2,  59.7,  122.7, 274.6, 55.1, 352.1, 87.7}));
//     Array pos_b6 = wrap2pi(deg2rad(     {62.2,  64.9,  142.5, 270.5, 53.1, 347.5,  88.7}));
//     Array pos_bb6 = wrap2pi(deg2rad(    {315.3, 33.9,  162.4, 264.0, 84.6, 306.9, 334.0}));
//
//     Array v0(6);
//     Array vf(6);
//     Array a0(6);
//     Array af(6);
//
//     // Task 1 Link6
//     task.add_angular_task(0, pos_0_link6, pos_w1_10cm, v0, vf, a0, af);
//     task.add_angular_task(0, pos_w1_10cm, pos_w1, v0, vf, a0, af);
//     gripper_toggle(1, 0) = 1;
//     task.add_angular_task(0, pos_w1, pos_w1_10cm, v0, vf, a0, af);
//     task.add_angular_task(0, pos_w1_10cm, pos_wb1, v0, vf, a0, af);
//     gripper_toggle(3, 0) = 1;
//
//     // Task 2 Link6
//     // task.add_angular_task(0, pos_wb1, pos_b2_10cm_link6, v0, vf, a0, af);
//     // task.add_angular_task(0, pos_b2_10cm_link6, pos_b2_link6, v0, vf, a0, af);
//     // gripper_toggle(5, 0) = 1;
//     // task.add_angular_task(0, pos_b2_link6, pos_b2_10cm_link6, v0, vf, a0, af);
//     // task.add_angular_task(0, pos_b2_10cm_link6, pos_bb2_link6, v0, vf, a0, af);
//     // gripper_toggle(7, 0) = 1;
//
//     // Task 3 Link6
//     // task.add_angular_task(0, pos_bb2, pos_w3_10cm, v0, vf, a0, af);
//     task.add_angular_task(0, pos_wb1, pos_w3_10cm, v0, vf, a0, af);
//     task.add_angular_task(0, pos_w3_10cm, pos_w3, v0, vf, a0, af);
//     gripper_toggle(5, 0) = 1;
//     task.add_angular_task(0, pos_w3, pos_w3_10cm, v0, vf, a0, af);
//     task.add_angular_task(0, pos_w3_10cm, pos_wb3, v0, vf, a0, af); // did not work after 20 restarts
//     gripper_toggle(7, 0) = 1;
//
//     // Task 4 Link6
//     task.add_angular_task(0, pos_wb3, pos_w4_10cm, v0, vf, a0, af);
//     task.add_angular_task(0, pos_w4_10cm, pos_w4, v0, vf, a0, af); // task is not valid (did not work after 5 restarts)
//     gripper_toggle(9, 0) = 1;
//     task.add_angular_task(0, pos_w4, pos_w4_10cm, v0, vf, a0, af);
//     task.add_angular_task(0, pos_w4_10cm, pos_wb4, v0, vf, a0, af);
//     gripper_toggle(11, 0) = 1;
//
//     // // Back to home
//     task.add_angular_task(0, pos_wb4, pos_0_link6, v0, vf, a0, af);
//     // task.add_angular_task(0, pos_w1_10cm, get_Link6_home(), v0, vf, a0, af);
//
//
//     // Gen3 tasks
//     v0.resize(7);
//     vf.resize(7);
//     a0.resize(7);
//     af.resize(7);
//
//     // Task 1 Gen3
//     task.add_angular_task(1, pos_0_gen3, pos_b2_10cm_gen3, v0, vf, a0, af);
//     task.add_angular_task(1, pos_b2_10cm_gen3, pos_b2_gen3, v0, vf, a0, af); // task does not work (pos_b2 is not right)
//     gripper_toggle(1, 1) = 1;
//     task.add_angular_task(1, pos_b2_gen3, pos_b2_10cm_gen3, v0, vf, a0, af);
//     task.add_angular_task(1, pos_b2_10cm_gen3, pos_bb2_gen3, v0, vf, a0, af);
//     gripper_toggle(3, 1) = 1;
//
//     // Task 2 Gen3
//     task.add_angular_task(1, pos_bb2_gen3, pos_b5_10cm, v0, vf, a0, af);
//     task.add_angular_task(1, pos_b5_10cm, pos_b5, v0, vf, a0, af); // task does not work (pos_b5 is not right)
//     gripper_toggle(5, 1) = 1;
//     task.add_angular_task(1, pos_b5, pos_b5_10cm, v0, vf, a0, af);
//     task.add_angular_task(1, pos_b5_10cm, pos_bb5, v0, vf, a0, af);
//     gripper_toggle(7, 1) = 1;
//
//     // Task 3 Gen3
//     task.add_angular_task(1, pos_bb5, pos_b6_10cm, v0, vf, a0, af);
//     task.add_angular_task(1, pos_b6_10cm, pos_b6, v0, vf, a0, af); // task does not work (pos_b6 is not right)
//     gripper_toggle(9, 1) = 1;
//     task.add_angular_task(1, pos_b6, pos_b6_10cm, v0, vf, a0, af);
//     task.add_angular_task(1, pos_b6_10cm, pos_bb6, v0, vf, a0, af);
//     gripper_toggle(11, 1) = 1;
//
//     // Back to home
//     task.add_angular_task(1, pos_bb6, pos_0_gen3, v0, vf, a0, af);
//     // task.add_angular_task(1, pos_b2_10cm, get_gen3_home(), v0, vf, a0, af);
//
//     task.gripper_toggle = gripper_toggle;
//
//     // B-spline
//     Bspline bspline_link6(10, 50, 5, 6);
//     Bspline bspline_gen3(8, 200, 5, 7);
//     std::vector<Bspline> robot_bspline;
//     robot_bspline.reserve(2);
//     robot_bspline.push_back(bspline_link6);
//     robot_bspline.push_back(bspline_gen3);
//
//     // Guess
//     Guess guess; // using default guess
//     guess.n_shot = 300;
//     std::vector<Guess> robot_guess;
//     robot_guess.push_back(guess);
//     robot_guess.push_back(guess);
//
//     // Constraints
//     Constraints<Manipulator> constraints;
//     constraints.position = true;
//     constraints.velocity = true;
//     constraints.acceleration = true;
//     constraints.torque = true;
//     constraints.tcp_speed = true;
//     constraints.self_collisions = true;
//     constraints.external_collisions = true;
//     constraints.n_collision_constraints = 100;
//     constraints.n_collision_skip = 2;
//
//     std::vector<Constraints<Manipulator>> robot_constraints;
//     robot_constraints.push_back(constraints);
//     robot_constraints.push_back(constraints);
//
//     // Objective
//     Objective<Manipulator> objective;
//     objective.K_time = 1;
//
//     std::vector<Objective<Manipulator>> robot_objective;
//     robot_objective.push_back(objective);
//     robot_objective.push_back(objective);
//
//     // --- Setup Optimization ---
//     Array hierarchy = {0, 1};
//     auto hierarchy_optim = create_hierarchical_object(task, hierarchy, robot_constraints, robot_objective, robot_guess, robot_bspline);
//     return hierarchy_optim;
// }

inline host_fn Optimization get_generic_link6_opt() {
  // Manip
  Manipulator manip = get_generic_Link6();

  // Task
  auto task = get_link6_task();

  // Create optimization
  Optimization opt(manip, task);

  // world
  World world = get_lab_world();
  opt.set_world(world);

  // Constraints
  ConstraintSelection constraints;
  constraints.position            = true;
  constraints.velocity            = true;
  constraints.acceleration        = true;
  constraints.torque              = true;
  constraints.tcp_speed           = true;
  constraints.self_collisions     = true;
  constraints.external_collisions = true;
  opt.set_constraints(constraints);

  // Objective
  Objective objective;
  objective.K_time = 1;
  opt.set_objective(objective);

  // B-spline
  Bspline bspline(manip.n_joints);
  opt.set_bspline(bspline);

  // Guess
  opt.guess.type   = Guess::random;
  opt.guess.n_shot = 100;

  return opt;
}

// Optimization<Link6> get_hardcoded_link6_opt() {
//     // Manip
//     Link6 manip;
//
//     // Task
//     auto task = get_link6_task();
//
//     // Create optimization
//     Optimization<Link6> opt(manip, task);
//
//     // world
//     World world = get_lab_world();
//     opt.set_world(world);
//
//     // Constraints
//     Constraints<Link6> constraints;
//     constraints.position = true;
//     constraints.velocity = true;
//     constraints.acceleration = true;
//     constraints.torque = true;
//     constraints.tcp_speed = true;
//     constraints.self_collisions = true;
//     constraints.external_collisions = true;
//     opt.set_constraints(constraints);
//
//     // Objective
//     Objective<Link6> objective;
//     objective.K_time = 1;
//     opt.set_objective(objective);
//
//     // B-spline
//     Bspline bspline(manip.joints);
//     opt.set_bspline(bspline);
//
//     // Guess
//     opt.guess.type = Guess::random;
//     opt.guess.n_shot = 100;
//
//     return opt;
// }

inline host_fn Optimization get_generic_gen3_opt() {
  // Manip
  Manipulator manip = get_generic_gen3();

  // Task
  auto task = get_gen3_task();

  // Create optimization
  Optimization opt(manip, task);

  // world
  World world = get_lab_world();
  opt.set_world(world);

  // Constraints
  ConstraintSelection constraints;
  constraints.position            = true;
  constraints.velocity            = true;
  constraints.acceleration        = true;
  constraints.torque              = true;
  constraints.tcp_speed           = true;
  constraints.self_collisions     = true;
  constraints.external_collisions = true;
  opt.set_constraints(constraints);

  // Objective
  Objective objective;
  objective.K_time = 1;
  opt.set_objective(objective);

  // B-spline
  Bspline bspline(manip.n_joints);
  opt.set_bspline(bspline);

  // Guess
  opt.guess.type   = Guess::random;
  opt.guess.n_shot = 100;

  return opt;
}

// Optimization get_gen3_gen3_opt() {
//     // Manip
//     Gen3 manip;
//     manip.p_base  = {1.4, 0, 0.1564f};
//     manip.Q_base = {-1, 0, 0, 0, -1, 0, 0, 0, 1};
//
//     // Task
//     auto task = get_gen3_task();
//
//     // Create optimization
//     Optimization opt(manip, task);
//
//     // world
//     World world = get_lab_world();
//     opt.set_world(world);
//
//     // Constraints
//     ConstraintSelection constraints;
//     constraints.position = true;
//     constraints.velocity = true;
//     constraints.acceleration = true;
//     constraints.torque = true;
//     constraints.tcp_speed = true;
//     constraints.self_collisions = true;
//     constraints.external_collisions = true;
//     opt.set_constraints(constraints);
//
//     // Objective
//     Objective objective;
//     objective.K_time = 1;
//     opt.set_objective(objective);
//
//     // B-spline
//     Bspline bspline(manip.joints);
//     opt.set_bspline(bspline);
//
//     // Guess
//     opt.guess.type = Guess::random;
//     opt.guess.n_shot = 100;
//
//     return opt;
// }

inline host_fn Array get_Link6_home() {
  return {0, -0.349, 1.92, 0, 0.698, 0};
}

// Case D1 from paper "Singularity Analysis of Kinova’s Link 6 Robot Arm via Grassmann Line Geometry"
// (https://espace2.etsmtl.ca/id/eprint/28453/1/Bonev-I-2024-28453.pdf)
inline host_fn Array get_Link6_singularity() {
  return deg2rad({0, 17.19, 22.25, 62.10, 160, 0});
}

inline host_fn Array get_gen3_home() {
  return wrap2pi(deg2rad(Array({0, 15, 180, -130, 0, 55, 90})));
}


} // namespace blast
