#pragma once
// #include <world/scenes.hpp>
#include "blast"

namespace blast {

inline host_fn Manipulator get_generic_Link6() {
  // Manipulator
  constexpr u32 joints = 6;
  // limits
  ManipulatorLimits limits;
  limits.position_max = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad
  limits.position_min = -limits.position_max;

  limits.velocity_max = {3.4907f, 3.4907f, 3.4907f, 5.5851f, 5.5851f, 5.5851f};                                   // rad/s
                                                                                                                  //   limits.vmin = -limits.velocity_max;

  limits.acceleration_max = {deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600)}; // rad/s^2
                                                                                                                  //   limits.amin = -limits.acceleration_max;

  limits.torque_max = {210, 210, 210, 100, 100, 100};                                                             // Nm
                                                                                                                  //   limits.tau_min = -limits.torque_max;

  limits.tool_speed_max = 2.0;

  // kinematic properties
  ManipulatorKinematics kinematics; // using default Q_base
  kinematics.joint_offsets = {
          Vec3{0.11024, -0.06926, -0.1375},
          Vec3{0.0, 0.4850, 0.0},
          Vec3{0.0, -0.15216, -0.0917},
          Vec3{0.0, -0.06296, -0.22275},
          Vec3{0.08703, 0.0860, -0.07692},
          Vec3{0.0, 0.0, -0.0920}}; // vector to next joint
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
  kinematics.static_rotations[6]  = {1, 0, 0, 0, -1, 0, 0, 0, -1};

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

  // capsules & internal collision data
  // Collision model
  ManipulatorCapsules collisions;
  Sphere              base;
  base.center               = {0, 0, 0.0}; // because this is relative to p_base and p_base is {0, 0, 0.053}
  base.radius               = 0.2375;
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
  model_caps.radius      = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 2
  model_caps.joint_frame = 1;
  // model_caps.p1 = {0, 0, -0.065};
  model_caps.p1 = {0, 0, -0.08};
  // model_caps.p2 = {0, 0.485, -0.065};
  model_caps.p2     = {0, 0.485, -0.08};
  model_caps.radius = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 3
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0, -0.065};
  model_caps.p2          = {0, 0, 0.085};
  model_caps.radius      = 0.065;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 4
  model_caps.joint_frame = 2;
  model_caps.p1          = {0, 0.00695, -0.0917};
  model_caps.p2          = {0, -0.36805, -0.0917};
  model_caps.radius      = 0.061;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 5
  model_caps.joint_frame = 4;
  model_caps.p1          = {0, 0, 0};
  model_caps.p2          = {0, 0, -0.08};
  model_caps.radius      = 0.060;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 6
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0, 0.08583};
  model_caps.p2          = {0, 0, -0.06417};
  model_caps.radius      = 0.060;
  collisions.capsule_list.push_back(model_caps);

  // Capsule 7
  model_caps.joint_frame = 5;
  model_caps.p1          = {0, 0.02125, -0.007};
  // model_caps.p2 = {0, 0.02125, -0.013};
  model_caps.p2     = {0, 0.02125, 0.143};
  model_caps.radius = 0.085;
  collisions.capsule_list.push_back(model_caps); // todo: change camera capsule

  // create manipulator link6
  Manipulator link6(joints, limits, kinematics, &dynamics, &collisions);

  return link6;
}

inline host_fn Task get_link6_task() {
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
  return Task::stop_to_stop(pi, pf);
}

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
  constraints.tool_speed          = true;
  constraints.self_collisions     = true;
  constraints.external_collisions = true;
  opt.set_constraints(constraints);

  // Objective
  Objective objective;
  objective.time_weight = 1;
  opt.set_objective(objective);

  // B-spline
  Bspline bspline(manip.n_joints);
  opt.set_bspline(bspline);

  // Guess
  opt.guess.type           = Guess::random;
  opt.guess.n_random_shots = 100;

  return opt;
}

inline host_fn std::vector<Task> get_UR5e_tasks() {
  std::vector<Array> start_pos;
  start_pos.push_back({1.94822, 0.473555, -0.0255247, -0.448375, 0.370356, -3.12883});
  start_pos.push_back({2.4316, 0.0965, 1.1824, 1.8626, -0.8567, 0.0250});
  start_pos.push_back({-0.704868, 2.17507, 0.834305, 0.132228, 2.27799, 0.00328674});
  start_pos.push_back({-0.819005, 2.50686, 0.344732, -2.85083, -2.39379, -2.93992});
  start_pos.push_back({2.04549, 0.675615, 0.137326, -0.814844, 0.449157, -3.08683});
  start_pos.push_back({2.3499, 1.44305, -1.11122, -0.333013, 0.751039, -3.09482});
  start_pos.push_back({2.36232, 1.06656, -0.132829, -0.933749, 0.788944, -3.13519});
  start_pos.push_back({2.22452, 1.15771, -0.467326, -0.690436, 0.641695, -3.12205});

  std::vector<Array> end_pos;
  end_pos.push_back({2.5825, 0.0700, -0.3892, 0.3196, 0.9927, -3.17328});
  end_pos.push_back({0.1706, -2.5340, -1.2649, -2.4857, -1.4118, 2.9011});
  end_pos.push_back({0.3572, 2.9674, 0.4629, 2.8527, -1.2170, 3.1408});
  end_pos.push_back({-0.0436368, -2.65335, 0.593745, 2.05947, -1.61882, -3.12196});
  end_pos.push_back({-2.6558, -0.6703, 0.0599, -2.5311, -2.0486, -0.0127});
  end_pos.push_back({3.12394, -0.553571, -0.509807, 1.06314, 1.52902, -3.10058});
  end_pos.push_back({2.4552, -0.1792, -1.0268, -1.9354, -0.8723, -0.0198});
  end_pos.push_back({2.14686, -0.283198, -0.264971, -2.59327, -0.566506, -0.0169222});

  std::vector<Task> task_list;
  task_list.reserve(64);
  for (auto& s: start_pos)
    for (auto& e: end_pos)
      task_list.push_back(Task::stop_to_stop(s, e));
  return task_list;
}

inline host_fn World get_kitchen_open_doors() {
  World world;
  world.add_box(Vec3{0.427, -1.0847, 1.1974}, Vec3{0.36029, 0.06, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, 0.48791, 1.1974}, Vec3{0.36029, 0.06, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 1.5298}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 1.2074}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.427, -0.29838, 0.86498}, Vec3{0.36029, 0.84629, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.72729, -0.29838, 1.2074}, Vec3{0.06, 0.84629, 0.36242}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{-0.20911, 0.86968, 1.2009}, Vec3{0.032523, 0.42256, 0.36242}, Mat3{-0.70706, -0.70716, 0, -0.70716, -0.70706, 0, 0, 0, 1});
  world.add_box(Vec3{-0.20549, -1.4701, 1.2009}, Vec3{0.037648, 0.42257, 0.36244}, Mat3{-0.70706, 0.70716, 0, 0.70716, -0.70706, 0, 0, 0, 1});
  world.add_box(Vec3{-0.47232, -0.77127, -0.48346}, Vec3{0.023735, 0.41378, 0.41552}, Mat3{4.6327e-05, 0, 1, 0, 1, 0, 1, 0, 4.6327e-05});
  world.add_box(Vec3{0.35698, 0.425, -0.066059}, Vec3{0.41378, 0.8425, 0.43736}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -1.1546, -0.066059}, Vec3{0.41378, 0.03, 0.43736}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, 0.31136}, Vec3{0.41378, 0.38356, 0.06}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, -0.055}, Vec3{0.41378, 0.38356, 0.025}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  world.add_box(Vec3{0.35698, -0.80106, -0.47336}, Vec3{0.41378, 0.38356, 0.03}, Mat3{1, 0, 0, 0, 1, 0, 0, 0, 1});
  return world;
}

inline Manipulator get_generic_ur5e() {
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

  // mirrored just for consistency
  // collisions.collision_matrix(0, 3) = 1;
  // collisions.collision_matrix(0, 4) = 1;
  // collisions.collision_matrix(0, 5) = 1;
  // collisions.collision_matrix(0, 6) = 1;
  //
  // collisions.collision_matrix(1, 4) = 1;
  // collisions.collision_matrix(1, 5) = 1;
  // collisions.collision_matrix(1, 6) = 1;
  //
  // collisions.collision_matrix(3, 6) = 1;

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
