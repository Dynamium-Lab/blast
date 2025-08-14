#include <blast>
#include <iostream>

#include "../tests/test_helper/test_helper.hpp"
#include "manipulator/Link6.hpp"

using blast::Array;
using blast::Matrix;

inline void add_task(Matrix& task, Array pos1, Array pos2, Array vel1 = {}, Array vel2 = {}, Array acc1 = {}, Array acc2 = {}) {
  Assert(pos1.size == pos2.size);
  task.resize(pos1.size, 6);
  for (int j = 0; j < task.rows; j++) {
    task(j, 0) = pos1[j];
    task(j, 3) = pos2[j];
    if (vel1.size != 0 && vel2.size != 0) {
      Assert(vel1.size == vel2.size == task.rows);
      task(j, 1) = vel1[j];
      task(j, 4) = vel2[j];
    }
    if (acc1.size != 0 && acc2.size != 0) {
      Assert(acc1.size == acc2.size == task.rows);
      task(j, 2) = acc1[j];
      task(j, 5) = acc2[j];
    }
  }
}

int main() {
  using namespace blast;

  // --- Test characteristics ---
  bool position            = true;
  bool velocity            = true;
  bool acceleration        = true;
  bool torque              = true;
  bool tcp_speed           = true;
  bool self_collisions     = true;
  bool external_collisions = true;




  // setup manip
  // setup world
  //...












  // // Manipulator
  // u32 joints = 6;
  // // limits
  // ManipulatorLimits limits;
  // limits.pmax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad
  // limits.pmin = -limits.pmax;
  //
  // limits.vmax = {3.4907f, 3.4907f, 3.4907f, 5.5851f, 5.5851f, 5.5851f}; // rad/s
  // limits.vmin = -limits.vmax;
  //
  // limits.amax = {deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600)}; // rad/s^2
  // limits.amin = -limits.amax;
  //
  // limits.tau_max = {210, 210, 210, 100, 100, 100}; // Nm
  // limits.tau_min = -limits.tau_max;
  //
  // limits.tcp_max = 2.0;
  //
  // // kinematic properties
  // ManipulatorKinematics kinematics; // using default Q_base
  // kinematics.dv = {
  //         {0.11024, -0.06926, -0.1375},
  //         {0.0, 0.4850, 0.0},
  //         {0.0, -0.15216, -0.0917},
  //         {0.0, -0.06296, -0.22275},
  //         {0.08703, 0.0860, -0.07692},
  //         {0.0, 0.0, -0.0920}}; // vector to next joint
  // kinematics.ev = {
  //         {0, 0, 1},
  //         {0, 0, 1},
  //         {0, 0, 1},
  //         {0, 0, 1},
  //         {0, 0, 1},
  //         {0, 0, 1}}; // direction vectors of joint
  // kinematics.p_j0   = {0.0, 0.0, 0.0530f};
  // kinematics.p_base = {0.0, 0.0, 0.0};
  // kinematics.Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  // kinematics.Q_static.resize(7);
  // kinematics.Q_static[0] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  // kinematics.Q_static[1] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  // kinematics.Q_static[2] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  // kinematics.Q_static[3] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  // kinematics.Q_static[4] = {1, 0, 0, 0, 0, -1, 0, 1, 0};
  // kinematics.Q_static[5] = {0, 0, 1, 0, 1, 0, -1, 0, 0};
  // kinematics.Q_static[6] = {1, 0, 0, 0, -1, 0, 0, 0, -1};
  //
  // // dynamic properties
  // ManipulatorDynamics dynamics;
  // dynamics.m = {
  //         4.8257f,
  //         5.9860f,
  //         3.4159f,
  //         2.0849f,
  //         2.0076f,
  //         1.5193f}; // link masses
  // dynamics.I = {
  //         {0.0192746f, -0.00239802f, -0.00896331f, -0.00239802f, 0.03087806f, 0.0016298f, -0.00896331f, 0.0016298f, 0.02134949f},
  //         {0.25899206f, -2.89E-05f, -1.23E-06f, -2.89E-05f, 0.01755445f, -0.02128064f, -1.23E-06f, -0.02128064f, 0.25291674f},
  //         {0.01742043f, -3.55E-06f, 8.4E-07f, -3.55E-06f, 0.01119175f, 0.00518163f, 8.4E-07f, 0.00518163f, 0.01212876f},
  //         {0.02454276f, 2.61E-06f, 1.799E-05f, 2.61E-06f, 0.02385702f, 0.00315758f, 1.799E-05f, 0.00315758f, 0.00294903f},
  //         {0.00734684f, 0.00124927f, -0.00090156f, 0.00124927f, 0.00464684f, -0.00236128f, -0.00090156f, -0.00236128f, 0.00589508f},
  //         {0.00390762f, -1.13E-06f, 1.16E-06f, -1.13E-06f, 0.00390722f, -2.21E-05f, 1.16E-06f, -2.21E-05f, 0.0013928f}}; // Inertial tensors
  // dynamics.av = {
  //         {0.03930119f, -0.00705889f, -0.08462154f},
  //         {2.53E-06f, 0.18829586f, -0.03988382f},
  //         {4.64E-06f, -0.02451414f, -0.02997969f},
  //         {-0.00010793f, -0.01056422f, -0.08091102f},
  //         {0.01243595f, 0.03284165f, -0.04091434f},
  //         {0.0f, 0.00050624f, -0.00388589f}}; // centers of mass
  //
  // // create manipulator that is a model of the Kinova Link6
  // Manipulator link6(joints, limits, kinematics, &dynamics);

  auto link6 = get_generic_Link6();

  // Collision model
  ManipulatorCapsules collisions;
  Sphere              base;
  base.c                    = {0, 0, 0};
  base.r                    = 0.2375;
  collisions.base_sphere    = base;
  collisions.collision_base = {0, 0, 0, 1, 1, 1, 1};

  collisions.collision_matrix.resize(7, 7);
  collisions.collision_matrix(5, 0) = 1;
  collisions.collision_matrix(6, 0) = 1;
  collisions.collision_matrix(4, 1) = 1;
  collisions.collision_matrix(5, 1) = 1;
  collisions.collision_matrix(6, 1) = 1;

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
  model_caps.p1          = {0, 0, -0.065};
  model_caps.p2          = {0, 0.485, -0.065};
  model_caps.r           = 0.065;
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
  model_caps.p1          = {0, 0.02125, -0.085};
  model_caps.p2          = {0, 0.02125, 0.065};
  model_caps.r           = 0.085;
  collisions.capsule_list.push_back(model_caps);

  link6.set_capsules(collisions);

  // Task
  Matrix task(link6.joints, 6);
  // Array pos_zero = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // Array pos_home = {0, -0.349, 1.92, 0, 0.698, 0};

  // for (u32 j = 0; j < link6.joints; j++) {
  //     task(j, 0) = pos_zero[j];
  //     task(j, 3) = pos_home[j];
  // }

  Array pos_wb1     = deg2rad({51.851436614990234,
                               -13.578636169433594,
                               107.87167358398438,
                               3.6194305419921875,
                               33.133209228515625,
                               51.21833801269531});
  Array pos_w1_10cm = deg2rad({-40.445762634277344,
                               -26.876392364501953,
                               83.60868835449219,
                               1.49383544921875,
                               19.951095581054688,
                               -42.22943115234375});
  for (int i = 0; i < link6.joints; i++) {
    task(i, 0) = pos_wb1[i];
    task(i, 3) = pos_w1_10cm[i];
  }

  // B-spline
  Bspline bspline(link6.joints);

  // Guess
  Guess guess; // using default guess
  guess.type = Guess::custom;
  Array x0(bspline.x_len(task), 1.0);
  guess.x0 = x0;

  // Constraints
  ConstraintSelection constraints;
  constraints.position            = position;
  constraints.velocity            = velocity;
  constraints.acceleration        = acceleration;
  constraints.torque              = torque;
  constraints.tcp_speed           = tcp_speed;
  constraints.self_collisions     = self_collisions;
  constraints.external_collisions = external_collisions;

  // Objective
  Objective objective;
  objective.K_time = 1;

  // World
  World world;
  add_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world);  // vertical plate (no coll)
  add_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // table

  // TestManip Link6
  Optimization opt(link6, task);
  // Guess guess_link6(result_hc.x0);

  opt.set_bspline(bspline);
  opt.set_guess(guess);
  opt.set_world(world);
  opt.set_objective(objective);
  opt.set_constraints(constraints);

  auto result = optimize(&opt);

  cout << "Compute time: " << result.compute_time << endl;
  cout << "Trajectory time: " << result.x[result.x.size - 1] << endl;
  cout << "Trajectory success: " << result.success << endl;
  cout << "Number of restarts: " << result.num_restart << endl;

  return 0;
}
