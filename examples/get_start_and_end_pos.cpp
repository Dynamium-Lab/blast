//
// Created by nikos on 2025-09-15.
//
#include <blast>

using namespace blast;

inline Manipulator get_generic_ur5e() {
  // Manipulator
  blast::u32 joints = 6;

  // limits
  blast::ManipulatorLimits limits;
  limits.pmax = {6.283200, 6.283200, 3.141600, 6.283200, 6.283200, 6.283200};
  limits.pmin = {-6.283200, -6.283200, -3.141600, -6.283200, -6.283200, -6.283200};
  limits.vmax = {3.141600, 3.141600, 3.141600, 3.141600, 3.141600, 3.141600};
  // limits.vmin = {-3.141600, -3.141600, -3.141600, -3.141600, -3.141600, -3.141600};
  limits.amax = {13.96, 13.96, 13.96, 13.96, 13.96, 13.96};
  // limits.amin = {-13.96, -13.96, -13.96, -13.96, -13.96, -13.96};
  limits.tau_max = {150.000000, 150.000000, 150.000000, 28.000000, 28.000000, 28.000000};
  // limits.tau_min = {-150.000000, -150.000000, -150.000000, -28.000000, -28.000000, -28.000000};
  limits.tcp_max = 2.0; // todo: verify this value

  blast::ManipulatorKinematics kinematics;
  kinematics.dv = {
          Vec3{0, 0, 0},
          {-0.425, 0, 0},
          {-0.3922, 0, 0.1333},
          {0, -0.0997, -0},
          {0, 0.0996, -0},
          {0, 0, 0} // to end effector
  };
  kinematics.ev = {
          Vec3{0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}};
  // kinematics.Q_static.resize(6);
  kinematics.Q_static[0] = {-1.000000, 0.000000, 0.000000, -0.000000, -1.000000, 0.000000, 0.000000, -0.000000, 1.000000}; // modified
  kinematics.Q_static[1] = {1.000000, 0.000000, 0.000000, -0.000000, -0.000000, 1.000000, 0.000000, -1.000000, -0.000000};
  kinematics.Q_static[2] = {1.000000, 0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000000, 1.000000};
  kinematics.Q_static[3] = {1.000000, 0.000000, 0.000000, -0.000000, 1.000000, 0.000000, 0.000000, -0.000000, 1.000000};
  kinematics.Q_static[4] = {1.000000, 0.000000, 0.000000, -0.000000, -0.000000, 1.000000, 0.000000, -1.000000, -0.000000};
  kinematics.Q_static[5] = {1.000000, -0.000000, 0.000000, 0.000000, -0.000000, -1.000000, 0.000000, 1.000000, -0.000000};
  // kinematics.Q_static[6] = {1.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 1.000000};
  kinematics.p_j0 = {0.000000, 0.000000, 0.162500};

  blast::ManipulatorDynamics dynamics;
  dynamics.m = {3.700000, 8.393000, 2.275000, 1.219000, 1.219000, 0.187900};
  dynamics.I = {
          Mat3{0.0103, 0, 0, 0, 0.0103, 0, 0, 0, 0.0067},
          {0.1339, 0, 0, 0, 0.1339, 0, 0, 0, 0.0151},
          {0.0312, 0, 0, 0, 0.0312, 0, 0, 0, 0.0041},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0026, 0, 0, 0, 0.0026, 0, 0, 0, 0.0022},
          {0.0001, 0, 0, 0, 0.0001, 0, 0, 0, 0.0001}};
  dynamics.av = {
          Vec3{0, 0, 0},
          {-0.2125, 0, 0.138},
          {-0.1961, 0, 0.007},
          {0, 0, 0},
          {0, 0, 0},
          {0, 0, -0.0229}};

  blast::ManipulatorCapsules collisions;

  blast::Sphere base;
  base.c                    = {0, 0, 0.0325};
  base.r                    = 0.09;
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
  capsule.r           = 0.09;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 1;
  capsule.p1          = {-0.42, 0, 0.1375};
  capsule.p2          = {0, 0, 0.1375};
  capsule.r           = 0.06;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {0, 0, 0.02};
  capsule.p2          = {0, 0, 0.18};
  capsule.r           = 0.065;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 2;
  capsule.p1          = {-0.373156, 0, 0.00850418};
  capsule.p2          = {0.000440611, 0.000121443, 0.00850418};
  capsule.r           = 0.05;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, 0, 0};
  capsule.p2          = {0, 0, -0.155};
  capsule.r           = 0.0425;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 3;
  capsule.p1          = {0, 0.03, 0};
  capsule.p2          = {0, -0.1, 0};
  capsule.r           = 0.0425;
  collisions.capsule_list.push_back(capsule);

  capsule.joint_frame = 5;
  capsule.p1          = {0, 0, -0.14};
  capsule.p2          = {0, 0, -0.01};
  capsule.r           = 0.038;
  collisions.capsule_list.push_back(capsule);

  blast::Manipulator generic_manip(joints, limits, kinematics, &dynamics, &collisions);
  return generic_manip;
}

int main() {
  auto manip   = get_generic_ur5e();
  manip.p_base = {-0.5, -0.3, 0.35};
  manip.Q_base = {-1, 0, 0, 0, -1, 0, 0, 0, 1};

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

  ManipulatorTempData manip_data;

  Matrix p_j_start(start_pos.size(), 3);
  Matrix p_j_end(start_pos.size(), 3);
  for (int i = 0; i < start_pos.size(); ++i) {
    auto joint_pos_start = start_pos[i];
    forward_kinematics(manip, manip_data, joint_pos_start);
    p_j_start(i, 0) = manip_data.p_j[6].x;
    p_j_start(i, 1) = manip_data.p_j[6].y;
    p_j_start(i, 2) = manip_data.p_j[6].z;

    auto joint_pos_end = end_pos[i];
    forward_kinematics(manip, manip_data, joint_pos_end);
    p_j_end(i, 0) = manip_data.p_j[6].x;
    p_j_end(i, 1) = manip_data.p_j[6].y;
    p_j_end(i, 2) = manip_data.p_j[6].z;
  }

  print_to_csv(p_j_start, "../../../examples/p_ee_start.csv");
  print_to_csv(p_j_end, "../../../examples/p_ee_end.csv");
}
