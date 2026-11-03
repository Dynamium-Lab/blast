#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <blast>

using namespace blast;

// user guide: https://www.kinovarobotics.com/uploads/User-Guide-Gen3-R07.pdf
struct Gen3 {
  // basic manipulator properties
  int joints = 7;

  // All limits are from webapp
  Array pmax    = {INF_REAL, 2.25, INF_REAL, (real) 2.58, INF_REAL, (real) 2.1, INF_REAL};                            // rad
  Array pmin    = -pmax;                                                                                              // rad
  Array vmax    = {(real) 1.745, (real) 1.745, (real) 1.745, (real) 1.745, (real) 2.443, (real) 2.443, (real) 2.443}; // rad/s
  Array vmin    = -vmax;                                                                                              // rad/s
  Array amax    = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL};
  Array amin    = -amax;                                                                                              // rad/s^2
  Array tau_max = {52, 52, 52, 52, 17, 17, 17};                                                                       // Nm
  Array tau_min = -tau_max;                                                                                           // Nm

  // tcp & elbow speed limits
  real tcp_max = 0.5; // from user guide

  // kinematic properties
  Vec3 p_base = {0, 0, 0};
  Mat3 Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  Vec3 p_j0 = {0, 0, (real) 0.1564};

  Vec3 dv[7] = {
          {0.0, 0.0054, -0.1284},
          {0.0, -0.2104, -0.0064},
          {0.0, -0.0064, -0.2104},
          {0.0, -0.2084, -0.0064},
          {0.0, 0.0, -0.1059},
          {0.0, -0.1059, 0.0},
          {0.0, 0.0, -0.0615 /*- 0.164*/} // todo: add gripper capsule
  }; // vector to next joint

  Vec3 ev[7] = {
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint

  // dynamic properties
  real m[7] = {
          (real) 1.377,
          (real) 1.1636,
          (real) 1.1636,
          (real) 0.93,
          (real) 0.678,
          (real) 0.678,
          (real) 0.364}; // link masses

  Mat3 I[7] = {
          {(real) 0.004570, (real) 0.000001, (real) 0.000002, (real) 0.000001, (real) 0.004831, (real) 0.000448, (real) 0.000002, (real) 0.000448, (real) 0.001409},
          {(real) 0.011088, (real) 0.000005, (real) 0, (real) 0.000005, (real) 0.001072, -(real) 0.000691, (real) 0, -(real) 0.000691, (real) 0.011255},
          {(real) 0.010932, (real) 0, -(real) 0.000007, (real) 0, (real) 0.011127, (real) 0.000606, -(real) 0.000007, (real) 0.000606, (real) 0.001043},
          {(real) 0.008147, -(real) 0.000001, (real) 0, -(real) 0.000001, (real) 0.000631, -(real) 0.000500, (real) 0, -(real) 0.000500, (real) 0.008316},
          {(real) 0.001596, (real) 0, (real) 0, (real) 0, (real) 0.001607, (real) 0.000256, (real) 0, (real) 0.000256, (real) 0.000399},
          {(real) 0.001641, (real) 0, (real) 0, (real) 0, (real) 0.000410, -(real) 0.000278, (real) 0, -(real) 0.000278, (real) 0.001641},
          {(real) 0.000214, (real) 0, (real) 0.000001, (real) 0, (real) 0.000223, -(real) 0.000002, (real) 0.000001, -(real) 0.000002, (real) 0.000240}}; // Inertial tensors

  Vec3 av[7] = {
          {-(real) 0.000023, -(real) 0.010364, -(real) 0.073360},
          {-(real) 0.000044, -(real) 0.099580, -(real) 0.013278},
          {-(real) 0.000044, -(real) 0.006641, -(real) 0.117892},
          {-(real) 0.000018, -(real) 0.075478, -(real) 0.015006},
          {(real) 0.000001, -(real) 0.009432, -(real) 0.063883},
          {(real) 0.000001, -(real) 0.045483, -(real) 0.009650},
          {-(real) 0.000093, (real) 0.000132, -(real) 0.022905}}; // centers of mass

  Vec3 sv[7] = {
          {-dv[0] + av[0]},
          {-dv[1] + av[1]},
          {-dv[2] + av[2]},
          {-dv[3] + av[3]},
          {-dv[4] + av[4]},
          {-dv[5] + av[5]},
          {-dv[6] + av[6]}}; // centers of mass from next joint

  Array _efforts;            // put the efforts temporarily when computing the constraints

  // Internal variables
  std::vector<Mat3>    Q_static               = {};                      // static rotation to next joint todo: remove this
  std::vector<Mat3>    _rotations             = std::vector<Mat3>(7);    // put the rotation matrices temporarily when computing the constraints
  std::vector<Mat3>    _rotations_mult        = std::vector<Mat3>(7);    // put the rotation matrices multiplications temporarily when computing the constraints
  std::vector<Vec3>    _p_j                   = std::vector<Vec3>(8);    // put the joint coordinates temporarily when computing the constraints
  std::vector<Capsule> _capsule_list          = std::vector<Capsule>(3); // put the capsules temporarily when computing the constraints
  u32                  _n_caps                = 3;
  u32                  _n_internal_collisions = 2;

  inline void compute_rotation_matrices(const Array& joint_position) {
    Array s(7);
    Array c(7);
    sincos(joint_position, s, c);

    // note: these are stored column-wise
    _rotations[0] = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    _rotations[1] = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    _rotations[2] = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
    _rotations[3] = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    _rotations[4] = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    _rotations[5] = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
    _rotations[6] = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};
  }
};

TEST_CASE("Gen3 compute_rotation_matrices() test", "[Generic]") {
  int  n_tests = 100;
  real epsilon = 1e-6;
  Gen3 expected_manip;
  auto generic_manip          = make_Kinova_Gen3();
  generic_manip.base_position = {0.0, 0.0, 0.0};
  generic_manip.base_rotation = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  ManipulatorTempData data_expected;
  ManipulatorTempData data_generic;

  for (int i = 0; i < n_tests; i++) {
    Array test_position(7);
    fill_random(test_position, PI);

    forward_kinematics(generic_manip, data_generic, test_position);
    expected_manip.compute_rotation_matrices(test_position);

    for (int j = 0; j < expected_manip._rotations.size(); j++) {
      CHECK(is_close(data_generic.rotations[j], expected_manip._rotations[j], epsilon));
    }
  }
}
