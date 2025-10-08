#pragma once

#include <blast>

namespace blast {

// user guide: https://artifactory.kinovaapps.com/ui/native/generic-documentation-public/Documentation/Link%206/Technical%20documentation/User%20Guide/
struct Link6 {
  // basic manipulator properties
  u32 joints = 6;

  // All limits are from webapp
  Array pmax    = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL};                         // rad
  Array pmin    = -pmax;                                                                                // rad
  Array vmax    = {3.4907f, 3.4907f, 3.4907f, 5.5851f, 5.5851f, 5.5851f};                               // rad/s
  Array vmin    = -vmax;                                                                                // rad/s
  Array amax    = {deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600)}; // rad/s^2
  Array amin    = -amax;                                                                                // rad/s^2
  Array tau_max = {210, 210, 210, 100, 100, 100};                                                       // Nm
  Array tau_min = -tau_max;                                                                             // Nm

  // tcp & elbow speed limits
  real tcp_max = 2.0;
  // real elbow_max = 2.0;

  // kinematic properties
  Vec3 p_base = {0, 0, 0};
  Mat3 Q_base = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  Vec3 p_j0 = {0, 0, 0.0530f};

  Vec3 dv[6] = {
          {0.11024, -0.06926, -0.1375},
          {0.0, 0.4850, 0.0},
          {0.0, -0.15216, -0.0917},
          {0.0, -0.06296, -0.22275},
          {0.08703, 0.0860, -0.07692},
          {0.0, 0.0, -0.0920}}; // vector to next joint

  Vec3 ev[6] = {
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1},
          {0, 0, 1}}; // direction vectors of joint

  // dynamic properties
  real m[6] = {
          4.8257f,
          5.9860f,
          3.4159f,
          2.0849f,
          2.0076f,
          1.5193f}; // link masses
  Mat3 I[6] = {
          {0.0192746f, -0.00239802f, -0.00896331f, -0.00239802f, 0.03087806f, 0.0016298f, -0.00896331f, 0.0016298f, 0.02134949f},
          {0.25899206f, -2.89E-05f, -1.23E-06f, -2.89E-05f, 0.01755445f, -0.02128064f, -1.23E-06f, -0.02128064f, 0.25291674f},
          {0.01742043f, -3.55E-06f, 8.4E-07f, -3.55E-06f, 0.01119175f, 0.00518163f, 8.4E-07f, 0.00518163f, 0.01212876f},
          {0.02454276f, 2.61E-06f, 1.799E-05f, 2.61E-06f, 0.02385702f, 0.00315758f, 1.799E-05f, 0.00315758f, 0.00294903f},
          {0.00734684f, 0.00124927f, -0.00090156f, 0.00124927f, 0.00464684f, -0.00236128f, -0.00090156f, -0.00236128f, 0.00589508f},
          {0.00390762f, -1.13E-06f, 1.16E-06f, -1.13E-06f, 0.00390722f, -2.21E-05f, 1.16E-06f, -2.21E-05f, 0.0013928f}}; // Inertial tensors
  Vec3 av[6] = {
          {0.03930119f, -0.00705889f, -0.08462154f},
          {2.53E-06f, 0.18829586f, -0.03988382f},
          {4.64E-06f, -0.02451414f, -0.02997969f},
          {-0.00010793f, -0.01056422f, -0.08091102f},
          {0.01243595f, 0.03284165f, -0.04091434f},
          {0.0f, 0.00050624f, -0.00388589f}}; // centers of mass

  Vec3 sv[6] = {
          {-dv[0] + av[0]},
          {-dv[1] + av[1]},
          {-dv[2] + av[2]},
          {-dv[3] + av[3]},
          {-dv[4] + av[4]},
          {-dv[5] + av[5]}}; // centers of mass from next joint

  // Internal variables
  Array _efforts;                  // put the efforts temporarily when computing the constraints

  std::vector<Mat3> Q_static = {}; // static rotation to next joint todo: remove this
  std::vector<Mat3> _Q;            // put the rotation matrices temporarily when computing the constraints
  std::vector<Mat3> _Q_mult;       // put the rotation matrices multiplications temporarily when computing the constraints
  std::vector<Vec3> _p_j;          // put the joint coordinates temporarily when computing the constraints

  Sphere sph_base = {p_base, 0.2375};

  std::vector<Capsule> _capsule_list; // put the capsules temporarily when computing the constraints
  u32                  _n_caps                = 7;
  u32                  _n_internal_collisions = 9;

  Link6() {
    _Q.resize(6);
    _Q_mult.resize(6);
    _p_j.resize(7); // with end effector
    _capsule_list.resize(7);
  }

  void compute_rotation_matrices(const Array& joint_position) {
    Array s(6);
    Array c(6);
    blast::sincos(joint_position, s, c);

    // note: these are stored column-wise
    _Q[0] = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1}; // base -> 0
    _Q[1] = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};  // 0    -> 1
    _Q[2] = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1}; // 1    -> 2
    _Q[3] = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};  // 2    -> 3
    _Q[4] = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};  // 3    -> 4
    _Q[5] = {0, s[5], c[5], 0, c[5], -s[5], -1, 0, 0};   // 4    -> 5
  }

  void compute_capsules() {
    Vec3 z2  = {_Q_mult[1].col_copy(2)};
    Vec3 y2  = {_Q_mult[1].col_copy(1)};
    Vec3 z3  = {_Q_mult[2].col_copy(2)};
    Vec3 y3  = {_Q_mult[2].col_copy(1)};
    Vec3 z4  = {_Q_mult[3].col_copy(2)};
    Vec3 z5  = {_Q_mult[4].col_copy(2)};
    Vec3 z6  = {_Q_mult[5].col_copy(2)};
    Vec3 zee = -z6;
    Vec3 yee = {-_Q_mult[5].col_copy(1)};

    Vec3 p1;
    Vec3 p2;
    real r;

    // Capsule 1
    p1               = _p_j[1] - 0.065 * z2;
    p2               = p1 + 0.11 * z2;
    r                = 0.065;
    _capsule_list[0] = {p1, p2, r};

    // Capsule 2
    p1               = _p_j[1] - 0.08 * z2;
    p2               = p1 + 0.485 * y2;
    r                = 0.065;
    _capsule_list[1] = {p1, p2, r};

    // Capsule 3
    p1               = _p_j[2] - 0.065 * z3;
    p2               = p1 + 0.150 * z3;
    r                = 0.065;
    _capsule_list[2] = {p1, p2, r};

    // Capsule 4
    p1               = _p_j[2] - 0.0917 * z3 + 0.00695 * y3;
    p2               = p1 - 0.375 * z4;
    r                = 0.061;
    _capsule_list[3] = {p1, p2, r};

    // Capsule 5
    p1               = _p_j[4];
    p2               = p1 - 0.08 * z5;
    r                = 0.060;
    _capsule_list[4] = {p1, p2, r};

    // Capsule 6
    p1               = _p_j[5] + 0.08583 * z6;
    p2               = p1 - 0.15 * z6;
    r                = 0.060;
    _capsule_list[5] = {p1, p2, r};

    // Capsule 7 todo: Refine capsule
    p1               = _p_j[6] - 0.085 * zee - 0.02125 * yee;
    p2               = p1 - 0.150 * zee;
    r                = 0.085;
    _capsule_list[6] = {p1, p2, r};
  }

  Array internal_collisions() const {
    return {distance(_capsule_list[0], _capsule_list[5]),
            distance(_capsule_list[0], _capsule_list[6]),
            distance(_capsule_list[1], _capsule_list[4]),
            distance(_capsule_list[1], _capsule_list[5]),
            distance(_capsule_list[1], _capsule_list[6]),
            distance(_capsule_list[3], sph_base),
            distance(_capsule_list[4], sph_base),
            distance(_capsule_list[5], sph_base),
            distance(_capsule_list[6], sph_base)};
  }
};

} // namespace blast
