#pragma once

#include <blast>

namespace blast {

// user guide: https://www.kinovarobotics.com/uploads/User-Guide-Gen3-R07.pdf
// todo: Validate limits
struct Gen3 {
    // basic manipulator properties
    int     joints = 7;

    // All limits are from webapp
    Array   pmax = {INF_REAL, 2.25, INF_REAL, 2.58f, INF_REAL, 2.1f, INF_REAL}; // rad
    Array   pmin = - pmax; // rad
    Array   vmax = {1.745f, 1.745f, 1.745f, 1.745f, 2.443f, 2.443f, 2.443f}; // rad/s
    Array   vmin = - vmax; // rad/s
    Array   amax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL};
    Array   amin = - amax; // rad/s^2
    Array   tau_max = {52, 52, 52, 52, 17, 17, 17}; // Nm
    Array   tau_min = - tau_max; // Nm

    // tcp & elbow speed limits
    real tcp_max = 0.5; // from user guide

    // kinematic properties
    Vec3    p_base;
    Mat3    Q_base;

    Vec3    p_j0 = {0, 0, 0.1564f};

    Vec3    dv[7] = {
        {0.0, 0.0054, -0.1284},
        {0.0, -0.2104, -0.0064},
        {0.0, -0.0064, -0.2104},
        {0.0, -0.2084, -0.0064},
        {0.0, 0.0, -0.1059},
        {0.0, -0.1059, 0.0},
        {0.0, 0.0, -0.0615 /*- 0.164*/} // todo: add gripper capsule
    }; // vector to next joint

    Vec3    ev[7] = {
        {0, 0, 1},
        {0, 0, 1},
        {0, 0, 1},
        {0, 0, 1},
        {0, 0, 1},
        {0, 0, 1},
        {0, 0, 1}
    }; // direction vectors of joint

    // dynamic properties
    real    m[7] = {
        1.377f,
        1.1636f,
        1.1636f,
        0.93f,
        0.678f,
        0.678f,
        0.364f
    }; // link masses

    Mat3    I[7] = {
        {0.004570f, 0.000001f, 0.000002f, 0.000001f, 0.004831f, 0.000448f, 0.000002f, 0.000448f, 0.001409f},
        {0.011088f, 0.000005f, 0.000000f, 0.000005f, 0.001072f, -0.000691f, 0.000000f, - 0.000691f, 0.011255f},
        {0.010932f, 0.000000f, -0.000007f, 0.000000f, 0.011127f, 0.000606f, -0.000007f, 0.000606f, 0.001043f},
        {0.008147f, -0.000001f, 0.000000f, -0.000001f, 0.000631f, -0.000500f, 0.000000f, -0.000500f, 0.008316f},
        {0.001596f, 0.000000f, 0.000000f, 0.000000f, 0.001607f, 0.000256f, 0.000000f, 0.000256f, 0.000399f},
        {0.001641f, 0.000000f, 0.000000f, 0.000000f, 0.000410f, -0.000278f, 0.000000f, -0.000278f, 0.001641f},
        {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, -0.000002f, 0.000001f, -0.000002f, 0.000240f}
    };  // Inertial tensors

    Vec3    av[7] = {
        {-0.000023f, -0.010364f, -0.073360f},
        {-0.000044f, -0.099580f, -0.013278f},
        {-0.000044f, -0.006641f, -0.117892f},
        {-0.000018f, -0.075478f, -0.015006f},
        {0.000001f, -0.009432f, -0.063883f},
        {0.000001f, -0.045483f, -0.009650f},
        {-0.000093f, 0.000132f, -0.022905f}
    }; // centers of mass

    Vec3    sv[7] = {
        {-dv[0] + av[0]},
        {-dv[1] + av[1]},
        {-dv[2] + av[2]},
        {-dv[3] + av[3]},
        {-dv[4] + av[4]},
        {-dv[5] + av[5]},
        {-dv[6] + av[6]}
    }; // centers of mass from next joint

    Array  _efforts; // put the efforts temporarily when computing the constraints

    // Internal variables
    std::vector<Mat3> Q_static = {};  // static rotation to next joint todo: remove this
    std::vector<Mat3>    _Q;          // put the rotation matrices temporarily when computing the constraints
    std::vector<Mat3>    _Q_mult;     // put the rotation matrices multiplications temporarily when computing the constraints
    std::vector<Vec3>    _p_j;        // put the joint coordinates temporarily when computing the constraints

    std::vector<Capsule> _capsule_list; // put the capsules temporarily when computing the constraints
    u32                  _n_caps = 3;
    u32                  _n_internal_collisions = 2;

    Gen3() {
        _Q.resize(7);
        _Q_mult.resize(7);
        _p_j.resize(8); // with end effector
        _capsule_list.resize(3);
    }

    inline void compute_rotation_matrices(const Array &joint_position) {
        Array s(7);
        Array c(7);
        sincos(joint_position, s, c);

        // note: these are stored column-wise
        _Q[0] = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        _Q[1] = {c[1], 0,  s[1], -s[1], 0,  c[1], 0, -1, 0};
        _Q[2] = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
        _Q[3] = {c[3], 0,  s[3], -s[3], 0,  c[3], 0, -1, 0};
        _Q[4] = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        _Q[5] = {c[5], 0,  s[5], -s[5], 0,  c[5], 0, -1, 0};
        _Q[6] = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};
    }

    inline void compute_capsules() {
        Vec3 p1;
        Vec3 p2;
        real r;

        // Capsule 1
        p1 = _p_j[1];
        p2 = _p_j[3];
        r = 0.055;
        _capsule_list[0] = {p1, p2, r};

        // Capsule 2
        p1 = _p_j[3];
        p2 = _p_j[5];
        r = 0.055;
        _capsule_list[1] = {p1, p2, r};

        // Capsule 3
        p1 = _p_j[5];
        p2 = _p_j[7];
        r = 0.055;
        _capsule_list[2] = {p1, p2, r};
    }

    inline Array internal_collisions() {
        // capsule covering the base
        // Capsule caps_base;
        // caps_base.p1 = p_base + {0.0, 0.0, 0.0};
        // caps_base.p2 = _p_j[1];
        // caps_base.r = 0.0875;

        Sphere sphere;
        Vec3 sphere_center_relative = {0, 0, -0.075};
        sphere.c = p_base + sphere_center_relative; // todo: Verify this sphere
        // sphere.c = {0, 0, -0.1564f};
        sphere.r = 0.15;

        real dist1 = distance(_capsule_list[0], _capsule_list[2]);
        real dist2 = distance(_capsule_list[2], sphere);

        return {dist1, dist2};
    }

};

} // namespace blast