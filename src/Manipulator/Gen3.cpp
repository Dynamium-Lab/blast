#include "blast.h"
#include "blast_error.h"

namespace blast {

real clamped_root(real slope, real h0, real h1);
void compute_intersection(real* sValue, i32* classify, real b, real f00, real f10, i32* edge, real end[][2]);
void compute_minimum_parameters(i32* edge, real end[][2], real b, real c, real e, real g00, real g10, real g01, real g11, real* parameter);
real two_segment_distance_sqr(Vec3 P0, Vec3 P1, Vec3 Q0, Vec3 Q1);


Gen3::Gen3() {
    p_base = {0, 0, 0.1564f};
    pmax.resize(7);
    vmax.resize(7);
    tau_max.resize(7);
    pmax = {INF_REAL, INF_REAL, INF_REAL, 2.58f, INF_REAL, 2.1f, INF_REAL};
    vmax = {1.745f, 1.745f, 1.745f, 1.745f, 2.443f, 2.443f, 2.443f};
    tau_max = {52, 52, 52, 52, 17, 17, 17};
    set_payload(0);
}

void Gen3::set_payload(const real m_payload, const Vec3 cg_payload, const Mat3 I_payload) {
    m[0] = 1.377f;
    m[1] = 1.1636f;
    m[2] = 1.1636f;
    m[3] = 0.93f;
    m[4] = 0.678f;
    m[5] = 0.678f;
    m[6] = 0.364f;

    I[0] = {0.004570f, 0.000001f, 0.000002f, 0.000001f, 0.004831f, 0.000448f, 0.000002f, 0.000448f, 0.001409f};
    I[1] = {0.011088f, 0.000005f, 0.000000f, 0.000005f, 0.001072f, 0.000691f, 0.000000f, 0.000691f, 0.011255f};
    I[2] = {0.010932f, 0.000000f, 0.000007f, 0.000000f, 0.011127f, 0.000606f, 0.000007f, 0.000606f, 0.001043f};
    I[3] = {0.008147f, 0.000001f, 0.000000f, 0.000001f, 0.000631f, 0.000500f, 0.000000f, 0.000500f, 0.008316f};
    I[4] = {0.001596f, 0.000000f, 0.000000f, 0.000000f, 0.001607f, 0.000256f, 0.000000f, 0.000256f, 0.000399f};
    I[5] = {0.001641f, 0.000000f, 0.000000f, 0.000000f, 0.000410f, 0.000278f, 0.000000f, 0.000278f, 0.001641f};
    I[6] = {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, 0.000002f, 0.000001f, 0.000002f, 0.000240f};

    // center of mass
    av[0] = {-0.000023f, -0.010364f, -0.073360f};
    av[1] = {-0.000044f, -0.099580f, -0.013278f};
    av[2] = {-0.000044f, -0.006641f, -0.117892f};
    av[3] = {-0.000018f, -0.075478f, -0.015006f};
    av[4] = {0.000001f, -0.009432f, -0.063883f};
    av[5] = {0.000001f, -0.045483f, -0.009650f};
    av[6] = {-0.000093f, 0.000132f, -0.022905f};

    // vector to next joint
    dv[0] = {0.0, 0.0054, -0.1284};
    dv[1] = {0.0, -0.2104, -0.0064};
    dv[2] = {0.0, -0.0064, -0.2104};
    dv[3] = {0.0, -0.2084, -0.0064};
    dv[4] = {0.0, 0.0, -0.1059};
    dv[5] = {0.0, -0.1059, 0.0};
    dv[6] = {0.0, 0.0, -0.0615 - 0.164};
    // todo: add option to know if tool is closed or opened (difference of 0.0135 in z)

    // unit joint direction
    ev[0] = {0, 0, 1};
    ev[1] = {0, 0, 1};
    ev[2] = {0, 0, 1};
    ev[3] = {0, 0, 1};
    ev[4] = {0, 0, 1};
    ev[5] = {0, 0, 1};
    ev[6] = {0, 0, 1};

    // Adjust end effector link with payload
    m[6] += 0.921f;
    av[6] = {-0.000093f, 0.000132f, -0.022905f};
    Vec3 av_tool(0, 0, -0.06f - 0.0615f);
    av[6] = (0.364f * av[6] + 0.921f * av_tool) * (1 / (0.364f + 0.921f));
    sv[6] = dv[6] - av[6];
    I[6] = {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, 0.000002f, 0.000001f, 0.000002f, 0.000240f};

    // Modify payload
    Vec3 av_payload = {0.0f, 0.0f, -0.115f};
    av_payload += cg_payload;

    auto av_old = av[6];
    auto m_old = m[6];

    auto m_new = m_old + m_payload;
    auto av_new = (m_old*av_old + m_payload*av_payload) / m_new;
    auto delta_av = av_new - av_old; // shift in center of mass
    auto av_to_mass = av_payload - av_new; // vector from payload to new center of mass

    av[6] = av_new;
    m[6] = m_new;
    I[6](0, 0) += m_old*delta_av.x*delta_av.x + m_payload*av_to_mass.x*av_to_mass.x; // todo: Validate
    I[6](1, 1) += m_old*delta_av.y*delta_av.y + m_payload*av_to_mass.y*av_to_mass.y;
    I[6](2, 2) += m_old*delta_av.z*delta_av.z + m_payload*av_to_mass.z*av_to_mass.z;
    I[6] += I_payload;

    // center of mass from next joint
    sv[0] = -dv[0] + av[0];
    sv[1] = -dv[1] + av[1];
    sv[2] = -dv[2] + av[2];
    sv[3] = -dv[3] + av[3];
    sv[4] = -dv[4] + av[4];
    sv[5] = -dv[5] + av[5];
    sv[6] = -dv[6] + av[6];

}

void Gen3::set_payload_without_gripper(const real m_payload, const Vec3 cg_payload, const Mat3 I_payload) {
    m[0] = 1.377f;
    m[1] = 1.1636f;
    m[2] = 1.1636f;
    m[3] = 0.93f;
    m[4] = 0.678f;
    m[5] = 0.678f;
    m[6] = 0.364f;

    I[0] = {0.004570f, 0.000001f, 0.000002f, 0.000001f, 0.004831f, 0.000448f, 0.000002f, 0.000448f, 0.001409f};
    I[1] = {0.011088f, 0.000005f, 0.000000f, 0.000005f, 0.001072f, 0.000691f, 0.000000f, 0.000691f, 0.011255f};
    I[2] = {0.010932f, 0.000000f, 0.000007f, 0.000000f, 0.011127f, 0.000606f, 0.000007f, 0.000606f, 0.001043f};
    I[3] = {0.008147f, 0.000001f, 0.000000f, 0.000001f, 0.000631f, 0.000500f, 0.000000f, 0.000500f, 0.008316f};
    I[4] = {0.001596f, 0.000000f, 0.000000f, 0.000000f, 0.001607f, 0.000256f, 0.000000f, 0.000256f, 0.000399f};
    I[5] = {0.001641f, 0.000000f, 0.000000f, 0.000000f, 0.000410f, 0.000278f, 0.000000f, 0.000278f, 0.001641f};
    I[6] = {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, 0.000002f, 0.000001f, 0.000002f, 0.000240f};

    // center of mass
    av[0] = {-0.000023f, -0.010364f, -0.073360f};
    av[1] = {-0.000044f, -0.099580f, -0.013278f};
    av[2] = {-0.000044f, -0.006641f, -0.117892f};
    av[3] = {-0.000018f, -0.075478f, -0.015006f};
    av[4] = {0.000001f, -0.009432f, -0.063883f};
    av[5] = {0.000001f, -0.045483f, -0.009650f};
    av[6] = {-0.000093f, 0.000132f, -0.022905f};

    // vector to next joint
    dv[0] = {0.0, 0.0054, -0.1284};
    dv[1] = {0.0, -0.2104, -0.0064};
    dv[2] = {0.0, -0.0064, -0.2104};
    dv[3] = {0.0, -0.2084, -0.0064};
    dv[4] = {0.0, 0.0, -0.1059};
    dv[5] = {0.0, -0.1059, 0.0};
    dv[6] = {0.0, 0.0, -0.0615 - 0.164};
    // todo: add option to know if tool is closed or opened (difference of 0.0135 in z)

    // unit joint direction
    ev[0] = {0, 0, 1};
    ev[1] = {0, 0, 1};
    ev[2] = {0, 0, 1};
    ev[3] = {0, 0, 1};
    ev[4] = {0, 0, 1};
    ev[5] = {0, 0, 1};
    ev[6] = {0, 0, 1};


    // Modify payload
    Vec3 av_payload = {0.0f, 0.0f, -0.115f};
    av_payload += cg_payload;

    auto av_old = av[6];
    auto m_old = m[6];

    auto m_new = m_old + m_payload;
    auto av_new = (m_old*av_old + m_payload*av_payload) /m_new;
    auto sv_new = dv[6] - av_new;
    auto delta_av = av_new - av_old; // shift in center of mass
    auto av_to_mass = av_payload - av_new; // vector from payload to new center of mass

    av[6] = av_new;
    sv[6] = sv_new;
    m[6] = m_new;
    I[6](0, 0) += m_old*delta_av.x*delta_av.x + m_payload*av_to_mass.x*av_to_mass.x; // todo: Validate
    I[6](1, 1) += m_old*delta_av.y*delta_av.y + m_payload*av_to_mass.y*av_to_mass.y;
    I[6](2, 2) += m_old*delta_av.z*delta_av.z + m_payload*av_to_mass.z*av_to_mass.z;
    I[6] += I_payload;

    // center of mass from next joint
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];
    sv[6] = dv[6] - av[6];
}

void Gen3::internal_constraints(const Trajectory& traj, real* dst) {
    const auto points = traj.pos.cols;
    dynamics(traj);

    for (u32 i = 0; i < points; i++) {
        // 5 self collisions
        auto p = traj.pos.col(i); Assert(p.is_alias);
        auto tmp_coll = internal_collisions(p);
        dst[0] = -tmp_coll[0]; // dist1sqr
        dst[1] = -tmp_coll[1]; // dist2sqr
        dst[2] = -tmp_coll[2]; // distTJ4
        dst[3] = -tmp_coll[3]; // distTJ6
        dst[4] = -tmp_coll[4]; // distTEE
        dst += 5;

        // 2 position limits
        dst[0] = (std::abs(traj.pos(3, i)) - pmax[3]) / pmax[3];
        dst[1] = (std::abs(traj.pos(5, i)) - pmax[5]) / pmax[5];
        dst += 2;

        // 7 velocity limits
        for (int j = 0; j < 7; j++)
            dst[j] = (std::abs(traj.vel(j, i)) - vmax[j]) / vmax[j];
        dst += 7;

        // 7 torque limits
        auto f = _efforts.col(i); Assert(f.is_alias);
        for (int j = 0; j < 7; j++)
            dst[j] = (std::abs(f[j]) - tau_max[j]) / tau_max[j];
        dst += 7;
    }
}

bool Gen3::validate_task(const Matrix &task, World *world) {
    Trajectory traj(2, 7);
    traj.pos.col(0) = task.col(0);
    traj.pos.col(1) = task.col(3);

    traj.vel.col(0) = task.col(1);
    traj.vel.col(1) = task.col(4);

    traj.acc.col(0) = task.col(2);
    traj.acc.col(1) = task.col(5);

    Array con(ncon(2));
    internal_constraints(traj, con.data);

    auto max_con = max(con);

    if (world) {
        Gen3 manip;
        std::vector<Capsule> capsules = manip.robot_capsules(traj.pos, 1);
        auto worst_collision = - test_collision(&capsules, world, 1);
        max_con = (worst_collision[0] > max_con) ? worst_collision[0] : max_con;
    }

    return max_con <= 0;
}

int Gen3::ncon(int points) {
    return (5 + 2 + 7 * 2) * points;
}

void Gen3::dynamics(const Trajectory& traj) {
    const auto points = traj.pos.cols;
    const auto joints = traj.pos.rows;
    if (_efforts.cols != points || _efforts.rows != joints)
        _efforts.resize(joints, points);


    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t, Q7t;
    Vec3 w1, w2, w3, w4, w5, w6, w7;
    Vec3 wd1, wd2, wd3, wd4, wd5, wd6, wd7;
    Vec3 cdd0 = {0, 0, 9.81f};
    Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6, cdd7;
    Vec3 f1, f2, f3, f4, f5, f6, f7;
    Vec3 n1, n2, n3, n4, n5, n6, n7;

    // loop all points
    for (u32 i = 0; i < points; i++) {
        auto p = traj.pos.col(i);
        auto v = traj.vel.col(i);
        auto a = traj.acc.col(i);

        Array s(7);
        Array c(7);
        blast::sincos(p, s, c);

        // note: these are stored column-wise
        Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
        Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
        Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
        Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
        Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};
        Q1t = transpose(Q1);
        Q2t = transpose(Q2);
        Q3t = transpose(Q3);
        Q4t = transpose(Q4);
        Q5t = transpose(Q5);
        Q6t = transpose(Q6);
        Q7t = transpose(Q7);

        // note: This is the Newton algorithm in 'Element de robotique' course notes.
        //       Careful because some variables are named differently and uses a slightly different conventions.
        //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
        //       joint turns with respect to the coordinate frame.
        //-- kinematics
        w1 = v[0] * ev[0];
        w2 = Q2t * w1 + v[1] * ev[1];
        w3 = Q3t * w2 + v[2] * ev[2];
        w4 = Q4t * w3 + v[3] * ev[3];
        w5 = Q5t * w4 + v[4] * ev[4];
        w6 = Q6t * w5 + v[5] * ev[5];
        w7 = Q7t * w6 + v[6] * ev[6];

        wd1 = a[0] * ev[0];
        cdd1 = Q1t * cdd0 + cross(wd1, av[0]) + cross(w1, cross(w1, av[0]));

        wd2 = Q2t * wd1 + a[1] * ev[1] + v[1] * cross(Q2t * w1, ev[1]);
        cdd2 = Q2t * cdd1 + cross(wd2, av[1]) + cross(w2, cross(w2, av[1])) - Q2t * cross(wd1, sv[0]) - Q2t * cross(w1, cross(w1, sv[0]));

        wd3 = Q3t * wd2 + a[2] * ev[2] + v[2] * cross(Q3t * w2, ev[2]);
        cdd3 = Q3t * cdd2 + cross(wd3, av[2]) + cross(w3, cross(w3, av[2])) - Q3t * cross(wd2, sv[1]) - Q3t * cross(w2, cross(w2, sv[1]));

        wd4 = Q4t * wd3 + a[3] * ev[3] + v[3] * cross(Q4t * w3, ev[3]);
        cdd4 = Q4t * cdd3 + cross(wd4, av[3]) + cross(w4, cross(w4, av[3])) - Q4t * cross(wd3, sv[2]) - Q4t * cross(w3, cross(w3, sv[2]));

        wd5 = Q5t * wd4 + a[4] * ev[4] + v[4] * cross(Q5t * w4, ev[4]);
        cdd5 = Q5t * cdd4 + cross(wd5, av[4]) + cross(w5, cross(w5, av[4])) - Q5t * cross(wd4, sv[3]) - Q5t * cross(w4, cross(w4, sv[3]));

        wd6 = Q6t * wd5 + a[5] * ev[5] + v[5] * cross(Q6t * w5, ev[5]);
        cdd6 = Q6t * cdd5 + cross(wd6, av[5]) + cross(w6, cross(w6, av[5])) - Q6t * cross(wd5, sv[4]) - Q6t * cross(w5, cross(w5, sv[4]));

        wd7 = Q7t * wd6 + a[6] * ev[6] + v[6] * cross(Q7t * w6, ev[6]);
        cdd7 = Q7t * cdd6 + cross(wd7, av[6]) + cross(w7, cross(w7, av[6])) - Q7t * cross(wd6, sv[5]) - Q7t * cross(w6, cross(w6, sv[5]));

        //-- dynamics
        f7 = m[6] * cdd7;
        n7 = I[6] * wd7 + cross(w7, I[6] * w7) + cross(av[6], f7);

        f6 = m[5] * cdd6 + Q7 * f7;
        n6 = I[5] * wd6 + cross(w6, I[5] * w6) + Q7 * n7 + cross(av[5], f6) - cross(sv[5], (Q7 * f7));

        f5 = m[4] * cdd5 + Q6 * f6;
        n5 = I[4] * wd5 + cross(w5, I[4] * w5) + Q6 * n6 + cross(av[4], f5) - cross(sv[4], (Q6 * f6));

        f4 = m[3] * cdd4 + Q5 * f5;
        n4 = I[3] * wd4 + cross(w4, I[3] * w4) + Q5 * n5 + cross(av[3], f4) - cross(sv[3], (Q5 * f5));

        f3 = m[2] * cdd3 + Q4 * f4;
        n3 = I[2] * wd3 + cross(w3, I[2] * w3) + Q4 * n4 + cross(av[2], f3) - cross(sv[2], (Q4 * f4));

        f2 = m[1] * cdd2 + Q3 * f3;
        n2 = I[1] * wd2 + cross(w2, I[1] * w2) + Q3 * n3 + cross(av[1], f2) - cross(sv[1], (Q3 * f3));

        f1 = m[0] * cdd1 + Q2 * f2;
        n1 = I[0] * wd1 + cross(w1, I[0] * w1) + Q2 * n2 + cross(av[0], f1) - cross(sv[0], (Q2 * f2));

        //-- extract torques (last element of each moment vector)
        _efforts(0, i) = n1.z;
        _efforts(1, i) = n2.z;
        _efforts(2, i) = n3.z;
        _efforts(3, i) = n4.z;
        _efforts(4, i) = n5.z;
        _efforts(5, i) = n6.z;
        _efforts(6, i) = n7.z;
    }
}

Array Gen3::dynamics(const Array& pos, const Array& vel, const Array& acc) {

    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t, Q7t;
    Vec3 w1, w2, w3, w4, w5, w6, w7;
    Vec3 wd1, wd2, wd3, wd4, wd5, wd6, wd7;
    Vec3 cdd0 = {0, 0, 9.81f};
    Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6, cdd7;
    Vec3 f1, f2, f3, f4, f5, f6, f7;
    Vec3 n1, n2, n3, n4, n5, n6, n7;

    // loop all points
    auto p = pos;
    auto v = vel;
    auto a = acc;

    Array s(7);
    Array c(7);
    blast::sincos(p, s, c);

    // note: these are stored column-wise
    Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
    Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
    Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};
    Q1t = transpose(Q1);
    Q2t = transpose(Q2);
    Q3t = transpose(Q3);
    Q4t = transpose(Q4);
    Q5t = transpose(Q5);
    Q6t = transpose(Q6);
    Q7t = transpose(Q7);

    // note: This is the Newton algorithm in 'Element de robotique' course notes.
    //       Careful because some variables are named differently and uses a slightly different conventions.
    //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
    //       joint turns with respect to the coordinate frame.
    //-- kinematics
    w1 = v[0] * ev[0];
    w2 = Q2t * w1 + v[1] * ev[1];
    w3 = Q3t * w2 + v[2] * ev[2];
    w4 = Q4t * w3 + v[3] * ev[3];
    w5 = Q5t * w4 + v[4] * ev[4];
    w6 = Q6t * w5 + v[5] * ev[5];
    w7 = Q7t * w6 + v[6] * ev[6];

    wd1 = a[0] * ev[0];
    cdd1 = Q1t * cdd0 + cross(wd1, av[0]) + cross(w1, cross(w1, av[0]));

    wd2 = Q2t * wd1 + a[1] * ev[1] + v[1] * cross(Q2t * w1, ev[1]);
    cdd2 = Q2t * cdd1 + cross(wd2, av[1]) + cross(w2, cross(w2, av[1])) - Q2t * cross(wd1, sv[0]) - Q2t * cross(w1, cross(w1, sv[0]));

    wd3 = Q3t * wd2 + a[2] * ev[2] + v[2] * cross(Q3t * w2, ev[2]);
    cdd3 = Q3t * cdd2 + cross(wd3, av[2]) + cross(w3, cross(w3, av[2])) - Q3t * cross(wd2, sv[1]) - Q3t * cross(w2, cross(w2, sv[1]));

    wd4 = Q4t * wd3 + a[3] * ev[3] + v[3] * cross(Q4t * w3, ev[3]);
    cdd4 = Q4t * cdd3 + cross(wd4, av[3]) + cross(w4, cross(w4, av[3])) - Q4t * cross(wd3, sv[2]) - Q4t * cross(w3, cross(w3, sv[2]));

    wd5 = Q5t * wd4 + a[4] * ev[4] + v[4] * cross(Q5t * w4, ev[4]);
    cdd5 = Q5t * cdd4 + cross(wd5, av[4]) + cross(w5, cross(w5, av[4])) - Q5t * cross(wd4, sv[3]) - Q5t * cross(w4, cross(w4, sv[3]));

    wd6 = Q6t * wd5 + a[5] * ev[5] + v[5] * cross(Q6t * w5, ev[5]);
    cdd6 = Q6t * cdd5 + cross(wd6, av[5]) + cross(w6, cross(w6, av[5])) - Q6t * cross(wd5, sv[4]) - Q6t * cross(w5, cross(w5, sv[4]));

    wd7 = Q7t * wd6 + a[6] * ev[6] + v[6] * cross(Q7t * w6, ev[6]);
    cdd7 = Q7t * cdd6 + cross(wd7, av[6]) + cross(w7, cross(w7, av[6])) - Q7t * cross(wd6, sv[5]) - Q7t * cross(w6, cross(w6, sv[5]));

    //-- dynamics
    f7 = m[6] * cdd7;
    n7 = I[6] * wd7 + cross(w7, I[6] * w7) + cross(av[6], f7);

    f6 = m[5] * cdd6 + Q7 * f7;
    n6 = I[5] * wd6 + cross(w6, I[5] * w6) + Q7 * n7 + cross(av[5], f6) - cross(sv[5], (Q7 * f7));

    f5 = m[4] * cdd5 + Q6 * f6;
    n5 = I[4] * wd5 + cross(w5, I[4] * w5) + Q6 * n6 + cross(av[4], f5) - cross(sv[4], (Q6 * f6));

    f4 = m[3] * cdd4 + Q5 * f5;
    n4 = I[3] * wd4 + cross(w4, I[3] * w4) + Q5 * n5 + cross(av[3], f4) - cross(sv[3], (Q5 * f5));

    f3 = m[2] * cdd3 + Q4 * f4;
    n3 = I[2] * wd3 + cross(w3, I[2] * w3) + Q4 * n4 + cross(av[2], f3) - cross(sv[2], (Q4 * f4));

    f2 = m[1] * cdd2 + Q3 * f3;
    n2 = I[1] * wd2 + cross(w2, I[1] * w2) + Q3 * n3 + cross(av[1], f2) - cross(sv[1], (Q3 * f3));

    f1 = m[0] * cdd1 + Q2 * f2;
    n1 = I[0] * wd1 + cross(w1, I[0] * w1) + Q2 * n2 + cross(av[0], f1) - cross(sv[0], (Q2 * f2));

    //-- extract torques (last element of each moment vector)
    Array result(7);

    result[0] = n1.z;
    result[1] = n2.z;
    result[2] = n3.z;
    result[3] = n4.z;
    result[4] = n5.z;
    result[5] = n6.z;
    result[6] = n7.z;

    return result;
}

Array Gen3::forward_kinematics(const Array &pos) {
    Array s(7);
    Array c(7);
    blast::sincos(pos, s, c);
    Mat3 Q1(c[0], -s[0],   0,    -s[0],  -c[0],   0,    0,   0, -1);
    Mat3 Q2(c[1],   0,    s[1],  -s[1],    0,    c[1],  0,  -1,  0);
    Mat3 Q3(c[2],   0,   -s[2],  -s[2],    0,   -c[2],  0,   1,  0);
    Mat3 Q4(c[3],   0,    s[3],  -s[3],    0,    c[3],  0,  -1,  0);
    Mat3 Q5(c[4],   0,   -s[4],  -s[4],    0,   -c[4],  0,   1,  0);
    Mat3 Q6(c[5],   0,    s[5],  -s[5],    0,    c[5],  0,  -1,  0);
    Mat3 Q7(c[6],   0,   -s[6],  -s[6],    0,   -c[6],  0,   1,  0);

    auto Q = Q1;
    auto p_ee = p_base
                + Q * dv[0]
                + (Q *= Q2) * dv[1]
                + (Q *= Q3) * dv[2]
                + (Q *= Q4) * dv[3]
                + (Q *= Q5) * dv[4]
                + (Q *= Q6) * dv[5]
                + (Q *= Q7) * dv[6];

    Array pose(6);
    pose[0] = p_ee.x;
    pose[1] = p_ee.y;
    pose[2] = p_ee.z;
    // todo: Q(0, 0) and Q(2, 2) must not be zero!!
    // pose[3] = atan2(Q(1, 0), Q(0, 0));
    // pose[4] = atan2(-Q(2, 0), sqrt(Q(2, 1)*Q(2, 1) + Q(2, 2)*Q(2, 2)));
    // pose[5] = atan2(Q(2, 1), Q(2, 2));
    return pose;
}

Matrix Gen3::forward_kinematics(const Matrix &pos) {
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Matrix poses(12, pos.cols);
    Vec3 p_ee;
    for (u32 point = 0; point < pos.cols; point++) {
        Array s(8);
        Array c(8);
        blast::sincos(pos.col(point), s, c);
        Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
        Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
        Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
        Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
        Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};
        auto Q = Q1;
        p_ee = p_base
               + Q * dv[0]
               + (Q *= Q2) * dv[1]
               + (Q *= Q3) * dv[2]
               + (Q *= Q4) * dv[3]
               + (Q *= Q5) * dv[4]
               + (Q *= Q6) * dv[5]
               + (Q *= Q7) * dv[6];
        poses(0, point) = p_ee.x;
        poses(1, point) = p_ee.y;
        poses(2, point) = p_ee.z;
        poses(3, point) = Q[0];
        poses(4, point) = Q[1];
        poses(5, point) = Q[2];
        poses(6, point) = Q[3];
        poses(7, point) = Q[4];
        poses(8, point) = Q[5];
        poses(9, point) = Q[6];
        poses(10, point) = Q[7];
        poses(11, point) = Q[8];
    }
    return poses;
}

Matrix Gen3::jacobian(const Array &pos) {
    Array s(8);
    Array c(8);
    blast::sincos(pos, s, c);
    Mat3 Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Mat3 Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Mat3 Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
    Mat3 Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    Mat3 Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Mat3 Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
    Mat3 Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};

    // unit vectors in 1st reference
    Vec3 e[7];
    auto Q_tmp = Q1;
    e[0] = Q_tmp * ev[0];
    e[1] = (Q_tmp *= Q2) * ev[1];
    e[2] = (Q_tmp *= Q3) * ev[2];
    e[3] = (Q_tmp *= Q4) * ev[3];
    e[4] = (Q_tmp *= Q5) * ev[4];
    e[5] = (Q_tmp *= Q6) * ev[5];
    e[6] = (Q_tmp *= Q7) * ev[6];

    Vec3 r[7];
    r[6] = dv[6];
    r[5] = dv[5] + Q7 * r[6];
    r[4] = dv[4] + Q6 * r[5];
    r[3] = dv[3] + Q5 * r[4];
    r[2] = dv[2] + Q4 * r[3];
    r[1] = dv[1] + Q3 * r[2];
    r[0] = dv[0] + Q2 * r[1];

    Q_tmp = Q1;
    r[0] = (Q_tmp) * r[0];
    r[1] = (Q_tmp *= Q2) * r[1];
    r[2] = (Q_tmp *= Q3) * r[2];
    r[3] = (Q_tmp *= Q4) * r[3];
    r[4] = (Q_tmp *= Q5) * r[4];
    r[5] = (Q_tmp *= Q6) * r[5];
    r[6] = (Q_tmp *= Q7) * r[6];

    auto cr0 = cross(e[0], r[0]);
    auto cr1 = cross(e[1], r[1]);
    auto cr2 = cross(e[2], r[2]);
    auto cr3 = cross(e[3], r[3]);
    auto cr4 = cross(e[4], r[4]);
    auto cr5 = cross(e[5], r[5]);
    auto cr6 = cross(e[6], r[6]);

    // jacobian matrix
    Matrix J(6, 7);
    J(0, 0) = e[0].x;
    J(1, 0) = e[0].y;
    J(2, 0) = e[0].z;
    J(0, 1) = e[1].x;
    J(1, 1) = e[1].y;
    J(2, 1) = e[1].z;
    J(0, 2) = e[2].x;
    J(1, 2) = e[2].y;
    J(2, 2) = e[2].z;
    J(0, 3) = e[3].x;
    J(1, 3) = e[3].y;
    J(2, 3) = e[3].z;
    J(0, 4) = e[4].x;
    J(1, 4) = e[4].y;
    J(2, 4) = e[4].z;
    J(0, 5) = e[5].x;
    J(1, 5) = e[5].y;
    J(2, 5) = e[5].z;
    J(0, 6) = e[6].x;
    J(1, 6) = e[6].y;
    J(2, 6) = e[6].z;

    J(3, 0) = cr0.x;
    J(4, 0) = cr0.y;
    J(5, 0) = cr0.z;
    J(3, 1) = cr1.x;
    J(4, 1) = cr1.y;
    J(5, 1) = cr1.z;
    J(3, 2) = cr2.x;
    J(4, 2) = cr2.y;
    J(5, 2) = cr2.z;
    J(3, 3) = cr3.x;
    J(4, 3) = cr3.y;
    J(5, 3) = cr3.z;
    J(3, 4) = cr4.x;
    J(4, 4) = cr4.y;
    J(5, 4) = cr4.z;
    J(3, 5) = cr5.x;
    J(4, 5) = cr5.y;
    J(5, 5) = cr5.z;
    J(3, 6) = cr6.x;
    J(4, 6) = cr6.y;
    J(5, 6) = cr6.z;

    return J;
}

Array Gen3::internal_collisions(const Array &joint_position) {
    Array s(7);
    Array c(7);
    blast::sincos(joint_position, s, c);
    Mat3 Q1(c[0], -s[0],   0,    -s[0],  -c[0],   0,    0,   0, -1);
    Mat3 Q2(c[1],   0,    s[1],  -s[1],    0,    c[1],  0,  -1,  0);
    Mat3 Q3(c[2],   0,   -s[2],  -s[2],    0,   -c[2],  0,   1,  0);
    Mat3 Q4(c[3],   0,    s[3],  -s[3],    0,    c[3],  0,  -1,  0);
    Mat3 Q5(c[4],   0,   -s[4],  -s[4],    0,   -c[4],  0,   1,  0);
    Mat3 Q6(c[5],   0,    s[5],  -s[5],    0,    c[5],  0,  -1,  0);
    Mat3 Q7(c[6],   0,   -s[6],  -s[6],    0,   -c[6],  0,   1,  0);
    auto Q = Q1;
    auto p_j2 = p_base + Q * dv[0];
    auto p_j3 = p_j2 + (Q *= Q2) * dv[1];
    auto p_j4 = p_j3 + (Q *= Q3) * dv[2];
    auto p_j5 = p_j4 + (Q *= Q4) * dv[3];
    auto p_j6 = p_j5 + (Q *= Q5) * dv[4];
    auto p_j7 = p_j6 + (Q *= Q6) * dv[5];
    auto p_ee = p_j7 + (Q *= Q7) * dv[6];
    const real r1sqr = 0.0875 * 0.0875;
    const real r2sqr = 0.0875 * 0.0875;
    // Self collisions sqr
    real dist1sqr = two_segment_distance_sqr({0, 0, 0}, p_j2, p_j6, p_ee) - r1sqr;
    // todo: bad way to check a sphere
    real dist2J2sqr = two_segment_distance_sqr(p_j2, p_j2, p_j6, p_ee) - r2sqr; // d^2 from j2 to J6-EE capsule
    real dist2Msqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r1sqr;  // d^2 from j2-j3 capsule to j6-EE capsule
    real dist2sqr = std::min(dist2J2sqr, dist2Msqr);
    // Collision with table sqr
    const real r4table = 0.05; // todo: validate dimensions
    const real r6table = 0.04; // todo: validate dimensions
    const Vec3 p_table(0, 0, -0.0025); // todo: Correct coords (z or y) ??
    real distTJ4 = p_j4.z - p_table.z - r4table;
    real distTJ6 = p_j6.z - p_table.z - r6table;
    real distTEE = p_ee.z - p_table.z - r6table;
    return {dist1sqr, dist2sqr, distTJ4, distTJ6, distTEE};
}

std::vector<Capsule> Gen3::robot_capsules(const Matrix& pos, int n_skip) {
    const int points = pos.cols;
    Matrix result_capsules(21, points/n_skip);
    Mat3 Q, Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Vec3 p_j2, p_j3, p_j4, p_j5, p_j6, p_j7, p_ee;
    Vec3 p_orig(0, 0, 0);
    Array s(7);
    Array c(7);
    int counter = 0;
    for (int i = 0; i < points; i += n_skip) {
        blast::sincos(pos.col(i), s, c);

        Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
        Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
        Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
        Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
        Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};

        Q = Q1;
        p_j2 = p_base + Q * dv[0];
        p_j3 = p_j2 + (Q *= Q2) * dv[1];
        p_j4 = p_j3 + (Q *= Q3) * dv[2];
        p_j5 = p_j4 + (Q *= Q4) * dv[3];
        p_j6 = p_j5 + (Q *= Q5) * dv[4];
        p_j7 = p_j6 + (Q *= Q6) * dv[5];
        p_ee = p_j7 + (Q *= Q7) * dv[6];

        u32 idx = 0;
        // Capsule 1
        result_capsules(0, counter) = p_j2.x;
        result_capsules(1, counter) = p_j2.y;
        result_capsules(2, counter) = p_j2.z;
        result_capsules(3, counter) = p_j4.x;
        result_capsules(4, counter) = p_j4.y;
        result_capsules(5, counter) = p_j4.z;
        result_capsules(6, counter) = 0.055;
        idx += 7;

        // Capsule 2
        result_capsules(idx + 0, counter) = p_j4.x;
        result_capsules(idx + 1, counter) = p_j4.y;
        result_capsules(idx + 2, counter) = p_j4.z;
        result_capsules(idx + 3, counter) = p_j6.x;
        result_capsules(idx + 4, counter) = p_j6.y;
        result_capsules(idx + 5, counter) = p_j6.z;
        result_capsules(idx + 6, counter) = 0.055;
        idx += 7;

        // Capsule 3
        result_capsules(idx + 0, counter) = p_j6.x;
        result_capsules(idx + 1, counter) = p_j6.y;
        result_capsules(idx + 2, counter) = p_j6.z;
        result_capsules(idx + 3, counter) = p_ee.x;
        result_capsules(idx + 4, counter) = p_ee.y;
        result_capsules(idx + 5, counter) = p_ee.z;
        result_capsules(idx + 6, counter) = 0.055;
        counter++;
    }

    auto caps_size = result_capsules.cols;
    std::vector<Capsule> capsules;
    capsules.resize(caps_size * 3); // 3 capsules for each point along the trajectory
    real radius = 0.055; // Hard coded radius of all robot capsules
    for (u32 i = 0; i < caps_size; i++) {
        auto caps_tmp = result_capsules.col(i);
        for (u32 j = 0; j < 3 ; j++) {
            capsules[i*3 + j].p1 = {caps_tmp[0 + 7*j], caps_tmp[1 + 7*j], caps_tmp[2 + 7*j]};
            capsules[i*3 + j].p2 = {caps_tmp[3 + 7*j], caps_tmp[4 + 7*j], caps_tmp[5 + 7*j]};
            capsules[i*3 + j].r = caps_tmp[6 + 7*j];
        }
    }

    return capsules;
}


// Support functions ---------------------------------------------------------------------------------------

real clamped_root(real slope, real h0, real h1) {
//note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    real r;
    if (h0 < 0) {
        if (h1 > 0) {
            r = -h0 / slope;
            if (r > 1)
                r = 0.5;
        }
        else
            r = 1;
    }
    else
        r = 0;
    return r;
}

void compute_intersection(real* sValue, i32* classify, real b, real f00, real f10, i32* edge, real end[][2]) {
// note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    real const zero = 0;
    real const half = (real)0.5;
    real const one = 1;
    if (classify[0] < 0) {
        edge[0] = 0;
        end[0][0] = zero;
        end[0][1] = f00 / b;
        if (end[0][1] < zero || end[0][1] > one)
            end[0][1] = half;
        if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
    }
    else if (classify[0] == 0) {
        edge[0] = 2;
        end[0][0] = sValue[0];
        end[0][1] = zero;
        if (classify[1] < 0) {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
        else if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else {
            edge[1] = 1;
            end[1][0] = one;
            end[1][1] = f10 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
    }
    else {
        edge[0] = 1;
        end[0][0] = one;
        end[0][1] = f10 / b;
        if (end[0][1] < zero || end[0][1] > one)
            end[0][1] = half;
        if (classify[1] == 0) {
            edge[1] = 3;
            end[1][0] = sValue[1];
            end[1][1] = one;
        }
        else {
            edge[1] = 0;
            end[1][0] = zero;
            end[1][1] = f00 / b;
            if (end[1][1] < zero || end[1][1] > one)
                end[1][1] = half;
        }
    }
}

void compute_minimum_parameters(i32* edge, real end[][2], real b, real c, real e, real g00, real g10, real g01, real g11, real* parameter) {
// note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    real const zero = 0;
    real const one = 1;
    real const delta = end[1][1] - end[0][1];
    real h0 = delta * (-b * end[0][0] + c * end[0][1] - e);
    if (h0 >= zero) {
        if (edge[0] == 0) {
            parameter[0] = zero;
            parameter[1] = clamped_root(c, g00, g01);
        }
        else if (edge[0] == 1) {
            parameter[0] = one;
            parameter[1] = clamped_root(c, g10, g11);
        }
        else {
            parameter[0] = end[0][0];
            parameter[1] = end[0][1];
        }
    }
    else {
        real h1 = delta * (-b * end[1][0] + c * end[1][1] - e);
        if (h1 <= zero) {
            if (edge[1] == 0) {
                parameter[0] = zero;
                parameter[1] = clamped_root(c, g00, g01);
            }
            else if (edge[1] == 1) {
                parameter[0] = one;
                parameter[1] = clamped_root(c, g10, g11);
            }
            else {
                parameter[0] = end[1][0];
                parameter[1] = end[1][1];
            }
        }
        else {
            real z = clamp( h0/(h0 - h1), 0, 1 );
            real omz = one - z;
            parameter[0] = omz * end[0][0] + z * end[1][0];
            parameter[1] = omz * end[0][1] + z * end[1][1];
        }
    }
}

real two_segment_distance_sqr(Vec3 P0, Vec3 P1, Vec3 Q0, Vec3 Q1) {
// note: adapted from https://www.geometrictools.com/GTE/Mathematics/DistSegmentSegment.h
    auto const P1mP0 = P1 - P0;
    auto const Q1mQ0 = Q1 - Q0;
    auto const P0mQ0 = P0 - Q0;
    real a = dot(P1mP0, P1mP0);
    real b = dot(P1mP0, Q1mQ0);
    real c = dot(Q1mQ0, Q1mQ0);
    real d = dot(P1mP0, P0mQ0);
    real e = dot(Q1mQ0, P0mQ0);
    real f00 = d;
    real f10 = f00 + a;
    real f01 = f00 - b;
    real f11 = f10 - b;
    real g00 = -e;
    real g10 = g00 - b;
    real g01 = g00 + c;
    real g11 = g10 + c;
    real parameter[2] = {0, 0};
    if (a > 0 && c > 0) {
        real sValue[2] = {
            clamped_root(a, f00, f10),
            clamped_root(a, f01, f11)
        };
        i32 classify[2] = {0, 0};
        for (size_t i = 0; i < 2; ++i) {
            if (sValue[i] <= 0)
                classify[i] = -1;
            else if (sValue[i] >= 1)
                classify[i] = 1;
            else
                classify[i] = 0;
        }
        if (classify[0] == -1 && classify[1] == -1) {
            parameter[0] = 0;
            parameter[1] = clamped_root(c, g00, g01);
        }
        else if (classify[0] == +1 && classify[1] == +1) {
            parameter[0] = 1;
            parameter[1] = clamped_root(c, g10, g11);
        }
        else {
            i32 edge[2] = { 0, 0 };
            real end[2][2];
            compute_intersection(sValue, classify, b, f00, f10, edge, end);
            compute_minimum_parameters(edge, end, b, c, e, g00, g10, g01, g11, parameter);
        }
    }
    else    {
        if (a > 0) {
            parameter[0] = clamped_root(a, f00, f10);
            parameter[1] = 0;
        }
        else     if (c > 0) {
            parameter[0] = 0;
            parameter[1] = clamped_root(c, g00, g01);
        }
        else {
            parameter[0] = 0;
            parameter[1] = 0;
        }
    }
    Vec3 closest0 = P0 + parameter[0]*P1mP0;
    Vec3 closest1 = Q0 + parameter[1]*Q1mQ0;
    Vec3 diff = closest0 - closest1;

    // auto result = sqrt(dot(diff, diff));
    return dot(diff, diff);
}

}