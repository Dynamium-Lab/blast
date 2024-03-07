#pragma once

#include "blast_math.hpp"
#include "blast_trajectory.hpp"
#include "collisions/world.h"

namespace blast {

/*todo:
    - in set_payload:
        - Adjust av_payload
    - in set_payload_without gripper :
        - Adjust av_payload
    - in validate_task;
        - Test that no external objects are inside base shpere
    - in forward_kinematics:
        - Validate pose
    - in inverse_kinematics:
        - test pinv vs pinv_svd
    - in robot_capsule:
        - Create struct for Link6 capsules + sphere
    - in BLAST_ENABLE_TESTS:
        - Add tests

*/

struct Link6 {
    // basic manipulator properties
    int     joints = 6;

    // actuator limits
    Array   pmax; // rad
    Array   vmax; // rad/s
    Array   amax; // rad/s^2
    Array   tau_max; // Nm

    // kinematic properties
    Vec3    p_base;
    Vec3    dv[6]; // vector to next joint
    Vec3    ev[6]; // direction vectors of joint

    // dynamic properties
    real    m[6];  // link masses
    Mat3    I[6];  // Inertial tensors
    Vec3    av[6]; // centers of mass
    Vec3    sv[6]; // centers of mass from next joint

    Matrix  _efforts; // put the efforts temporarily when computing the constraints

    // initialization
    Link6();
    void    set_payload(const real m_payload, const Vec3 cg_payload = {}, const Mat3 I_payload = {});
    void    set_payload_without_gripper(const real m_payload, const Vec3 cg_payload = {}, const Mat3 I_payload = {});

    // optimization
    void    internal_constraints(const Trajectory& traj, real* dst);
    bool    validate_task(const Matrix &task);
    int     ncon(int points);

    // dynamics
    void    dynamics(const Trajectory& traj); // note: results stored in _efforts

    // kinematics
    Array   forward_kinematics(const Array &pos);
    Matrix  forward_kinematics(const Matrix &pos);
    Array   inverse_kinematics(const Array& pose, const Array& initial_joint_position);
    Matrix  jacobian(const Array &joint_position);

    // collisions
    Array   internal_collisions(const Array &joint_position);
    Matrix  robot_capsules(const Array &joint_position);
};


// initialization -------------------

inline Link6::Link6() {
    p_base = {0, 0, 0.0530f};
    pmax.resize(6);
    vmax.resize(6);
    amax.resize(6);
    tau_max.resize(6);
    pmax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad todo: make sure this is true
    vmax = {3.4907f, 3.4907f, 3.4907f, 5.5851f, 5.5851f, 5.5851f}; // rad/s
    amax = {deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600)}; // rad/s^2
    tau_max = {210, 210, 210, 100, 100, 100}; // Nm
    set_payload(0);
}

inline void Link6::set_payload(const real m_payload, const Vec3 cg_payload, const Mat3 I_payload) {
    // link mass
    m[0] = 4.8257f;
    m[1] = 5.9860f;
    m[2] = 3.4159f;
    m[3] = 2.0849f;
    m[4] = 2.0076f;
    m[5] = 1.5193f;

    // inertial tensors
    I[0] = {0.0192746f, -0.00239802f, -0.00896331f, -0.00239802f, 0.03087806f, 0.0016298f, -0.00896331f, 0.0016298f, 0.02134949f};
    I[1] = {0.25899206f, -2.89E-05f, -1.23E-06f, -2.89E-05f, 0.01755445f, -0.02128064f, -1.23E-06f, -0.02128064f, 0.25291674f};
    I[2] = {0.01742043f, -3.55E-06f, 8.4E-07f, -3.55E-06f, 0.01119175f, 0.00518163f, 8.4E-07f, 0.00518163f, 0.01212876f};
    I[3] = {0.02454276f, 2.61E-06f, 1.799E-05f, 2.61E-06f, 0.02385702f, 0.00315758f, 1.799E-05f, 0.00315758f, 0.00294903f};
    I[4] = {0.00734684f, 0.00124927f, -0.00090156f, 0.00124927f, 0.00464684f, -0.00236128f, -0.00090156f, -0.00236128f, 0.00589508f};
    I[5] = {0.00390762f, -1.13E-06f, 1.16E-06f, -1.13E-06f, 0.00390722f, -2.21E-05f, 1.16E-06f, -2.21E-05f, 0.0013928f};

    // center of mass
    av[0] = {0.03930119f, -0.00705889f, -0.08462154f};
    av[1] = {2.53E-06f, 0.18829586f, -0.03988382f};
    av[2] = {4.64E-06f, -0.02451414f, -0.02997969f};
    av[3] = {-0.00010793f, -0.01056422f, -0.08091102f};
    av[4] = {0.01243595f, 0.03284165f, -0.04091434f};
    av[5] = {0.0f, 0.00050624f, -0.00388589f};

    // vector to next joint
    dv[0] = {0.11024, -0.06926, -0.1375};   // 0 -> 1
    dv[1] = {0.0, 0.4850, 0.0};             // 1 -> 2
    dv[2] = {0.0, -0.15216, -0.0917};       // 2 -> 3
    dv[3] = {0.0, -0.06296, -0.22275};      // 3 -> 4
    dv[4] = {0.08703, 0.0860, -0.07692};    // 4 -> 5
    dv[5] = {0.0, 0.0, -0.0920};             // 5 -> endeffector
    Vec3 dv_tool = {0.0, 0.0, - 0.0185 - 0.0185 - 0.163}; // endeffector -> vision + adapter + gripper
    // todo: add option to know if tool is closed or opened (difference of 0.0135 in z)

    // unit joint direction
    ev[0] = {0, 0, 1};
    ev[1] = {0, 0, 1};
    ev[2] = {0, 0, 1};
    ev[3] = {0, 0, 1};
    ev[4] = {0, 0, 1};
    ev[5] = {0, 0, 1};

    // Adjust end effector link with payload
    auto m_ee = m[5];
    real m_vision = 0.4192f;
    real m_adapter = 0.2101f;
    real m_gripper = 0.831f;

    m[5] += m_vision + m_adapter + m_gripper;
    Vec3 av_vision(-0.0094, -0.033, 0.0137);
    Vec3 av_adapter(0.0035, 0.0004, 0.0185 + 0.0055);
    Vec3 av_gripper(0.0, 0.0, 0.0185 + 0.0185 + 0.0473);
    Vec3 av_ee = av[5];
    auto av_new = (1.5193f * av_ee + 0.4192f * av_vision + 0.2101f * av_adapter + 0.831f * av_gripper) * (1 / (m[5]));
    av[5] = av_new;

    Mat3 I_vision   = {0.000502f, 0.000013f, -0.000003f, 0.000013f, 0.000409f, -0.00002f, -0.000003f, -0.00002f, 0.000823f};
    Mat3 I_adapter  = {0.000101f, 0.000000f, -0.000005f, 0.000000f, 0.000014f, 0.000000f, -0.000005f, 0.000000f, 0.000218f};
    Mat3 I_gripper  = {0.001050f, 0.000000f,  0.000000f, 0.000000f, 0.001730f, 0.000000f,  0.000000f, 0.000000f, 0.000990f};

    auto delta_av = av_new - av_ee;
    auto av_to_vision = av_new - av_vision;
    auto av_to_adapter = av_new - av_adapter;
    auto av_to_gripper = av_new - av_gripper;


    I[5](0, 0) += m_ee*delta_av.x*delta_av.x + m_vision*av_to_vision.x*av_to_vision.x + m_adapter*av_to_adapter.x*av_to_adapter.x + m_gripper*av_to_gripper.x*av_to_gripper.x; // todo: Validate
    I[5](1, 1) += m_ee*delta_av.y*delta_av.y + m_vision*av_to_vision.y*av_to_vision.y + m_adapter*av_to_adapter.y*av_to_adapter.y + m_gripper*av_to_gripper.y*av_to_gripper.y;
    I[5](2, 2) += m_ee*delta_av.z*delta_av.z + m_vision*av_to_vision.z*av_to_vision.z + m_adapter*av_to_adapter.z*av_to_adapter.z + m_gripper*av_to_gripper.z*av_to_gripper.z;
    I[5] += I_vision + I_adapter + I_gripper;

    // Modify payload
    Vec3 av_payload = {0.0f, 0.0f, 0.115f}; // todo: good ?
    av_payload += cg_payload;

    auto av_old = av[5];
    auto m_old = m[5];

    auto m_new = m_old + m_payload;
    auto av_new = (m_old*av_old + m_payload*av_payload) / m_new;
    auto delta_av = av_new - av_old; // shift in center of mass
    auto av_to_mass = av_payload - av_new; // vector from payload to new center of mass

    av[5] = av_new;
    m[5] = m_new;
    I[5](0, 0) += m_old*delta_av.x*delta_av.x + m_payload*av_to_mass.x*av_to_mass.x; // todo: Validate
    I[5](1, 1) += m_old*delta_av.y*delta_av.y + m_payload*av_to_mass.y*av_to_mass.y;
    I[5](2, 2) += m_old*delta_av.z*delta_av.z + m_payload*av_to_mass.z*av_to_mass.z;
    I[5] += I_payload;

    // center of mass from next joint
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] + dv_tool - av[5];
}

inline void Link6::set_payload_without_gripper(const real m_payload, const Vec3 cg_payload, const Mat3 I_payload) {
    // link mass
    m[0] = 4.8257f;
    m[1] = 5.9860f;
    m[2] = 3.4159f;
    m[3] = 2.0849f;
    m[4] = 2.0076f;
    m[5] = 1.5193f; // todo: modify with gripper (as of now, no gripper)

    // inertial tensors
    I[0] = {0.0192746f, -0.00239802f, -0.00896331f, -0.00239802f, 0.03087806f, 0.0016298f, -0.00896331f, 0.0016298f, 0.02134949f};
    I[1] = {0.25899206f, -2.89E-05f, -1.23E-06f, -2.89E-05f, 0.01755445f, -0.02128064f, -1.23E-06f, -0.02128064f, 0.25291674f};
    I[2] = {0.01742043f, -3.55E-06f, 8.4E-07f, -3.55E-06f, 0.01119175f, 0.00518163f, 8.4E-07f, 0.00518163f, 0.01212876f};
    I[3] = {0.02454276f, 2.61E-06f, 1.799E-05f, 2.61E-06f, 0.02385702f, 0.00315758f, 1.799E-05f, 0.00315758f, 0.00294903f};
    I[4] = {0.00734684f, 0.00124927f, -0.00090156f, 0.00124927f, 0.00464684f, -0.00236128f, -0.00090156f, -0.00236128f, 0.00589508f};
    I[5] = {0.00390762f, -1.13E-06f, 1.16E-06f, -1.13E-06f, 0.00390722f, -2.21E-05f, 1.16E-06f, -2.21E-05f, 0.0013928f}; // todo: modify with gripper (as of now, no gripper)

    // center of mass
    av[0] = {0.03930119f, -0.00705889f, -0.08462154f};
    av[1] = {2.53E-06f, 0.18829586f, -0.03988382f};
    av[2] = {4.64E-06f, -0.02451414f, -0.02997969f};
    av[3] = {-0.00010793f, -0.01056422f, -0.08091102f};
    av[4] = {0.01243595f, 0.03284165f, -0.04091434f};
    av[5] = {0.0f, 0.00050624f, -0.00388589f};

    // vector to next joint
    dv[0] = {0.11024, -0.06926, -0.1375};   // 0 -> 1
    dv[1] = {0.0, 0.4850, 0.0};             // 1 -> 2
    dv[2] = {0.0, -0.15216, -0.0917};       // 2 -> 3
    dv[3] = {0.0, -0.06296, -0.22275};      // 3 -> 4
    dv[4] = {0.08703, 0.0860, -0.07692};    // 4 -> 5
    dv[5] = {0.0, 0.0, -0.0920};             // 5 -> endeffector

    // unit joint direction
    ev[0] = {0, 0, 1};
    ev[1] = {0, 0, 1};
    ev[2] = {0, 0, 1};
    ev[3] = {0, 0, 1};
    ev[4] = {0, 0, 1};
    ev[5] = {0, 0, 1};

    // Modify payload
    // todo: FIX Vec3 av_payload = {0.0f, 0.0f, -0.115f};
    Vec3 av_payload = {0, 0, 0}; // temporary
    av_payload += cg_payload;

    auto av_old = av[5];
    auto m_old = m[5];

    auto m_new = m_old + m_payload;
    auto av_new = (m_old*av_old + m_payload*av_payload) / m_new;
    auto delta_av = av_new - av_old; // shift in center of mass
    auto av_to_mass = av_payload - av_new; // vector from payload to new center of mass

    av[5] = av_new;
    m[5] = m_new;
    I[5](0, 0) += m_old*delta_av.x*delta_av.x + m_payload*av_to_mass.x*av_to_mass.x; // todo: Validate
    I[5](1, 1) += m_old*delta_av.y*delta_av.y + m_payload*av_to_mass.y*av_to_mass.y;
    I[5](2, 2) += m_old*delta_av.z*delta_av.z + m_payload*av_to_mass.z*av_to_mass.z;
    I[5] += I_payload;

    // center of mass from next joint
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];
}


// optimization -------------------

inline void Link6::internal_constraints(const Trajectory& traj, real* dst) {
    const auto points = traj.pos.cols;
    dynamics(traj);

    for (u32 i = 0; i < points; i++) {
        // todo: self collisions
        auto p = traj.pos.col(i); Assert(p.is_alias);
        auto tmp_coll = internal_collisions(p);
        dst[0] = -tmp_coll[0]; // dist1
        dst[1] = -tmp_coll[1]; // dist2
        dst[2] = -tmp_coll[2]; // dist3
        dst[3] = -tmp_coll[3]; // dist4
        dst[4] = -tmp_coll[4]; // dist5
        dst[5] = -tmp_coll[4]; // dist6
        dst[6] = -tmp_coll[4]; // dist7
        dst += 7;

        // 6 velocity limits
        for (int j = 0; j < (int)joints; j++)
            dst[j] = (abs(traj.vel(j, i)) - vmax[j]) / vmax[j];
        dst += joints;

        // 6 acceleration limits
        for (int j = 0; j < (int)joints; j++)
            dst[j] = (abs(traj.acc(j, i)) - amax[j]) / amax[j];
        dst += joints;

        // 6 torque limits
        auto f = _efforts.col(i); Assert(f.is_alias);
        for (int j = 0; j < (int)joints; j++)
            dst[j] = (abs(f[j]) - tau_max[j]) / tau_max[j];
        dst += joints;
    }
}

inline bool Link6::validate_task(const Matrix &task) {
    Trajectory traj(2, 6);
    traj.pos.col(0) = task.col(0);
    traj.pos.col(1) = task.col(3);

    traj.vel.col(0) = task.col(1);
    traj.vel.col(1) = task.col(4);

    traj.acc.col(0) = task.col(2);
    traj.acc.col(1) = task.col(5);

    Array con(ncon(2));
    internal_constraints(traj, con.data);

    return array_max(con) <= 0;
}

inline int Link6::ncon(int points) {
    return (7 + 6 * 3) * points;
}


// dynamics -------------------

inline void Link6::dynamics(const Trajectory& traj) {
    const auto points = traj.pos.cols;
    const auto joints = traj.pos.rows;
    if (_efforts.cols != points || _efforts.rows != joints)
        _efforts.resize(joints, points);


    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t;
    Vec3 w1, w2, w3, w4, w5, w6;
    Vec3 wd1, wd2, wd3, wd4, wd5, wd6;
    Vec3 cdd0 = {0, 0, 9.81f};
    Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6;
    Vec3 f1, f2, f3, f4, f5, f6;
    Vec3 n1, n2, n3, n4, n5, n6;

    // loop all points
    for (int i = 0; i < points; i++) {
        auto p = traj.pos.col(i);
        auto v = traj.vel.col(i);
        auto a = traj.acc.col(i);

        Array s(6);
        Array c(6);
        blast::sincos(p, s, c);

        // note: these are stored column-wise
        Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
        Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
        Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
        Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0}; // todo: double check
        Q1t = transpose(Q1);
        Q2t = transpose(Q2);
        Q3t = transpose(Q3);
        Q4t = transpose(Q4);
        Q5t = transpose(Q5);
        Q6t = transpose(Q6);

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

        //-- dynamics
        f6 = m[5] * cdd6;
        n6 = I[5] * wd6 + cross(w6, I[5] * w6) + cross(av[5], f6);

        f5 = m[4] * cdd5 + Q6 * f6;
        n5 = I[4] * wd5 + cross(w5, I[4] * w5) + Q6 * n6 + cross(av[4], f5) + cross(sv[4], (Q6 * f6));

        f4 = m[3] * cdd4 + Q5 * f5;
        n4 = I[3] * wd4 + cross(w4, I[3] * w4) + Q5 * n5 + cross(av[3], f4) + cross(sv[3], (Q5 * f5));

        f3 = m[2] * cdd3 + Q4 * f4;
        n3 = I[2] * wd3 + cross(w3, I[2] * w3) + Q4 * n4 + cross(av[2], f3) + cross(sv[2], (Q4 * f4));

        f2 = m[1] * cdd2 + Q3 * f3;
        n2 = I[1] * wd2 + cross(w2, I[1] * w2) + Q3 * n3 + cross(av[1], f2) + cross(sv[1], (Q3 * f3));

        f1 = m[0] * cdd1 + Q2 * f2;
        n1 = I[0] * wd1 + cross(w1, I[0] * w1) + Q2 * n2 + cross(av[0], f1) + cross(sv[0], (Q2 * f2));

        //-- extract torques (last element of each moment vector)
        _efforts(0, i) = n1.z;
        _efforts(1, i) = n2.z;
        _efforts(2, i) = n3.z;
        _efforts(3, i) = n4.z;
        _efforts(4, i) = n5.z;
        _efforts(5, i) = n6.z;
    }
}


// kinematics -------------------

inline Array Link6::forward_kinematics(const Array &pos) {
    Array s(6);
    Array c(6);
    blast::sincos(pos, s, c);
    Mat3 Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Mat3 Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
    Mat3 Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
    Mat3 Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
    Mat3 Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Mat3 Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0};

    auto Q = Q1;
    auto p_ee = p_base
                + Q * dv[0]
                + (Q *= Q2) * dv[1]
                + (Q *= Q3) * dv[2]
                + (Q *= Q4) * dv[3]
                + (Q *= Q5) * dv[4]
                + (Q *= Q6) * dv[5];

    Array pose(6);
    pose[0] = p_ee.x;
    pose[1] = p_ee.y;
    pose[2] = p_ee.z;
    pose[3] = atan2(Q(2, 1), Q(2, 2)) + PI;
    pose[4] = atan2(-Q(2, 0), sqrt(Q(2, 1)*Q(2, 1) + Q(2, 2)*Q(2, 2)));
    pose[5] = atan2(Q(1, 0), Q(0, 0));
    return pose;
}

inline Matrix Link6::forward_kinematics(const Matrix &pos) {
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    Matrix poses(12, pos.cols);
    Vec3 p_ee;
    for (u32 point = 0; point < pos.cols; point++) {
        Array s(8);
        Array c(8);
        blast::sincos(pos.col(point), s, c);
        Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
        Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
        Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
        Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0};
        auto Q = Q1;
        p_ee = p_base
               + Q * dv[0]
               + (Q *= Q2) * dv[1]
               + (Q *= Q3) * dv[2]
               + (Q *= Q4) * dv[3]
               + (Q *= Q5) * dv[4]
               + (Q *= Q6) * dv[5];
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

inline Array Link6::inverse_kinematics(const Array& pose, const Array& initial_joint_position) {

    const double tolerance = 0.001; // Tolerance for convergence
    const int max_iter = 100; // Maximum number of iterations

    Array current_joint_angles = initial_joint_position;

    // Iterate until convergence or maximum iterations reached
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate the current end effector position using forward kinematics
        Array current_pose = forward_kinematics(current_joint_angles);
        Array delta_pose = pose - current_pose;

        // Check if the end effector is close enough to the desired position
        if (is_close(pose, current_pose, tolerance))
            break;

        // Calculate the Jacobian matrix
        Matrix jacobian_matrix = jacobian(current_joint_angles);

        Matrix jacobian_pinv = pinv(jacobian_matrix);

        //current_joint_angles = current_joint_angles + jacobian_pinv * delta_pose;
    }

    return current_joint_angles;
}

inline Matrix Link6::jacobian(const Array &pos) {
    Array s(8);
    Array c(8);
    blast::sincos(pos, s, c);
    Mat3 Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Mat3 Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
    Mat3 Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
    Mat3 Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
    Mat3 Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Mat3 Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0};

    // unit vectors in 1st reference
    Vec3 e[6];
    auto Q_tmp = Q1;
    e[0] = Q_tmp * ev[0];
    e[1] = (Q_tmp *= Q2) * ev[1];
    e[2] = (Q_tmp *= Q3) * ev[2];
    e[3] = (Q_tmp *= Q4) * ev[3];
    e[4] = (Q_tmp *= Q5) * ev[4];
    e[5] = (Q_tmp *= Q6) * ev[5];

    Vec3 r[6];
    r[5] = dv[5];
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


    auto cr0 = cross(e[0], r[0]);
    auto cr1 = cross(e[1], r[1]);
    auto cr2 = cross(e[2], r[2]);
    auto cr3 = cross(e[3], r[3]);
    auto cr4 = cross(e[4], r[4]);
    auto cr5 = cross(e[5], r[5]);

    // jacobian matrix
    Matrix J(6, 6);
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

    return J;
}


// collisions -------------------

inline Array Link6::internal_collisions(const Array &joint_position) {
    // Link6 : size 7x5
    Matrix capsules = robot_capsules(joint_position); // Format : [p1x p1y p1z p2x p2y p2z r; ... (for each capsule)]

    // shpere covering base link
    sphere sph_base;
    sph_base.c = p_base;
    sph_base.r = 0.2375;

    // capsule covering first link
    capsule caps1;
    caps1.p1 = {capsules(0, 0), capsules(1, 0), capsules(2, 0)};
    caps1.p2 = {capsules(3, 0), capsules(4, 0), capsules(5, 0)};
    caps1.r = capsules(6, 0);

    // capsule covering second link
    capsule caps2;
    caps2.p1 = {capsules(0, 1), capsules(1, 1), capsules(2, 1)};
    caps2.p2 = {capsules(3, 1), capsules(4, 1), capsules(5, 1)};
    caps2.r = capsules(6, 1);

    // capsule covering joint 5
    capsule caps3;
    caps3.p1 = {capsules(0, 2), capsules(1, 2), capsules(2, 2)};
    caps3.p2 = {capsules(3, 2), capsules(4, 2), capsules(5, 2)};
    caps3.r = capsules(6, 2);

    // capsule covering last link
    capsule caps4;
    caps4.p1 = {capsules(0, 3), capsules(1, 3), capsules(2, 3)};
    caps4.p2 = {capsules(3, 3), capsules(4, 3), capsules(5, 3)};
    caps4.r = capsules(6, 3);

    // capsule covering gripper & vision module
    capsule caps5;
    caps5.p1 = {capsules(0, 4), capsules(1, 4), capsules(2, 4)};
    caps5.p2 = {capsules(3, 4), capsules(4, 4), capsules(5, 4)};
    caps5.r = capsules(6, 4);

    real dist1 = distmin(caps1, caps3);
    real dist2 = distmin(caps1, caps4);
    real dist3 = distmin(caps1, caps5);

    real dist4 = distmin(caps2, sph_base);
    real dist5 = distmin(caps3, sph_base);
    real dist6 = distmin(caps4, sph_base);
    real dist7 = distmin(caps5, sph_base);


    // Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    // Array s(6);
    // Array c(6);
    // blast::sincos(joint_position, s, c);
    // // note: these are stored column-wise
    // Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    // Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
    // Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
    // Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
    // Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    // Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0}; // todo: double check
    // Vec3 p_orig(0, 0, 0);
    // Vec3 p_j2;
    // Vec3 p_j3;
    // Vec3 p_j4;
    // Vec3 p_j5;
    // Vec3 p_j6;
    // Vec3 p_ee;
    // auto p_tmp = p_base;
    // auto Q_tmp = Q1;
    // p_tmp += Q_tmp * dv[0];
    // p_j2 = p_tmp;
    // p_tmp += (Q_tmp *= Q2) * dv[1];
    // p_j3 = p_tmp;
    // p_tmp += (Q_tmp *= Q3) * dv[2];
    // p_j4 = p_tmp;
    // p_tmp += (Q_tmp *= Q4) * dv[3];
    // p_j5 = p_tmp;
    // p_tmp += (Q_tmp *= Q5) * dv[4];
    // p_j6 = p_tmp;
    // p_tmp += (Q_tmp *= Q6) * dv[5];
    // p_ee = p_tmp;
    // const real r2sqr = 0.140 * 0.140; // size of capsule 2
    // const real r3sqr = 0.100 * 0.100; // size of capsule 3


    return {dist1, dist2, dist3, dist4, dist5, dist6, dist7};
}

inline Matrix Link6::robot_capsules(const Array &joint_position) {
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;

    Array s(6);
    Array c(6);
    blast::sincos(joint_position, s, c);

    // note: these are stored column-wise
    Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
    Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
    Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
    Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0}; // todo: double check

    Vec3 p_orig(0, 0, 0);
    Vec3 p_j2;
    Vec3 p_j3;
    Vec3 p_j4;
    Vec3 p_j5;
    Vec3 p_j6;
    Vec3 p_ee;

    Vec3 z2;
    Vec3 z3;
    Vec3 z4;
    Vec3 z5;
    Vec3 z6;
    Vec3 zee;
    Vec3 yee;

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_j2 = p_tmp;
    z2 = Q_tmp.col(2);
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_j3 = p_tmp;
    z3 = Q_tmp.col(2);
    p_tmp += (Q_tmp *= Q3) * dv[2];
    // p_j4 = p_tmp;
    z4 = Q_tmp.col(2);
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_j5 = p_tmp;
    z5 = Q_tmp.col(2);
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_j6 = p_tmp;
    z6 = Q_tmp.col(2);
    p_tmp += (Q_tmp *= Q6) * dv[5];
    zee = Q_tmp.col(2);
    yee = Q_tmp.col(1);
    p_ee = p_tmp;

    Matrix capsules(7, 5);

    Vec3 p1;
    Vec3 p2;

    // Capsule 1
    p1 = p_j2 - 0.02218*z2 + 0.00010104*(p_j3-p_j2);
    p2 = p1 + 0.8845*(p_j3-p_j2);
    capsules(0, 0) = p1.x;
    capsules(1, 0) = p1.y;
    capsules(2, 0) = p1.z;
    capsules(3, 0) = p2.x;
    capsules(4, 0) = p2.y;
    capsules(5, 0) = p2.z;
    capsules(6, 0) = 0.110; // radius

    // Capsule 2
    p1 = p_j2 - 0.11388*z3;
    p2 = p1 - 0.375*z4;
    capsules(0, 1) = p1.x;
    capsules(1, 1) = p1.y;
    capsules(2, 1) = p1.z;
    capsules(3, 1) = p2.x;
    capsules(4, 1) = p2.y;
    capsules(5, 1) = p2.z;
    capsules(6, 1) = 0.061; // radius

    // Capsule 3
    p1 = p_j5;
    p2 = p1 - 0.08*z5;
    capsules(0, 2) = p1.x;
    capsules(1, 2) = p1.y;
    capsules(2, 2) = p1.z;
    capsules(3, 2) = p2.x;
    capsules(4, 2) = p2.y;
    capsules(5, 2) = p2.z;
    capsules(6, 2) = 0.060; // radius

    // Capsule 4
    p1 = p_j6 + 0.08583*z6;
    p2 = p1 - 0.15*z6;
    capsules(0, 3) = p1.x;
    capsules(1, 3) = p1.y;
    capsules(2, 3) = p1.z;
    capsules(3, 3) = p2.x;
    capsules(4, 3) = p2.y;
    capsules(5, 3) = p2.z;
    capsules(6, 3) = 0.060; // radius

    // Capsule 5
    p1 = p_ee - 0.01289*zee - 0.2125*yee;
    p2 = p1 - 0.150*zee;
    capsules(0, 3) = p1.x;
    capsules(1, 3) = p1.y;
    capsules(2, 3) = p1.z;
    capsules(3, 3) = p2.x;
    capsules(4, 3) = p2.y;
    capsules(5, 3) = p2.z;
    capsules(6, 3) = 0.085; // radius

    return capsules;
}



#ifdef BLAST_ENABLE_TESTS
TEST_CASE("SelfCollisionGen3", "[Manipulator]") {
    Gen3 manip;
    Array theta1(7);
    theta1 = {0, 15, 180, 230, 360, 55, 90};
    theta1 = deg2rad(theta1);
    Array theta2(7);
    theta2 = {0, 0, 0, 0, 0, 0, 0};
    theta2 = deg2rad(theta2);
    Array theta3(7);
    theta3 = {0, 16, 180, 221, 358, 284, 88};
    theta3 = deg2rad(theta3);
    Array theta4(7);
    theta4 = {347, 47, 158, 212, 341, 300, 8};
    theta4 = deg2rad(theta4);

    auto dist_sqr_min_1 = manip.internal_collisions(theta1);
    auto dist_sqr_min_2 = manip.internal_collisions(theta2);
    auto dist_sqr_min_3 = manip.internal_collisions(theta3);
    auto dist_sqr_min_4 = manip.internal_collisions(theta4);

    CHECK(dist_sqr_min_1[0] > 0);
    CHECK(dist_sqr_min_1[1] > 0);
    CHECK(dist_sqr_min_2[0] > 0);
    CHECK(dist_sqr_min_2[1] > 0);
    CHECK(dist_sqr_min_3[0] < 0);
    CHECK(dist_sqr_min_3[1] < 0);
    CHECK(dist_sqr_min_4[0] < 0);
    CHECK(dist_sqr_min_4[1] > 0);
}

#endif // tests
};
