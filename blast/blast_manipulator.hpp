#pragma once

#include "blast_math.hpp"
// #include "dynamics/blast_dynamics_MDA.hpp"
#include <vector>
#include <cmath>

namespace blast {
using std::vector;

struct Manipulator {
    u32 joints;
    Array pmax;
    Array pmin;
    Array vmax;
    Array vmin;
    Array amax;
    Array amin;
    Array tau_max;
    Array tau_min;

    Manipulator(u32 njoints) : joints(njoints), vmax(njoints), vmin(njoints), pmax(njoints), pmin(njoints), amax(njoints), amin(njoints), tau_max(njoints), tau_min(njoints) {}

    // todo: remove?
    virtual Array constraints(const Trajectory &traj) {
        return Array();
    }

    virtual void constraints(const Trajectory &traj, real *dst) {}

    // todo: remove?
    virtual bool validate_task(const Matrix &task) {
        return false;
    }

    // return the number of constraints with the given 'npoints'
    virtual u32 ncon(u32 npoints) {
        return 0;
    }
};

struct ManipulatorEigen {
    int joints;
    VectorXd pmax;
    VectorXd pmin;
    VectorXd vmax;
    VectorXd vmin;
    VectorXd amax;
    VectorXd amin;
    VectorXd tau_max;
    VectorXd tau_min;
    ManipulatorEigen(int njoints) : joints(njoints), vmax(njoints), vmin(njoints), pmax(njoints), pmin(njoints), amax(njoints), amin(njoints), tau_max(njoints), tau_min(njoints) {}
    // todo: remove?
    virtual VectorXd constraints(const TrajectoryEigen &traj) {
        return VectorXd();
    }
    virtual void constraints(const TrajectoryEigen &traj, real *dst) {}
    // todo: remove?
    virtual bool validate_task(const MatrixXd &task) {
        return false;
    }
    // return the number of constraints with the given 'npoints'
    virtual int ncon(int npoints) {
        return 0;
    }
};

struct ManipulatorUR5 : public Manipulator {
    bool is_init = false;

    // DH parameters
    real a2, a3;
    real d4, d5, d6;
    Vec3 base_p;
    real L2, L3, L5;

    // dynamics parameters
    real m1, m2, m3, m4, m5, m6;
    Mat3 I1, I2, I3, I4, I5, I6;
    Vec3 a1v, a2v, a3v, a4v, a5v, a6v;
    Vec3 s1, s2, s3, s4, s5, s6;
    Vec3 e1, e2, e3, e4, e5, e6;

    ManipulatorUR5() : Manipulator(6) {}

    void dynamics(const Trajectory &traj, Matrix &efforts);
    void dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc, Matrix &efforts);
    void init_dynamics(real mass = 0);
};

struct Gen3Lite : public Manipulator {
    Vec3 p_base;
    real m[6];
    Mat3 I[6];
    Vec3 av[6];
    Vec3 dv[6];
    Vec3 sv[6];
    Vec3 ev[6];
    // default constructor
    Gen3Lite();
    // compute joint torque as a function of trajector (traj)
    void dynamics(const Trajectory &traj, Matrix &efforts);
    void dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc, Matrix &efforts);
    // compute forward kinematics for 1 point
    Array forward_kinematics(Array& joint_position);
    Matrix jacobian(const Array& joint_position);
};

struct Gen3_7DOF : public Manipulator {
    Vec3 p_base;
    real m[7];
    Mat3 I[7];
    Vec3 av[7];
    Vec3 dv[7];
    Vec3 sv[7];
    Vec3 ev[7];

    Matrix _efforts; // put the efforts temporarily when computing the constraints

    // default constructor
    Gen3_7DOF();

    // compute new center of mass and robot mass with payload
    void set_payload(const real mass);
    void set_payload_without_gripper(const real mass);

    // compute joint torque as a function of trajector (traj)
    void dynamics(const Trajectory &traj, Matrix &efforts);
    void dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc, Matrix &efforts);
    void dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc); // note: cache the results internally

    // compute forward kinematics for 1 point
    Array forward_kinematics(const Array &joint_position);
    Matrix forward_kinematics(const Matrix &joint_positions);

    Array inverse_kinematics(const Array& pose, const Array& initial_joint_position);

    // compute jacobian matrix
    Matrix jacobian(const Array &joint_position);

    // check collision
    Array collision_check(const Array &joint_position);

    // check all constraints on the manipulator for 1 point
    Array constraints(const Array &pos, const Array &vel, const Array &acc);

    // check all constraints on the manipulator for a trajectory
    Array constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc);

    // check all constraints and put the result in 'dst'
    void constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc, real *dst);

    // check all constraints on the manipulator for a trajectory.
    // - note: any contraint that is positive is not respected.
    virtual Array constraints(const Trajectory &traj) override;

    // check all constraints and put the result in 'dst'
    // - note: any contraint that is positive is not respected.
    virtual void constraints(const Trajectory &traj, real *dst) override;

    // check that the task given is feasible (collision and joint limits)
    virtual bool validate_task(const Matrix &task) override;
    virtual Matrix robot_capsules(const Matrix &joint_position, const int n_skip);

    virtual u32 ncon(u32 npoints) override {
        // 5 collision results
        // 2 position constraints
        // 7 velocity constraints
        // 7 torque constraints
        //** (for each point in the trajectory) **
        return (5 + 2 + 7 * 2) * npoints;
    }
};

struct Gen3Eigen : public ManipulatorEigen {
    Vector3d p_base;
    double   m[7];
    Matrix3d I[7];
    Vector3d av[7];
    Vector3d dv[7];
    Vector3d sv[7];
    Vector3d ev[7];

    MatrixXd _efforts; // put the efforts temporarily when computing the constraints

    // default constructor
    Gen3Eigen();

    // compute new center of mass and robot mass with payload
    void set_payload(double mass = 0, bool add_gripper = false);

    // compute joint torque as a function of trajector (traj)
    void dynamics(const TrajectoryEigen &traj, MatrixXd &efforts);
    void dynamics(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc, MatrixXd &efforts);
    void dynamics(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc); // note: cache the results internally

    // compute forward kinematics for 1 point
    VectorXd forward_kinematics(const VectorXd &joint_position);

    // todo: this is not implemented!!
    VectorXd inverse_kinematics(const VectorXd& pose, const VectorXd& initial_joint_position);

    // compute jacobian matrix
    MatrixXd jacobian(const VectorXd &joint_position);

    // check collision
    VectorXd collision_check(const VectorXd &joint_position);

    // check all constraints on the manipulator for 1 point
    VectorXd constraints(const VectorXd &pos, const VectorXd &vel, const VectorXd &acc);

    // check all constraints on the manipulator for a trajectory
    VectorXd constraints(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc);

    // check all constraints and put the result in 'dst'
    void constraints(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc, double *dst);

    // check all constraints on the manipulator for a trajectory.
    // - note: any contraint that is positive is not respected.
    virtual VectorXd constraints(const TrajectoryEigen &traj) override;

    // check all constraints and put the result in 'dst'
    // - note: any contraint that is positive is not respected.
    virtual void constraints(const TrajectoryEigen &traj, double *dst) override;

    // check that the task given is feasible (collision and joint limits)
    virtual bool validate_task(const MatrixXd &task) override;

    virtual int ncon(int npoints) override {
        // 5 collision results
        // 2 position constraints
        // 7 velocity constraints
        // 7 torque constraints
        //** (for each point in the trajectory) **
        return (5 + 2 + 7 * 2) * npoints;
    }
};

struct Link6 : public Manipulator {
    Vec3 p_base;
    real m[6];
    Mat3 I[6];
    Vec3 av[6];
    Vec3 dv[6];
    Vec3 sv[6];
    Vec3 ev[6];

    Matrix _efforts; // put the efforts temporarily when computing the constraints

    // default constructor
    Link6();

    // // compute new center of mass and robot mass with payload
    // void set_payload(const real mass);
    // void set_payload_without_gripper(const real mass);

    // // compute joint torque as a function of trajector (traj)
    void dynamics(const Trajectory &traj, Matrix &efforts);
    void dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc, Matrix &efforts);
    void dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc); // note: cache the results internally

    // compute forward kinematics for 1 point
    Array forward_kinematics(const Array &joint_position);
    Matrix forward_kinematics(const Matrix &joint_positions);

    // Array inverse_kinematics(const Array& pose, const Array& initial_joint_position);

    // // compute jacobian matrix
    // Matrix jacobian(const Array &joint_position);

    // check collision
    Array collision_check(const Array &joint_position);

    // check all constraints on the manipulator for 1 point
    Array constraints(const Array &pos, const Array &vel, const Array &acc);

    // check all constraints on the manipulator for a trajectory
    Array constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc);

    // check all constraints and put the result in 'dst'
    void constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc, real *dst);

    // check all constraints on the manipulator for a trajectory.
    // - note: any contraint that is positive is not respected.
    virtual Array constraints(const Trajectory &traj) override;

    // check all constraints and put the result in 'dst'
    // - note: any contraint that is positive is not respected.
    virtual void constraints(const Trajectory &traj, real *dst) override;

    // check that the task given is feasible (collision and joint limits)
    virtual bool validate_task(const Matrix &task) override;

    virtual u32 ncon(u32 npoints) override {
        // 5 collision results
        // 0 position constraints
        // 6 velocity constraints
        // 6 torque constraints
        //** (for each point in the trajectory) **
        return (5 + 6 * 2) * npoints;
    }
};

//------ Universal Robots UR5e manipulator functions ---------------------------------

host_fn void ManipulatorUR5::dynamics(const Trajectory &traj, Matrix &efforts) {
    const auto points = traj.t.size;
    Assert(is_init);

    real vel1, vel2, vel3, vel4, vel5, vel6;
    real acc1, acc2, acc3, acc4, acc5, acc6;

    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t;
    Vec3 w12, w23, w34, w45, w56, w67;
    Vec3 wd12, wd23, wd34, wd45, wd56, wd67;
    Vec3 cdd01, cdd12, cdd23, cdd34, cdd45, cdd56, cdd67;
    cdd01.z = 9.81f;
    Vec3 f1, f2, f3, f4, f5, f6;
    Vec3 n1, n2, n3, n4, n5, n6;

    // loop all points
    for (u32 i = 0; i < points; i++) {
        vel1 = traj.vel(0, i);
        vel2 = traj.vel(1, i);
        vel3 = traj.vel(2, i);
        vel4 = traj.vel(3, i);
        vel5 = traj.vel(4, i);
        vel6 = traj.vel(5, i);
        acc1 = traj.acc(0, i);
        acc2 = traj.acc(1, i);
        acc3 = traj.acc(2, i);
        acc4 = traj.acc(3, i);
        acc5 = traj.acc(4, i);
        acc6 = traj.acc(5, i);

        // SIMD compute sines and cosines note: approx 10% faster
        // const auto p = &traj.pos.data[i * joints];

        const auto p = traj.pos.col(i);
        Array s(p.size);
        Array c(p.size);
        blast::sincos(p, s, c);

// #if BLAST_SIZEOF_REAL == 8
//         real s[6];
//         real c[6];
//         // first 4
//         {
//             __m256d s_tmp;
//             __m256d c_tmp;
//             __m256d angle_v = _mm256_load_pd(p);
//             s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//             _mm256_storeu_pd(s, s_tmp);
//             _mm256_storeu_pd(c, c_tmp);
//         }
//         // 5 and 6
//         {
//             __m128d angle_v = _mm_load_pd(p + 4);
//             __m128d s_tmp;
//             __m128d c_tmp;
//             s_tmp = _mm_sincos_pd(&c_tmp, angle_v);
//             _mm_storeu_pd(s + 4, s_tmp);
//             _mm_storeu_pd(c + 4, c_tmp);
//         }
// #else
//         real s[8];
//         real c[8];
//         __m256 s_tmp;
//         __m256 c_tmp;
//         __m256 angle_v = _mm256_load_ps(p);
//         s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//         _mm256_storeu_ps(s, s_tmp);
//         _mm256_storeu_ps(c, c_tmp);
// #endif

        Q1(0, 0) = c[0];
        Q1(1, 0) = s[0];
        Q1(2, 1) = 1.0;
        Q1(0, 2) = s[0];
        Q1(1, 2) = -c[0];

        Q2(0, 0) = -s[1];
        Q2(1, 0) = c[1];
        Q2(0, 1) = -c[1];
        Q2(1, 1) = -s[1];
        Q2(2, 2) = 1.0;

        Q3(0, 0) = c[2];
        Q3(1, 0) = s[2];
        Q3(0, 1) = -s[2];
        Q3(1, 1) = c[2];
        Q3(2, 2) = 1.0;

        Q4(0, 0) = -s[3];
        Q4(1, 0) = c[3];
        Q4(2, 1) = 1.0;
        Q4(0, 2) = c[3];
        Q4(1, 2) = s[3];

        Q5(0, 0) = c[4];
        Q5(1, 0) = s[4];
        Q5(2, 1) = 1.0;
        Q5(0, 2) = s[4];
        Q5(1, 2) = -c[4];

        Q6(0, 0) = c[5];
        Q6(1, 0) = s[5];
        Q6(0, 1) = -s[5];
        Q6(1, 1) = c[5];
        Q6(2, 2) = 1.0;

        Q1t = transpose(Q1);
        Q2t = transpose(Q2);
        Q3t = transpose(Q3);
        Q4t = transpose(Q4);
        Q5t = transpose(Q5);
        Q6t = transpose(Q6);

        // note: this is the Newton algorithm in 'Element de robotique' course notes.
        //-- kinematics
        w12 = vel1 * e1;
        w23 = Q2t * w12 + vel2 * e2;
        w34 = Q3t * w23 + vel3 * e3;
        w45 = Q4t * w34 + vel4 * e4;
        w56 = Q5t * w45 + vel5 * e5;
        w67 = Q6t * w56 + vel6 * e6;

        wd12 = acc1 * e1;
        cdd12 = Q1t * cdd01 + cross(wd12, (a1v + s1)) + cross(w12, cross(w12, (a1v + s1)));

        wd23 = Q2t * wd12 + acc2 * e2 + vel2 * cross((Q2t * w12), (e2));
        cdd23 = Q2t * cdd12 + cross(wd23, a2v + s2) + cross(w23, cross(w23, (a2v + s2))) - Q2t * cross(wd12, s1) - Q2t * cross(w12, cross(w12, s1));

        wd34 = Q3t * wd23 + acc3 * e3 + vel3 * cross((Q3t * w23), (e3));
        cdd34 = Q3t * cdd23 + cross(wd34, a3v + s3) + cross(w34, cross(w34, (a3v + s3))) - Q3t * cross(wd23, s2) - Q3t * cross(w23, cross(w23, s2));

        wd45 = Q4t * wd34 + acc4 * e4 + vel4 * cross((Q4t * w34), (e4));
        cdd45 = Q4t * cdd34 + cross(wd45, a4v + s4) + cross(w45, cross(w45, (a4v + s4))) - Q4t * cross(wd34, s3) - Q4t * cross(w34, cross(w34, s3));

        wd56 = Q5t * wd45 + acc5 * e5 + vel5 * cross((Q5t * w45), (e5));
        cdd56 = Q5t * cdd45 + cross(wd56, a5v + s5) + cross(w56, cross(w56, (a5v + s5))) - Q5t * cross(wd45, s4) - Q5t * cross(w45, cross(w45, s4));

        wd67 = Q6t * wd56 + acc6 * e6 + vel6 * cross((Q6t * w56), (e6));
        cdd67 = Q6t * cdd56 + cross(wd67, a6v + s6) + cross(w67, cross(w67, (a6v + s6))) - Q6t * cross(wd56, s5) - Q6t * cross(w56, cross(w56, s5));

        //-- dynamics
        f6 = Q6 * (m6 * cdd67);
        n6 = Q6 * (I6 * wd67 + cross(w67, I6 * w67) + cross(a6v + s6, Q6t * f6));

        f5 = Q5 * (m5 * cdd56 + f6);
        n5 = Q5 * (I5 * wd56 + cross(w56, I5 * w56) + n6 + cross(a5v + s5, Q5t * f5) - cross(s5, f6));

        f4 = Q4 * (m4 * cdd45 + f5);
        n4 = Q4 * (I4 * wd45 + cross(w45, I4 * w45) + n5 + cross(a4v + s4, Q4t * f4) - cross(s4, f5));

        f3 = Q3 * (m3 * cdd34 + f4);
        n3 = Q3 * (I3 * wd34 + cross(w34, I3 * w34) + n4 + cross(a3v + s3, Q3t * f3) - cross(s3, f4));

        f2 = Q2 * (m2 * cdd23 + f3);
        n2 = Q2 * (I2 * wd23 + cross(w23, I2 * w23) + n3 + cross(a2v + s2, Q2t * f2) - cross(s2, f3));

        f1 = Q1 * (m1 * cdd12 + f2);
        n1 = Q1 * (I1 * wd12 + cross(w12, I1 * w12) + n2 + cross(a1v + s1, Q1t * f1) - cross(s1, f2));

        //-- extract torques (last element of each moment vector)
        efforts(0, i) = n1.z;
        efforts(1, i) = n2.z;
        efforts(2, i) = n3.z;
        efforts(3, i) = n4.z;
        efforts(4, i) = n5.z;
        efforts(5, i) = n6.z;
    }
}

host_fn void ManipulatorUR5::init_dynamics(real mass) {
    is_init = true;

    base_p.z = 0.089159f;

    tau_max = {
        1,
        1,
        1,
        1,
        1,
        1
    };
    tau_min = -tau_max;

    // kinematic parameters
    a2 = 0.425f;
    a3 = 0.39225f;
    d4 = -0.10915f;
    d5 = 0.09465f;
    d6 = 0.0823f;

    // link mass
    m1 = 3.7f;
    m2 = 8.393f;
    m3 = 2.33f;
    m4 = 1.219f;
    m5 = 1.219f;
    m6 = 0.1879f;

    // vectors from one joint towards the next
    a2v.x = a2;
    a3v.x = a3;
    a4v.y = d4;
    a5v.y = d5;
    a6v.z = d6;

    // unit vector of joint orientation (elements de robotique ch.5, p.176)
    e1.y = 1;
    e2.z = 1;
    e3.z = 1;
    e4.y = 1;
    e5.y = 1;
    e6.z = 1;

    // Vector from next joint to center of mass (elements de robotique ch.5, p.174)
    s1.y = -0.02561f;
    s1.z = -0.00193f;
    s2.x = -0.2125f;
    s2.z = -0.11336f;
    s3.x = -0.15f;
    s3.z = -0.0265f;
    s4.y = 0.0018f;
    s4.z = 0.01634f;
    s5.y = -0.0018f;
    s5.z = 0.01634f;
    s6.z = -0.001159f;

    // interial tensors
    I1(0, 0) = 0.007f;
    I1(1, 1) = 0.007f;
    I1(2, 2) = 0.007f;

    I2(0, 0) = 0.015f;
    I2(1, 1) = 0.253f;
    I2(2, 2) = 0.253f;

    I3(0, 0) = 0.0033f;
    I3(1, 1) = 0.06f;
    I3(2, 2) = 0.06f;

    I4(0, 0) = 0.015f;
    I4(1, 1) = 0.015f;
    I4(2, 2) = 0.015f;

    I5(0, 0) = 0.015f;
    I5(1, 1) = 0.015f;
    I5(2, 2) = 0.015f;

    I6(0, 0) = 0.000082f;
    I6(1, 1) = 0.000082f;
    I6(2, 2) = 0.000082f;

    // define new center of mass of the last link (with added mass)
    if (mass) {
        Vec3 s6_old = s6;
        real m6_old = m6;
        m6 += mass;
        s6 *= m6_old / m6;

        I6(0, 0) += m6_old * (s6.x - s6_old.x) * (s6.x - s6_old.x) + mass * s6.x * s6.x;
        I6(1, 1) += m6_old * (s6.y - s6_old.y) * (s6.y - s6_old.y) + mass * s6.y * s6.y;
        I6(2, 2) += m6_old * (s6.z - s6_old.z) * (s6.z - s6_old.z) + mass * s6.z * s6.z;
    }
}

//------ Kinova Gen3 Lite manipulator functions --------------------------------------

host_fn Gen3Lite::Gen3Lite() : Manipulator(6) {
    // position of the first joint with respect to the table in the center of the base
    p_base = {0, 0, 0.1283f};

    // mass
    m[0] = 0.959744f;
    m[1] = 1.177561f;
    m[2] = 0.597676f;
    m[3] = 0.526934f;
    m[4] = 0.580973f;
    m[5] = 0.2018f;

    // inertial tensor
    I[0] = {0.0016595f, 2e-08f, 3.6e-07f, 2e-08f, 0.00140355f, 0.00034927f, 3.6e-07f, 0.00034927f, 0.00089493f};
    I[1] = {0.0114928f, 1e-06f, 1.6e-07f, 1e-06f, 0.00102851f, 0.00140765f, 1.6e-07f, 0.00140765f, 0.01133492f};
    I[2] = {0.0016326f, 7.11E-06f, 1.54e-06f, 7.11E-06f, 0.00029798f, 9.587E-05f, 1.54E-06f, 9.587E-05f, 0.00169091f};
    I[3] = {0.0006910f, 2.4E-07f, 0.0001648f, 2.4E-07f, 0.00078519f, 7.4E-07f, 0.00016483f, 7.4E-07f, 0.00034115f};
    I[4] = {0.0002127f, 5.21E-06f, 2.91E-06f, 5.21E-06f, 0.00106371f, 1.1E-07f, 2.91E-06f, 1.1E-07f, 0.00108465f};
    I[5] = {0.0003428f, 1.9e-7f, 1e-7f, 1.9e-7f, 0.00028915f, 0.00000027f, 1e-7f, 0.00000027f, 0.00013076f};

    // center of mass
    av[0] = {0.000025f, 0.022135f, 0.099377f};
    av[1] = {0.029983f, 0.211548f, 0.045303f};
    av[2] = {0.030156f, -0.095022f, 0.007356f};
    av[3] = {0.005752f, 0.010004f, 0.087192f};
    av[4] = {0.080565f, 0.009804f, 0.018728f};
    av[5] = {0.009930f, 0.009950f, 0.061360f};

    // vector to next joint
    // dv[0] = {0.000f,   0.000f,   0.1283f}; // todo: remove and add 0.1300f to dv[5].z
    dv[0] = {0.000f,    -0.030f,   0.115f};
    dv[1] = {0.000f,    0.280f,   0.000f};
    dv[2] = {0.000f,   -0.140f,   0.020f};
    dv[3] = {0.0285f,   0.000f,   0.105f};
    dv[4] = {-0.105f,    0.000f,   0.0285f};
    dv[5] = {0.000f,     0.000f,  0.1300f};

    // center of mass (from next joint)
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];

    // unit joint direction
    ev[0] = { 0,  0,  1};
    ev[1] = { 0,  0,  1};
    ev[2] = { 0,  0,  1};
    ev[3] = { 0,  0,  1};
    ev[4] = { 0,  0,  1};
    ev[5] = { 0,  0,  1};

    // kinematic and dynamic constraints
    pmax = {2.69f, 2.69f, 2.69f, 2.59f, 2.57f, 2.59f}; // rad
    vmax = {1.6f, 1.6f, 1.6f, 1.6f, 1.6f, 3.2f};       // rad/s
    tau_max = {10, 14, 10, 7, 7, 7};                   // Nm
    vmin = -vmax;
    pmin = -pmax;
    tau_min = -tau_max;
}

host_fn void Gen3Lite::dynamics(const Trajectory &traj, Matrix &efforts) {
    dynamics(traj.pos, traj.vel, traj.acc, efforts);
}

host_fn void Gen3Lite::dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc, Matrix &efforts) {

    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t;
    Vec3 w12, w23, w34, w45, w56, w67;
    Vec3 wd12, wd23, wd34, wd45, wd56, wd67;
    Vec3 cdd01, cdd12, cdd23, cdd34, cdd45, cdd56, cdd67;
    cdd01.z = 9.81f;
    Vec3 f1, f2, f3, f4, f5, f6;
    Vec3 n1, n2, n3, n4, n5, n6;

    const u32 joints = pos.rows;

    // loop all points
    for (u32 i = 0; i < pos.cols; i++) {
        auto v = &vel.data[i * joints];
        auto a = &acc.data[i * joints];

        // SIMD compute sines and cosines note: approx 10% faster
        const auto p = pos.col(i);
        Array s(8);
        Array c(8);
        blast::sincos(p, s, c);
//         auto p = &pos.data[i * joints];
// #if BLAST_SIZEOF_REAL == 8 // double precision
//         real s[6];
//         real c[6];
//         // first 4
//         {
//             __m256d s_tmp;
//             __m256d c_tmp;
//             __m256d angle_v = _mm256_load_pd(p);
//             s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//             _mm256_storeu_pd(s, s_tmp);
//             _mm256_storeu_pd(c, c_tmp);
//         }
//         // 5 and 6
//         {
//             __m128d angle_v = _mm_load_pd(p + 4);
//             __m128d s_tmp;
//             __m128d c_tmp;
//             s_tmp = _mm_sincos_pd(&c_tmp, angle_v);
//             _mm_storeu_pd(s + 4, s_tmp);
//             _mm_storeu_pd(c + 4, c_tmp);
//         }
// #else
//         real s[8];
//         real c[8];
//         __m256 s_tmp;
//         __m256 c_tmp;
//         __m256 angle_v = _mm256_load_ps(p);
//         s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//         _mm256_storeu_ps(s, s_tmp);
//         _mm256_storeu_ps(c, c_tmp);
// #endif

        Q1 = {c[0],   s[0],     0,        -s[0],     c[0],      0,              0,      0,    1};
        Q2 = {c[1],      0,  s[1],        -s[1],        0,   c[1],              0,     -1,    0};
        Q3 = {c[2],  -s[2],     0,        -s[2],    -c[2],      0,              0,      0,   -1};
        Q4 = {c[3],      0,  s[3],        -s[3],        0,   c[3],              0,     -1,    0};
        Q5 = {   0,   s[4], -c[4],           0,      c[4],     s[4],            1,      0,    0};
        Q6 = {   0,   s[5],  c[5],            0,     c[5],  -s[5],             -1,      0,    0};
        Q1t = transpose(Q1);
        Q2t = transpose(Q2);
        Q3t = transpose(Q3);
        Q4t = transpose(Q4);
        Q5t = transpose(Q5);
        Q6t = transpose(Q6);

        // note: this is the Newton algorithm in 'Element de robotique' course notes.
        //-- kinematics

        w12 =  v[0] * ev[0];
        w23 = Q2t * w12 + v[1] * ev[1];
        w34 = Q3t * w23 + v[2] * ev[2];
        w45 = Q4t * w34 + v[3] * ev[3];
        w56 = Q5t * w45 + v[4] * ev[4];
        w67 = Q6t * w56 + v[5] * ev[5];

        wd12 = a[0] * ev[0];
        cdd12 = Q1t * cdd01 + cross(wd12, av[0]) + cross(w12, cross(w12, av[0]));

        wd23 = Q2t * wd12 + a[1] * ev[1] + v[1] * cross(Q2t * w12, ev[1]);
        cdd23 = Q2t * cdd12 + cross(wd23, av[1]) + cross(w23, cross(w23, av[1])) - Q2t * cross(wd12, sv[0]) - Q2t * cross(w12, cross(w12, sv[0]));

        wd34 = Q3t * wd23 + a[2] * ev[2] + v[2] * cross(Q3t * w23, ev[2]);
        cdd34 = Q3t * cdd23 + cross(wd34, av[2]) + cross(w34, cross(w34, av[2])) - Q3t * cross(wd23, sv[1]) - Q3t * cross(w23, cross(w23, sv[1]));

        wd45 = Q4t * wd34 + a[3] * ev[3] + v[3] * cross(Q4t * w34, ev[3]);
        cdd45 = Q4t * cdd34 + cross(wd45, av[3]) + cross(w45, cross(w45, av[3])) - Q4t * cross(wd34, sv[2]) - Q4t * cross(w34, cross(w34, sv[2]));

        wd56 = Q5t * wd45 + a[4] * ev[4] + v[4] * cross(Q5t * w45, ev[4]);
        cdd56 = Q5t * cdd45 + cross(wd56, av[4]) + cross(w56, cross(w56, av[4])) - Q5t * cross(wd45, sv[3]) - Q5t * cross(w45, cross(w45, sv[3]));

        wd67 = Q6t * wd56 + a[5] * ev[5] + v[5] * cross(Q6t * w56, ev[5]);
        cdd67 = Q6t * cdd56 + cross(wd67, av[5]) + cross(w67, cross(w67, av[5])) - Q6t * cross(wd56, sv[4]) - Q6t * cross(w56, cross(w56, sv[4]));

        //-- dynamics
        f6 = Q6 * (m[5] * cdd67);
        n6 = Q6 * (I[5] * wd67 + cross(w67, I[5] * w67) + cross(av[5], Q6t * f6));

        f5 = Q5 * (m[4] * cdd56 + f6);
        n5 = Q5 * (I[4] * wd56 + cross(w56, I[4] * w56) + n6 + cross(av[4], Q5t * f5) - cross(sv[4], f6));

        f4 = Q4 * (m[3] * cdd45 + f5);
        n4 = Q4 * (I[3] * wd45 + cross(w45, I[3] * w45) + n5 + cross(av[3], Q4t * f4) - cross(sv[3], f5));

        f3 = Q3 * (m[2] * cdd34 + f4);
        n3 = Q3 * (I[2] * wd34 + cross(w34, I[2] * w34) + n4 + cross(av[2], Q3t * f3) - cross(sv[2], f4));

        f2 = Q2 * (m[1] * cdd23 + f3);
        n2 = Q2 * (I[1] * wd23 + cross(w23, I[1] * w23) + n3 + cross(av[1], Q2t * f2) - cross(sv[1], f3));

        f1 = Q1 * (m[0] * cdd12 + f2);
        n1 = Q1 * (I[0] * wd12 + cross(w12, I[0] * w12) + n2 + cross(av[0], Q1t * f1) - cross(sv[0], f2));

        //-- extract torques (last element of each moment vector)
        efforts(0, i) = n1.z;
        efforts(1, i) = n2.z;
        efforts(2, i) = n3.z;
        efforts(3, i) = n4.z;
        efforts(4, i) = n5.z;
        efforts(5, i) = n6.z;
    }
}

host_fn Array Gen3Lite::forward_kinematics(Array &joint_position) {
    // auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;

    Array s(8);
    Array c(8);
    blast::sincos(joint_position, s, c);

// #if BLAST_SIZEOF_REAL == 8 // double precision
//     real s[6];
//     real c[6];
//     // first 4
//     {
//         __m256d s_tmp;
//         __m256d c_tmp;
//         __m256d angle_v = _mm256_load_pd(p);
//         s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//         _mm256_storeu_pd(s, s_tmp);
//         _mm256_storeu_pd(c, c_tmp);
//     }
//     // 5 and 6
//     {
//         __m128d angle_v = _mm_load_pd(p + 4);
//         __m128d s_tmp;
//         __m128d c_tmp;
//         s_tmp = _mm_sincos_pd(&c_tmp, angle_v);
//         _mm_storeu_pd(s + 4, s_tmp);
//         _mm_storeu_pd(c + 4, c_tmp);
//     }
// #else // single precision
//     real s[8];
//     real c[8];
//     __m256 s_tmp;
//     __m256 c_tmp;
//     __m256 angle_v = _mm256_load_ps(p);
//     s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//     _mm256_storeu_ps(s, s_tmp);
//     _mm256_storeu_ps(c, c_tmp);
// #endif

    // note: these are stored column-wise
    Q1 = {c[0], s[0], 0, -s[0], c[0], 0, 0, 0, 1};
    Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
    Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    Q5 = {0, s[4], -c[4], 0, c[4], s[4], 1, 0, 0};
    Q6 = {0, s[5], c[5], 0, c[5], -s[5], -1, 0, 0};

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_tmp += (Q_tmp *= Q3) * dv[2];
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_tmp += (Q_tmp *= Q6) * dv[5];

    Array pose(6);
    pose[0] = p_tmp.x;
    pose[1] = p_tmp.y;
    pose[2] = p_tmp.z;
    // todo: Q_tmp(0, 0) and Q_tmp(2, 2) must not be zero!!
    // pose[3] = atan2(Q_tmp(1, 0), Q_tmp(0, 0));
    // pose[4] = atan2(-Q_tmp(2, 0), sqrt(Q_tmp(2, 1) * Q_tmp(2, 1) + Q_tmp(2, 2) * Q_tmp(2, 2)));
    // pose[5] = atan2(Q_tmp(2, 1), Q_tmp(2, 2));
    return pose;
}

host_fn Matrix Gen3Lite::jacobian(const Array& joint_position) {

    // auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;

    Array s(8);
    Array c(8);
    blast::sincos(joint_position, s, c);

// #if BLAST_USE_DOUBLES
//     real s[8];
//     real c[8];
//     __m256d s_tmp;
//     __m256d c_tmp;
//     for (u32 i = 0; i < 8; i += 4) {
//         __m256d angle_v = _mm256_load_pd(p + i);
//         s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//         _mm256_storeu_pd(s+i, s_tmp);
//         _mm256_storeu_pd(c+i, c_tmp);
//     }
// #else
//     real s[8];
//     real c[8];
//     __m256 s_tmp;
//     __m256 c_tmp;
//     __m256 angle_v = _mm256_load_ps(p);
//     s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//     _mm256_storeu_ps(s, s_tmp);
//     _mm256_storeu_ps(c, c_tmp);
// #endif

    // note: these are stored column-wise
    Q1 = {c[0], s[0], 0, -s[0], c[0], 0, 0, 0, 1};
    Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
    Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    Q5 = {0, s[4], -c[4], 0, c[4], s[4], 1, 0, 0};
    Q6 = {0, s[5], c[5], 0, c[5], -s[5], -1, 0, 0};

    // unit vectors in 1st reference
    Vec3 e_1[6];
    auto Q_tmp = Q1;
    e_1[0] = Q_tmp * ev[0];
    e_1[1] = (Q_tmp *= Q2) * ev[1];
    e_1[2] = (Q_tmp *= Q3) * ev[2];
    e_1[3] = (Q_tmp *= Q4) * ev[3];
    e_1[4] = (Q_tmp *= Q5) * ev[4];
    e_1[5] = (Q_tmp *= Q6) * ev[5];

    Vec3 r[7];
    r[5] = dv[5];
    r[4] = dv[4] + Q6 * r[6];
    r[3] = dv[3] + Q5 * r[5];
    r[2] = dv[2] + Q4 * r[4];
    r[1] = dv[1] + Q3 * r[3];
    r[0] = dv[0] + Q2 * r[2];

    Q_tmp = Q1;
    r[0] = (Q_tmp) * r[1];
    r[1] = (Q_tmp *= Q2) * r[2];
    r[2] = (Q_tmp *= Q3) * r[3];
    r[3] = (Q_tmp *= Q4) * r[4];
    r[4] = (Q_tmp *= Q5) * r[5];
    r[5] = (Q_tmp *= Q6) * r[6];

    auto cr0 = cross(e_1[0], r[0]);
    auto cr1 = cross(e_1[1], r[1]);
    auto cr2 = cross(e_1[2], r[2]);
    auto cr3 = cross(e_1[3], r[3]);
    auto cr4 = cross(e_1[4], r[4]);
    auto cr5 = cross(e_1[5], r[5]);

    // jacobian matrix
    Matrix J(6, 6);
    J(0, 0) = e_1[0].x;
    J(1, 0) = e_1[0].y;
    J(2, 0) = e_1[0].z;
    J(0, 1) = e_1[1].x;
    J(1, 1) = e_1[1].y;
    J(2, 1) = e_1[1].z;
    J(0, 2) = e_1[2].x;
    J(1, 2) = e_1[2].y;
    J(2, 2) = e_1[2].z;
    J(0, 3) = e_1[3].x;
    J(1, 3) = e_1[3].y;
    J(2, 3) = e_1[3].z;
    J(0, 4) = e_1[4].x;
    J(1, 4) = e_1[4].y;
    J(2, 4) = e_1[4].z;
    J(0, 5) = e_1[5].x;
    J(1, 5) = e_1[5].y;
    J(2, 5) = e_1[5].z;

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

//------ Kinova Gen3 7DOF manipulator functions ---------------------------------------

host_fn Gen3_7DOF::Gen3_7DOF() : Manipulator(7) {
    // position of the first joint with respect to the table in the center of the base
    p_base = {0, 0, 0.1564f};

    // mass
    m[0] = 1.377f;
    m[1] = 1.1636f;
    m[2] = 1.1636f;
    m[3] = 0.93f;
    m[4] = 0.678f;
    m[5] = 0.678f;
    m[6] = 0.364f + 0.921f; // todo: modify with/without gripper ?

    // inertial tensor
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
    Vec3 av_tool(0, 0, -0.06f - 0.0615f); // todo: modify with/without gripper ?
    av[6] = (0.364f * av[6] + 0.921f * av_tool) * (1 / (0.364f + 0.921f));

    // vector to next joint
    dv[0] = {0.0, 0.0054, -0.1284};
    dv[1] = {0.0, -0.2104, -0.0064};
    dv[2] = {0.0, -0.0064, -0.2104};
    dv[3] = {0.0, -0.2084, -0.0064};
    dv[4] = {0.0, 0.0, -0.1059};
    dv[5] = {0.0, -0.1059, 0.0};
    dv[6] = {0.0, 0.0, -0.0615 - 0.164};

    // todo: add option to know if tool is closed or opened (difference of 0.0135 in z)

    // center of mass (from next joint)
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];
    sv[6] = dv[6] - av[6];

    // unit joint direction
    ev[0] = {0, 0, 1};
    ev[1] = {0, 0, 1};
    ev[2] = {0, 0, 1};
    ev[3] = {0, 0, 1};
    ev[4] = {0, 0, 1};
    ev[5] = {0, 0, 1};
    ev[6] = {0, 0, 1};

    // kinematic and dynamic constraints
    pmax = {INF_REAL, INF_REAL, INF_REAL, 2.58f, INF_REAL, 2.1f, INF_REAL}; // rad
    // pmin = -pmax;
    vmax = {1.745f, 1.745f, 1.745f, 1.745f, 2.443f, 2.443f, 2.443f}; // rad/s
    // vmin = -vmax;
    tau_max = {52, 52, 52, 52, 17, 17, 17}; // Nm
    // tau_min = -tau_max;
}

host_fn void Gen3_7DOF::set_payload(const real mass) {
    // Set to default
    m[6] = 0.364f + 0.921f;
    av[6] = {-0.000093f, 0.000132f, -0.022905f};
    Vec3 av_tool(0, 0, -0.06f - 0.0615f); // todo: add as global in manip Gen3_7DOF?
    av[6] = (0.364f * av[6] + 0.921f * av_tool) * (1 / (0.364f + 0.921f));
    sv[6] = dv[6] - av[6];
    I[6] = {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, 0.000002f, 0.000001f, 0.000002f, 0.000240f};

    // Modify payload
    Vec3 av_payload = {0.0f, 0.0f, -0.115f}; // todo: modify to real center of mass length (design prototype) todo: add to manip ?
    auto av_old = av[6];
    auto m_old = m[6];

    auto m_new = m_old + mass;
    auto av_new = (m_old*av_old + mass*av_payload) * (1/m_new);
    auto sv_new = dv[6] - av_new;
    auto delta_av = av_new - av_old; // shift in center of mass
    auto av_to_mass = av_payload - av_new; // vector from payload to new center of mass

    av[6] = av_new;
    sv[6] = sv_new;
    m[6] = m_new;
    I[6](0, 0) += m_old*delta_av.x*delta_av.x + mass*av_to_mass.x*av_to_mass.x; // todo: Validate
    I[6](1, 1) += m_old*delta_av.y*delta_av.y + mass*av_to_mass.y*av_to_mass.y;
    I[6](2, 2) += m_old*delta_av.z*delta_av.z + mass*av_to_mass.z*av_to_mass.z;
}

host_fn void Gen3_7DOF::set_payload_without_gripper(const real mass) {
    // Set to default
    m[6] = 0.364f;
    av[6] = {-0.000093f, 0.000132f, -0.022905f};
    sv[6] = dv[6] - av[6];
    I[6] = {0.000214f, 0.000000f, 0.000001f, 0.000000f, 0.000223f, 0.000002f, 0.000001f, 0.000002f, 0.000240f};

    // Modify payload
    Vec3 av_payload = {0.0f, 0.0f, -0.115f}; // todo: modify to real center of mass length (design prototype) todo: add to manip ?
    auto av_old = av[6];
    auto m_old = m[6];

    auto m_new = m_old + mass;
    auto av_new = (m_old*av_old + mass*av_payload) * (1/m_new);
    auto sv_new = dv[6] - av_new;
    auto delta_av = av_new - av_old; // shift in center of mass
    auto av_to_mass = av_payload - av_new; // vector from payload to new center of mass

    av[6] = av_new;
    sv[6] = sv_new;
    m[6] = m_new;
    I[6](0, 0) += m_old*delta_av.x*delta_av.x + mass*av_to_mass.x*av_to_mass.x; // todo: Validate
    I[6](1, 1) += m_old*delta_av.y*delta_av.y + mass*av_to_mass.y*av_to_mass.y;
    I[6](2, 2) += m_old*delta_av.z*delta_av.z + mass*av_to_mass.z*av_to_mass.z;
}

host_fn void Gen3_7DOF::dynamics(const Trajectory &traj, Matrix &efforts) {
    dynamics(traj.pos, traj.vel, traj.acc, efforts);
}

host_fn void Gen3_7DOF::dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc, Matrix &efforts) {

    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t, Q7t;
    Vec3 w1, w2, w3, w4, w5, w6, w7;
    Vec3 wd1, wd2, wd3, wd4, wd5, wd6, wd7;
    Vec3 cdd0 = {0, 0, 9.81f};
    Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6, cdd7;
    Vec3 f1, f2, f3, f4, f5, f6, f7;
    Vec3 n1, n2, n3, n4, n5, n6, n7;

    const u32 joints = pos.rows;
    const u32 points = pos.cols;

    // loop all points
    for (u32 i = 0; i < points; i++) {
        auto v = &vel.data[i * joints];
        auto a = &acc.data[i * joints];

        // SIMD compute sines and cosines note: approx 10% faster
        auto p = pos.col(i);
        Array s(7);
        Array c(7);
        blast::sincos(p, s, c);
//         auto p = &pos.data[i * joints];

// #if BLAST_SIZEOF_REAL == 8
//         real s[8];
//         real c[8];
//         __m256d s_tmp;
//         __m256d c_tmp;
//         for (u32 j = 0; j < 8; j += 4) {
//             __m256d angle_v = _mm256_load_pd(p + j);
//             s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//             _mm256_storeu_pd(s + j, s_tmp);
//             _mm256_storeu_pd(c + j, c_tmp);
//         }
// #else
//         real s[8];
//         real c[8];
//         __m256 s_tmp;
//         __m256 c_tmp;
//         __m256 angle_v = _mm256_load_ps(p);
//         s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//         _mm256_storeu_ps(s, s_tmp);
//         _mm256_storeu_ps(c, c_tmp);
// #endif

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
        n6 = I[5] * wd6 + cross(w6, I[5] * w6) + Q7 * n7 + cross(av[5], f6) + cross(sv[5], (Q7 * f7));

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
        efforts(0, i) = n1.z;
        efforts(1, i) = n2.z;
        efforts(2, i) = n3.z;
        efforts(3, i) = n4.z;
        efforts(4, i) = n5.z;
        efforts(5, i) = n6.z;
        efforts(6, i) = n7.z;
    }
}

host_fn void Gen3_7DOF::dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc) {
    const auto points = pos.cols;
    const auto joints = pos.rows;
    if (_efforts.cols != points || _efforts.rows != joints)
        _efforts.resize(joints, points);

    dynamics(pos, vel, acc, _efforts);
}

host_fn Array Gen3_7DOF::forward_kinematics(const Array &joint_position) {

    //  - note: manual SIMD (10% better performance than using sincos function on arrays like commented below)
    // real s[8];
    // real c[8];
    // auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;

    Array s(8);
    Array c(8);
    blast::sincos(joint_position, s, c);

// #if BLAST_SIZEOF_REAL == 8
//     __m256d s_tmp;
//     __m256d c_tmp;
//     for (u32 i = 0; i < 8; i += 4) {
//         __m256d angle_v = _mm256_load_pd(p + i);
//         s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//         _mm256_storeu_pd(s + i, s_tmp);
//         _mm256_storeu_pd(c + i, c_tmp);
//     }
// #else
//     __m256 s_tmp;
//     __m256 c_tmp;
//     __m256 angle_v = _mm256_load_ps(p);
//     s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//     _mm256_storeu_ps(s, s_tmp);
//     _mm256_storeu_ps(c, c_tmp);
// #endif

    // todo: cleanup
    // create aliases instead of actual arrays to allocate the memory on the stack (1.7x performance)
    // real s_data[8]; // nearest factor of 4
    // real c_data[8]; // nearest factor of 4
    // Array s(s_data, 8);
    // Array c(c_data, 8);
    // Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    // sincos(joint_position, s, c);

    // note: these are stored column-wise
    Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
    Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
    Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_tmp += (Q_tmp *= Q3) * dv[2];
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_tmp += (Q_tmp *= Q6) * dv[5];
    p_tmp += (Q_tmp *= Q7) * dv[6];

    Array pose(6);
    pose[0] = p_tmp.x;
    pose[1] = p_tmp.y;
    pose[2] = p_tmp.z;
    // todo: Q_tmp(0, 0) and Q_tmp(2, 2) must not be zero!!
    // pose[3] = atan2(Q_tmp(1, 0), Q_tmp(0, 0));
    // pose[4] = atan2(-Q_tmp(2, 0), sqrt(Q_tmp(2, 1)*Q_tmp(2, 1) + Q_tmp(2, 2)*Q_tmp(2, 2)));
    // pose[5] = atan2(Q_tmp(2, 1), Q_tmp(2, 2));
    return pose;
}

host_fn Matrix Gen3_7DOF::forward_kinematics(const Matrix &joint_positions) {
    auto p = joint_positions.data;
    Matrix pose(12, joint_positions.cols);

    for (u32 point = 0; point < joint_positions.cols; point++) {

        Array s(8);
        Array c(8);
        blast::sincos(joint_positions.col(point), s, c);

// #if BLAST_USE_DOUBLES
//         real s[8];
//         real c[8];
//         __m256d s_tmp;
//         __m256d c_tmp;
//         for (u32 i = 0; i < 8; i += 4) {
//             __m256d angle_v = _mm256_load_pd(p + i);
//             s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//             _mm256_storeu_pd(s + i, s_tmp);
//             _mm256_storeu_pd(c + i, c_tmp);
//         }
// #else
//         real s[8];
//         real c[8];
//         __m256 s_tmp;
//         __m256 c_tmp;
//         __m256 angle_v = _mm256_load_ps(p);
//         s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//         _mm256_storeu_ps(s, s_tmp);
//         _mm256_storeu_ps(c, c_tmp);
// #endif

#if 1
        Mat3 Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Mat3 Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
        Mat3 Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
        Mat3 Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
        Mat3 Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Mat3 Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
        Mat3 Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};
        auto p_tmp = p_base;
        auto Q_tmp = Q1;
        p_tmp += Q_tmp * dv[0];
        p_tmp += (Q_tmp *= Q2) * dv[1];
        p_tmp += (Q_tmp *= Q3) * dv[2];
        p_tmp += (Q_tmp *= Q4) * dv[3];
        p_tmp += (Q_tmp *= Q5) * dv[4];
        p_tmp += (Q_tmp *= Q6) * dv[5];
        p_tmp += (Q_tmp *= Q7) * dv[6];
        pose(0, point) = p_tmp.x;
        pose(1, point) = p_tmp.y;
        pose(2, point) = p_tmp.z;
        pose(3, point) = Q_tmp[0];
        pose(4, point) = Q_tmp[1];
        pose(5, point) = Q_tmp[2];
        pose(6, point) = Q_tmp[3];
        pose(7, point) = Q_tmp[4];
        pose(8, point) = Q_tmp[5];
        pose(9, point) = Q_tmp[6];
        pose(10, point) = Q_tmp[7];
        pose(11, point) = Q_tmp[8];
#else
        // mat4 implementation (not as fast as rotation matrix + vector)
        const auto p1 = p_base + dv[0];
        const Mat4 T1 = {c[0], -s[0], 0, 0, -s[0], -c[0], 0, 0, 0, 0, -1, 0, p1.x, p1.y, p1.z, 1};
        const Mat4 T2 = {c[1], 0, s[1], 0, -s[1], 0, c[1], 0, 0, -1, 0, 0, dv[1].x, dv[1].y, dv[1].z, 1};
        const Mat4 T3 = {c[2], 0, -s[2], 0, -s[2], 0, -c[2], 0, 0, 1, 0, 0, dv[2].x, dv[2].y, dv[2].z, 1};
        const Mat4 T4 = {c[3], 0, s[3], 0, -s[3], 0, c[3], 0, 0, -1, 0, 0, dv[3].x, dv[3].y, dv[3].z, 1};
        const Mat4 T5 = {c[4], 0, -s[4], 0, -s[4], 0, -c[4], 0, 0, 1, 0, 0, dv[4].x, dv[4].y, dv[4].z, 1};
        const Mat4 T6 = {c[5], 0, s[5], 0, -s[5], 0, c[5], 0, 0, -1, 0, 0, dv[5].x, dv[5].y, dv[5].z, 1};
        const Mat4 T7 = {c[6], 0, -s[6], 0, -s[6], 0, -c[6], 0, 0, 1, 0, 0, dv[6].x, dv[6].y, dv[6].z, 1};
        const auto T = T1 * T2 * T3 * T4 * T5 * T6 * T7;
        pose(0, point) = T(0, 3);
        pose(1, point) = T(1, 3);
        pose(2, point) = T(2, 3);
        pose(3, point) = T[0];
        pose(4, point) = T[1];
        pose(5, point) = T[2];
        pose(6, point) = T[4];
        pose(7, point) = T[5];
        pose(8, point) = T[6];
        pose(9, point) = T[8];
        pose(10, point) = T[9];
        pose(11, point) = T[10];
#endif
        p += joints;
    }
    return pose;
}

host_fn Array Gen3_7DOF::inverse_kinematics(const Array& pose, const Array& initial_joint_position) {

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

host_fn Matrix Gen3_7DOF::jacobian(const Array &joint_position) {

    // auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;

    Array s(8);
    Array c(8);
    blast::sincos(joint_position, s, c);

// #if BLAST_USE_DOUBLES
//     real s[8];
//     real c[8];
//     __m256d s_tmp;
//     __m256d c_tmp;
//     for (u32 i = 0; i < 8; i += 4) {
//         __m256d angle_v = _mm256_load_pd(p + i);
//         s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//         _mm256_storeu_pd(s + i, s_tmp);
//         _mm256_storeu_pd(c + i, c_tmp);
//     }
// #else
//     real s[8];
//     real c[8];
//     __m256 s_tmp;
//     __m256 c_tmp;
//     __m256 angle_v = _mm256_load_ps(p);
//     s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//     _mm256_storeu_ps(s, s_tmp);
//     _mm256_storeu_ps(c, c_tmp);
// #endif

    // note: these are stored column-wise
    Q1 = {c[0], -s[0],  0, -s[0], -c[0],   0,   0,  0, -1};
    Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
    Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
    Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};

    // unit vectors in 1st reference
    Vec3 e_1[7];
    auto Q_tmp = Q1;
    e_1[0] = Q_tmp * ev[0];
    e_1[1] = (Q_tmp *= Q2) * ev[1];
    e_1[2] = (Q_tmp *= Q3) * ev[2];
    e_1[3] = (Q_tmp *= Q4) * ev[3];
    e_1[4] = (Q_tmp *= Q5) * ev[4];
    e_1[5] = (Q_tmp *= Q6) * ev[5];
    e_1[6] = (Q_tmp *= Q7) * ev[6];

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

    auto cr0 = cross(e_1[0], r[0]);
    auto cr1 = cross(e_1[1], r[1]);
    auto cr2 = cross(e_1[2], r[2]);
    auto cr3 = cross(e_1[3], r[3]);
    auto cr4 = cross(e_1[4], r[4]);
    auto cr5 = cross(e_1[5], r[5]);
    auto cr6 = cross(e_1[6], r[6]);

    // jacobian matrix
    Matrix J(6, 7);
    J(0, 0) = e_1[0].x;
    J(1, 0) = e_1[0].y;
    J(2, 0) = e_1[0].z;
    J(0, 1) = e_1[1].x;
    J(1, 1) = e_1[1].y;
    J(2, 1) = e_1[1].z;
    J(0, 2) = e_1[2].x;
    J(1, 2) = e_1[2].y;
    J(2, 2) = e_1[2].z;
    J(0, 3) = e_1[3].x;
    J(1, 3) = e_1[3].y;
    J(2, 3) = e_1[3].z;
    J(0, 4) = e_1[4].x;
    J(1, 4) = e_1[4].y;
    J(2, 4) = e_1[4].z;
    J(0, 5) = e_1[5].x;
    J(1, 5) = e_1[5].y;
    J(2, 5) = e_1[5].z;
    J(0, 6) = e_1[6].x;
    J(1, 6) = e_1[6].y;
    J(2, 6) = e_1[6].z;

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

host_fn Array Gen3_7DOF::collision_check(const Array &joint_position) {
    // real s[8];
    // real c[8];
    // auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;

    Array s(7);
    Array c(7);
    blast::sincos(joint_position, s, c);

    // // todo: make it work with floats!
    // __m256d s_tmp;
    // __m256d c_tmp;
    // for (u32 i = 0; i < 8; i += 4) {
    //     __m256d angle_v = _mm256_load_pd(p + i);
    //     s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
    //     _mm256_storeu_pd(s + i, s_tmp);
    //     _mm256_storeu_pd(c + i, c_tmp);
    // }

    // note: these are stored column-wise
    Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
    Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
    Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
    Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
    Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};

    Vec3 p_orig(0, 0, 0);
    Vec3 p_j2;
    Vec3 p_j3;
    Vec3 p_j4;
    Vec3 p_j6;
    Vec3 p_ee;

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_j2 = p_tmp;
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_j3 = p_tmp;
    p_tmp += (Q_tmp *= Q3) * dv[2];
    p_j4 = p_tmp;
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_j6 = p_tmp;
    p_tmp += (Q_tmp *= Q6) * dv[5];
    p_tmp += (Q_tmp *= Q7) * dv[6];
    p_ee = p_tmp;

    const real r1sqr = 0.0875 * 0.0875;
    const real r2sqr = 0.0875 * 0.0875;

    // Self collisions sqr
    real dist1sqr = two_segment_distance_sqr(p_orig, p_j2, p_j6, p_ee) - r1sqr;
    real dist2J2sqr = two_segment_distance_sqr(p_j2, p_j2, p_j6, p_ee) - r2sqr; // distance sqr from J2 to J6 EE line (sphere)
    real dist2Msqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r1sqr;  // distance sqr from J2 J3 line to J6 EE line (capsule)
    real dist2sqr = dist2J2sqr <= dist2Msqr ? dist2J2sqr : dist2Msqr;

    // todo: remove collision detection from here

    // Collision with table sqr
    const real r4table = 0.05; // todo: validate dimensions
    const real r6table = 0.04; // todo: validate dimensions

    const Vec3 p_table(0, 0, -0.0025); // todo: Correct coords (z or y) ??
    real distTJ4 = p_j4.z - p_table.z - r4table;
    real distTJ6 = p_j6.z - p_table.z - r6table;
    real distTEE = p_ee.z - p_table.z - r6table;

    // Array of distance min sqr and distance from table sqr
    Array distMin(5);
    distMin = {dist1sqr, dist2sqr, distTJ4, distTJ6, distTEE};
    // todo: Add collision_check function for more obstacles; Array distMin(5) hard coded size

    return distMin;
}

host_fn Array Gen3_7DOF::constraints(const Array &pos, const Array &vel, const Array &acc) {
    Matrix p(pos);
    Matrix v(vel);
    Matrix a(acc);

    dynamics(p, v, a);
    // 5 collision results
    // 2 position constraints
    // 7 velocity constraints
    // 7 torque constraints
    Array result(5 + 2 + 7 * 2);

    // distance to collision >= 0
    auto tmp_coll = collision_check(pos);
    result[0] = -tmp_coll[0]; // dist1sqr
    result[1] = -tmp_coll[1]; // dist2sqr
    result[2] = -tmp_coll[2]; // distTJ4
    result[3] = -tmp_coll[3]; // distTJ6
    result[4] = -tmp_coll[4]; // distTEE

    // position
    result[5] = (abs(pos[3]) - pmax[3]) / pmax[3];
    result[6] = (abs(pos[5]) - pmax[5]) / pmax[5];

    // velocity
    for (int j = 0; j < 7; j++)
        result[j+7] = (abs(vel[j]) - vmax[j]) / vmax[j];

    auto f = _efforts.col(0); Assert(f.is_alias);
    for (int j = 0; j < 7; j++)
        result[j+14] = (abs(f[j]) - tau_max[j]) / tau_max[j];

    return result;
}

host_fn Array Gen3_7DOF::constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc) {
    const auto points = pos.cols;
    Array result(ncon(points));
    constraints(pos, vel, acc, result.data);
    return result;
}

host_fn Array Gen3_7DOF::constraints(const Trajectory &traj) {
    Array result(ncon(traj.t.size));
    constraints(traj.pos, traj.vel, traj.acc, result.data);
    return result;
}

host_fn void Gen3_7DOF::constraints(const Trajectory &traj, real *dst) {
    constraints(traj.pos, traj.vel, traj.acc, dst);
}

host_fn void Gen3_7DOF::constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc, real *dst) {
    const auto points = pos.cols;

    dynamics(pos, vel, acc);

    for (u32 i = 0; i < points; i++) {
        // collision
        auto p = pos.col(i); Assert(p.is_alias);
        auto tmp_coll = collision_check(p);
        dst[0] = -tmp_coll[0]; // dist1sqr
        dst[1] = -tmp_coll[1]; // dist2sqr
        dst[2] = -tmp_coll[2]; // distTJ4
        dst[3] = -tmp_coll[3]; // distTJ6
        dst[4] = -tmp_coll[4]; // distTEE
        dst += 5;

        // position
        dst[0] = (abs(pos(3, i)) - pmax[3]) / pmax[3];
        dst[1] = (abs(pos(5, i)) - pmax[5]) / pmax[5];
        dst += 2;

        // velocity
        for (int j = 0; j < 7; j++)
            dst[j] = (abs(vel(j, i)) - vmax[j]) / vmax[j];
        dst += 7;

        auto f = _efforts.col(i); Assert(f.is_alias);
        for (int j = 0; j < 7; j++)
            dst[j] = (abs(f[j]) - tau_max[j]) / tau_max[j];
        dst += 7;
    }
}

host_fn bool Gen3_7DOF::validate_task(const Matrix &task) {
    auto ci = constraints(task.col(0), task.col(1), task.col(2));
    auto cf = constraints(task.col(3), task.col(4), task.col(5));
    return array_max(ci) > 0 && array_max(cf) > 0 ? false: true;
}

host_fn Matrix Gen3_7DOF::robot_capsules(const Matrix &joint_positions, const int n_skip) {
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    int n_col = (joint_positions.cols - joint_positions.cols%n_skip)/n_skip + 1;
    Matrix result_capsules(12, n_col);

    Array s(7);
    Array c(7);
    for (u32 col = 0; col < n_col; col++) {
        int idx = col*n_skip < joint_positions.cols ? col*n_skip : joint_positions.cols - 1;
        auto joint_position_tmp = joint_positions.col(idx); //todo: double check??
        blast::sincos(joint_position_tmp, s, c);

        // note: these are stored column-wise
        Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Q2 = {c[1], 0, s[1], -s[1], 0, c[1], 0, -1, 0};
        Q3 = {c[2], 0, -s[2], -s[2], 0, -c[2], 0, 1, 0};
        Q4 = {c[3], 0, s[3], -s[3], 0, c[3], 0, -1, 0};
        Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Q6 = {c[5], 0, s[5], -s[5], 0, c[5], 0, -1, 0};
        Q7 = {c[6], 0, -s[6], -s[6], 0, -c[6], 0, 1, 0};

        Vec3 p_orig(0, 0, 0);
        Vec3 p_j2;
        Vec3 p_j4;
        Vec3 p_j6;
        Vec3 p_ee;

        auto p_tmp = p_base;
        auto Q_tmp = Q1;
        p_tmp += Q_tmp * dv[0];
        p_j2 = p_tmp;
        p_tmp += (Q_tmp *= Q2) * dv[1];
        p_tmp += (Q_tmp *= Q3) * dv[2];
        p_j4 = p_tmp;
        p_tmp += (Q_tmp *= Q4) * dv[3];
        p_tmp += (Q_tmp *= Q5) * dv[4];
        p_j6 = p_tmp;
        p_tmp += (Q_tmp *= Q6) * dv[5];
        p_tmp += (Q_tmp *= Q7) * dv[6];
        p_ee = p_tmp;

        result_capsules(0, col) = p_j2.x;
        result_capsules(1, col) = p_j2.y;
        result_capsules(2, col) = p_j2.z;

        result_capsules(3, col) = p_j4.x;
        result_capsules(4, col) = p_j4.y;
        result_capsules(5, col) = p_j4.z;

        result_capsules(6, col) = p_j6.x;
        result_capsules(7, col) = p_j6.y;
        result_capsules(8, col) = p_j6.z;

        result_capsules(9, col) = p_ee.x;
        result_capsules(10, col) = p_ee.y;
        result_capsules(11, col) = p_ee.z;
    }

    return result_capsules;
}

//------ Kinova Gen3 7DOF manipulator functions (using Eigen) ---------------------------------------

host_fn Gen3Eigen::Gen3Eigen() : ManipulatorEigen(7) {
    // position of the first joint with respect to the table in the center of the base
    p_base = {0, 0, 0.1564f};

    set_payload(0, true);

    // todo: add option to know if tool is closed or opened (difference of 0.0135 in z)?

    // unit joint direction
    ev[0] = {0, 0, 1};
    ev[1] = {0, 0, 1};
    ev[2] = {0, 0, 1};
    ev[3] = {0, 0, 1};
    ev[4] = {0, 0, 1};
    ev[5] = {0, 0, 1};
    ev[6] = {0, 0, 1};

    // kinematic and dynamic constraints
    pmax << INF_REAL, INF_REAL, INF_REAL, 2.58, INF_REAL, 2.1, INF_REAL; // rad
    // pmin = -pmax;
    vmax << 1.745, 1.745, 1.745, 1.745, 2.443, 2.443, 2.443; // rad/s
    // vmin = -vmax;
    tau_max << 52, 52, 52, 52, 17, 17, 17; // Nm
    // tau_min = -tau_max;
}

host_fn void Gen3Eigen::set_payload(double mass, bool add_gripper) {
    // mass
    m[0] = 1.377;
    m[1] = 1.1636;
    m[2] = 1.1636;
    m[3] = 0.93;
    m[4] = 0.678;
    m[5] = 0.678;
    m[6] = 0.364;

    // inertial tensor
    I[0] << 0.004570, 0.000001, 0.000002,
    0.000001, 0.004831, 0.000448,
    0.000002, 0.000448, 0.001409;
    I[1] << 0.011088, 0.000005, 0.000000,
    0.000005, 0.001072, 0.000691,
    0.000000, 0.000691, 0.011255;
    I[2] << 0.010932, 0.000000, 0.000007,
    0.000000, 0.011127, 0.000606,
    0.000007, 0.000606, 0.001043;
    I[3] << 0.008147, 0.000001, 0.000000,
    0.000001, 0.000631, 0.000500,
    0.000000, 0.000500, 0.008316;
    I[4] << 0.001596, 0.000000, 0.000000,
    0.000000, 0.001607, 0.000256,
    0.000000, 0.000256, 0.000399;
    I[5] << 0.001641, 0.000000, 0.000000,
    0.000000, 0.000410, 0.000278,
    0.000000, 0.000278, 0.001641;
    I[6] << 0.000214, 0.000000, 0.000001,
    0.000000, 0.000223, 0.000002,
    0.000001, 0.000002, 0.000240;

    // center of mass
    av[0] = {-0.000023, -0.010364, -0.073360};
    av[1] = {-0.000044, -0.099580, -0.013278};
    av[2] = {-0.000044, -0.006641, -0.117892};
    av[3] = {-0.000018, -0.075478, -0.015006};
    av[4] = {0.000001, -0.009432, -0.063883};
    av[5] = {0.000001, -0.045483, -0.009650};
    av[6] = {-0.000093, 0.000132, -0.022905};

    // vector to next joint
    dv[0] = {0.0, 0.0054, -0.1284};
    dv[1] = {0.0, -0.2104, -0.0064};
    dv[2] = {0.0, -0.0064, -0.2104};
    dv[3] = {0.0, -0.2084, -0.0064};
    dv[4] = {0.0, 0.0, -0.1059};
    dv[5] = {0.0, -0.1059, 0.0};
    dv[6] = {0.0, 0.0, -0.0615};

    // Modify payload
    if (add_gripper) {
        m[6] += 0.921;
        dv[6][2] -= 0.164;

        Vector3d av_tool(0, 0, -0.06 - 0.0615);
        av[6] = (0.364 * av[6] + 0.921 * av_tool) * (1 / (0.364 + 0.921));

        Vector3d av_payload = {0.0f, 0.0f, -0.115f}; // todo: modify to real center of mass length (design prototype) todo: add to manip ?

        Vector3d av_old = av[6];
        double m_old = m[6];

        double m_new = m_old + mass;
        Vector3d av_new = (m_old*av_old + mass*av_payload) * (1/m_new);
        Vector3d sv_new = dv[6] - av_new;
        Vector3d delta_av = av_new - av_old;        // shift in center of mass
        Vector3d av_to_mass = av_payload - av_new;  // vector from payload to new center of mass

        av[6] = av_new;
        sv[6] = sv_new;
        m[6] = m_new;
        I[6](0, 0) += m_old*delta_av[0]*delta_av[0] + mass*av_to_mass[0]*av_to_mass[0];
        I[6](1, 1) += m_old*delta_av[1]*delta_av[1] + mass*av_to_mass[1]*av_to_mass[1];
        I[6](2, 2) += m_old*delta_av[2]*delta_av[2] + mass*av_to_mass[2]*av_to_mass[2];
    }
    else {
        // Set to default
        m[6] = 0.364f;
        av[6] = {-0.000093f, 0.000132f, -0.022905f};
        sv[6] = dv[6] - av[6];
        // Modify payload
        Vector3d av_payload = {0.0f, 0.0f, -0.115f}; // todo: modify to real center of mass length (design prototype) todo: add to manip ?
        auto av_old = av[6];
        auto m_old = m[6];
        m[6] = m_old + mass;
        av[6] = (m_old*av_old + mass*av_payload) * (1/m[6]);
        sv[6] = dv[6] - av[6];
        auto delta_av = av[6] - av_old; // shift in center of mass
        auto av_to_mass = av_payload - av[6]; // vector from payload to new center of mass

        I[6](0, 0) += m_old*delta_av[0]*delta_av[0] + mass*av_to_mass[0]*av_to_mass[0]; // todo: Validate
        I[6](1, 1) += m_old*delta_av[1]*delta_av[1] + mass*av_to_mass[1]*av_to_mass[1];
        I[6](2, 2) += m_old*delta_av[2]*delta_av[2] + mass*av_to_mass[2]*av_to_mass[2];
    }

// center of mass (from next joint)
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];
    sv[6] = dv[6] - av[6];
}

host_fn void Gen3Eigen::dynamics(const TrajectoryEigen &traj, MatrixXd &efforts) {
    dynamics(traj.pos, traj.vel, traj.acc, efforts);
}

host_fn void Gen3Eigen::dynamics(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc, MatrixXd &efforts) {

    Matrix3d Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Matrix3d Q1t, Q2t, Q3t, Q4t, Q5t, Q6t, Q7t;
    Vector3d w1, w2, w3, w4, w5, w6, w7;
    Vector3d wd1, wd2, wd3, wd4, wd5, wd6, wd7;
    Vector3d cdd0 = {0, 0, 9.81f};
    Vector3d cdd1, cdd2, cdd3, cdd4, cdd5, cdd6, cdd7;
    Vector3d f1, f2, f3, f4, f5, f6, f7;
    Vector3d n1, n2, n3, n4, n5, n6, n7;

    const int joints = pos.rows();
    const int points = pos.cols();

    // loop all points
    for (int i = 0; i < points; i++) {
        ArrayXd p = pos.col(i);
        auto v = &vel.data()[i * joints]; // todo: use blocks?
        auto a = &acc.data()[i * joints]; // todo: use blocks?

        ArrayXd s = p.sin();
        ArrayXd c = p.cos();

        Q1t <<
            c[0], -s[0], 0,
            -s[0], -c[0], 0,
            0, 0, -1;
        Q2t <<
            c[1], 0,  s[1],
            -s[1], 0,  c[1],
            0, -1, 0;
        Q3t <<
            c[2], 0, -s[2],
            -s[2], 0, -c[2],
            0,  1, 0;
        Q4t <<
            c[3], 0,  s[3],
            -s[3], 0,  c[3],
            0, -1, 0;
        Q5t <<
            c[4], 0, -s[4],
            -s[4], 0, -c[4],
            0,  1, 0;
        Q6t <<
            c[5], 0,  s[5],
            -s[5], 0,  c[5],
            0, -1, 0;
        Q7t <<
            c[6], 0, -s[6],
            -s[6], 0, -c[6],
            0,  1, 0;

        Q1 = Q1t.transpose();
        Q2 = Q2t.transpose();
        Q3 = Q3t.transpose();
        Q4 = Q4t.transpose();
        Q5 = Q5t.transpose();
        Q6 = Q6t.transpose();
        Q7 = Q7t.transpose();

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
        cdd1 = Q1t * cdd0 + wd1.cross(av[0]) + w1.cross(w1.cross(av[0]));

        wd2 = Q2t * wd1 + a[1] * ev[1] + v[1] * (Q2t * w1).cross(ev[1]);
        cdd2 = Q2t * cdd1 + wd2.cross(av[1]) + w2.cross(w2.cross(av[1])) - Q2t * wd1.cross(sv[0]) - Q2t * w1.cross(w1.cross(sv[0]));

        wd3 = Q3t * wd2 + a[2] * ev[2] + v[2] * (Q3t * w2).cross(ev[2]);
        cdd3 = Q3t * cdd2 + wd3.cross(av[2]) + w3.cross(w3.cross(av[2])) - Q3t * wd2.cross(sv[1]) - Q3t * w2.cross(w2.cross(sv[1]));

        wd4 = Q4t * wd3 + a[3] * ev[3] + v[3] * (Q4t * w3).cross(ev[3]);
        cdd4 = Q4t * cdd3 + wd4.cross(av[3]) + w4.cross(w4.cross(av[3])) - Q4t * wd3.cross(sv[2]) - Q4t * w3.cross(w3.cross(sv[2]));

        wd5 = Q5t * wd4 + a[4] * ev[4] + v[4] * (Q5t * w4).cross(ev[4]);
        cdd5 = Q5t * cdd4 + wd5.cross(av[4]) + w5.cross(w5.cross(av[4])) - Q5t * wd4.cross(sv[3]) - Q5t * w4.cross(w4.cross(sv[3]));

        wd6 = Q6t * wd5 + a[5] * ev[5] + v[5] * (Q6t * w5).cross(ev[5]);
        cdd6 = Q6t * cdd5 + wd6.cross(av[5]) + w6.cross(w6.cross(av[5])) - Q6t * wd5.cross(sv[4]) - Q6t * w5.cross(w5.cross(sv[4]));

        wd7 = Q7t * wd6 + a[6] * ev[6] + v[6] * (Q7t * w6).cross(ev[6]);
        cdd7 = Q7t * cdd6 + wd7.cross(av[6]) + w7.cross(w7.cross(av[6])) - Q7t * wd6.cross(sv[5]) - Q7t * w6.cross(w6.cross(sv[5]));

        //-- dynamics
        f7 = m[6] * cdd7;
        n7 = I[6] * wd7 + w7.cross(I[6] * w7) + av[6].cross(f7);

        f6 = m[5] * cdd6 + Q7 * f7;
        n6 = I[5] * wd6 + w6.cross(I[5] * w6) + Q7 * n7 + av[5].cross(f6) + sv[5].cross((Q7 * f7));

        f5 = m[4] * cdd5 + Q6 * f6;
        n5 = I[4] * wd5 + w5.cross(I[4] * w5) + Q6 * n6 + av[4].cross(f5) + sv[4].cross((Q6 * f6));

        f4 = m[3] * cdd4 + Q5 * f5;
        n4 = I[3] * wd4 + w4.cross(I[3] * w4) + Q5 * n5 + av[3].cross(f4) + sv[3].cross((Q5 * f5));

        f3 = m[2] * cdd3 + Q4 * f4;
        n3 = I[2] * wd3 + w3.cross(I[2] * w3) + Q4 * n4 + av[2].cross(f3) + sv[2].cross((Q4 * f4));

        f2 = m[1] * cdd2 + Q3 * f3;
        n2 = I[1] * wd2 + w2.cross(I[1] * w2) + Q3 * n3 + av[1].cross(f2) + sv[1].cross((Q3 * f3));

        f1 = m[0] * cdd1 + Q2 * f2;
        n1 = I[0] * wd1 + w1.cross(I[0] * w1) + Q2 * n2 + av[0].cross(f1) + sv[0].cross((Q2 * f2));

        //-- extract torques (last element of each moment vector)
        efforts(0, i) = n1[2];
        efforts(1, i) = n2[2];
        efforts(2, i) = n3[2];
        efforts(3, i) = n4[2];
        efforts(4, i) = n5[2];
        efforts(5, i) = n6[2];
        efforts(6, i) = n7[2];
    }
}

host_fn void Gen3Eigen::dynamics(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc) {
    const auto points = pos.cols();
    const auto joints = pos.rows();
    if (_efforts.cols() != points || _efforts.rows() != joints)
        _efforts.resize(joints, points);

    dynamics(pos, vel, acc, _efforts);
}

host_fn VectorXd Gen3Eigen::forward_kinematics(const VectorXd &joint_position) {

    //  - note: manual SIMD (10% better performance than using sincos function on arrays like commented below)
    // real s[8];
    // real c[8];
    // auto p = joint_position.data;
    Matrix3d Q1, Q2, Q3, Q4, Q5, Q6, Q7;

    ArrayXd p(joint_position);
    ArrayXd s = p.sin();
    ArrayXd c = p.cos();

    // note: these are stored column-wise
    Q1 <<
       c[0], -s[0], 0,
       -s[0], -c[0], 0,
       0, 0, -1;
    Q2 <<
       c[1], 0,  s[1],
       -s[1], 0,  c[1],
       0, -1, 0;
    Q3 <<
       c[2], 0, -s[2],
       -s[2], 0, -c[2],
       0,  1, 0;
    Q4 <<
       c[3], 0,  s[3],
       -s[3], 0,  c[3],
       0, -1, 0;
    Q5 <<
       c[4], 0, -s[4],
       -s[4], 0, -c[4],
       0,  1, 0;
    Q6 <<
       c[5], 0,  s[5],
       -s[5], 0,  c[5],
       0, -1, 0;
    Q7 <<
       c[6], 0, -s[6],
       -s[6], 0, -c[6],
       0,  1, 0;
    Q1.transposeInPlace();
    Q2.transposeInPlace();
    Q3.transposeInPlace();
    Q4.transposeInPlace();
    Q5.transposeInPlace();
    Q6.transposeInPlace();
    Q7.transposeInPlace();

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_tmp += (Q_tmp *= Q3) * dv[2];
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_tmp += (Q_tmp *= Q6) * dv[5];
    p_tmp += (Q_tmp *= Q7) * dv[6];

    VectorXd pose(6);
    pose[0] = p_tmp[0];
    pose[1] = p_tmp[1];
    pose[2] = p_tmp[2];
    // todo: Q_tmp(0, 0) and Q_tmp(2, 2) must not be zero!!
    // pose[3] = atan2(Q_tmp(1, 0), Q_tmp(0, 0));
    // pose[4] = atan2(-Q_tmp(2, 0), sqrt(Q_tmp(2, 1)*Q_tmp(2, 1) + Q_tmp(2, 2)*Q_tmp(2, 2)));
    // pose[5] = atan2(Q_tmp(2, 1), Q_tmp(2, 2));
    return pose;
}

host_fn VectorXd Gen3Eigen::inverse_kinematics(const VectorXd& pose, const VectorXd& initial_joint_position) {
    const double tolerance = 0.001; // Tolerance for convergence
    const int max_iter = 100; // Maximum number of iterations

    VectorXd current_joint_angles = initial_joint_position;

    // // Iterate until convergence or maximum iterations reached
    // for (int iter = 0; iter < max_iter; ++iter) {
    //     // Calculate the current end effector position using forward kinematics
    //     VectorXd current_pose = forward_kinematics(current_joint_angles);
    //     VectorXd delta_pose = pose - current_pose;
    //     // Check if the end effector is close enough to the desired position
    //     if (is_close(pose, current_pose, tolerance))
    //         break;
    //     // Calculate the Jacobian matrix
    //     MatrixXd jacobian_matrix = jacobian(current_joint_angles);
    //     MatrixXd jacobian_pinv = pinv(jacobian_matrix);
    //     //current_joint_angles = current_joint_angles + jacobian_pinv * delta_pose;
    // }

    return current_joint_angles;
}

host_fn MatrixXd Gen3Eigen::jacobian(const VectorXd &joint_position) {

    // auto p = joint_position.data;
    Matrix3d Q1, Q2, Q3, Q4, Q5, Q6, Q7;

    ArrayXd p(joint_position);
    ArrayXd s = p.sin();
    ArrayXd c = p.cos();

    // note: these are stored column-wise
    Q1 <<
       c[0], -s[0], 0,
       -s[0], -c[0], 0,
       0, 0, -1;
    Q2 <<
       c[1], 0,  s[1],
       -s[1], 0,  c[1],
       0, -1, 0;
    Q3 <<
       c[2], 0, -s[2],
       -s[2], 0, -c[2],
       0,  1, 0;
    Q4 <<
       c[3], 0,  s[3],
       -s[3], 0,  c[3],
       0, -1, 0;
    Q5 <<
       c[4], 0, -s[4],
       -s[4], 0, -c[4],
       0,  1, 0;
    Q6 <<
       c[5], 0,  s[5],
       -s[5], 0,  c[5],
       0, -1, 0;
    Q7 <<
       c[6], 0, -s[6],
       -s[6], 0, -c[6],
       0,  1, 0;
    Q1.transposeInPlace();
    Q2.transposeInPlace();
    Q3.transposeInPlace();
    Q4.transposeInPlace();
    Q5.transposeInPlace();
    Q6.transposeInPlace();
    Q7.transposeInPlace();

    // unit vectors in 1st reference
    Vector3d e_1[7];
    auto Q_tmp = Q1;
    e_1[0] = Q_tmp * ev[0];
    e_1[1] = (Q_tmp *= Q2) * ev[1];
    e_1[2] = (Q_tmp *= Q3) * ev[2];
    e_1[3] = (Q_tmp *= Q4) * ev[3];
    e_1[4] = (Q_tmp *= Q5) * ev[4];
    e_1[5] = (Q_tmp *= Q6) * ev[5];
    e_1[6] = (Q_tmp *= Q7) * ev[6];

    Vector3d r[7];
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

    auto cr0 = e_1[0].cross(r[0]);
    auto cr1 = e_1[1].cross(r[1]);
    auto cr2 = e_1[2].cross(r[2]);
    auto cr3 = e_1[3].cross(r[3]);
    auto cr4 = e_1[4].cross(r[4]);
    auto cr5 = e_1[5].cross(r[5]);
    auto cr6 = e_1[6].cross(r[6]);

    // jacobian matrix
    Eigen::Matrix<double, 6, 7> J;
    J(0, 0) = e_1[0][0];
    J(1, 0) = e_1[0][1];
    J(2, 0) = e_1[0][2];
    J(0, 1) = e_1[1][0];
    J(1, 1) = e_1[1][1];
    J(2, 1) = e_1[1][2];
    J(0, 2) = e_1[2][0];
    J(1, 2) = e_1[2][1];
    J(2, 2) = e_1[2][2];
    J(0, 3) = e_1[3][0];
    J(1, 3) = e_1[3][1];
    J(2, 3) = e_1[3][2];
    J(0, 4) = e_1[4][0];
    J(1, 4) = e_1[4][1];
    J(2, 4) = e_1[4][2];
    J(0, 5) = e_1[5][0];
    J(1, 5) = e_1[5][1];
    J(2, 5) = e_1[5][2];
    J(0, 6) = e_1[6][0];
    J(1, 6) = e_1[6][1];
    J(2, 6) = e_1[6][2];

    J(3, 0) = cr0[0];
    J(4, 0) = cr0[1];
    J(5, 0) = cr0[2];
    J(3, 1) = cr1[0];
    J(4, 1) = cr1[1];
    J(5, 1) = cr1[2];
    J(3, 2) = cr2[0];
    J(4, 2) = cr2[1];
    J(5, 2) = cr2[2];
    J(3, 3) = cr3[0];
    J(4, 3) = cr3[1];
    J(5, 3) = cr3[2];
    J(3, 4) = cr4[0];
    J(4, 4) = cr4[1];
    J(5, 4) = cr4[2];
    J(3, 5) = cr5[0];
    J(4, 5) = cr5[1];
    J(5, 5) = cr5[2];
    J(3, 6) = cr6[0];
    J(4, 6) = cr6[1];
    J(5, 6) = cr6[2];

    return J;
}

host_fn VectorXd Gen3Eigen::collision_check(const VectorXd &joint_position) {
    Matrix3d Q1, Q2, Q3, Q4, Q5, Q6, Q7;

    ArrayXd p(joint_position);
    ArrayXd s = p.sin();
    ArrayXd c = p.cos();

    // note: these are stored column-wise
    Q1 <<
       c[0], -s[0], 0,
       -s[0], -c[0], 0,
       0, 0, -1;
    Q2 <<
       c[1], 0,  s[1],
       -s[1], 0,  c[1],
       0, -1, 0;
    Q3 <<
       c[2], 0, -s[2],
       -s[2], 0, -c[2],
       0,  1, 0;
    Q4 <<
       c[3], 0,  s[3],
       -s[3], 0,  c[3],
       0, -1, 0;
    Q5 <<
       c[4], 0, -s[4],
       -s[4], 0, -c[4],
       0,  1, 0;
    Q6 <<
       c[5], 0,  s[5],
       -s[5], 0,  c[5],
       0, -1, 0;
    Q7 <<
       c[6], 0, -s[6],
       -s[6], 0, -c[6],
       0,  1, 0;
    Q1.transposeInPlace();
    Q2.transposeInPlace();
    Q3.transposeInPlace();
    Q4.transposeInPlace();
    Q5.transposeInPlace();
    Q6.transposeInPlace();
    Q7.transposeInPlace();

    Vector3d p_orig(0, 0, 0);
    Vector3d p_j2;
    Vector3d p_j3;
    Vector3d p_j4;
    Vector3d p_j6;
    Vector3d p_ee;

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_j2 = p_tmp;
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_j3 = p_tmp;
    p_tmp += (Q_tmp *= Q3) * dv[2];
    p_j4 = p_tmp;
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_j6 = p_tmp;
    p_tmp += (Q_tmp *= Q6) * dv[5];
    p_tmp += (Q_tmp *= Q7) * dv[6];
    p_ee = p_tmp;

    const double r1sqr = 0.0875 * 0.0875;
    const double r2sqr = 0.0875 * 0.0875;

    // Self collisions sqr
    double dist1sqr     = two_segment_distance_sqr(p_orig, p_j2, p_j6, p_ee) - r1sqr;
    double dist2J2sqr   = two_segment_distance_sqr(p_j2, p_j2, p_j6, p_ee) - r2sqr; // distance sqr from J2 to J6 EE line (sphere)
    double dist2Msqr    = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r1sqr;  // distance sqr from J2 J3 line to J6 EE line (capsule)
    double dist2sqr = std::min(dist2J2sqr, dist2Msqr);

    // todo: remove collision detection from here

    // Collision with table sqr
    const double r4table = 0.05; // todo: validate dimensions
    const double r6table = 0.04; // todo: validate dimensions

    const Vector3d p_table(0, 0, -0.0025); // todo: Correct coords (z or y) ??
    double distTJ4 = p_j4[2] - p_table[2] - r4table;
    double distTJ6 = p_j6[2] - p_table[2] - r6table;
    double distTEE = p_ee[2] - p_table[2] - r6table;

    // VectorXd of distance min sqr and distance from table sqr
    VectorXd distMin(5);
    distMin << dist1sqr, dist2sqr, distTJ4, distTJ6, distTEE;
    // todo: Add collision_check function for more obstacles; VectorXd distMin(5) hard coded size

    return distMin;
}

host_fn VectorXd Gen3Eigen::constraints(const VectorXd &pos, const VectorXd &vel, const VectorXd &acc) {
    MatrixXd p(pos);
    MatrixXd v(vel);
    MatrixXd a(acc);

    dynamics(p, v, a);
    // 5 collision results
    // 2 position constraints
    // 7 velocity constraints
    // 7 torque constraints
    VectorXd result(5 + 2 + 7 * 2);

    // distance to collision >= 0
    auto tmp_coll = collision_check(pos);
    result[0] = -tmp_coll[0]; // dist1sqr
    result[1] = -tmp_coll[1]; // dist2sqr
    result[2] = -tmp_coll[2]; // distTJ4
    result[3] = -tmp_coll[3]; // distTJ6
    result[4] = -tmp_coll[4]; // distTEE

    // position
    result[5] = (abs(pos[3]) - pmax[3]) / pmax[3];
    result[6] = (abs(pos[5]) - pmax[5]) / pmax[5];

    // velocity
    for (int j = 0; j < 7; j++)
        result[j+7] = (abs(vel[j]) - vmax[j]) / vmax[j];

    auto f = _efforts.col(0);
    for (int j = 0; j < 7; j++)
        result[j+14] = (abs(f[j]) - tau_max[j]) / tau_max[j];

    return result;
}

host_fn VectorXd Gen3Eigen::constraints(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc) {
    const auto points = pos.cols();
    VectorXd result(ncon(points));
    constraints(pos, vel, acc, result.data());
    return result;
}

host_fn VectorXd Gen3Eigen::constraints(const TrajectoryEigen &traj) {
    VectorXd result(ncon(traj.t.size()));
    constraints(traj.pos, traj.vel, traj.acc, result.data());
    return result;
}

host_fn void Gen3Eigen::constraints(const TrajectoryEigen &traj, double *dst) {
    constraints(traj.pos, traj.vel, traj.acc, dst);
}

host_fn void Gen3Eigen::constraints(const MatrixXd &pos, const MatrixXd &vel, const MatrixXd &acc, double *dst) {
    const auto points = pos.cols();

    dynamics(pos, vel, acc);

    for (int i = 0; i < points; i++) {
        // collision
        VectorXd p = pos.col(i);
        VectorXd tmp_coll = collision_check(p);
        dst[0] = -tmp_coll[0]; // dist1sqr
        dst[1] = -tmp_coll[1]; // dist2sqr
        dst[2] = -tmp_coll[2]; // distTJ4
        dst[3] = -tmp_coll[3]; // distTJ6
        dst[4] = -tmp_coll[4]; // distTEE
        dst += 5;

        // position
        dst[0] = (abs(pos(3, i)) - pmax[3]) / pmax[3];
        dst[1] = (abs(pos(5, i)) - pmax[5]) / pmax[5];
        dst += 2;

        // velocity
        for (int j = 0; j < 7; j++)
            dst[j] = (abs(vel(j, i)) - vmax[j]) / vmax[j];
        dst += 7;

        auto f = _efforts.col(i);
        for (int j = 0; j < 7; j++)
            dst[j] = (abs(f[j]) - tau_max[j]) / tau_max[j];
        dst += 7;
    }
}

host_fn bool Gen3Eigen::validate_task(const MatrixXd &task) {
    auto ci = constraints(VectorXd(task.col(0)), task.col(1), task.col(2));
    auto cf = constraints(VectorXd(task.col(3)), task.col(4), task.col(5));
    return (ArrayXd(ci) <= 0).all() && (ArrayXd(cf) <= 0).all(); // todo: verify
}

//------ Kinova Link6 6DOF manipulator functions ---------------------------------------

host_fn Link6::Link6() : Manipulator(6) {
    p_base = {0, 0, 0.0530f};

    // mass
    m[0] = 4.8257f;
    m[1] = 5.9860f;
    m[2] = 3.4159f;
    m[3] = 2.0849f;
    m[4] = 2.0076f;
    m[5] = 1.5193f; // todo: modify with gripper (as of now, no gripper)

    // inertial tensor
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
    av[5] = {0.0f, 0.00050624f, -0.00388589f}; // todo: modify with gripper (as of now, no gripper)

    // vector to next joint
    dv[0] = {0.11024, -0.06926, -0.1375};   // 0 -> 1
    dv[1] = {0.0, 0.4850, 0.0};             // 1 -> 2
    dv[2] = {0.0, -0.15216, -0.0917};       // 2 -> 3
    dv[3] = {0.0, -0.06296, -0.22275};      // 3 -> 4
    dv[4] = {0.08703, 0.0860, -0.07692};    // 4 -> 5
    dv[5] = {0.0, 0.0, -0.920};             // 5 -> endeffector (todo: add gripper)

    // todo: add option to know if tool is closed or opened (difference of 0.0135 in z)

    // center of mass (from next joint)
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];

    // unit joint direction
    ev[0] = {0, 0, 1};
    ev[1] = {0, 0, 1};
    ev[2] = {0, 0, 1};
    ev[3] = {0, 0, 1};
    ev[4] = {0, 0, 1};
    ev[5] = {0, 0, 1};

    // kinematic and dynamic constraints
    pmax = {INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL, INF_REAL}; // rad todo: make sure this is true
    vmax = {3.4907f, 3.4907f, 3.4907f, 5.5851f, 5.5851f, 5.5851f}; // rad/s
    amax = {deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600), deg2rad(600)}; // rad/s^2
    tau_max = {210, 210, 210, 100, 100, 100}; // Nm
}

host_fn void Link6::dynamics(const Trajectory &traj, Matrix &efforts) {
    dynamics(traj.pos, traj.vel, traj.acc, efforts);
}

host_fn void Link6::dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc, Matrix &efforts) {

    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t;
    Vec3 w1, w2, w3, w4, w5, w6;
    Vec3 wd1, wd2, wd3, wd4, wd5, wd6;
    Vec3 cdd0 = {0, 0, 9.81f};
    Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6;
    Vec3 f1, f2, f3, f4, f5, f6;
    Vec3 n1, n2, n3, n4, n5, n6;

    const u32 joints = pos.rows;
    const u32 points = pos.cols;

    // loop all points
    for (u32 i = 0; i < points; i++) {
        auto v = &vel.data[i * joints];
        auto a = &acc.data[i * joints];

        // SIMD compute sines and cosines note: approx 10% faster
        auto p = pos.col(i);
        Array s(6);
        Array c(6);
        blast::sincos(p, s, c);
//         auto p = &pos.data[i * joints];

// #if BLAST_SIZEOF_REAL == 8
//         real s[8];
//         real c[8];
//         __m256d s_tmp;
//         __m256d c_tmp;
//         for (u32 j = 0; j < 8; j += 4) {
//             __m256d angle_v = _mm256_load_pd(p + j);
//             s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//             _mm256_storeu_pd(s + j, s_tmp);
//             _mm256_storeu_pd(c + j, c_tmp);
//         }
// #else
//         real s[8];
//         real c[8];
//         __m256 s_tmp;
//         __m256 c_tmp;
//         __m256 angle_v = _mm256_load_ps(p);
//         s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//         _mm256_storeu_ps(s, s_tmp);
//         _mm256_storeu_ps(c, c_tmp);
// #endif

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
        efforts(0, i) = n1.z;
        efforts(1, i) = n2.z;
        efforts(2, i) = n3.z;
        efforts(3, i) = n4.z;
        efforts(4, i) = n5.z;
        efforts(5, i) = n6.z;
    }
}

host_fn void Link6::dynamics(const Matrix &pos, const Matrix &vel, const Matrix &acc) {
    const auto points = pos.cols;
    const auto joints = pos.rows;
    if (_efforts.cols != points || _efforts.rows != joints)
        _efforts.resize(joints, points);

    dynamics(pos, vel, acc, _efforts);
}

host_fn Array Link6::forward_kinematics(const Array &joint_position) {

    //  - note: manual SIMD (10% better performance than using sincos function on arrays like commented below)
    // real s[8];
    // real c[8];
    // auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;

    Array s(8);
    Array c(8);
    blast::sincos(joint_position, s, c);

// #if BLAST_SIZEOF_REAL == 8
//     __m256d s_tmp;
//     __m256d c_tmp;
//     for (u32 i = 0; i < 8; i += 4) {
//         __m256d angle_v = _mm256_load_pd(p + i);
//         s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
//         _mm256_storeu_pd(s + i, s_tmp);
//         _mm256_storeu_pd(c + i, c_tmp);
//     }
// #else
//     __m256 s_tmp;
//     __m256 c_tmp;
//     __m256 angle_v = _mm256_load_ps(p);
//     s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
//     _mm256_storeu_ps(s, s_tmp);
//     _mm256_storeu_ps(c, c_tmp);
// #endif

    // todo: cleanup
    // create aliases instead of actual arrays to allocate the memory on the stack (1.7x performance)
    // real s_data[8]; // nearest factor of 4
    // real c_data[8]; // nearest factor of 4
    // Array s(s_data, 8);
    // Array c(c_data, 8);
    // Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    // sincos(joint_position, s, c);

    // note: these are stored column-wise
    Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
    Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
    Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
    Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
    Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
    Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0}; // todo: double check

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_tmp += (Q_tmp *= Q3) * dv[2];
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_tmp += (Q_tmp *= Q6) * dv[5];

    Array pose(6);
    pose[0] = p_tmp.x;
    pose[1] = p_tmp.y;
    pose[2] = p_tmp.z;
    // todo: Q_tmp(0, 0) and Q_tmp(2, 2) must not be zero!!
    // pose[3] = atan2(Q_tmp(1, 0), Q_tmp(0, 0));
    // pose[4] = atan2(-Q_tmp(2, 0), sqrt(Q_tmp(2, 1)*Q_tmp(2, 1) + Q_tmp(2, 2)*Q_tmp(2, 2)));
    // pose[5] = atan2(Q_tmp(2, 1), Q_tmp(2, 2));
    return pose;
}

host_fn Matrix Link6::forward_kinematics(const Matrix &joint_positions) {
    Matrix pose(12, joint_positions.cols);

    for (u32 point = 0; point < joint_positions.cols; point++) {
        Array s(8);
        Array c(8);
        blast::sincos(joint_positions.col(point), s, c);

        Mat3 Q1 = {c[0], -s[0], 0, -s[0], -c[0], 0, 0, 0, -1};
        Mat3 Q2 = {c[1], 0, -s[1], -s[1], 0, -c[1], 0, 1, 0};
        Mat3 Q3 = {c[2], -s[2], 0, -s[2], -c[2], 0, 0, 0, -1};
        Mat3 Q4 = {c[3], 0, -s[3], -s[3], 0, -c[3], 0, 1, 0};
        Mat3 Q5 = {c[4], 0, -s[4], -s[4], 0, -c[4], 0, 1, 0};
        Mat3 Q6 = {0, s[5], c[5], 0, c[5], s[5], -1, 0, 0};
        auto p_tmp = p_base;
        auto Q_tmp = Q1;
        p_tmp += Q_tmp * dv[0];
        p_tmp += (Q_tmp *= Q2) * dv[1];
        p_tmp += (Q_tmp *= Q3) * dv[2];
        p_tmp += (Q_tmp *= Q4) * dv[3];
        p_tmp += (Q_tmp *= Q5) * dv[4];
        p_tmp += (Q_tmp *= Q6) * dv[5];
        pose(0, point) = p_tmp.x;
        pose(1, point) = p_tmp.y;
        pose(2, point) = p_tmp.z;
        pose(3, point) = Q_tmp[0];
        pose(4, point) = Q_tmp[1];
        pose(5, point) = Q_tmp[2];
        pose(6, point) = Q_tmp[3];
        pose(7, point) = Q_tmp[4];
        pose(8, point) = Q_tmp[5];
        pose(9, point) = Q_tmp[6];
        pose(10, point) = Q_tmp[7];
        pose(11, point) = Q_tmp[8];
    }
    return pose;
}

host_fn Array Link6::collision_check(const Array &joint_position) {
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
    Vec3 p_j6;
    Vec3 p_ee;

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp * dv[0];
    p_j2 = p_tmp;
    p_tmp += (Q_tmp *= Q2) * dv[1];
    p_j3 = p_tmp;
    p_tmp += (Q_tmp *= Q3) * dv[2];
    p_j4 = p_tmp;
    p_tmp += (Q_tmp *= Q4) * dv[3];
    p_tmp += (Q_tmp *= Q5) * dv[4];
    p_j6 = p_tmp;
    p_tmp += (Q_tmp *= Q6) * dv[5];
    p_ee = p_tmp;

    const real r1sqr = 0.140 * 0.140; // size of capsule 1
    const real r2sqr = 0.140 * 0.140; // size of capsule 2
    const real r3sqr = 0.100 * 0.100; // size of capsule 3

    // Self collisions sqr
    real dist1sqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r1sqr; // todo: Fix and finish

    // real dist1sqr = two_segment_distance_sqr(p_orig, p_j2, p_j6, p_ee) - r1sqr;
    // // real dist2J2sqr = ((p_ee.x - p_j2.x) + (p_j6.x - p_j2.x)) * ((p_ee.x - p_j2.x) + (p_j6.x - p_j2.x)) + ((p_ee.y - p_j2.y) + (p_j6.y - p_j2.y)) * ((p_ee.y - p_j2.y) + (p_j6.y - p_j2.y)) + ((p_ee.z - p_j2.z) + (p_j6.z - p_j2.z)) * ((p_ee.z - p_j2.z) + (p_j6.z - p_j2.z)) - r2sqr;
    // real dist2J2sqr = two_segment_distance_sqr(p_j2, p_j2, p_j6, p_ee) - r2sqr; // distance sqr from J2 to J6 EE line (sphere)
    // real dist2Msqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r1sqr;  // distance sqr from J2 J3 line to J6 EE line (capsule)
    // real dist2sqr = dist2J2sqr <= dist2Msqr ? dist2J2sqr : dist2Msqr;

    // todo: remove collision detection from here

    // Collision with table sqr
    // const real r4table = 0.05; // todo: validate dimensions
    // const real r6table = 0.04; // todo: validate dimensions

    // const Vec3 p_table(0, 0, -0.0025); // todo: Correct coords (z or y) ??
    // real distTJ4 = p_j4.z - p_table.z - r4table;
    // real distTJ6 = p_j6.z - p_table.z - r6table;
    // real distTEE = p_ee.z - p_table.z - r6table;

    // Array of distance min sqr and distance from table sqr
    Array distMin(2);
    distMin = {dist1sqr};
    // todo: Add collision_check function for more obstacles; Array distMin(5) hard coded size

    return distMin;
}

host_fn Array Link6::constraints(const Array &pos, const Array &vel, const Array &acc) {
    Matrix p(pos);
    Matrix v(vel);
    Matrix a(acc);

    dynamics(p, v, a);
    // 5 collision results
    // 0 position constraints
    // 6 velocity constraints
    // 6 torque constraints
    Array result(ncon(1));

    // distance to collision >= 0
    auto tmp_coll = collision_check(pos);
    result[0] = -tmp_coll[0]; // dist1sqr
    result[1] = -tmp_coll[1]; // dist2sqr
    result[2] = -tmp_coll[2]; // distTJ4
    result[3] = -tmp_coll[3]; // distTJ6
    result[4] = -tmp_coll[4]; // distTEE

    // velocity
    for (int j = 0; j < joints; j++)
        result[j+5] = (abs(vel[j]) - vmax[j]) / vmax[j];

    auto f = _efforts.col(0); Assert(f.is_alias);
    for (int j = 0; j < joints; j++)
        result[j+5+joints] = (abs(f[j]) - tau_max[j]) / tau_max[j];

    return result;
}

host_fn Array Link6::constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc) {
    const auto points = pos.cols;
    Array result(ncon(points));
    constraints(pos, vel, acc, result.data);
    return result;
}

host_fn Array Link6::constraints(const Trajectory &traj) {
    Array result(ncon(traj.t.size));
    constraints(traj.pos, traj.vel, traj.acc, result.data);
    return result;
}

host_fn void Link6::constraints(const Trajectory &traj, real *dst) {
    constraints(traj.pos, traj.vel, traj.acc, dst);
}

host_fn void Link6::constraints(const Matrix &pos, const Matrix &vel, const Matrix &acc, real *dst) {
    const auto points = pos.cols;

    dynamics(pos, vel, acc);

    for (u32 i = 0; i < points; i++) {
        // collision
        Array p = pos.col(i); Assert(p.is_alias);
        auto tmp_coll = collision_check(p);
        dst[0] = -tmp_coll[0]; // dist1sqr
        dst[1] = -tmp_coll[1]; // dist2sqr
        dst[2] = -tmp_coll[2]; // distTJ4
        dst[3] = -tmp_coll[3]; // distTJ6
        dst[4] = -tmp_coll[4]; // distTEE
        dst += 5;

        // velocity
        for (int j = 0; j < joints; j++)
            dst[j] = (abs(vel(j, i)) - vmax[j]) / vmax[j];
        dst += joints;

        auto f = _efforts.col(i); Assert(f.is_alias);
        for (int j = 0; j < joints; j++)
            dst[j] = (abs(f[j]) - tau_max[j]) / tau_max[j];
        dst += joints;
    }
}

host_fn bool Link6::validate_task(const Matrix &task) {
    auto ci = constraints(task.col(0), task.col(1), task.col(2));
    auto cf = constraints(task.col(3), task.col(4), task.col(5));
    return array_max(ci) > 0 && array_max(cf) > 0 ? false: true;
}


// note: CUDA stuff, only enabled if compiling for Nvidia GPUs
#ifdef BLAST_ENABLE_TESTS
TEST_CASE("SelfCollision", "[Manipulator]") {
    using namespace blast;
    Gen3_7DOF manip;
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

    auto dist_sqr_min_1 = manip.collision_check(theta1);
    auto dist_sqr_min_2 = manip.collision_check(theta2);
    auto dist_sqr_min_3 = manip.collision_check(theta3);
    auto dist_sqr_min_4 = manip.collision_check(theta4);

    CHECK(dist_sqr_min_1[0] > 0);
    CHECK(dist_sqr_min_1[1] > 0);
    CHECK(dist_sqr_min_2[0] > 0);
    CHECK(dist_sqr_min_2[1] > 0);
    CHECK(dist_sqr_min_3[0] < 0);
    CHECK(dist_sqr_min_3[1] < 0);
    CHECK(dist_sqr_min_4[0] < 0);
    CHECK(dist_sqr_min_4[1] > 0);
}

TEST_CASE("Link6Checks", "[Manipulator]") {
    using namespace blast;
    Link6 manip;

    const u32 points = 10;
    const u32 joints = 6;
    const u32 p = 5;
    const u32 ncontrol = 12;

    real amp = 10;
    Matrix task(joints, 6);
    for (u32 i = 0; i < task.rows; i++) {
        for (u32 j = 0; j < task.cols; j++) {
            task(i, j) = amp * get_random();
        }
    }
    // random optimization vector
    Array x(joints * (ncontrol - 6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x[x.size - 1] = 3.f;
    Bspline bspline(ncontrol, points, p, joints);
    bspline.compute_trajectory(x, task);
    Array r = manip.constraints(bspline.traj);
}

TEST_CASE("EigenCorrectness", "[Manipulator]") {
    using namespace blast;
    Gen3_7DOF manip;
    Gen3Eigen manip_eigen;

    const u32 points = 10;
    const u32 joints = 7;
    const u32 p = 5;
    const u32 ncontrol = 12;

    // random task
    real amp = 10;
    Matrix task(joints, 6);
    for (u32 i = 0; i < task.rows; i++) {
        for (u32 j = 0; j < task.cols; j++) {
            task(i, j) = amp * get_random();
        }
    }
    MatrixXd task_eigen(joints, 6);
    for (u32 i = 0; i < task.size; i++)
        task_eigen.data()[i] = task.data[i];
    // random optimization vector
    Array x(joints * (ncontrol - 6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x[x.size - 1] = 3.f;
    VectorXd x_eigen(joints * (ncontrol - 6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x_eigen[i] = x[i];
    Bspline bspline(ncontrol, points, p, joints);
    BsplineEigen bspline_eigen(ncontrol, points, p, joints);
    bspline.compute_trajectory(x, task);
    bspline_eigen.compute_trajectory(x_eigen, task_eigen);

    Array r = manip.constraints(bspline.traj);
    VectorXd r_eigen = manip_eigen.constraints(bspline_eigen.traj);

    double max_err = 0;
    for (int i = 0; i < r.size; i++) {
        max_err = std::max(max_err, std::abs(r[i] - r_eigen[i]));
        if (max_err > 0.1) {
            break;
        }
    }

    REQUIRE(max_err < 0.01);
}

TEST_CASE("ManipSpeedTest", "[Manipulator]") {
    using namespace blast;
    Gen3_7DOF manip;
    Gen3Eigen manip_eigen;

    const u32 points = 256;
    const u32 joints = 7;
    const u32 p = 5;
    const u32 ncontrol = 24;

    // random task
    real amp = 10;
    Matrix task(joints, 6);
    for (u32 i = 0; i < task.rows; i++) {
        for (u32 j = 0; j < task.cols; j++) {
            task(i, j) = amp * get_random();
        }
    }
    MatrixXd task_eigen(joints, 6);
    for (u32 i = 0; i < task.size; i++)
        task_eigen.data()[i] = task.data[i];
    // random optimization vector
    Array x(joints * (ncontrol - 6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x[x.size - 1] = 3.f;
    VectorXd x_eigen(joints * (ncontrol - 6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x_eigen[i] = x[i];
    Bspline bspline(ncontrol, points, p, joints);
    BsplineEigen bspline_eigen(ncontrol, points, p, joints);

    Array r(manip.ncon(points));
    BENCHMARK("Blast constraints speed") {
        bspline.compute_trajectory(x, task);
        r = manip.constraints(bspline.traj);
        return r[122];
    };

    VectorXd r_eigen(manip_eigen.ncon(points));
    BENCHMARK("Blast + Eigen constraints speed") {
        bspline_eigen.compute_trajectory(x_eigen, task_eigen);
        r_eigen = manip_eigen.constraints(bspline_eigen.traj);
        return r_eigen[122];
    };

}
#endif // tests
} // namespace blast
