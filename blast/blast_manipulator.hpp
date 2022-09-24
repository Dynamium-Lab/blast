#pragma once

#include "blast_math.hpp"
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

    virtual void dynamics(const Pva& pva, Matrix& efforts) = 0;
    virtual void dynamics(const Matrix& pos, const Matrix& vel, const Matrix& acc, Matrix& efforts) = 0;

    virtual Array validate(const Pva& pva) {
        return Array();
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

    virtual void dynamics(const Pva& pva, Matrix& efforts) override;
    virtual void dynamics(const Matrix& pos, const Matrix& vel, const Matrix& acc, Matrix& efforts) override;
    void init_dynamics(real mass=0);
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

    // compute joint torque as a function of trajector (pva)
    virtual void dynamics(const Pva& pva, Matrix& efforts) override;
    virtual void dynamics(const Matrix& pos, const Matrix& vel, const Matrix& acc, Matrix& efforts) override;

    // compute forward kinematics for 1 point
    Array forward_kinematics(Array& joint_position);
};

struct Gen3_7DOF : public Manipulator {
    Vec3 p_base;
    real m[7];
    Mat3 I[7];
    Vec3 av[7];
    Vec3 dv[7];
    Vec3 sv[7];
    Vec3 ev[7];

    // default constructor
    Gen3_7DOF();

    // compute joint torque as a function of trajector (pva)
    virtual void dynamics(const Pva& pva, Matrix& efforts) override;
    virtual void dynamics(const Matrix& pos, const Matrix& vel, const Matrix& acc, Matrix& efforts) override;

    // compute forward kinematics for 1 point
    Array forward_kinematics(const Array& joint_position);
    Matrix forward_kinematics(const Matrix& joint_positions);

    // compute jacobian matrix
    Matrix Gen3_7DOF::jacobian_matrix(const Array& joint_position);

    // check collision
    Array collision_dist_sqr(const Array& joint_position);

    // check all constraints on the manipulator for 1 point
    Array validate(const Array& pos, const Array& vel, const Array& acc);

    // check all constraints on the manipulator for a trajectory
    Array validate(const Matrix& pos, const Matrix& vel, const Matrix& acc);

    // check all constraints on the manipulator for a trajectory
    virtual Array validate(const Pva& pva) override;
};


// NOTE: the following structure is only usefull when using CUDA for Nvidia GPUs
#ifdef __NVCC__
// note: must be a global variable because it is __constant__
// note: much faster if it's __constant__ because all threads access the same location (heavy use of broadcast operations)
const int MANIP_BROADCAST_ARENA_SIZE = 256;
__constant__ real manip_broadcast_arena[MANIP_BROADCAST_ARENA_SIZE];

struct cuGen3_7DOF {
    //note: for the mapping to work, all the data must be on the stack!!
    Vec3 p_base;
    real m[7];
    Mat3 I[7];
    Vec3 av[7];
    Vec3 dv[7];
    Vec3 sv[7];
    Vec3 ev[7];

    real pmax[7];
    real vmax[7];
    real tau_max[7];

    real* device_constraints; // 21 * npoints

    real* host_constraints; // 21 * npoints

    // compute the parameters according to the mass
    __host__ void init(real mass, u32 npoints);

    // fetch the latest computed constraints from the device to the host
    __host__ void fetch_constraints(u32 npoints);

    // compute the constraints based on the current trajectory point
    __device__ void compute_constraints(const real pos[7], const real vel[7], const real acc[7], real* con);
};
#endif




//------ Universal Robots UR5e manipulator functions ---------------------------------

inline void ManipulatorUR5::dynamics(const Pva& pva, Matrix& efforts) {
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
    for (u32 i = 0; i < pva.points; i++) {
        vel1 = pva.vel(0, i);
        vel2 = pva.vel(1, i);
        vel3 = pva.vel(2, i);
        vel4 = pva.vel(3, i);
        vel5 = pva.vel(4, i);
        vel6 = pva.vel(5, i);
        acc1 = pva.acc(0, i);
        acc2 = pva.acc(1, i);
        acc3 = pva.acc(2, i);
        acc4 = pva.acc(3, i);
        acc5 = pva.acc(4, i);
        acc6 = pva.acc(5, i);

        // SIMD compute sines and cosines note: approx 10% faster
        const auto p = &pva.pos.data[i * pva.joints];

#if BLAST_SIZEOF_REAL == 8
        real s[6];
        real c[6];
        // first 4
        {
            __m256d s_tmp;
            __m256d c_tmp;
            __m256d angle_v = _mm256_load_pd(p);
            s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
            _mm256_storeu_pd(s, s_tmp);
            _mm256_storeu_pd(c, c_tmp);
        }
        // 5 and 6
        {
            __m128d angle_v = _mm_load_pd(p+4);
            __m128d s_tmp;
            __m128d c_tmp;
            s_tmp = _mm_sincos_pd(&c_tmp, angle_v);
            _mm_storeu_pd(s+4, s_tmp);
            _mm_storeu_pd(c+4, c_tmp);
        }
#else
        real s[8];
        real c[8];
        __m256 s_tmp;
        __m256 c_tmp;
        __m256 angle_v = _mm256_load_ps(p);
        s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
        _mm256_storeu_ps(s, s_tmp);
        _mm256_storeu_ps(c, c_tmp);
#endif

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
        w12 = vel1*e1;
        w23 = Q2t*w12 + vel2*e2;
        w34 = Q3t*w23 + vel3*e3;
        w45 = Q4t*w34 + vel4*e4;
        w56 = Q5t*w45 + vel5*e5;
        w67 = Q6t*w56 + vel6*e6;

        wd12 = acc1*e1;
        cdd12 = Q1t*cdd01
                + cross(wd12, (a1v + s1))
                + cross(w12, cross(w12, (a1v+s1)));

        wd23 = Q2t*wd12 + acc2*e2 + vel2*cross((Q2t*w12), (e2));
        cdd23 = Q2t*cdd12
                + cross(wd23, a2v + s2)
                + cross(w23, cross(w23, (a2v+s2)))
                - Q2t*cross(wd12, s1)
                - Q2t*cross(w12, cross(w12, s1));

        wd34 = Q3t*wd23 + acc3*e3 + vel3*cross((Q3t*w23), (e3));
        cdd34 = Q3t*cdd23
                + cross(wd34, a3v + s3)
                + cross(w34, cross(w34, (a3v+s3)))
                - Q3t*cross(wd23, s2)
                - Q3t*cross(w23, cross(w23, s2));

        wd45 = Q4t*wd34 + acc4*e4 + vel4*cross((Q4t*w34), (e4));
        cdd45 = Q4t*cdd34
                + cross(wd45, a4v + s4)
                + cross(w45, cross(w45, (a4v+s4)))
                - Q4t*cross(wd34, s3)
                - Q4t*cross(w34, cross(w34, s3));

        wd56 = Q5t*wd45 + acc5*e5 + vel5*cross((Q5t*w45), (e5));
        cdd56 = Q5t*cdd45
                + cross(wd56, a5v + s5)
                + cross(w56, cross(w56, (a5v+s5)))
                - Q5t*cross(wd45, s4)
                - Q5t*cross(w45, cross(w45, s4));

        wd67 = Q6t*wd56 + acc6*e6 + vel6*cross((Q6t*w56), (e6));
        cdd67 = Q6t*cdd56
                + cross(wd67, a6v + s6)
                + cross(w67, cross(w67, (a6v+s6)))
                - Q6t*cross(wd56, s5)
                - Q6t*cross(w56, cross(w56, s5));

        //-- dynamics
        f6 = Q6*(m6*cdd67);
        n6 = Q6*(I6*wd67 + cross(w67, I6*w67) + cross(a6v+s6, Q6t*f6));

        f5 = Q5*(m5*cdd56 + f6);
        n5 = Q5*(I5*wd56 + cross(w56, I5*w56) + n6 + cross(a5v+s5, Q5t*f5) - cross(s5, f6));

        f4 = Q4*(m4*cdd45 + f5);
        n4 = Q4*(I4*wd45 + cross(w45, I4*w45) + n5 + cross(a4v+s4, Q4t*f4) - cross(s4, f5));

        f3 = Q3*(m3*cdd34 + f4);
        n3 = Q3*(I3*wd34 + cross(w34, I3*w34) + n4 + cross(a3v+s3, Q3t*f3) - cross(s3, f4));

        f2 = Q2*(m2*cdd23 + f3);
        n2 = Q2*(I2*wd23 + cross(w23, I2*w23) + n3 + cross(a2v+s2, Q2t*f2) - cross(s2, f3));

        f1 = Q1*(m1*cdd12 + f2);
        n1 = Q1*(I1*wd12 + cross(w12, I1*w12) + n2 + cross(a1v+s1, Q1t*f1) - cross(s1, f2));

        //-- extract torques (last element of each moment vector)
        efforts(0, i) = n1.z;
        efforts(1, i) = n2.z;
        efforts(2, i) = n3.z;
        efforts(3, i) = n4.z;
        efforts(4, i) = n5.z;
        efforts(5, i) = n6.z;
    }
}

inline void ManipulatorUR5::init_dynamics(real mass) {
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

        I6(0, 0) += m6_old*(s6.x-s6_old.x)*(s6.x-s6_old.x) + mass*s6.x*s6.x;
        I6(1, 1) += m6_old*(s6.y-s6_old.y)*(s6.y-s6_old.y) + mass*s6.y*s6.y;
        I6(2, 2) += m6_old*(s6.z-s6_old.z)*(s6.z-s6_old.z) + mass*s6.z*s6.z;
    }
}


//------ Kinova Gen3 Lite manipulator functions --------------------------------------

inline Gen3Lite::Gen3Lite() : Manipulator(6) {
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
    I[0] = {0.0016595f, 2e-08f,    3.6e-07f,   2e-08f,    0.00140355f, 0.00034927f, 3.6e-07f,    0.00034927f, 0.00089493f};
    I[1] = {0.0114928f, 1e-06f,    1.6e-07f,   1e-06f,    0.00102851f, 0.00140765f, 1.6e-07f,    0.00140765f, 0.01133492f};
    I[2] = {0.0016326f, 7.11E-06f, 1.54e-06f,  7.11E-06f, 0.00029798f, 9.587E-05f,  1.54E-06f,   9.587E-05f,  0.00169091f};
    I[3] = {0.0006910f, 2.4E-07f,  0.0001648f, 2.4E-07f,  0.00078519f, 7.4E-07f,    0.00016483f, 7.4E-07f,    0.00034115f};
    I[4] = {0.0002127f, 5.21E-06f, 2.91E-06f,  5.21E-06f, 0.00106371f, 1.1E-07f,    2.91E-06f,   1.1E-07f,    0.00108465f};
    I[5] = {0.0003428f, 1.9e-7f,   1e-7f,      1.9e-7f,   0.00028915f, 0.00000027f, 1e-7f,       0.00000027f, 0.00013076f};

    // center of mass
    av[0] = {0.000025f,  0.022135f, 0.099377f};
    av[1] = {0.029983f,  0.211548f, 0.045303f};
    av[2] = {0.030156f, -0.095022f, 0.007356f};
    av[3] = {0.005752f,  0.010004f, 0.087192f};
    av[4] = {0.080565f,  0.009804f, 0.018728f};
    av[5] = {0.009930f,  0.009950f, 0.061360f};

    // vector to next joint
    dv[0] = {0.000f,   -0.030f,   0.115f};
    dv[1] = {0.000f,    0.280f,   0.000f};
    dv[2] = {0.000f,   -0.140f,   0.020f};
    dv[3] = {0.0285f,   0.000f,   0.105f};
    dv[4] = {-0.105f,   0.000f,   0.0285f};
    dv[5] = {0.000f,    0.000f,   0.130f};

    // center of mass (from next joint)
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];

    // unit joint direction
    ev[0] = { 0,  0,  1};
    ev[1] = { 0, -1,  0};
    ev[2] = { 0,  0, -1};
    ev[3] = { 0, -1,  0};
    ev[4] = { 1,  0,  0};
    ev[5] = {-1,  0,  0};

    // kinematic and dynamic constraints
    pmax = {2.69f, 2.69f, 2.69f, 2.59f, 2.57f, 2.59f}; // rad
    vmax = {1.6f,  1.6f,  1.6f,  1.6f,  1.6f,  3.2f};  // rad/s
    tau_max = {10,   14,   10,   7,    7,    7};    // Nm
    vmin = -vmax;
    pmin = -pmax;
    tau_min = -tau_max;
}

inline void Gen3Lite::dynamics(const Pva& pva, Matrix& efforts) {
    dynamics(pva.pos, pva.vel, pva.acc, efforts);
}

inline void Gen3Lite::dynamics(const Matrix& pos, const Matrix& vel, const Matrix& acc, Matrix& efforts) {

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
        auto p = &pos.data[i * joints];
#if BLAST_SIZEOF_REAL == 8 // double precision
        real s[6];
        real c[6];
        // first 4
        {
            __m256d s_tmp;
            __m256d c_tmp;
            __m256d angle_v = _mm256_load_pd(p);
            s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
            _mm256_storeu_pd(s, s_tmp);
            _mm256_storeu_pd(c, c_tmp);
        }
        // 5 and 6
        {
            __m128d angle_v = _mm_load_pd(p+4);
            __m128d s_tmp;
            __m128d c_tmp;
            s_tmp = _mm_sincos_pd(&c_tmp, angle_v);
            _mm_storeu_pd(s+4, s_tmp);
            _mm_storeu_pd(c+4, c_tmp);
        }
#else
        real s[8];
        real c[8];
        __m256 s_tmp;
        __m256 c_tmp;
        __m256 angle_v = _mm256_load_ps(p);
        s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
        _mm256_storeu_ps(s, s_tmp);
        _mm256_storeu_ps(c, c_tmp);
#endif

        Q1 = {c[0],  s[0],  0,        -s[0],  c[0],  0,         0,  0,  1};
        Q2 = {c[1],  0,     s[1],     -s[1],  0,     c[1],      0, -1,  0};
        Q3 = {c[2], -s[2],  0,        -s[2], -c[2],  0,         0,  0, -1};
        Q4 = {c[3],  0,     s[3],     -s[3],  0,     c[3],      0, -1,  0};
        Q5 = {0,     s[4], -c[4],      0,     c[4],  s[4],      1,  0,  0};
        Q6 = {0,     s[5],  c[5],      0,     c[5], -s[5],     -1,  0,  0};
        Q1t = transpose(Q1);
        Q2t = transpose(Q2);
        Q3t = transpose(Q3);
        Q4t = transpose(Q4);
        Q5t = transpose(Q5);
        Q6t = transpose(Q6);

        // note: this is the Newton algorithm in 'Element de robotique' course notes.
        //-- kinematics
        w12 = v[0]*ev[0];
        w23 = Q2t*w12 + v[1]*ev[1];
        w34 = Q3t*w23 + v[2]*ev[2];
        w45 = Q4t*w34 + v[3]*ev[3];
        w56 = Q5t*w45 + v[4]*ev[4];
        w67 = Q6t*w56 + v[5]*ev[5];

        wd12 = a[0]*ev[0];
        cdd12 = Q1t*cdd01
                + cross(wd12, av[0])
                + cross(w12, cross(w12, av[0]));

        wd23 = Q2t*wd12 + a[1]*ev[1] + v[1]*cross(Q2t*w12, ev[1]);
        cdd23 = Q2t*cdd12
                + cross(wd23, av[1])
                + cross(w23, cross(w23, av[1]))
                - Q2t*cross(wd12, sv[0])
                - Q2t*cross(w12, cross(w12, sv[0]));

        wd34 = Q3t*wd23 + a[2]*ev[2] + v[2]*cross(Q3t*w23, ev[2]);
        cdd34 = Q3t*cdd23
                + cross(wd34, av[2])
                + cross(w34, cross(w34, av[2]))
                - Q3t*cross(wd23, sv[1])
                - Q3t*cross(w23, cross(w23, sv[1]));

        wd45 = Q4t*wd34 + a[3]*ev[3] + v[3]*cross(Q4t*w34, ev[3]);
        cdd45 = Q4t*cdd34
                + cross(wd45, av[3])
                + cross(w45, cross(w45, av[3]))
                - Q4t*cross(wd34, sv[2])
                - Q4t*cross(w34, cross(w34, sv[2]));

        wd56 = Q5t*wd45 + a[4]*ev[4] + v[4]*cross(Q5t*w45, ev[4]);
        cdd56 = Q5t*cdd45
                + cross(wd56, av[4])
                + cross(w56, cross(w56, av[4]))
                - Q5t*cross(wd45, sv[3])
                - Q5t*cross(w45, cross(w45, sv[3]));

        wd67 = Q6t*wd56 + a[5]*ev[5] + v[5]*cross(Q6t*w56, ev[5]);
        cdd67 = Q6t*cdd56
                + cross(wd67, av[5])
                + cross(w67, cross(w67, av[5]))
                - Q6t*cross(wd56, sv[4])
                - Q6t*cross(w56, cross(w56, sv[4]));

        //-- dynamics
        f6 = Q6*(m[5]*cdd67);
        n6 = Q6*(I[5]*wd67 + cross(w67, I[5]*w67) + cross(av[5], Q6t*f6));

        f5 = Q5*(m[4]*cdd56 + f6);
        n5 = Q5*(I[4]*wd56 + cross(w56, I[4]*w56) + n6 + cross(av[4], Q5t*f5) - cross(sv[4], f6));

        f4 = Q4*(m[3]*cdd45 + f5);
        n4 = Q4*(I[3]*wd45 + cross(w45, I[3]*w45) + n5 + cross(av[3], Q4t*f4) - cross(sv[3], f5));

        f3 = Q3*(m[2]*cdd34 + f4);
        n3 = Q3*(I[2]*wd34 + cross(w34, I[2]*w34) + n4 + cross(av[2], Q3t*f3) - cross(sv[2], f4));

        f2 = Q2*(m[1]*cdd23 + f3);
        n2 = Q2*(I[1]*wd23 + cross(w23, I[1]*w23) + n3 + cross(av[1], Q2t*f2) - cross(sv[1], f3));

        f1 = Q1*(m[0]*cdd12 + f2);
        n1 = Q1*(I[0]*wd12 + cross(w12, I[0]*w12) + n2 + cross(av[0], Q1t*f1) - cross(sv[0], f2));

        //-- extract torques (last element of each moment vector)
        efforts(0, i) = n1.z;
        efforts(1, i) = n2.z;
        efforts(2, i) = n3.z;
        efforts(3, i) = n4.z;
        efforts(4, i) = n5.z;
        efforts(5, i) = n6.z;
    }
}

inline Array Gen3Lite::forward_kinematics(Array& joint_position) {
    auto p = joint_position.data;


#if BLAST_SIZEOF_REAL == 8 // double precision
    real s[6];
    real c[6];
    // first 4
    {
        __m256d s_tmp;
        __m256d c_tmp;
        __m256d angle_v = _mm256_load_pd(p);
        s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
        _mm256_storeu_pd(s, s_tmp);
        _mm256_storeu_pd(c, c_tmp);
    }
    // 5 and 6
    {
        __m128d angle_v = _mm_load_pd(p+4);
        __m128d s_tmp;
        __m128d c_tmp;
        s_tmp = _mm_sincos_pd(&c_tmp, angle_v);
        _mm_storeu_pd(s+4, s_tmp);
        _mm_storeu_pd(c+4, c_tmp);
    }
#else // single precision
    real s[8];
    real c[8];
    __m256 s_tmp;
    __m256 c_tmp;
    __m256 angle_v = _mm256_load_ps(p);
    s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
    _mm256_storeu_ps(s, s_tmp);
    _mm256_storeu_ps(c, c_tmp);
#endif

    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    Q1 = {c[0],  s[0],  0,        -s[0],  c[0],  0,         0,  0,  1};
    Q2 = {c[1],  0,     s[1],     -s[1],  0,     c[1],      0, -1,  0};
    Q3 = {c[2], -s[2],  0,        -s[2], -c[2],  0,         0,  0, -1};
    Q4 = {c[3],  0,     s[3],     -s[3],  0,     c[3],      0, -1,  0};
    Q5 = {0,     s[4], -c[4],      0,     c[4],  s[4],      1,  0,  0};
    Q6 = {0,     s[5],  c[5],      0,     c[5], -s[5],     -1,  0,  0};

    auto Q_tmp = Q1;
    auto p_tmp = Q_tmp*dv[0];
    p_tmp += (Q_tmp*=Q2)*dv[1];
    p_tmp += (Q_tmp*=Q3)*dv[2];
    p_tmp += (Q_tmp*=Q4)*dv[3];
    p_tmp += (Q_tmp*=Q5)*dv[4];
    p_tmp += (Q_tmp*=Q6)*dv[5];
    p_tmp += p_base;

    Array pose(6);
    pose[0] = p_tmp.x;
    pose[1] = p_tmp.y;
    pose[2] = p_tmp.z;
    // todo: Q_tmp(0, 0) and Q_tmp(2, 2) must not be zero!!
    pose[3] = atan2(Q_tmp(1, 0), Q_tmp(0, 0));
    pose[4] = atan2(-Q_tmp(2, 0), sqrt(Q_tmp(2, 1)*Q_tmp(2, 1) + Q_tmp(2, 2)*Q_tmp(2, 2)));
    pose[5] = atan2(Q_tmp(2, 1), Q_tmp(2, 2));
    return pose;
}


//------ Kinova Gen3 7DOF manipulator functions ---------------------------------------

inline Gen3_7DOF::Gen3_7DOF() : Manipulator(7) {
    // position of the first joint with respect to the table in the center of the base
    p_base = {0, 0, 0.1564f};

    // mass
    m[0] = 1.377f;
    m[1] = 1.1636f;
    m[2] = 1.1636f;
    m[3] = 0.93f;
    m[4] = 0.678f;
    m[5] = 0.678f;
    m[6] = 0.364f + 0.921f;

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
    av[4] = { 0.000001f, -0.009432f, -0.063883f};
    av[5] = { 0.000001f, -0.045483f, -0.009650f};
    av[6] = {-0.000093f,  0.000132f, -0.022905f};
    Vec3 av_tool(0, 0, -0.06f-0.0615f);
    av[6] = (0.364f*av[6] + 0.921f*av_tool) * (1 / (0.364f+0.921f));

    // vector to next joint
    dv[0] = { 0.0f,  0.0054f, -0.1284f};
    dv[1] = { 0.0f, -0.2104f, -0.0064f};
    dv[2] = { 0.0f, -0.0064f, -0.2104f};
    dv[3] = { 0.0f, -0.2084f, -0.0064f};
    dv[4] = { 0.0f,  0.0f,    -0.1059f};
    dv[5] = { 0.0f, -0.1059f,  0.0f};
    dv[6] = { 0.0f,  0.0f,    -0.0615f-0.164f};

    // center of mass (from next joint)
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];
    sv[6] = dv[6] - av[6];

    // unit joint direction
    ev[0] = { 0,  0,  1};
    ev[1] = { 0,  0,  1};
    ev[2] = { 0,  0,  1};
    ev[3] = { 0,  0,  1};
    ev[4] = { 0,  0,  1};
    ev[5] = { 0,  0,  1};
    ev[6] = { 0,  0,  1};

    // kinematic and dynamic constraints
    pmax = {inf, inf, inf, 2.58f, inf, 2.1f, inf}; // rad
    // pmin = -pmax;
    vmax = {1.745f, 1.745f, 1.745f, 1.745f, 2.443f, 2.443f, 2.443f};  // rad/s
    // vmin = -vmax;
    tau_max = {52, 52, 52, 52, 17, 17, 17};    // Nm
    // tau_min = -tau_max;
}

inline void Gen3_7DOF::dynamics(const Pva& pva, Matrix& efforts) {
    dynamics(pva.pos, pva.vel, pva.acc, efforts);
}

inline void Gen3_7DOF::dynamics(const Matrix& pos, const Matrix& vel, const Matrix& acc, Matrix& efforts) {

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
        auto p = &pos.data[i * joints];

#if BLAST_SIZEOF_REAL == 8
        real s[8];
        real c[8];
        __m256d s_tmp;
        __m256d c_tmp;
        for (u32 j = 0; j < 8; j += 4) {
            __m256d angle_v = _mm256_load_pd(p + j);
            s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
            _mm256_storeu_pd(s+j, s_tmp);
            _mm256_storeu_pd(c+j, c_tmp);
        }
#else
        real s[8];
        real c[8];
        __m256 s_tmp;
        __m256 c_tmp;
        __m256 angle_v = _mm256_load_ps(p);
        s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
        _mm256_storeu_ps(s, s_tmp);
        _mm256_storeu_ps(c, c_tmp);
#endif

        // note: these are stored column-wise
        Q1 = {c[0], -s[0],  0,        -s[0], -c[0],   0,        0,  0, -1};
        Q2 = {c[1],   0,   s[1],      -s[1],   0,    c[1],      0, -1,  0};
        Q3 = {c[2],   0,  -s[2],      -s[2],   0,   -c[2],      0,  1,  0};
        Q4 = {c[3],   0,   s[3],      -s[3],   0,    c[3],      0, -1,  0};
        Q5 = {c[4],   0,  -s[4],      -s[4],   0,   -c[4],      0,  1,  0};
        Q6 = {c[5],   0,   s[5],      -s[5],   0,    c[5],      0, -1,  0};
        Q7 = {c[6],   0,  -s[6],      -s[6],   0,   -c[6],      0,  1,  0};
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
        w1 = v[0]*ev[0];
        w2 = Q2t*w1 + v[1]*ev[1];
        w3 = Q3t*w2 + v[2]*ev[2];
        w4 = Q4t*w3 + v[3]*ev[3];
        w5 = Q5t*w4 + v[4]*ev[4];
        w6 = Q6t*w5 + v[5]*ev[5];
        w7 = Q7t*w6 + v[6]*ev[6];

        wd1 = a[0]*ev[0];
        cdd1 = Q1t*cdd0 + cross(wd1, av[0]) + cross(w1, cross(w1, av[0]));

        wd2 = Q2t*wd1 + a[1]*ev[1] + v[1]*cross(Q2t*w1, ev[1]);
        cdd2 = Q2t*cdd1 + cross(wd2, av[1]) + cross(w2, cross(w2, av[1])) - Q2t*cross(wd1, sv[0]) - Q2t*cross(w1, cross(w1, sv[0]));

        wd3 = Q3t*wd2 + a[2]*ev[2] + v[2]*cross(Q3t*w2, ev[2]);
        cdd3 = Q3t*cdd2 + cross(wd3, av[2]) + cross(w3, cross(w3, av[2])) - Q3t*cross(wd2, sv[1]) - Q3t*cross(w2, cross(w2, sv[1]));

        wd4 = Q4t*wd3 + a[3]*ev[3] + v[3]*cross(Q4t*w3, ev[3]);
        cdd4 = Q4t*cdd3 + cross(wd4, av[3]) + cross(w4, cross(w4, av[3])) - Q4t*cross(wd3, sv[2]) - Q4t*cross(w3, cross(w3, sv[2]));

        wd5 = Q5t*wd4 + a[4]*ev[4] + v[4]*cross(Q5t*w4, ev[4]);
        cdd5 = Q5t*cdd4 + cross(wd5, av[4]) + cross(w5, cross(w5, av[4])) - Q5t*cross(wd4, sv[3]) - Q5t*cross(w4, cross(w4, sv[3]));

        wd6 = Q6t*wd5 + a[5]*ev[5] + v[5]*cross(Q6t*w5, ev[5]);
        cdd6 = Q6t*cdd5 + cross(wd6, av[5]) + cross(w6, cross(w6, av[5])) - Q6t*cross(wd5, sv[4]) - Q6t*cross(w5, cross(w5, sv[4]));

        wd7 = Q7t*wd6 + a[6]*ev[6] + v[6]*cross(Q7t*w6, ev[6]);
        cdd7 = Q7t*cdd6 + cross(wd7, av[6]) + cross(w7, cross(w7, av[6])) - Q7t*cross(wd6, sv[5]) - Q7t*cross(w6, cross(w6, sv[5]));

        //-- dynamics
        f7 = m[6]*cdd7;
        n7 = I[6]*wd7 + cross(w7, I[6]*w7) + cross(av[6], f7);

        f6 = m[5]*cdd6 + Q7*f7;
        n6 = I[5]*wd6 + cross(w6, I[5]*w6) + Q7*n7 + cross(av[5], f6) + cross(sv[5], (Q7*f7));

        f5 = m[4]*cdd5 + Q6*f6;
        n5 = I[4]*wd5 + cross(w5, I[4]*w5) + Q6*n6 + cross(av[4], f5) + cross(sv[4], (Q6*f6));

        f4 = m[3]*cdd4 + Q5*f5;
        n4 = I[3]*wd4 + cross(w4, I[3]*w4) + Q5*n5 + cross(av[3], f4) + cross(sv[3], (Q5*f5));

        f3 = m[2]*cdd3 + Q4*f4;
        n3 = I[2]*wd3 + cross(w3, I[2]*w3) + Q4*n4 + cross(av[2], f3) + cross(sv[2], (Q4*f4));

        f2 = m[1]*cdd2 + Q3*f3;
        n2 = I[1]*wd2 + cross(w2, I[1]*w2) + Q3*n3 + cross(av[1], f2) + cross(sv[1], (Q3*f3));

        f1 = m[0]*cdd1 + Q2*f2;
        n1 = I[0]*wd1 + cross(w1, I[0]*w1) + Q2*n2 + cross(av[0], f1) + cross(sv[0], (Q2*f2));

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

inline Array Gen3_7DOF::forward_kinematics(const Array& joint_position) {

    //  - note: manual SIMD (10% better performance than using sincos function on arrays like commented below)
    real s[8];
    real c[8];
    auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;

#if BLAST_SIZEOF_REAL == 8
    __m256d s_tmp;
    __m256d c_tmp;
    for (u32 i = 0; i < 8; i += 4) {
        __m256d angle_v = _mm256_load_pd(p + i);
        s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
        _mm256_storeu_pd(s+i, s_tmp);
        _mm256_storeu_pd(c+i, c_tmp);
    }
#else
    __m256 s_tmp;
    __m256 c_tmp;
    __m256 angle_v = _mm256_load_ps(p);
    s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
    _mm256_storeu_ps(s, s_tmp);
    _mm256_storeu_ps(c, c_tmp);
#endif

    // create aliases instead of actual arrays to allocate the memory on the stack (1.7x performance)
    // real s_data[8]; // nearest factor of 4
    // real c_data[8]; // nearest factor of 4
    // Array s(s_data, 8);
    // Array c(c_data, 8);
    // Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    // sincos(joint_position, s, c);

    // note: these are stored column-wise
    Q1 = {c[0], -s[0],  0,        -s[0], -c[0],   0,        0,  0, -1};
    Q2 = {c[1],   0,   s[1],      -s[1],   0,    c[1],      0, -1,  0};
    Q3 = {c[2],   0,  -s[2],      -s[2],   0,   -c[2],      0,  1,  0};
    Q4 = {c[3],   0,   s[3],      -s[3],   0,    c[3],      0, -1,  0};
    Q5 = {c[4],   0,  -s[4],      -s[4],   0,   -c[4],      0,  1,  0};
    Q6 = {c[5],   0,   s[5],      -s[5],   0,    c[5],      0, -1,  0};
    Q7 = {c[6],   0,  -s[6],      -s[6],   0,   -c[6],      0,  1,  0};

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp*dv[0];
    p_tmp += (Q_tmp*=Q2)*dv[1];
    p_tmp += (Q_tmp*=Q3)*dv[2];
    p_tmp += (Q_tmp*=Q4)*dv[3];
    p_tmp += (Q_tmp*=Q5)*dv[4];
    p_tmp += (Q_tmp*=Q6)*dv[5];
    p_tmp += (Q_tmp*=Q7)*dv[6];

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

inline Matrix Gen3_7DOF::forward_kinematics(const Matrix& joint_positions) {
    auto p = joint_positions.data;
    Matrix pose(12, joint_positions.cols);

    for (u32 point = 0; point < joint_positions.cols; point++) {

#if BLAST_SIZEOF_REAL == 8
        real s[8];
        real c[8];
        __m256d s_tmp;
        __m256d c_tmp;
        for (u32 i = 0; i < 8; i += 4) {
            __m256d angle_v = _mm256_load_pd(p + i);
            s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
            _mm256_storeu_pd(s+i, s_tmp);
            _mm256_storeu_pd(c+i, c_tmp);
        }
#else
        real s[8];
        real c[8];
        __m256 s_tmp;
        __m256 c_tmp;
        __m256 angle_v = _mm256_load_ps(p);
        s_tmp = _mm256_sincos_ps(&c_tmp, angle_v);
        _mm256_storeu_ps(s, s_tmp);
        _mm256_storeu_ps(c, c_tmp);
#endif


#if 1
        Mat3 Q1 = {c[0], -s[0],  0,        -s[0], -c[0],   0,        0,  0, -1};
        Mat3 Q2 = {c[1],   0,   s[1],      -s[1],   0,    c[1],      0, -1,  0};
        Mat3 Q3 = {c[2],   0,  -s[2],      -s[2],   0,   -c[2],      0,  1,  0};
        Mat3 Q4 = {c[3],   0,   s[3],      -s[3],   0,    c[3],      0, -1,  0};
        Mat3 Q5 = {c[4],   0,  -s[4],      -s[4],   0,   -c[4],      0,  1,  0};
        Mat3 Q6 = {c[5],   0,   s[5],      -s[5],   0,    c[5],      0, -1,  0};
        Mat3 Q7 = {c[6],   0,  -s[6],      -s[6],   0,   -c[6],      0,  1,  0};
        auto p_tmp = p_base;
        auto Q_tmp = Q1;
        p_tmp += Q_tmp*dv[0];
        p_tmp += (Q_tmp*=Q2)*dv[1];
        p_tmp += (Q_tmp*=Q3)*dv[2];
        p_tmp += (Q_tmp*=Q4)*dv[3];
        p_tmp += (Q_tmp*=Q5)*dv[4];
        p_tmp += (Q_tmp*=Q6)*dv[5];
        p_tmp += (Q_tmp*=Q7)*dv[6];
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
        const Mat4 T1 = {c[0], -s[0],  0,   0, -s[0], -c[0],   0,   0,   0,  0, -1,  0,  p1.x, p1.y, p1.z, 1};
        const Mat4 T2 = {c[1],   0,   s[1], 0, -s[1],   0,   c[1],  0,   0, -1,  0,  0,  dv[1].x, dv[1].y, dv[1].z, 1};
        const Mat4 T3 = {c[2],   0,  -s[2], 0, -s[2],   0,  -c[2],  0,   0,  1,  0,  0,  dv[2].x, dv[2].y, dv[2].z, 1};
        const Mat4 T4 = {c[3],   0,   s[3], 0, -s[3],   0,   c[3],  0,   0, -1,  0,  0,  dv[3].x, dv[3].y, dv[3].z, 1};
        const Mat4 T5 = {c[4],   0,  -s[4], 0, -s[4],   0,  -c[4],  0,   0,  1,  0,  0,  dv[4].x, dv[4].y, dv[4].z, 1};
        const Mat4 T6 = {c[5],   0,   s[5], 0, -s[5],   0,   c[5],  0,   0, -1,  0,  0,  dv[5].x, dv[5].y, dv[5].z, 1};
        const Mat4 T7 = {c[6],   0,  -s[6], 0, -s[6],   0,  -c[6],  0,   0,  1,  0,  0,  dv[6].x, dv[6].y, dv[6].z, 1};
        const auto T = T1*T2*T3*T4*T5*T6*T7;
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

inline Matrix Gen3_7DOF::jacobian_matrix(const Array& joint_position) {

    //  - note: manual SIMD (10% better performance than using sincos function on arrays like commented below)
    real s[8];
    real c[8];
    auto p = joint_position.data;
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;

    __m256d s_tmp;
    __m256d c_tmp;
    for (u32 i = 0; i < 8; i += 4) {
        __m256d angle_v = _mm256_load_pd(p + i);
        s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
        _mm256_storeu_pd(s+i, s_tmp);
        _mm256_storeu_pd(c+i, c_tmp);
    }

    // note: these are stored column-wise
    Q1 = {c[0], -s[0],  0,        -s[0], -c[0],   0,        0,  0, -1};
    Q2 = {c[1],   0,   s[1],      -s[1],   0,    c[1],      0, -1,  0};
    Q3 = {c[2],   0,  -s[2],      -s[2],   0,   -c[2],      0,  1,  0};
    Q4 = {c[3],   0,   s[3],      -s[3],   0,    c[3],      0, -1,  0};
    Q5 = {c[4],   0,  -s[4],      -s[4],   0,   -c[4],      0,  1,  0};
    Q6 = {c[5],   0,   s[5],      -s[5],   0,    c[5],      0, -1,  0};
    Q7 = {c[6],   0,  -s[6],      -s[6],   0,   -c[6],      0,  1,  0};

    // unit vectors in 1st reference
    Vec3 e_1[7];
    e_1[0] = ev[0];
    auto Q_tmp = Q2;
    e_1[1] = Q_tmp*ev[1];
    e_1[2] = (Q_tmp*=Q3)*ev[2];
    e_1[3] = (Q_tmp*=Q4)*ev[3];
    e_1[4] = (Q_tmp*=Q5)*ev[4];
    e_1[5] = (Q_tmp*=Q6)*ev[5];
    e_1[6] = (Q_tmp*=Q7)*ev[6];

    Vec3 r[7];
    r[6] = dv[6];
    r[5] = dv[5] + Q7*r[6];
    r[4] = dv[4] + Q6*r[5];
    r[3] = dv[3] + Q5*r[4];
    r[2] = dv[2] + Q4*r[3];
    r[1] = dv[1] + Q3*r[2];
    r[0] = dv[0] + Q2*r[1];

    Q_tmp = Q2;
    r[1] = (Q_tmp)*r[1];
    r[2] = (Q_tmp*=Q3)*r[2];
    r[3] = (Q_tmp*=Q4)*r[3];
    r[4] = (Q_tmp*=Q5)*r[4];
    r[5] = (Q_tmp*=Q6)*r[5];
    r[6] = (Q_tmp*=Q7)*r[6];


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

inline Array Gen3_7DOF::collision_dist_sqr(const Array& joint_position) {

    // create aliases instead of actual arrays to allocate the memory on the stack (performance)
    real s_data[8]; // nearest factor of 4
    real c_data[8]; // nearest factor of 4
    Array s(s_data, 8);
    Array c(c_data, 8);
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    sincos(joint_position, s, c);

    // note: these are stored column-wise
    Q1 = {c[0], -s[0],  0,        -s[0], -c[0],   0,        0,  0, -1};
    Q2 = {c[1],   0,   s[1],      -s[1],   0,    c[1],      0, -1,  0};
    Q3 = {c[2],   0,  -s[2],      -s[2],   0,   -c[2],      0,  1,  0};
    Q4 = {c[3],   0,   s[3],      -s[3],   0,    c[3],      0, -1,  0};
    Q5 = {c[4],   0,  -s[4],      -s[4],   0,   -c[4],      0,  1,  0};
    Q6 = {c[5],   0,   s[5],      -s[5],   0,    c[5],      0, -1,  0};
    Q7 = {c[6],   0,  -s[6],      -s[6],   0,   -c[6],      0,  1,  0};

    Vec3 p_orig(0, 0, 0);
    Vec3 p_j2;
    Vec3 p_j3;
    Vec3 p_j4;
    Vec3 p_j6;
    Vec3 p_ee;

    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp*dv[0];
    p_j2 = p_tmp;
    p_tmp += (Q_tmp*=Q2)*dv[1];
    p_j3 = p_tmp;
    p_tmp += (Q_tmp*=Q3)*dv[2];
    p_j4 = p_tmp;
    p_tmp += (Q_tmp*=Q4)*dv[3];
    p_tmp += (Q_tmp*=Q5)*dv[4];
    p_j6 = p_tmp;
    p_tmp += (Q_tmp*=Q6)*dv[5];
    p_tmp += (Q_tmp*=Q7)*dv[6];
    p_ee = p_tmp;

    const real r1_sqr = 0.09f*0.09f;
    const real r2_sqr = 0.09f*0.09f;

    // Self collisions sqr
    real dist1sqr = two_segment_distance_sqr(p_orig, p_j2, p_j6, p_ee) - r1_sqr;
    real dist2sqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r2_sqr;

    // Collision with table sqr
    const real r4table_sqr = 0.05f*0.05f;
    const real r6table_sqr = 0.04f*0.04f;

    const Vec3 p_table(0, 0, -0.0025f); // todo: Correct coords (z or y) ??
    real distTJ4sqr = p_j4.z - p_table.z;
    distTJ4sqr *= distTJ4sqr;
    real distTJ6sqr = p_j6.z - p_table.z;
    distTJ6sqr *= distTJ6sqr;
    real distTEEsqr = p_ee.z - p_table.z;
    distTEEsqr *= distTEEsqr;

    // Array of distance min sqr and distance from table sqr
    Array distSqrMin(5);
    distSqrMin = {dist1sqr, dist2sqr, distTJ4sqr-r4table_sqr, distTJ6sqr-r6table_sqr, distTEEsqr-r6table_sqr};

    return distSqrMin;
}

inline Array Gen3_7DOF::validate(const Array& pos, const Array& vel, const Array& acc) {
    Matrix p(pos);
    Matrix v(vel);
    Matrix a(acc);

    // 5 collision results
    // 2 position constraints
    // 7 velocity constraints
    // 7 torque constraints
    Array result(5 + 2 + 7*2);

    // distance to collision >= 0
    auto tmp_coll = collision_dist_sqr(pos);
    result[0] = tmp_coll[0]; // dist1sqr
    result[1] = tmp_coll[1]; // dist2sqr
    result[2] = tmp_coll[2]; // distTJ4sqr - r1_sqr
    result[3] = tmp_coll[3]; // distTJ6sqr - r1_sqr
    result[4] = tmp_coll[4]; // distTEEsqr - r1_sqr

    Matrix efforts(joints, 1);
    dynamics(p, v, a, efforts);

    // position constraints for joints 4 and 6
    result[5] = pmax[3] - abs(pos[3]);
    result[6] = pmax[5] - abs(pos[5]);

    // velocity constraints for all joints
    auto current_result = &result[7];
    for (u32 i = 0; i < 7; i++)
        current_result[i] = vmax[i] - abs(vel[i]);

    // torque constraints for all joints
    current_result += 7;
    for (u32 i = 0; i < 7; i++)
        current_result[i] = tau_max[i] - abs(efforts(i, 0));
    return result;
}

inline Array Gen3_7DOF::validate(const Matrix& pos, const Matrix& vel, const Matrix& acc) {
    const auto points = pos.cols;
    // 5 collision results
    // 2 position constraints
    // 7 velocity constraints
    // 7 torque constraints
    //** (for each point in the trajectory) **
    Array result((5 + 2 + 7*2)*points);

    auto current_result = result.data;
    Matrix efforts(joints, pos.cols); // todo: perf hit by constructing every time?
    dynamics(pos, vel, acc, efforts);

    for (u32 i = 0; i < points; i++) {

        Array p(&pos.data[i * joints], joints);
        Assert(p.is_alias);

        auto tmp_coll = collision_dist_sqr(p);
        current_result[0] = tmp_coll[0]; // dist1sqr
        current_result[1] = tmp_coll[1]; // dist2sqr
        current_result[2] = tmp_coll[2]; // distTJ4sqr - r1_sqr
        current_result[3] = tmp_coll[3]; // distTJ6sqr - r1_sqr
        current_result[4] = tmp_coll[4]; // distTEEsqr - r1_sqr
        current_result += 5;

        current_result[0] = (pmax[3] - abs(pos(3, i))) / pmax[3];
        current_result[1] = (pmax[5] - abs(pos(5, i))) / pmax[5];
        current_result += 2;

        for (int j = 0; j < 7; j++) { // todo: check performance impact of loop
            current_result[j] = (vmax[j] - abs(vel(j, i))) / vmax[j];
        }
        current_result += 7;

        auto f = efforts.col(i); // note: alias
        for (int j = 0; j < 7; j++) { // todo: check performance impact of loop
            current_result[j] = (tau_max[j] - abs(f[j])) / tau_max[j];
        }
        current_result += 7;
    }

    return result;
}

inline Array Gen3_7DOF::validate(const Pva& pva) {
    return validate(pva.pos, pva.vel, pva.acc);
}






//------ FOR GPU COMPUTATION ONLY ------------------------------------------------------------------------------------
#ifdef __NVCC__
host_fn void cuGen3_7DOF::init(real mass, u32 npoints) {
    // position of the first joint with respect to the table in the center of the base
    p_base = {0, 0, 0.1564f};

    // mass
    m[0] = 1.377f;
    m[1] = 1.1636f;
    m[2] = 1.1636f;
    m[3] = 0.93f;
    m[4] = 0.678f;
    m[5] = 0.678f;
    m[6] = 0.364f + 0.921f;

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
    av[4] = { 0.000001f, -0.009432f, -0.063883f};
    av[5] = { 0.000001f, -0.045483f, -0.009650f};
    av[6] = {-0.000093f,  0.000132f, -0.022905f};
    Vec3 av_tool(0, 0, -0.06f-0.0615f);
    av[6] = (0.364f*av[6] + 0.921f*av_tool) * (1 / (0.364f+0.921f));

    // vector to next joint
    dv[0] = { 0.0f,  0.0054f, -0.1284f};
    dv[1] = { 0.0f, -0.2104f, -0.0064f};
    dv[2] = { 0.0f, -0.0064f, -0.2104f};
    dv[3] = { 0.0f, -0.2084f, -0.0064f};
    dv[4] = { 0.0f,  0.0f,    -0.1059f};
    dv[5] = { 0.0f, -0.1059f,  0.0f};
    dv[6] = { 0.0f,  0.0f,    -0.0615f-0.164f};

    // center of mass (from next joint)
    sv[0] = dv[0] - av[0];
    sv[1] = dv[1] - av[1];
    sv[2] = dv[2] - av[2];
    sv[3] = dv[3] - av[3];
    sv[4] = dv[4] - av[4];
    sv[5] = dv[5] - av[5];
    sv[6] = dv[6] - av[6];

    // unit joint direction
    ev[0] = { 0,  0,  1};
    ev[1] = { 0,  0,  1};
    ev[2] = { 0,  0,  1};
    ev[3] = { 0,  0,  1};
    ev[4] = { 0,  0,  1};
    ev[5] = { 0,  0,  1};
    ev[6] = { 0,  0,  1};

    // kinematic and dynamic constraints
    // position (rad)
    pmax[0] = inf;
    pmax[1] = inf;
    pmax[2] = inf;
    pmax[3] = 2.58f;
    pmax[4] = inf;
    pmax[5] = 2.1f;
    pmax[6] = inf;

    // velocity (rad/s)
    vmax[0] = 1.745f;
    vmax[1] = 1.745f;
    vmax[2] = 1.745f;
    vmax[3] = 1.745f;
    vmax[4] = 2.443f;
    vmax[5] = 2.443f;
    vmax[6] = 2.443f;

    // torque (Nm)
    tau_max[0] = 52;
    tau_max[1] = 52;
    tau_max[2] = 52;
    tau_max[3] = 52;
    tau_max[4] = 17;
    tau_max[5] = 17;
    tau_max[6] = 17;

    assert_buffer_size<cuGen3_7DOF, sizeof(manip_broadcast_arena)>(); // compile time check
    cuda_check( cudaMalloc(&device_constraints, 21*npoints*sizeof(real)) );
    cuda_check( cudaMallocHost(&host_constraints, 21*npoints*sizeof(real)) );

    // note: this must be done last
    cuda_check( cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0) );
}

host_fn void cuGen3_7DOF::fetch_constraints(u32 npoints) {
    cuda_check( cudaMemcpy(host_constraints, device_constraints, 21*npoints*sizeof(real), cudaMemcpyDeviceToHost) );
}

dev_fn void cuGen3_7DOF::compute_constraints(const real pos[7], const real vel[7], const real acc[7], real* con) {
#if BLAST_USE_DOUBLES
    real s[7];
    real c[7];
    ::sincos(pos[0], &s[0], &c[0]);
    ::sincos(pos[1], &s[1], &c[1]);
    ::sincos(pos[2], &s[2], &c[2]);
    ::sincos(pos[3], &s[3], &c[3]);
    ::sincos(pos[4], &s[4], &c[4]);
    ::sincos(pos[5], &s[5], &c[5]);
    ::sincos(pos[6], &s[6], &c[6]);
#else
    real s[7];
    real c[7];
    sincosf(pos[0], &s[0], &c[0]);
    sincosf(pos[1], &s[1], &c[1]);
    sincosf(pos[2], &s[2], &c[2]);
    sincosf(pos[3], &s[3], &c[3]);
    sincosf(pos[4], &s[4], &c[4]);
    sincosf(pos[5], &s[5], &c[5]);
    sincosf(pos[6], &s[6], &c[6]);
#endif

    const Mat3 Q1( c[0], -s[0],  0,        -s[0], -c[0],   0,        0,  0, -1 );
    const Mat3 Q2( c[1],   0,   s[1],      -s[1],   0,    c[1],      0, -1,  0 );
    const Mat3 Q3( c[2],   0,  -s[2],      -s[2],   0,   -c[2],      0,  1,  0 );
    const Mat3 Q4( c[3],   0,   s[3],      -s[3],   0,    c[3],      0, -1,  0 );
    const Mat3 Q5( c[4],   0,  -s[4],      -s[4],   0,   -c[4],      0,  1,  0 );
    const Mat3 Q6( c[5],   0,   s[5],      -s[5],   0,    c[5],      0, -1,  0 );
    const Mat3 Q7( c[6],   0,  -s[6],      -s[6],   0,   -c[6],      0,  1,  0 );
    const Mat3 Q1t( transpose(Q1) );
    const Mat3 Q2t( transpose(Q2) );
    const Mat3 Q3t( transpose(Q3) );
    const Mat3 Q4t( transpose(Q4) );
    const Mat3 Q5t( transpose(Q5) );
    const Mat3 Q6t( transpose(Q6) );
    const Mat3 Q7t( transpose(Q7) );

    //-- Collision constraints
    Vec3 p_orig(0, 0, 0);
    Vec3 p_j2;
    Vec3 p_j3;
    Vec3 p_j4;
    Vec3 p_j6;
    Vec3 p_ee;
    auto p_tmp = p_base;
    auto Q_tmp = Q1;
    p_tmp += Q_tmp*dv[0];
    p_j2 = p_tmp;
    p_tmp += (Q_tmp*=Q2)*dv[1];
    p_j3 = p_tmp;
    p_tmp += (Q_tmp*=Q3)*dv[2];
    p_j4 = p_tmp;
    p_tmp += (Q_tmp*=Q4)*dv[3];
    p_tmp += (Q_tmp*=Q5)*dv[4];
    p_j6 = p_tmp;
    p_tmp += (Q_tmp*=Q6)*dv[5];
    p_tmp += (Q_tmp*=Q7)*dv[6];
    p_ee = p_tmp;
    const real r1_sqr = 0.09f*0.09f;
    const real r2_sqr = 0.09f*0.09f;
    // Self collisions sqr
    real dist1sqr = two_segment_distance_sqr(p_orig, p_j2, p_j6, p_ee) - r1_sqr;
    real dist2sqr = two_segment_distance_sqr(p_j2, p_j3, p_j6, p_ee) - r2_sqr;
    // Collision with table sqr
    const real r4table_sqr = 0.05f*0.05f;
    const real r6table_sqr = 0.04f*0.04f;
    const Vec3 p_table(0, 0, -0.0025f);
    real distTJ4sqr = p_j4.z - p_table.z;
    distTJ4sqr *= distTJ4sqr;
    real distTJ6sqr = p_j6.z - p_table.z;
    distTJ6sqr *= distTJ6sqr;
    real distTEEsqr = p_ee.z - p_table.z;
    distTEEsqr *= distTEEsqr;
    con[0] = dist1sqr;
    con[1] = dist2sqr;
    con[2] = distTJ4sqr - r4table_sqr;
    con[3] = distTJ6sqr - r6table_sqr;
    con[4] = distTEEsqr - r6table_sqr;


    //-- position constraints
    con[5] = (pmax[3] - abs(pos[3])) / pmax[3];
    con[6] = (pmax[5] - abs(pos[5])) / pmax[5];


    //-- velocity constraints
    auto current_result = &con[7];
    for (u32 i = 0; i < 7; i++)
        current_result[i] = (vmax[i] - abs(vel[i])) / vmax[i];
    current_result += 7;


    //-- Dynamic constraints
    // note: This is the Newton algorithm in 'Element de robotique' course notes.
    //       Careful because some variables are named differently and uses a slightly different conventions.
    //       For example, the ith coordinate frame turns with the ith joint, where in the course notes, the
    //       joint turns with respect to the coordinate frame.
    Vec3 w1, w2, w3, w4, w5, w6, w7;
    Vec3 wd1, wd2, wd3, wd4, wd5, wd6, wd7;
    Vec3 cdd0 = {0, 0, 9.81f};
    Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6, cdd7;
    Vec3 f1, f2, f3, f4, f5, f6, f7;
    Vec3 n1, n2, n3, n4, n5, n6, n7;

    w1 = vel[0]*ev[0];
    w2 = Q2t*w1 + vel[1]*ev[1];
    w3 = Q3t*w2 + vel[2]*ev[2];
    w4 = Q4t*w3 + vel[3]*ev[3];
    w5 = Q5t*w4 + vel[4]*ev[4];
    w6 = Q6t*w5 + vel[5]*ev[5];
    w7 = Q7t*w6 + vel[6]*ev[6];

    wd1 = acc[0]*ev[0];
    cdd1 = Q1t*cdd0 + cross(wd1, av[0]) + cross(w1, cross(w1, av[0]));
    wd2 = Q2t*wd1 + acc[1]*ev[1] + vel[1]*cross(Q2t*w1, ev[1]);
    cdd2 = Q2t*cdd1 + cross(wd2, av[1]) + cross(w2, cross(w2, av[1])) - Q2t*cross(wd1, sv[0]) - Q2t*cross(w1, cross(w1, sv[0]));
    wd3 = Q3t*wd2 + acc[2]*ev[2] + vel[2]*cross(Q3t*w2, ev[2]);
    cdd3 = Q3t*cdd2 + cross(wd3, av[2]) + cross(w3, cross(w3, av[2])) - Q3t*cross(wd2, sv[1]) - Q3t*cross(w2, cross(w2, sv[1]));
    wd4 = Q4t*wd3 + acc[3]*ev[3] + vel[3]*cross(Q4t*w3, ev[3]);
    cdd4 = Q4t*cdd3 + cross(wd4, av[3]) + cross(w4, cross(w4, av[3])) - Q4t*cross(wd3, sv[2]) - Q4t*cross(w3, cross(w3, sv[2]));
    wd5 = Q5t*wd4 + acc[4]*ev[4] + vel[4]*cross(Q5t*w4, ev[4]);
    cdd5 = Q5t*cdd4 + cross(wd5, av[4]) + cross(w5, cross(w5, av[4])) - Q5t*cross(wd4, sv[3]) - Q5t*cross(w4, cross(w4, sv[3]));
    wd6 = Q6t*wd5 + acc[5]*ev[5] + vel[5]*cross(Q6t*w5, ev[5]);
    cdd6 = Q6t*cdd5 + cross(wd6, av[5]) + cross(w6, cross(w6, av[5])) - Q6t*cross(wd5, sv[4]) - Q6t*cross(w5, cross(w5, sv[4]));
    wd7 = Q7t*wd6 + acc[6]*ev[6] + vel[6]*cross(Q7t*w6, ev[6]);
    cdd7 = Q7t*cdd6 + cross(wd7, av[6]) + cross(w7, cross(w7, av[6])) - Q7t*cross(wd6, sv[5]) - Q7t*cross(w6, cross(w6, sv[5]));

    //-- dynamics
    f7 = m[6]*cdd7;
    n7 = I[6]*wd7 + cross(w7, I[6]*w7) + cross(av[6], f7);
    f6 = m[5]*cdd6 + Q7*f7;
    n6 = I[5]*wd6 + cross(w6, I[5]*w6) + Q7*n7 + cross(av[5], f6) + cross(sv[5], (Q7*f7));
    f5 = m[4]*cdd5 + Q6*f6;
    n5 = I[4]*wd5 + cross(w5, I[4]*w5) + Q6*n6 + cross(av[4], f5) + cross(sv[4], (Q6*f6));
    f4 = m[3]*cdd4 + Q5*f5;
    n4 = I[3]*wd4 + cross(w4, I[3]*w4) + Q5*n5 + cross(av[3], f4) + cross(sv[3], (Q5*f5));
    f3 = m[2]*cdd3 + Q4*f4;
    n3 = I[2]*wd3 + cross(w3, I[2]*w3) + Q4*n4 + cross(av[2], f3) + cross(sv[2], (Q4*f4));
    f2 = m[1]*cdd2 + Q3*f3;
    n2 = I[1]*wd2 + cross(w2, I[1]*w2) + Q3*n3 + cross(av[1], f2) + cross(sv[1], (Q3*f3));
    f1 = m[0]*cdd1 + Q2*f2;
    n1 = I[0]*wd1 + cross(w1, I[0]*w1) + Q2*n2 + cross(av[0], f1) + cross(sv[0], (Q2*f2));

    //-- extract torques (last element of each moment vector)
    current_result[0] = (tau_max[0] - abs(n1.z)) / tau_max[0];;
    current_result[1] = (tau_max[1] - abs(n2.z)) / tau_max[1];;
    current_result[2] = (tau_max[2] - abs(n3.z)) / tau_max[2];;
    current_result[3] = (tau_max[3] - abs(n4.z)) / tau_max[3];;
    current_result[4] = (tau_max[4] - abs(n5.z)) / tau_max[4];;
    current_result[5] = (tau_max[5] - abs(n6.z)) / tau_max[5];;
    current_result[6] = (tau_max[6] - abs(n7.z)) / tau_max[6];;
}
#endif

} // namespace blast
