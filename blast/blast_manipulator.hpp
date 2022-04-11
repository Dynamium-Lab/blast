#pragma once

#include "blast_math.hpp"
#include <vector>
#include <cmath>

namespace blast {
using std::vector;


struct Manipulator {
    u32 joints;
    Array vmax;
    Array vmin;
    Array amax;
    Array amin;
    Array umax;
    Array umin;

    Manipulator(u32 njoints) : joints(njoints), vmax(njoints), vmin(njoints), amax(njoints), amin(njoints), umax(njoints), umin(njoints) {}

    virtual void dynamics(Pva& pva, Matrix& efforts) = 0;
};

struct ManipulatorGeneric : public Manipulator {
    // DH parameters
    vector<Vec3> av;    // vector from joint to next joint
    vector<Vec3> e;     // unit vectors of joints

    // inertial parameters
    vector<real> m;     // mass
    vector<Vec3> sv;    // vector from next joint to center of mass
    vector<Mat3> I;     // axis aligned inertial matrices

    ManipulatorGeneric(u32 njoints);

    virtual void dynamics(Pva& pva, Matrix& efforts) override;
    void forward_kinematics(Matrix& joint_position, Matrix& cartesian_position);
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

    virtual void dynamics(Pva& pva, Matrix& efforts) override;
    void init_dynamics(real mass=0);
};

struct Gen3Lite : public Manipulator {

    Gen3Lite();
    virtual void dynamics(Pva& pva, Matrix& efforts) override;
};




//------ Generic manipulator functions -----------------------------------------------------------------------------------

inline ManipulatorGeneric::ManipulatorGeneric(u32 njoints) :
    Manipulator(njoints),
    av(njoints),
    e(njoints),
    m(njoints),
    sv(njoints),
    I(njoints) {
}

inline void ManipulatorGeneric::dynamics(Pva& pva, Matrix& efforts) {
    Assert(pva.joints == joints);
    Assert(efforts.rows == joints);
    Assert(efforts.cols == pva.points);

    const auto joints_4 = joints%4? joints + (4-joints%4) : joints; // padding for better SIMD
    Array s(joints_4); // sines
    Array c(joints_4); // cosines

    vector<Mat3> Q(joints); // rotation matrices
    vector<Mat3> Qt(joints); // transpose rotation matrices
    vector<Vec3> w(joints); // angular velocity
    vector<Vec3> wd(joints); // angular acceleration
    vector<Vec3> cdd(joints); // cartesian acceleration of the center of mass
    vector<Vec3> f(joints); // forces on joint
    vector<Vec3> n(joints); // axis aligned torque (z is the joint torque)

    for (u32 i = 0; i < pva.points; i++) {
        auto p = &pva.pos(0, i);
        auto v = &pva.vel(0, i);
        auto a = &pva.acc(0, i);

        // SIMD compute sines and cosines
        for (u32 ji = 0; ji < joints_4/4; ji++) {
            const auto j = 4*ji;
            __m256d s_tmp;
            __m256d c_tmp;
            __m256d angle_v = _mm256_load_pd(&pva.pos(j, i));
            s_tmp = _mm256_sincos_pd(&c_tmp, angle_v);
            _mm256_storeu_pd(s.data+j, s_tmp);
            _mm256_storeu_pd(c.data+j, c_tmp);
        } // note: SIMD saved approx 10-12% compute time here

        for (u32 j = 0; j < joints; j++) {

            // todo: compute rotation matrices!!

            Qt[j] = transpose(Q[j]);
        }

        // forward kinematics
        {
            // first joint
            w[0] = v[0] * e[0];
            wd[0] = a[0] * e[0];
            cdd[0] = Qt[0] * Vec3{0, 0, 9.81f}; // simulate gravity by accelerating the base upward

            // all but first joint
            for (u32 j = 1; j < joints; j++) {
                w[j] = Qt[j]*w[j-1] + v[j]*e[j];
                wd[j] = Qt[j]*wd[j-1] + a[j]*e[j] + v[j]*cross(Qt[j]*wd[j-1], e[j]);
                cdd[j] = Qt[j]*cdd[j-1]
                         + cross(wd[j], av[j]+sv[j])
                         + cross(w[j], cross(w[j], av[j]+sv[j]))
                         - Qt[j]*cross(wd[j-1], sv[j-1])
                         - Qt[j]*cross(w[j-1], cross(w[j-1], sv[j-1]));
            }
        }

        // reverse dynamics
        {
            // last joint
            u32 idx = joints-1;
            f[idx] = Q[idx]*(m[idx]*cdd[idx]);
            n[idx] = Q[idx]*(I[idx]*wd[idx] + cross(w[idx], I[idx]*w[idx]) + cross(av[idx]+sv[idx], Qt[idx]*f[idx]));

            // all but last joint
            for (int j = (int)joints-2; j>=0; j--) {
                f[j] = Q[j]*(
                           m[j]*cdd[j]
                           + f[j+1]);
                n[j] = Q[j]*(
                           I[j]*wd[j]
                           + cross(w[j], I[j]*w[j])
                           + n[j+1]
                           + cross(av[j]+sv[j], Qt[j]*f[j])
                           - cross(sv[j], f[j+1]));
            }
        }

        // note: assuming all revolute joints here
        // assign joint torques to output
        for (u32 j = 0; j < joints; j++) {
            efforts(j, i) = n[j].z;
        }

    } // for points

} // dynamics()

inline void ManipulatorGeneric::forward_kinematics(Matrix& joint_position, Matrix& cartesian_position) {
    Assert(joint_position.rows == joints);
    Assert(joint_position.cols == cartesian_position.cols);
    Assert(cartesian_position.rows == 6);
    // todo: Implement me!
}

//------ Universal Robots UR5e manipulator functions ---------------------------------------------------------------------

inline void ManipulatorUR5::dynamics(Pva& pva, Matrix& efforts) {
    Assert(is_init);

    real vel1, vel2, vel3, vel4, vel5, vel6;
    real acc1, acc2, acc3, acc4, acc5, acc6;

    Mat3 Q1, Q2, Q3, Q4, Q5, Q6;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t;
    Vec3 w12, w23, w34, w45, w56, w67;
    Vec3 wd12, wd23, wd34, wd45, wd56, wd67;
    Vec3 cdd01, cdd12, cdd23, cdd34, cdd45, cdd56, cdd67;
    cdd01.z = 9.81;
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
        auto p = &pva.pos(0, i);
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

    base_p.z = 0.089159;

    umax = {
        1,
        1,
        1,
        1,
        1,
        1
    };
    umin = -umax;

    // kinematic parameters
    a2 = 0.425;
    a3 = 0.39225;
    d4 = -0.10915;
    d5 = 0.09465;
    d6 = 0.0823;

    // link mass
    m1 = 3.7;
    m2 = 8.393;
    m3 = 2.33;
    m4 = 1.219;
    m5 = 1.219;
    m6 = 0.1879;

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
    s1.y = -0.02561;
    s1.z = -0.00193;
    s2.x = -0.2125;
    s2.z = -0.11336;
    s3.x = -0.15;
    s3.z = -0.0265;
    s4.y = 0.0018;
    s4.z = 0.01634;
    s5.y = -0.0018;
    s5.z = 0.01634;
    s6.z = -0.001159;

    // interial tensors
    I1(0, 0) = 0.007;
    I1(1, 1) = 0.007;
    I1(2, 2) = 0.007;

    I2(0, 0) = 0.015;
    I2(1, 1) = 0.253;
    I2(2, 2) = 0.253;

    I3(0, 0) = 0.0033;
    I3(1, 1) = 0.06;
    I3(2, 2) = 0.06;

    I4(0, 0) = 0.015;
    I4(1, 1) = 0.015;
    I4(2, 2) = 0.015;

    I5(0, 0) = 0.015;
    I5(1, 1) = 0.015;
    I5(2, 2) = 0.015;

    I6(0, 0) = 0.000082;
    I6(1, 1) = 0.000082;
    I6(2, 2) = 0.000082;

    // define new center of mass of the last link (with added mass)
    if (mass) {
        Vec3 s6_old = s6;
        double m6_old = m6;
        m6 += mass;
        s6 *= m6_old / m6;

        I6(0, 0) += m6_old*(s6.x-s6_old.x)*(s6.x-s6_old.x) + mass*s6.x*s6.x;
        I6(1, 1) += m6_old*(s6.y-s6_old.y)*(s6.y-s6_old.y) + mass*s6.y*s6.y;
        I6(2, 2) += m6_old*(s6.z-s6_old.z)*(s6.z-s6_old.z) + mass*s6.z*s6.z;
    }
}

//------ Kinova Gen3 Lite manipulator functions --------------------------------------------------------------------------

inline Gen3Lite::Gen3Lite() :
    Manipulator(6) {}

inline void Gen3Lite::dynamics(Pva& pva, Matrix& efforts) {

}

} // namespace blast
