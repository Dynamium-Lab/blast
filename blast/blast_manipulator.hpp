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

    Manipulator(u32 njoints) : joints(njoints) {}

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



inline void ManipulatorGeneric::forward_kinematics(Matrix& joint_position, Matrix& cartesian_position) {
    Assert(joint_position.rows == joints);
    Assert(joint_position.cols == cartesian_position.cols);
    Assert(cartesian_position.rows == 6);

    // todo: Implement me!
}





//------ FUNCTIONS ------------------------------------------------------------------------------------

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

    // todo: construct these every run??
    Array s(joints); // sines
    Array c(joints); // cosines

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

        // all joints
        for (u32 j = 0; j < joints; j++) {
            // todo: SIMD??
            s[j] = sin(p[j]);
            c[j] = cos(p[j]);

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
                f[j] = Q[j]*(m[j]*cdd[j] + f[j+1]);
                n[j] = Q[j]*(I[j]*wd[j] + cross(w[j], I[j]*w[j]) + n[j+1] + cross(av[j]+sv[j], Qt[j]*f[j]) - cross(sv[j], f[j+1]));
            }
        }

        // note: assuming all revolute joints here
        // assign joint torques to output
        for (u32 j = 0; j < joints; j++) {
            efforts(j, i) = n[j].z;
        }

    } // for points

} // dynamics()



} // namespace blast
