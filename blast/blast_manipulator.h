#pragma once

#include "blast.h"
#include <vector>

namespace blast {

struct Gen3;
struct Link6;

struct Gen3 {
    // basic manipulator properties
    int     joints = 7;

    // actuator limits
    Array   pmax; // rad
    Array   vmax; // rad/s
    Array   tau_max; // Nm

    // kinematic properties
    Vec3    p_base;
    Vec3    dv[7]; // vector to next joint
    Vec3    ev[7]; // direction vectors of joint

    // dynamic properties
    real    m[7];  // link masses
    Mat3    I[7];  // Inertial tensors
    Vec3    av[7]; // centers of mass
    Vec3    sv[7]; // centers of mass from next joint

    Matrix  _efforts; // put the efforts temporarily when computing the constraints

    // initialization
    Gen3();
    void    set_payload(const real m_payload, const Vec3 cg_payload = {}, const Mat3 I_payload = {});
    void    set_payload_without_gripper(const real m_payload, const Vec3 cg_payload = {}, const Mat3 I_payload = {});

    // optimization
    void    internal_constraints(const Trajectory& traj, real* dst);
    bool    validate_task(const Matrix &task, World *world = nullptr);
    int     ncon(int points);

    // dynamics
    void    dynamics(const Trajectory& traj); // note: results stored in _efforts
    Array   dynamics(const Array& pos, const Array& vel, const Array& acc);

    // kinematics
    Array   forward_kinematics(const Array &pos);
    Matrix  forward_kinematics(const Matrix &pos);
    Matrix  jacobian(const Array &joint_position);

    // collisions
    Array   internal_collisions(const Array &joint_position);
    Matrix  cartesian_positions(const Matrix& pos);
    std::vector<Capsule>  robot_capsules(const Matrix& pos, int n_collision_skip);
};

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
    bool    validate_task(const Matrix &task, World *world = nullptr);
    int     ncon(int points);

    // dynamics
    void    dynamics(const Trajectory& traj); // note: results stored in _efforts
    Array   dynamics(const Array& pos, const Array& vel, const Array& acc);

    // kinematics
    Array   forward_kinematics(const Array &pos);
    Matrix  forward_kinematics(const Matrix &pos);
    Matrix  jacobian(const Array &joint_position);

    // collisions
    Array   internal_collisions(const Array &joint_position);
    Matrix  robot_capsules(const Array &joint_position);
    std::vector<Capsule>  robot_capsules(const Matrix &pos, const int n_skip);
};

struct R2 {
    const int     joints = 2;

    Vec3    p_base;
    Vec3    dv[2]; // vector to next joint
    Vec3    ev[2]; // direction vectors of joint
    real    m[2];  // link masses
    Mat3    I[2];  // Inertial tensors
    Vec3    av[2]; // centers of mass
    Vec3    sv[2]; // centers of mass from next joint

    R2();
    Matrix    dynamics(const Trajectory& traj);
    void      dynamics(const Trajectory& traj, Matrix& result);
};


blast_fn inline real clamped_root(real slope, real h0, real h1) {
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


blast_fn inline void compute_intersection(real* sValue, i32* classify, real b, real f00, real f10, i32* edge, real end[][2]) {
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

blast_fn inline void compute_minimum_parameters(i32* edge, real end[][2], real b, real c, real e, real g00, real g10, real g01, real g11, real* parameter) {
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


blast_fn inline real two_segment_distance_sqr(Vec3 P0, Vec3 P1, Vec3 Q0, Vec3 Q1) {
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