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

}