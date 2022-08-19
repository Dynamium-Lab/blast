#pragma once

#include "cuda/blast_cuda.cuh"
#include "blast_math.hpp"

namespace blast {

// note: must be a global variable because it is __constant__
// note: much faster if it's __constant__ because all threads access the same location (heavy use of broadcast operations)
const int MANIP_BROADCAST_ARENA_SIZE = 256;
__constant__ double manip_broadcast_arena[MANIP_BROADCAST_ARENA_SIZE];

struct cuGen3_7DOF {
    //note: for the mapping to work, all the data must be on the stack!!
    double      pmax[7];
    double      vmax[7];
    double      tau_max[7];
    Vec3        p_base;
    double      m[7];
    Mat3        I[7];
    Vec3        av[7];
    Vec3        dv[7];
    Vec3        sv[7];
    Vec3        ev[7];

    // compute the parameters according to the mass
    __host__ void init(double mass);

    // map to GPU using cudaMemcpyToSymbol
    __host__ void map_to_gpu();

    // compute the dynamics for the given point
    __device__ void dynamics(const double pos[7], const double vel[7], const double acc[7], double tor[7]);

    // compute the constraints based on the current trajectory point
    __device__ void validate(const double pos[7], const double vel[7], const double acc[7], double constraints[]);
};





//------ Kinova Gen3 7DOF manipulator functions ---------------------------------------

inline __host__ void cuGen3_7DOF::init(double mass) {
    // position of the first joint with respect to the table in the center of the base
    p_base = {0, 0, 0.1564};

    // mass
    m[0] = 1.377;
    m[1] = 1.1636;
    m[2] = 1.1636;
    m[3] = 0.93;
    m[4] = 0.678;
    m[5] = 0.678;
    m[6] = 0.364 + 0.921;

    // inertial tensor
    I[0] = {0.004570, 0.000001, 0.000002, 0.000001, 0.004831, 0.000448, 0.000002, 0.000448, 0.001409};
    I[1] = {0.011088, 0.000005, 0.000000, 0.000005, 0.001072, 0.000691, 0.000000, 0.000691, 0.011255};
    I[2] = {0.010932, 0.000000, 0.000007, 0.000000, 0.011127, 0.000606, 0.000007, 0.000606, 0.001043};
    I[3] = {0.008147, 0.000001, 0.000000, 0.000001, 0.000631, 0.000500, 0.000000, 0.000500, 0.008316};
    I[4] = {0.001596, 0.000000, 0.000000, 0.000000, 0.001607, 0.000256, 0.000000, 0.000256, 0.000399};
    I[5] = {0.001641, 0.000000, 0.000000, 0.000000, 0.000410, 0.000278, 0.000000, 0.000278, 0.001641};
    I[6] = {0.000214, 0.000000, 0.000001, 0.000000, 0.000223, 0.000002, 0.000001, 0.000002, 0.000240};

    // center of mass
    av[0] = {-0.000023, -0.010364, -0.073360};
    av[1] = {-0.000044, -0.099580, -0.013278};
    av[2] = {-0.000044, -0.006641, -0.117892};
    av[3] = {-0.000018, -0.075478, -0.015006};
    av[4] = { 0.000001, -0.009432, -0.063883};
    av[5] = { 0.000001, -0.045483, -0.009650};
    av[6] = {-0.000093,  0.000132, -0.022905};
    Vec3 av_tool(0, 0, -0.06-0.0615);
    av[6] = (0.364*av[6] + 0.921*av_tool) * (1 / (0.364+0.921));

    // vector to next joint
    dv[0] = { 0.0,  0.0054, -0.1284};
    dv[1] = { 0.0, -0.2104, -0.0064};
    dv[2] = { 0.0, -0.0064, -0.2104};
    dv[3] = { 0.0, -0.2084, -0.0064};
    dv[4] = { 0.0,  0.0,    -0.1059};
    dv[5] = { 0.0, -0.1059,  0.0};
    dv[6] = { 0.0,  0.0,    -0.0615-0.164};

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
    pmax[3] = 2.58;
    pmax[4] = inf;
    pmax[5] = 2.1;
    pmax[6] = inf;

    // velocity (rad/s)
    vmax[0] = 1.745;
    vmax[1] = 1.745;
    vmax[2] = 1.745;
    vmax[3] = 1.745;
    vmax[4] = 2.443;
    vmax[5] = 2.443;
    vmax[6] = 2.443;

    // torque (Nm)
    tau_max[0] = 52;
    tau_max[1] = 52;
    tau_max[2] = 52;
    tau_max[3] = 52;
    tau_max[4] = 17;
    tau_max[5] = 17;
    tau_max[6] = 17;
}

inline __host__ void cuGen3_7DOF::map_to_gpu() {
    assert_buffer_size<cuGen3_7DOF, sizeof(manip_broadcast_arena)>(); // compile time check
    cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0);
}

inline __device__ void cuGen3_7DOF::dynamics(const double* pos, const double* vel, const double* acc, double* tor) {
    Mat3 Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Mat3 Q1t, Q2t, Q3t, Q4t, Q5t, Q6t, Q7t;
    Vec3 w1, w2, w3, w4, w5, w6, w7;
    Vec3 wd1, wd2, wd3, wd4, wd5, wd6, wd7;
    Vec3 cdd0 = {0, 0, 9.81};
    Vec3 cdd1, cdd2, cdd3, cdd4, cdd5, cdd6, cdd7;
    Vec3 f1, f2, f3, f4, f5, f6, f7;
    Vec3 n1, n2, n3, n4, n5, n6, n7;

    float s[7];
    float c[7];
    sincosf(pos[0], &s[0], &c[0]);
    sincosf(pos[1], &s[1], &c[1]);
    sincosf(pos[2], &s[2], &c[2]);
    sincosf(pos[3], &s[3], &c[3]);
    sincosf(pos[4], &s[4], &c[4]);
    sincosf(pos[5], &s[5], &c[5]);
    sincosf(pos[6], &s[6], &c[6]);

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
    tor[0] = n1.z;
    tor[1] = n2.z;
    tor[2] = n3.z;
    tor[3] = n4.z;
    tor[4] = n5.z;
    tor[5] = n6.z;
    tor[6] = n7.z;
}

inline __device__ void cuGen3_7DOF::validate(const double* pos, const double* vel, const double* acc, double* constraints) {

}


} // namespace blast
