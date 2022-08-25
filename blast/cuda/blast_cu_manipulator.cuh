#pragma once

#include "cuda/blast_cuda.cuh"
#include "blast_math.hpp"

namespace blast {

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
    __host__ void fetch_constraints();

    // compute the constraints based on the current trajectory point
    __device__ void compute_constraints(const real pos[7], const real vel[7], const real acc[7], real con[21]);
};





//------ Kinova Gen3 7DOF manipulator functions ---------------------------------------

inline __host__ void cuGen3_7DOF::init(real mass, u32 npoints) {
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
    cuda_check( cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0) );
    cuda_check( cudaMalloc(&device_constraints, 21*npoints*sizeof(real)) );
    cuda_check( cudaMallocHost(&host_constraints, 21*npoints*sizeof(real)) );
}

inline __device__ void cuGen3_7DOF::compute_constraints(const real pos[7], const real vel[7], const real acc[7], real con[21]) {
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
    con[5] = pmax[3] - abs(pos[3]);
    con[6] = pmax[5] - abs(pos[5]);


    //-- velocity constraints
    auto current_result = &con[7];
    for (u32 i = 0; i < 7; i++)
        current_result[i] = vmax[i] - abs(vel[i]);
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
    current_result[0] = n1.z;
    current_result[1] = n2.z;
    current_result[2] = n3.z;
    current_result[3] = n4.z;
    current_result[4] = n5.z;
    current_result[5] = n6.z;
    current_result[6] = n7.z;
}


} // namespace blast
