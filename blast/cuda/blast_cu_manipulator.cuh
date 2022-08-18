#pragma once

#include "cuda/blast_cuda.cuh"
#include "blast_math.hpp"

namespace blast {

const int MANIP_BROADCAST_ARENA_SIZE = 256;
__constant__ double manip_broadcast_arena[MANIP_BROADCAST_ARENA_SIZE];

struct alignas(16) cuGen3_7DOF {
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
    __host__ void map_to_gpu(double* arena);

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
    pmax[0] = inf;
    pmax[1] = inf;
    pmax[2] = inf;
    pmax[3] = 2.58;
    pmax[4] = inf;
    pmax[5] = 2.1;
    pmax[6] = inf;

    vmax[0] = 1.745;
    vmax[1] = 1.745;
    vmax[2] = 1.745;
    vmax[3] = 1.745;
    vmax[4] = 2.443;
    vmax[5] = 2.443;
    vmax[6] = 2.443;

    tau_max[0] = 52;
    tau_max[1] = 52;
    tau_max[2] = 52;
    tau_max[3] = 52;
    tau_max[4] = 17;
    tau_max[5] = 17;
    tau_max[6] = 17;
}

inline __host__ void cuGen3_7DOF::map_to_gpu(double* arena) {
    assert_buffer_size<cuGen3_7DOF, sizeof(manip_broadcast_arena)>(); // compile time check

    cudaMemcpyToSymbol(manip_broadcast_arena, this, sizeof(*this), 0);
}

inline __device__ void cuGen3_7DOF::dynamics(const double* pos, const double* vel, const double* acc, double* tor) {

}

inline __device__ void cuGen3_7DOF::validate(const double* pos, const double* vel, const double* acc, double* constraints) {

}


} // namespace blast
