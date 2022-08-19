#pragma once

#include "cuda/blast_cuda.cuh"
#include "cuda/blast_cu_math.cuh"

#include "blast_math.hpp"
#include "blast_trajectory.hpp"


namespace blast {

struct cuPvaBspline {
    // data for both the host and the device
    unsigned points;
    unsigned joints;
    unsigned p;
    unsigned ncontrol;

    // data only for the device
    // note: not in __constant__ memory because each thread accesses a different slice (no broadcast benefit)
    double* device_pos;     // joints x points
    double* device_vel;     // joints x points
    double* device_acc;     // joints x points
    double* device_t;       // points
    double* device_basis_p; // nccontrol x points
    double* device_basis_v; // nccontrol x points
    double* device_basis_a; // nccontrol x points
    double* device_control; // nccontrol x joints
    double  dt;
    double  one_over_T;
    double  one_over_T2;

    // data only for the host
    PvaBspline* host_pva;



    // compute the basis functions
    __host__ void init(u32 points, u32 joints, u32 p, u32 ncontrol);

    // free all host and device memory
    __host__ void clear();

    // fetch the latest computed trajectory on the GPU to host_pva
    __host__ void fetch_pva();

    // compute the control points based on the optimization vector and the task
    __host__ void compute_control(const Array&x, Matrix& task);

    // send the last computed control points to device memory
    __host__ void send_control();

    // compute control points on host and send to device
    __host__ void compute_control_and_send(const Array&x, Matrix& task);


    // compute the trajectory for the given point
    __device__ void compute_trajectory(unsigned point);
};



//------ HOST FUNCTIONS ------------------------------------------------------------------------------------

__host__
inline void cuPvaBspline::init(u32 _points, u32 _joints, u32 _p, u32 _ncontrol) {
    static bool init = false;
    if (init) {
        clear();
    }
    init = true;

    points   = _points;
    joints   = _joints;
    p        = _p;
    ncontrol = _ncontrol;

    host_pva = new PvaBspline(ncontrol, points, p, joints);
    host_pva->compute_basis();

    // prepare and send basis function values
    const u32 basis_bytes = points*ncontrol*sizeof(double);
    cuda_check( cudaMalloc(&device_basis_p, basis_bytes) );
    cuda_check( cudaMalloc(&device_basis_v, basis_bytes) );
    cuda_check( cudaMalloc(&device_basis_a, basis_bytes) );
    cuda_check( cudaMemcpy(device_basis_p, host_pva->basis_p.data, basis_bytes, cudaMemcpyHostToDevice) );
    cuda_check( cudaMemcpy(device_basis_p, host_pva->basis_v.data, basis_bytes, cudaMemcpyHostToDevice) );
    cuda_check( cudaMemcpy(device_basis_p, host_pva->basis_a.data, basis_bytes, cudaMemcpyHostToDevice) );

    // prepare pos, vel, and acc matrices and t array
    cuda_check( cudaMalloc(&device_pos,     points*joints*sizeof(double)) );
    cuda_check( cudaMalloc(&device_vel,     points*joints*sizeof(double)) );
    cuda_check( cudaMalloc(&device_acc,     points*joints*sizeof(double)) );
    cuda_check( cudaMalloc(&device_t,       points*sizeof(double)) );
    cuda_check( cudaMalloc(&device_control, joints*ncontrol*sizeof(double)) );
}

__host__
inline void cuPvaBspline::clear() {
    if (device_basis_p) {
        cuda_check( cudaFree(device_basis_p) );
        device_basis_p = nullptr;
    }
    if (device_basis_v) {
        cuda_check( cudaFree(device_basis_v) );
        device_basis_v = nullptr;
    }
    if (device_basis_a) {
        cuda_check( cudaFree(device_basis_a) );
        device_basis_a = nullptr;
    }
    if (device_pos) {
        cuda_check( cudaFree(device_pos) );
        device_pos = nullptr;
    }
    if (device_vel) {
        cuda_check( cudaFree(device_vel) );
        device_vel = nullptr;
    }
    if (device_acc) {
        cuda_check( cudaFree(device_acc) );
        device_acc = nullptr;
    }
    if (device_t) {
        cuda_check( cudaFree(device_t) );
        device_t = nullptr;
    }

    delete(host_pva);
    host_pva = nullptr;
}

__host__
inline void cuPvaBspline::fetch_pva() {
    cuda_check( cudaMemcpy(host_pva->pos.data, device_pos, host_pva->pos.size * sizeof(double), cudaMemcpyDeviceToHost) );
    cuda_check( cudaMemcpy(host_pva->vel.data, device_vel, host_pva->vel.size * sizeof(double), cudaMemcpyDeviceToHost) );
    cuda_check( cudaMemcpy(host_pva->acc.data, device_acc, host_pva->acc.size * sizeof(double), cudaMemcpyDeviceToHost) );
}

__host__
inline void cuPvaBspline::compute_control(const Array&x, Matrix& task) {
    Assert(host_pva);

    dt = x.back() / (double)(points-1);
    one_over_T = 1.0 / x.back();
    one_over_T2 = one_over_T * one_over_T;

    host_pva->compute_control(x, task);
}

__host__
inline void cuPvaBspline::send_control() {
    cuda_check( cudaMemcpy(device_control, host_pva->control.data, joints*ncontrol*sizeof(double), cudaMemcpyHostToDevice) );
}

__host__
inline void cuPvaBspline::compute_control_and_send(const Array&x, Matrix& task) {
    compute_control(x, task);
    send_control();
}



//------ DEVICE FUNCTIONS ------------------------------------------------------------------------------------
__device__
inline void cuPvaBspline::compute_trajectory(unsigned point) {
    auto bp = &device_basis_p[point * ncontrol];
    auto bv = &device_basis_v[point * ncontrol];
    auto ba = &device_basis_a[point * ncontrol];

    for (int joint = 0; joint < joints; joint++) {
        auto c = &device_control[ncontrol * joint];
        device_pos[point*joints + joint] = dot(c, bp, ncontrol);
        device_vel[point*joints + joint] = dot(c, bv, ncontrol) * one_over_T;
        device_acc[point*joints + joint] = dot(c, ba, ncontrol) * one_over_T2;
    }
    device_t[point] = dt * point;
}

} // namespace blast
