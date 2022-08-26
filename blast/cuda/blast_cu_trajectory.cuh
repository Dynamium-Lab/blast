#pragma once

#include "blast.hpp"

#include "cuda/blast_cuda.cuh"



namespace blast {

struct cuPvaBspline {
    bool is_init = false;
    // data for both the host and the device
    unsigned points;
    unsigned joints;
    unsigned p;
    unsigned ncontrol;

    // data only for the device
    // note: not in __constant__ memory because each thread accesses a different slice (no broadcast benefit)
    real* device_pos     = nullptr; // joints x points
    real* device_vel     = nullptr; // joints x points
    real* device_acc     = nullptr; // joints x points
    real* device_t       = nullptr; // points
    real* device_basis_p = nullptr; // nccontrol x points
    real* device_basis_v = nullptr; // nccontrol x points
    real* device_basis_a = nullptr; // nccontrol x points
    real* device_control = nullptr; // nccontrol x joints todo: should this be in constant memory??
    real  dt;
    real  one_over_T;
    real  one_over_T2;

    // data only for the host
    PvaBspline* host = nullptr;



    // compute the basis functions
    host_fn void init(u32 points, u32 joints, u32 p, u32 ncontrol);

    // free all host and device memory
    host_fn void clear();

    // fetch the latest computed trajectory on the GPU to host_pva
    host_fn void fetch_pva();

    // compute the control points based on the optimization vector and the task
    host_fn void compute_control(const Array&x, const Matrix& task);

    // send the last computed control points to device memory
    host_fn void send_control();

    // compute control points on host and send to device
    host_fn void compute_control_and_send(const Array&x, const Matrix& task);


    // compute the trajectory for the given point
    dev_fn void compute_trajectory(unsigned point);
};



//------ HOST FUNCTIONS ------------------------------------------------------------------------------------

host_fn void cuPvaBspline::init(u32 _points, u32 _joints, u32 _p, u32 _ncontrol) {
    if (is_init)
        return;
    is_init = true;

    points   = _points;
    joints   = _joints;
    p        = _p;
    ncontrol = _ncontrol;

    host = new PvaBspline(ncontrol, points, p, joints);

    // prepare and send basis function values
    const u32 basis_bytes = points*ncontrol*sizeof(real);
    cuda_check( cudaMalloc(&device_basis_p, basis_bytes) );
    cuda_check( cudaMalloc(&device_basis_v, basis_bytes) );
    cuda_check( cudaMalloc(&device_basis_a, basis_bytes) );
    cuda_check( cudaMemcpy(device_basis_p, host->basis_p.data, basis_bytes, cudaMemcpyHostToDevice) );
    cuda_check( cudaMemcpy(device_basis_v, host->basis_v.data, basis_bytes, cudaMemcpyHostToDevice) );
    cuda_check( cudaMemcpy(device_basis_a, host->basis_a.data, basis_bytes, cudaMemcpyHostToDevice) );

    // prepare pos, vel, and acc matrices and t array
    cuda_check( cudaMalloc(&device_pos,     points*joints*sizeof(real)) );
    cuda_check( cudaMalloc(&device_vel,     points*joints*sizeof(real)) );
    cuda_check( cudaMalloc(&device_acc,     points*joints*sizeof(real)) );
    cuda_check( cudaMalloc(&device_t,       points*sizeof(real)) );
    cuda_check( cudaMalloc(&device_control, joints*ncontrol*sizeof(real)) );
}

host_fn void cuPvaBspline::clear() {
    printf("Clearing memory\n");
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
    delete(host);
    host = nullptr;
    is_init = false;
}

host_fn void cuPvaBspline::fetch_pva() {
    cuda_check( cudaMemcpy(host->pos.data, device_pos, host->pos.size * sizeof(real), cudaMemcpyDeviceToHost) );
    cuda_check( cudaMemcpy(host->vel.data, device_vel, host->vel.size * sizeof(real), cudaMemcpyDeviceToHost) );
    cuda_check( cudaMemcpy(host->acc.data, device_acc, host->acc.size * sizeof(real), cudaMemcpyDeviceToHost) );
    cuda_check( cudaMemcpy(host->t.data, device_t, host->t.size * sizeof(real), cudaMemcpyDeviceToHost) );
}

host_fn void cuPvaBspline::compute_control(const Array&x, const Matrix& task) {
    Assert(host);
    dt = x.back() / (real)(points-1);
    one_over_T = 1 / x.back();
    one_over_T2 = one_over_T * one_over_T;
    host->compute_control(x, task);
}

host_fn void cuPvaBspline::send_control() {
    cuda_check( cudaMemcpy(device_control, host->control.data, joints*ncontrol*sizeof(real), cudaMemcpyHostToDevice) );
}

host_fn void cuPvaBspline::compute_control_and_send(const Array&x, const Matrix& task) {
    compute_control(x, task);
    send_control();
}



//------ DEVICE FUNCTIONS ------------------------------------------------------------------------------------
dev_fn void cuPvaBspline::compute_trajectory(unsigned point) {
    device_t[point] = dt * point;
    auto bp = device_basis_p + point*ncontrol;
    auto bv = device_basis_v + point*ncontrol;
    auto ba = device_basis_a + point*ncontrol;
    for (int joint = 0; joint < joints; joint++) {
        auto c = device_control + joint*ncontrol;
        real position = 0;
        real velocity = 0;
        real acceleration = 0;
        for (int i = 0; i < ncontrol; i++) {
            position     += c[i] * bp[i];
            velocity     += c[i] * bv[i];
            acceleration += c[i] * ba[i];
        }
        device_pos[point*joints + joint] = position;
        device_vel[point*joints + joint] = velocity     * one_over_T;
        device_acc[point*joints + joint] = acceleration * one_over_T2;
    }
}

} // namespace blast
