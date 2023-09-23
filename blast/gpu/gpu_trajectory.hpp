#pragma once
#include "blast.hpp"

// note: CUDA stuff, only enabled if compiling for Nvidia GPUs
#ifdef __NVCC__
namespace blast {

struct cuManagedMultiBsplines;

__global__ void compute_trajectories_kernel(cuManagedMultiBsplines bspline);

// cuda gpu computed Bsplines with managed memory
// note: this computes N trajectories in parallel
struct cuManagedMultiBsplines {
    bool is_copy = false;
    // data for both the host and the device
    unsigned points;
    unsigned joints;
    unsigned p;
    unsigned ncontrol;
    unsigned ntrajectories;

    real *pos = nullptr;     // joints x points
    real *vel = nullptr;     // joints x points
    real *acc = nullptr;     // joints x points
    real *t   = nullptr;     // points
    real *control = nullptr; // nccontrol x joints todo: should this be in constant memory??
    real *dt;
    real *one_over_T;
    real *one_over_T2;

    Bspline *host;

    real *dev_basis_p = nullptr; // nccontrol x points (only on device)
    real *dev_basis_v = nullptr; // nccontrol x points (only on device)
    real *dev_basis_a = nullptr; // nccontrol x points (only on device)

    host_fn cuManagedMultiBsplines() = delete;

    //todo: check the order of the arguments, they are not the same as Bspline class
    host_fn cuManagedMultiBsplines(unsigned points, unsigned joints, unsigned p, unsigned ncontrol, unsigned ntrajectories)
        : points(points), joints(joints), p(p), ncontrol(ncontrol), ntrajectories(ntrajectories) {
        host = new Bspline(ncontrol, points, p, joints);

        // prepare and send basis function values
        // note: this is only done once, so we don't use managed memory
        const unsigned basis_bytes = points * ncontrol * sizeof(real);
        cuda_check(cudaMalloc(&dev_basis_p, basis_bytes));
        cuda_check(cudaMalloc(&dev_basis_v, basis_bytes));
        cuda_check(cudaMalloc(&dev_basis_a, basis_bytes));
        cuda_check(cudaMemcpy(dev_basis_p, host->basis_p.data, basis_bytes, cudaMemcpyHostToDevice));
        cuda_check(cudaMemcpy(dev_basis_v, host->basis_v.data, basis_bytes, cudaMemcpyHostToDevice));
        cuda_check(cudaMemcpy(dev_basis_a, host->basis_a.data, basis_bytes, cudaMemcpyHostToDevice));
        // Managed memory for trajectory data transfers
        cuda_check(cudaMallocManaged(&dt,           ntrajectories * sizeof(real)));
        cuda_check(cudaMallocManaged(&one_over_T,   ntrajectories * sizeof(real)));
        cuda_check(cudaMallocManaged(&one_over_T2,  ntrajectories * sizeof(real)));
        cuda_check(cudaMallocManaged(&pos,      ntrajectories * points * joints *   sizeof(real)));
        cuda_check(cudaMallocManaged(&vel,      ntrajectories * points * joints *   sizeof(real)));
        cuda_check(cudaMallocManaged(&acc,      ntrajectories * points * joints *   sizeof(real)));
        cuda_check(cudaMallocManaged(&t,        ntrajectories * points *            sizeof(real)));
        cuda_check(cudaMallocManaged(&control,  ntrajectories * joints * ncontrol * sizeof(real)));
    }

    host_fn cuManagedMultiBsplines(const cuManagedMultiBsplines& other) {
        *this = other;
        is_copy = true;
    }

    // free all device memory
    host_fn ~cuManagedMultiBsplines() {
        if (!is_copy) {
            cudaFree(dev_basis_p);
            cudaFree(dev_basis_v);
            cudaFree(dev_basis_a);
            cudaFree(pos);
            cudaFree(vel);
            cudaFree(acc);
            cudaFree(t);
            cudaFree(control);
            cudaFree(dt);
            cudaFree(one_over_T);
            cudaFree(one_over_T2);
            delete host;
        }
    }

    host_fn void compute_control(const Array &x, const Matrix &task) {
        const u32 xlen = host->xlen(task);
        Assert(x.size == xlen*ntrajectories);
        for (int i = 0; i < ntrajectories; i++) {
            const auto T = x[i * xlen + xlen - 1];
            dt[i] = T / (real)(points - 1);
            one_over_T[i] = 1.0 / T;
            one_over_T2[i] = 1.0 / (T*T);
            const Array x_i(x.data + i*xlen, xlen);         Assert(x_i.is_alias);
            auto dst = control + i * joints * ncontrol;
            host->compute_control(x_i, task, dst);
        }
    }

    //note: should only be used for debugging and testing (probably not optimal to copy trajectories back and forth)
    host_fn void compute_trajectories(const Array &x, const Matrix &task) {
        compute_control(x, task);
        cudaMemPrefetchAsync(control, ntrajectories * joints * ncontrol * sizeof(real), 0, 0);

        compute_trajectories_kernel<<<ntrajectories, points>>>(*this);
        cuda_check_kernel;
        cuda_check(cudaDeviceSynchronize());
    }

    // compute the trajectory for the given point on the GPU (only use from kernel)
    dev_fn void compute_trajectory_point(unsigned point, unsigned traj_id) {
        t[traj_id*points + point] = dt[traj_id] * point;
        const auto bp = &dev_basis_p[point * ncontrol];
        const auto bv = &dev_basis_v[point * ncontrol];
        const auto ba = &dev_basis_a[point * ncontrol];
        for (int joint = 0; joint < joints; joint++) {
            const auto c = &control[joint * ncontrol];
            real position = 0;
            real velocity = 0;
            real acceleration = 0;
            for (int i = 0; i < ncontrol; i++) {
                position += c[i] * bp[i];
                velocity += c[i] * bv[i];
                acceleration += c[i] * ba[i];
            }
            const auto id = traj_id*points*joints + point*joints + joint;
            pos[id] = position;
            vel[id] = velocity * one_over_T[traj_id];
            acc[id] = acceleration * one_over_T2[traj_id];
        }
    }
};

__global__ void compute_trajectories_kernel(cuManagedMultiBsplines bspline) {
    const auto point_id = threadIdx.x;
    const auto traj_id  = blockIdx.x;
    bspline.compute_trajectory_point(point_id, traj_id);
}

struct cuBspline {
    bool is_init = false;
    // data for both the host and the device
    unsigned points;
    unsigned joints;
    unsigned p;
    unsigned ncontrol;

    // data only for the device
    // note: not in __constant__ memory because each thread accesses a different slice (no broadcast benefit)
    real *device_pos = nullptr;     // joints x points
    real *device_vel = nullptr;     // joints x points
    real *device_acc = nullptr;     // joints x points
    real *device_t   = nullptr;     // points
    real *device_basis_p = nullptr; // nccontrol x points
    real *device_basis_v = nullptr; // nccontrol x points
    real *device_basis_a = nullptr; // nccontrol x points
    real *device_control = nullptr; // nccontrol x joints todo: should this be in constant memory??
    real dt;
    real one_over_T;
    real one_over_T2;

    // data only for the host
    Bspline *host;

    // compute the basis functions
    host_fn void init(u32 points, u32 joints, u32 p, u32 ncontrol);

    // free all host and device memory
    host_fn void clear();

    // fetch the latest computed trajectory on the GPU to host_pva
    host_fn void fetch_pva();

    // compute control points on host and send to device
    host_fn void compute_control_and_send(const Array &x, const Matrix &task);

    // compute the trajectory for the given point on the GPU (only use from kernel)
    dev_fn void compute_trajectory(unsigned point);
};

//------ HOST FUNCTIONS ------------------------------------------------------------------------------------
host_fn void cuBspline::init(u32 _points, u32 _joints, u32 _p, u32 _ncontrol) {
    if (is_init)
        return;
    is_init = true;

    points = _points;
    joints = _joints;
    p = _p;
    ncontrol = _ncontrol;

    host = new Bspline(points, joints, p, ncontrol);

    // prepare and send basis function values
    const u32 basis_bytes = points * ncontrol * sizeof(real);
    cuda_check(cudaMalloc(&device_basis_p, basis_bytes));
    cuda_check(cudaMalloc(&device_basis_v, basis_bytes));
    cuda_check(cudaMalloc(&device_basis_a, basis_bytes));
    cuda_check(cudaMemcpy(device_basis_p, host->basis_p.data, basis_bytes, cudaMemcpyHostToDevice));
    cuda_check(cudaMemcpy(device_basis_v, host->basis_v.data, basis_bytes, cudaMemcpyHostToDevice));
    cuda_check(cudaMemcpy(device_basis_a, host->basis_a.data, basis_bytes, cudaMemcpyHostToDevice));

    // prepare pos, vel, and acc matrices and t array
    cuda_check(cudaMalloc(&device_pos, points * joints * sizeof(real)));
    cuda_check(cudaMalloc(&device_vel, points * joints * sizeof(real)));
    cuda_check(cudaMalloc(&device_acc, points * joints * sizeof(real)));
    cuda_check(cudaMalloc(&device_t, points * sizeof(real)));
    cuda_check(cudaMalloc(&device_control, joints * ncontrol * sizeof(real)));
}

host_fn void cuBspline::clear() {
    printf("Clearing BSplines memory on the device\n");
    if (device_basis_p) {
        cuda_check(cudaFree(device_basis_p));
        device_basis_p = nullptr;
    }
    if (device_basis_v) {
        cuda_check(cudaFree(device_basis_v));
        device_basis_v = nullptr;
    }
    if (device_basis_a) {
        cuda_check(cudaFree(device_basis_a));
        device_basis_a = nullptr;
    }
    if (device_pos) {
        cuda_check(cudaFree(device_pos));
        device_pos = nullptr;
    }
    if (device_vel) {
        cuda_check(cudaFree(device_vel));
        device_vel = nullptr;
    }
    if (device_acc) {
        cuda_check(cudaFree(device_acc));
        device_acc = nullptr;
    }
    if (device_t) {
        cuda_check(cudaFree(device_t));
        device_t = nullptr;
    }
    is_init = false;
}

host_fn void cuBspline::fetch_pva() {
    cuda_check(cudaMemcpy(host->traj.pos.data, device_pos, host->traj.pos.size * sizeof(real), cudaMemcpyDeviceToHost));
    cuda_check(cudaMemcpy(host->traj.vel.data, device_vel, host->traj.vel.size * sizeof(real), cudaMemcpyDeviceToHost));
    cuda_check(cudaMemcpy(host->traj.acc.data, device_acc, host->traj.acc.size * sizeof(real), cudaMemcpyDeviceToHost));
    cuda_check(cudaMemcpy(host->traj.t.data, device_t, host->traj.t.size * sizeof(real), cudaMemcpyDeviceToHost));
}

host_fn void cuBspline::compute_control_and_send(const Array &x, const Matrix &task) {
    Assert(host);
    dt = x.back() / (real)(points - 1);
    one_over_T = 1 / x.back();
    one_over_T2 = one_over_T * one_over_T;
    host->compute_control(x, task);
    cuda_check(cudaMemcpy(device_control, host->control.data, joints * ncontrol * sizeof(real), cudaMemcpyHostToDevice));
}

//------ DEVICE FUNCTIONS ------------------------------------------------------------------------------------
dev_fn void cuBspline::compute_trajectory(const unsigned point) {
    device_t[point] = dt * point;
    const auto bp = device_basis_p + point * ncontrol;
    const auto bv = device_basis_v + point * ncontrol;
    const auto ba = device_basis_a + point * ncontrol;
    for (int joint = 0; joint < joints; joint++) {
        const auto c = device_control + joint * ncontrol;
        real position = 0;
        real velocity = 0;
        real acceleration = 0;
        for (int i = 0; i < ncontrol; i++) {
            position += c[i] * bp[i];
            velocity += c[i] * bv[i];
            acceleration += c[i] * ba[i];
        }
        device_pos[point * joints + joint] = position;
        device_vel[point * joints + joint] = velocity * one_over_T;
        device_acc[point * joints + joint] = acceleration * one_over_T2;
    }
}

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
#include "optional/blast_optional_utilities.hpp"
TEST_CASE("GpuCpuCorrectness", "[Trajectory]") {
    using namespace blast;

    printf("Testing on GPU with the following properties:\n");
    print_device_properties();

    const u32 points = 256;
    const u32 joints = 7;
    const u32 p = 5;
    const u32 ncontrol = 24;
    const u32 ntrajectories = 128;
    cuManagedMultiBsplines dev_bsplines(points, joints, p, ncontrol, ntrajectories);
    Bspline host_bspline(ncontrol, points, p, joints);

    // random task
    Matrix task(joints, 6);
    real amp = 10;
    for (u32 i = 0; i < task.rows; i++)
        for (u32 j = 0; j < task.cols; j++)
            task(i, j) = amp * get_random();

    // random optimization vector
    Array x(joints*(ncontrol-6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x.back() = abs(x.back());

    // copy the 'x' Array 'ntrajectories' times for the gpu version
    const u32 xlen = host_bspline.xlen(task);
    Array x_multi(ntrajectories * host_bspline.xlen(task));
    for (int i = 0; i < ntrajectories; i++) {
        for (int j = 0; j < xlen; j++) {
            x_multi[i*xlen + j] = x[j];
        }
    }

    // compute the baseline
    host_bspline.compute_trajectory(x, task);
    dev_bsplines.compute_trajectories(x_multi, task);

    // Test correctness of the trajectory
    for (int traj_id = 0; traj_id < ntrajectories; traj_id++) {
        for (int i = 0; i < (int)points; i++) {
            for (int j = 0; j < joints; j++) {
                const auto id = traj_id*points*joints + i*joints + j;
                CHECK((float)host_bspline.traj.pos(j, i) == (float)dev_bsplines.pos[id]);
                CHECK((float)host_bspline.traj.vel(j, i) == (float)dev_bsplines.vel[id]);
                CHECK((float)host_bspline.traj.acc(j, i) == (float)dev_bsplines.acc[id]);
            }
        }
    }
}
#endif

#endif