#pragma once
#include "blast.hpp"

#include "blast_optional_utilities.hpp"

namespace blast {

// Container for the position, velocity, and acceleration matrices and the time vector.
struct Pva {
    Matrix pos; // njoints x npoints
    Matrix vel; // njoints x npoints
    Matrix acc; // njoints x npoints
    Array t;    // npoints
    u32 joints = 0;
    u32 points = 0;

    Pva() = default;
    Pva(u32 njoints, u32 npoints);

    Pva &operator=(const Pva &);

    // compute the PVA using bsplines.
    virtual void compute_trajectory(const Array &x, Matrix &task) {
        (void)x;
        (void)task;
    }
    virtual u32 xlen(Matrix &task) {
        return 0;
    }
};

struct PvaBspline : public Pva {
    Matrix control; // nctrl x njoints
    Matrix basis_p; // nctrl x npoints
    Matrix basis_v; // nctrl x npoints
    Matrix basis_a; // nctrl x npoints
    u32 nctrl;
    u32 p;

    PvaBspline(u32 ncontrol, u32 npoints, u32 P, u32 dof);

    // Compute a trajectory from the given optimization vector
    //  - note: fastest when 'ncontrol' is a multiple of 4 (SIMD)
    void compute_trajectory(const Array &x, Matrix &task) override;
    u32 xlen(Matrix &task) override;

    void compute_basis();
    void compute_control(const Array &x, const Matrix &task);
};

// the following structure is only usefull when using CUDA for Nvidia GPUs
#ifdef __NVCC__
struct cuPvaBspline {
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
    real *device_t = nullptr;       // points
    real *device_basis_p = nullptr; // nccontrol x points
    real *device_basis_v = nullptr; // nccontrol x points
    real *device_basis_a = nullptr; // nccontrol x points
    real *device_control = nullptr; // nccontrol x joints todo: should this be in constant memory??
    real dt;
    real one_over_T;
    real one_over_T2;

    // data only for the host
    PvaBspline *host = nullptr;

    // compute the basis functions
    host_fn void init(u32 points, u32 joints, u32 p, u32 ncontrol);

    // free all host and device memory
    host_fn void clear();

    // fetch the latest computed trajectory on the GPU to host_pva
    host_fn void fetch_pva();

    // compute the control points based on the optimization vector and the task
    host_fn void compute_control(const Array &x, const Matrix &task);

    // send the last computed control points to device memory
    host_fn void send_control();

    // compute control points on host and send to device
    host_fn void compute_control_and_send(const Array &x, const Matrix &task);

    // compute the trajectory for the given point
    dev_fn void compute_trajectory(unsigned point);
};
#endif

//------ FUNCTIONS ------------------------------------------------------------------------------------
inline Pva::Pva(u32 njoints, u32 npoints) : pos(njoints, npoints),
    vel(njoints, npoints),
    acc(njoints, npoints),
    t(npoints),
    joints(njoints),
    points(npoints) {
}

inline Pva &Pva::operator=(const Pva &rhs) {
    pos = rhs.pos;
    vel = rhs.vel;
    acc = rhs.acc;
    t = rhs.t;
    joints = rhs.joints;
    points = rhs.points;
    return *this;
}

inline PvaBspline::PvaBspline(u32 ncontrol, u32 npoints, u32 P, u32 dof) : Pva(dof, npoints),
    nctrl(ncontrol),
    p(P),
    control(ncontrol, dof),
    basis_p(ncontrol, npoints),
    basis_v(ncontrol, npoints),
    basis_a(ncontrol, npoints) {
    compute_basis();
}

inline void PvaBspline::compute_basis() {
    u32 m = nctrl + p;
    Array knots(m + 1);
    {
        for (u32 i = m; i > m - p - 1; i--)
            knots[i] = 1.0f;
        const real du = 1.0f / (real)(m + 1 - 2 * (p + 1) + 1);
        for (u32 i = p + 1; i < m - p; i++)
            knots[i] = knots[i - 1] + du;
    }
    print(knots);

    Array N(m * (p + 1)); // triangle basis function
    const real du = 1.0f / (points - 1);
    zero(basis_p);
    zero(basis_v);
    zero(basis_a);
    real *basis_p_col = basis_p.data;
    real *basis_v_col = basis_v.data;
    real *basis_a_col = basis_a.data;
    for (u32 point = 0; point < points; point++) {
        const real u = point * du;

        for (u32 i = 0; i < m; i++)
            N[i] = u >= knots[i] && u < knots[i + 1] ? 1.0f : 0.0f;
        if (point == points - 1)
            N[nctrl - 1] = 1.0f;
        for (u32 pi = 1; pi <= p; pi++) {
            for (u32 i = 0; i < m - pi; i++) {
                if (knots[i + pi] != knots[i])
                    N[m * pi + i] = N[m * (pi - 1) + i] * (u - knots[i]) / (knots[i + pi] - knots[i]);
                else
                    N[m * pi + i] = 0.0f;
                if (knots[i + pi + 1] != knots[i + 1])
                    N[m * pi + i] += N[m * (pi - 1) + i + 1] * (knots[i + pi + 1] - u) / (knots[i + pi + 1] - knots[i + 1]);
            }
        }

        // position basis functions
        for (u32 i = 0; i < nctrl; i++)
            basis_p_col[i] = N[m * p + i];

        // velocity basis functions
        for (u32 i = 0; i < nctrl - 1; i++)
            basis_v_col[i] = -(real)p * N[m * (p - 1) + i + 1] / (knots[i + p + 1] - knots[i + 1]);
        for (u32 i = 1; i < nctrl; i++)
            basis_v_col[i] += (real)p * N[m * (p - 1) + i] / (knots[i + p] - knots[i]);

        // acceleration basis functions
        for (u32 i = 0; i < nctrl - 2; i++)
            basis_a_col[i] = (p * (p - 1)) * N[m * (p - 2) + i + 2] / ((knots[i + p + 1] - knots[i + 1]) * (knots[i + p + 1] - knots[i + 2]));
        for (u32 i = 1; i < nctrl - 1; i++)
            basis_a_col[i] -= (p * (p - 1)) * N[m * (p - 2) + i + 1] * (1.0f / (knots[i + p] - knots[i]) + 1.0f / (knots[i + p + 1] - knots[i + 1])) / (knots[i + p] - knots[i + 1]);
        for (u32 i = 2; i < nctrl; i++)
            basis_a_col[i] += (p * (p - 1)) * N[m * (p - 2) + i] / ((knots[i + p - 1] - knots[i]) * (knots[i + p] - knots[i]));

        // increment pointers
        basis_p_col += nctrl;
        basis_v_col += nctrl;
        basis_a_col += nctrl;
    }
}

inline void PvaBspline::compute_control(const Array &x, const Matrix &task) {
    Assert(nctrl >= 6);
    const real T = x[x.size - 1];
    const real du = 1.0f / (nctrl - p);
    const real T2 = T * T;
    const real one_over_T = 1 / T;
    const real one_over_T2 = one_over_T * one_over_T;


    u32 ctr_i = 0;
    u32 x_i = 0;
    auto ctr = control.data;

    const real Kv = T * du / p;
    const real Ka = 2 * T2 * du * du / (p * (p - 1));
    for (u32 joint = 0; joint < joints; joint++) {
        // Initial PVA
        const auto pi = task(joint, 0);
        const auto vi = task(joint, 1);
        const auto ai = task(joint, 2);
        ctr[ctr_i++] = isnan(pi) ? x[x_i++] : pi;
        ctr[ctr_i++] = isnan(vi) ? x[x_i++] : Kv * vi + ctr[ctr_i - 1];
        ctr[ctr_i++] = isnan(ai) ? x[x_i++] : Ka * ai + 3 * ctr[ctr_i -1] - 2 * ctr[ctr_i - 2];

        // From optimization vector
        for (u32 i = 3; i < nctrl - 3; i++)
            ctr[ctr_i++] = x[x_i++];

        // Final PVA
        const auto pf = task(joint, 3);
        const auto vf = task(joint, 4);
        const auto af = task(joint, 5);
        const real Pn = isnan(pf) ? x[x_i++] : pf;
        const real Pn_minus_1 = isnan(vf) ? x[x_i++] : Pn - Kv * vf;
        const real Pn_minus_2 = isnan(af) ? x[x_i++] : Ka * af - 2 * Pn + 3 * Pn_minus_1;
        ctr[ctr_i++] = Pn_minus_2;
        ctr[ctr_i++] = Pn_minus_1;
        ctr[ctr_i++] = Pn;
    }
}

inline u32 PvaBspline::xlen(Matrix &task) {
    Assert(task.rows == joints);
    Assert(task.cols == 6);
    auto results = joints * (nctrl - 6) + 1;
    for (u32 i = 0; i < task.size; i++)
        if (isnan(task.data[i]))
            results++;
    return results;
}

inline void PvaBspline::compute_trajectory(const Array &x, Matrix &task) {
    Assert(x.size == xlen(task));
    Assert(task.rows == joints);
    Assert(task.cols == 6);

    compute_control(x, task);

    const real T = x.back();
    const real dt = T / (points - 1);
    const real one_over_T = 1 / T;
    const real one_over_T2 = one_over_T * one_over_T;

    for (u32 point = 0; point < points; point++) {
        t[point] = dt * point;
        auto bp = basis_p.col(point);
        auto bv = basis_v.col(point);
        auto ba = basis_a.col(point);
        for (u32 joint = 0; joint < joints; joint++) {
            auto c = control.col(joint);
            pos(joint, point) = dot(c, bp);
            vel(joint, point) = dot(c, bv) * one_over_T;
            acc(joint, point) = dot(c, ba) * one_over_T2;
        }
    }
}

inline Pva compute_5order_trajectory(real T, Matrix &task) {
    const u32 joints = task.rows;
    const u32 points = (u32)ceil(T * 1000 + 1);

    Pva result(joints, points);

    const real deltaT = 0.001f;
    T = (points - 1) * deltaT;

    Matrix A(6, joints);

    for (u32 j = 0; j < joints; j++) {
        // todo: Fix for trajectory that has NaN values
        const auto p0 = task(j, 0);
        const auto v0 = task(j, 1) * T;
        const auto a0 = task(j, 2) * T * T;
        const auto pf = task(j, 3);
        const auto vf = task(j, 4) * T;
        const auto af = task(j, 5) * T * T;

        A(0, j) = p0;
        A(1, j) = v0;
        A(2, j) = a0 * 0.5f;
        A(3, j) = 0.5f * af - 1.5f * a0 - 10 * (p0 - pf) - 6 * v0 - 4 * vf;
        A(4, j) = 1.5f * a0 - af + 15 * (p0 - pf) + 8 * v0 + 7 * vf;
        A(5, j) = 0.5f * (af - a0) - 6 * (p0 - pf) - 3 * (v0 + vf);
    }

    for (u32 i = 0; i < points; i++) {
        const auto s = (real)i / (real)(points - 1);
        const auto s2 = s * s;
        const auto s3 = s2 * s;
        const auto s4 = s3 * s;
        const auto s5 = s4 * s;
        result.t[i] = i * deltaT;
        for (u32 j = 0; j < joints; j++) {
            const auto a = &A(0, j);

            const auto q = a[0] + a[1] * s + a[2] * s2 + a[3] * s3 + a[4] * s4 + a[5] * s5;
            const auto qd = a[1] + 2 * a[2] * s + 3 * a[3] * s2 + 4 * a[4] * s3 + 5 * a[5] * s4;
            const auto qdd = 2 * a[2] + 6 * a[3] * s + 12 * a[4] * s2 + 20 * a[5] * s3;

            result.pos(j, i) = q;
            result.vel(j, i) = qd / T;
            result.acc(j, i) = qdd / (T * T);
        }
    }

    return result;
}

//------ FOR GPU COMPUTATION ONLY ------------------------------------------------------------------------------------
#ifdef __NVCC__

host_fn void cuPvaBspline::init(u32 _points, u32 _joints, u32 _p, u32 _ncontrol) {
    if (is_init)
        return;
    is_init = true;

    points = _points;
    joints = _joints;
    p = _p;
    ncontrol = _ncontrol;

    host = new PvaBspline(ncontrol, points, p, joints);

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

host_fn void cuPvaBspline::clear() {
    printf("Clearing memory\n");
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
    delete (host);
    host = nullptr;
    is_init = false;
}

host_fn void cuPvaBspline::fetch_pva() {
    cuda_check(cudaMemcpy(host->pos.data, device_pos, host->pos.size * sizeof(real), cudaMemcpyDeviceToHost));
    cuda_check(cudaMemcpy(host->vel.data, device_vel, host->vel.size * sizeof(real), cudaMemcpyDeviceToHost));
    cuda_check(cudaMemcpy(host->acc.data, device_acc, host->acc.size * sizeof(real), cudaMemcpyDeviceToHost));
    cuda_check(cudaMemcpy(host->t.data, device_t, host->t.size * sizeof(real), cudaMemcpyDeviceToHost));
}

host_fn void cuPvaBspline::compute_control(const Array &x, const Matrix &task) {
    Assert(host);
    dt = x.back() / (real)(points - 1);
    one_over_T = 1 / x.back();
    one_over_T2 = one_over_T * one_over_T;
    host->compute_control(x, task);
}

host_fn void cuPvaBspline::send_control() {
    cuda_check(cudaMemcpy(device_control, host->control.data, joints * ncontrol * sizeof(real), cudaMemcpyHostToDevice));
}

host_fn void cuPvaBspline::compute_control_and_send(const Array &x, const Matrix &task) {
    compute_control(x, task);
    send_control();
}

//------ DEVICE FUNCTIONS ------------------------------------------------------------------------------------
dev_fn void cuPvaBspline::compute_trajectory(const unsigned point) {
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

#endif

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
TEST_CASE ("SplineTest", "[Trajectory]") {
    using namespace blast;

    const u32 nctrl = 8 * 3;
    const u32 npts = 4000;
    const u32 p = 5;
    const u32 njoints = 7;

    // random task
    real amp = 10;
    Matrix task(njoints, 6);
    for (u32 i = 0; i < task.rows; i++) {
        for (u32 j = 0; j < task.cols; j++) {
            task(i, j) = amp * get_random();
        }
    }

    // random optimization vector
    Array x(njoints * (nctrl - 6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x[x.size - 1] = 3.f;

    // Compute basis functions

    // Compute trajectory
    PvaBspline pva(nctrl, npts, p, njoints);
    pva.compute_trajectory(x, task);

    double init_max_pos_error = 0;
    double init_max_vel_error = 0;
    double init_max_acc_error = 0;
    double final_max_pos_error = 0;
    double final_max_vel_error = 0;
    double final_max_acc_error = 0;
    double max_vel_error = 0;
    double max_acc_error = 0;
    double dt = x[x.size - 1] / (double)(npts - 1);
    for (u32 i = 0; i < njoints; i++) {
        // boundary conditions
        init_max_pos_error = std::max(init_max_pos_error, std::abs((double)pva.pos(i, 0) - (double)task(i, 0)));
        init_max_vel_error = std::max(init_max_vel_error, std::abs((double)pva.vel(i, 0) - (double)task(i, 1)));
        init_max_acc_error = std::max(init_max_acc_error, std::abs((double)pva.acc(i, 0) - (double)task(i, 2)));
        final_max_pos_error = std::max(final_max_pos_error, std::abs((double)pva.pos(i, npts - 1) - (double)task(i, 3)));
        final_max_vel_error = std::max(final_max_vel_error, std::abs((double)pva.vel(i, npts - 1) - (double)task(i, 4)));
        final_max_acc_error = std::max(final_max_acc_error, std::abs((double)pva.acc(i, npts - 1) - (double)task(i, 5)));

        // derivatives
        for (u32 j = 1; j < npts - 1; j++) {
            double diff_p = ((double)pva.pos(i, j + 1) - (double)pva.pos(i, j - 1)) / (2.0 * dt);
            max_vel_error = std::max(max_vel_error, std::abs(diff_p - (double)pva.vel(i, j)));
            double diff_v = ((double)pva.vel(i, j + 1) - (double)pva.vel(i, j - 1)) / (2.0 * dt);
            max_acc_error = std::max(max_acc_error, std::abs(diff_v - (double)pva.acc(i, j)));
        }
    }

    REQUIRE(init_max_pos_error < 0.2);
    REQUIRE(init_max_vel_error < 0.2);
    REQUIRE(init_max_acc_error < 0.2);
    REQUIRE(final_max_pos_error < 0.2);
    REQUIRE(final_max_vel_error < 0.2);
    REQUIRE(final_max_acc_error < 0.2);
    REQUIRE(max_vel_error < 0.2);
    REQUIRE(max_acc_error < 0.9);
}
#endif
