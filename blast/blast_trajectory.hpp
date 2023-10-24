#pragma once
#include "blast.hpp"

#ifdef BLAST_DEBUG
#include "optional/blast_optional_utilities.hpp"
#endif

namespace blast {

struct Trajectory {
    Matrix pos; // njoints x npoints
    Matrix vel; // njoints x npoints
    Matrix acc; // njoints x npoints
    Array t;    // npoints

    Trajectory(u32 npoints, u32 njoints) :
        pos(njoints, npoints),
        vel(njoints, npoints),
        acc(njoints, npoints),
        t(npoints) {}

    Trajectory() = default;
};

struct Bspline {
    Trajectory traj;
    Matrix control; // nctrl x njoints
    Matrix basis_p; // nctrl x npoints
    Matrix basis_v; // nctrl x npoints
    Matrix basis_a; // nctrl x npoints
    u32 joints;
    u32 points;
    u32 nctrl;
    u32 p;

    host_fn Bspline(u32 ncontrol, u32 npoints, u32 P, u32 njoints);

    // Compute a trajectory from the given optimization vector
    //  - note: fastest when 'ncontrol' is a multiple of 4 (SIMD)
    host_fn void compute_trajectory(const Array &x, Matrix &task);
    host_fn u32 xlen(const Matrix &task) {
        Assert(task.rows == joints);
        Assert(task.cols == 6);
        auto results = joints * (nctrl - 6) + 1;
        for (u32 i = 0; i < task.size; i++)
            if (std::isnan(task.data[i]))
                results++;
        return results;
    }

    host_fn void compute_basis();
    host_fn void compute_basis_open();
    host_fn void compute_control(const Array &x, const Matrix &task);
    host_fn void compute_control(const Array &x, const Matrix &task, real *dst);
};

//------ FUNCTIONS ------------------------------------------------------------------------------------

host_fn Bspline::Bspline(u32 ncontrol, u32 npoints, u32 P, u32 njoints) :
    traj(npoints, njoints),
    control(ncontrol, njoints),
    basis_p(ncontrol, npoints),
    basis_v(ncontrol, npoints),
    basis_a(ncontrol, npoints),
    joints(njoints),
    points(npoints),
    nctrl(ncontrol),
    p(P) {
    compute_basis();
}

host_fn void Bspline::compute_basis() {
    u32 m = nctrl + p;
    Array knots(m + 1);
    {
        for (u32 i = m; i > m - p - 1; i--)
            knots[i] = 1.0f;
        const real du = 1.0f / (real)(m + 1 - 2 * (p + 1) + 1);
        for (u32 i = p + 1; i < m - p; i++)
            knots[i] = knots[i - 1] + du;
    }

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

host_fn void Bspline::compute_basis_open() {
    u32 m = nctrl + p;

    Array knots(m + 1);
    knots[0] = 0.0;
    knots.back() = 1.0;
    {
        const real du = 1.0f / (real)m;
        for (u32 i = 1; i < m; i++)
            knots[i] = i * du;
    }

#ifdef BLAST_DEBUG
    print(knots);
#endif

    Matrix N(m, (p + 1));  // triangle basis function
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
            N(i, 0) = u >= knots[i] && u <= knots[i + 1] ? 1.0f : 0.0f;
        for (u32 pi = 1; pi <= p; pi++) {
            for (u32 i = 0; i < m - pi; i++) {
                if (knots[i + pi] != knots[i])
                    N(i, pi) = N(i, pi-1) * (u - knots[i]) / (knots[i + pi] - knots[i]);
                else
                    N(i, pi) = 0.0f;
                if (knots[i + pi + 1] != knots[i + 1])
                    N(i, pi) += N(i+1, pi-1) * (knots[i + pi + 1] - u) / (knots[i + pi + 1] - knots[i + 1]);
            }
        }

        // position basis functions
        for (u32 i = 0; i < nctrl; i++)
            basis_p_col[i] = N(i, p);

        // velocity basis functions
        for (u32 i = 0; i < nctrl - 1; i++)
            basis_v_col[i] = -(real)p * N(i+1, p-1) / (knots[i + p + 1] - knots[i + 1]);
        for (u32 i = 1; i < nctrl; i++)
            basis_v_col[i] += (real)p * N(i, p-1) / (knots[i + p] - knots[i]);

        // acceleration basis functions
        for (u32 i = 0; i < nctrl - 2; i++)
            basis_a_col[i] = (real)(p * (p - 1)) * N(i+2, p-2) / ((knots[i + p + 1] - knots[i + 1]) * (knots[i + p + 1] - knots[i + 2]));
        for (u32 i = 1; i < nctrl - 1; i++)
            basis_a_col[i] -= (real)(p * (p - 1)) * N(i+1, p-2) * (1.0f / (knots[i + p] - knots[i]) + 1.0f / (knots[i + p + 1] - knots[i + 1])) / (knots[i + p] - knots[i + 1]);
        for (u32 i = 2; i < nctrl; i++)
            basis_a_col[i] += (p * (p - 1)) * N(i, p-2) / ((knots[i + p - 1] - knots[i]) * (knots[i + p] - knots[i]));

        // increment pointers
        basis_p_col += nctrl;
        basis_v_col += nctrl;
        basis_a_col += nctrl;
    }
}
host_fn void Bspline::compute_control(const Array &x, const Matrix &task) {
    compute_control(x, task, control.data);
}

host_fn void Bspline::compute_control(const Array &x, const Matrix &task, real *dst) {
    using std::isnan;
    Assert(nctrl >= 6);
    const real T = x[x.size - 1];
    const real du = 1.0f / (nctrl - p);
    const real T2 = T * T;
    const real one_over_T = 1 / T;
    const real one_over_T2 = one_over_T * one_over_T;

    u32 ctr_i = 0;
    u32 x_i = 0;
    auto ctr = dst;

    const real Kv = T * du / p;
    const real Ka = 2 * T2 * du * du / (p * (p - 1));
    for (u32 joint = 0; joint < joints; joint++) {
        // Initial PVA
        const auto pi = task(joint, 0);
        const auto vi = task(joint, 1);
        const auto ai = task(joint, 2);
        ctr[ctr_i++] = isnan(pi) ? x[x_i++] : pi;
        ctr[ctr_i++] = isnan(vi) ? x[x_i++] : Kv * vi + ctr[ctr_i - 1];
        ctr[ctr_i++] = isnan(ai) ? x[x_i++] : Ka * ai + 3 * ctr[ctr_i - 1] - 2 * ctr[ctr_i - 2];

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

/**
 * @brief Computes the trajectory of a B-spline curve.
 *
 * @param x An array of real numbers representing the knots of the B-spline curve.
 * @param task A matrix representing the control points of the B-spline curve.
 *
 * This function computes the trajectory of a B-spline curve given its knots and control points.
 * The knots are passed as an array `x` and the control points are passed as a matrix `task`.
 * The function first checks that the size of `x` is equal to the length of `task` along the x-axis,
 * and that the number of rows and columns in `task` are equal to the number of joints and 6, respectively.
 * Then, it calls the `compute_control` function to compute the control points of the B-spline curve.
 * Finally, it computes the position, velocity, and acceleration of each joint at each point along the trajectory
 * using the computed control points and the basis functions for position, velocity, and acceleration.
 */
host_fn void Bspline::compute_trajectory(const Array &x, Matrix &task) {
    Assert(x.size == xlen(task));
    Assert(task.rows == joints);
    Assert(task.cols == 6);

    compute_control(x, task);

    const real T = x.back();
    const real dt = T / (points - 1);
    const real one_over_T = 1 / T;
    const real one_over_T2 = one_over_T * one_over_T;

    for (u32 point = 0; point < points; point++) {
        traj.t[point] = dt * point;
        auto bp = basis_p.col(point);
        auto bv = basis_v.col(point);
        auto ba = basis_a.col(point);
        for (u32 joint = 0; joint < joints; joint++) {
            auto c = control.col(joint);
            traj.pos(joint, point) = dot(c, bp);
            traj.vel(joint, point) = dot(c, bv) * one_over_T;
            traj.acc(joint, point) = dot(c, ba) * one_over_T2;
        }
    }
}

host_fn Trajectory compute_5order_trajectory(real T, Matrix &task) {
    const u32 joints = task.rows;
    const u32 points = (u32)ceil(T * 1000 + 1);

    Trajectory result(points, joints);

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

} // namespace blast

#ifdef BLAST_ENABLE_TESTS
TEST_CASE("SplineTest", "[Trajectory]") {
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
    Bspline bspline(nctrl, npts, p, njoints);
    bspline.compute_trajectory(x, task);
    auto& traj = bspline.traj;

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
        init_max_pos_error = std::max(init_max_pos_error, std::abs((double)traj.pos(i, 0) - (double)task(i, 0)));
        init_max_vel_error = std::max(init_max_vel_error, std::abs((double)traj.vel(i, 0) - (double)task(i, 1)));
        init_max_acc_error = std::max(init_max_acc_error, std::abs((double)traj.acc(i, 0) - (double)task(i, 2)));
        final_max_pos_error = std::max(final_max_pos_error, std::abs((double)traj.pos(i, npts - 1) - (double)task(i, 3)));
        final_max_vel_error = std::max(final_max_vel_error, std::abs((double)traj.vel(i, npts - 1) - (double)task(i, 4)));
        final_max_acc_error = std::max(final_max_acc_error, std::abs((double)traj.acc(i, npts - 1) - (double)task(i, 5)));

        // derivatives
        for (u32 j = 1; j < npts - 1; j++) {
            double diff_p = ((double)traj.pos(i, j + 1) - (double)traj.pos(i, j - 1)) / (2.0 * dt);
            max_vel_error = std::max(max_vel_error, std::abs(diff_p - (double)traj.vel(i, j)));
            double diff_v = ((double)traj.vel(i, j + 1) - (double)traj.vel(i, j - 1)) / (2.0 * dt);
            max_acc_error = std::max(max_acc_error, std::abs(diff_v - (double)traj.acc(i, j)));
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
