#pragma once
#include "blast_math.hpp"

namespace blast {

// Container for the position, velocity, and acceleration matrices and the time vector.
struct Pva {
    Matrix pos; // njoints x npoints
    Matrix vel; // njoints x npoints
    Matrix acc; // njoints x npoints
    Array t;    // npoints
    u32 joints;
    u32 points;

    Pva(u32 npoints, u32 njoints);

    // compute the PVA using bsplines.
    virtual void compute_trajectory(Array&x, Matrix& task) = 0;
    virtual u32 xlen() = 0;
};

struct PvaBspline : public Pva {
    Matrix control; // nctrl x njoints
    Matrix basis_p; // nctrl x npoints
    Matrix basis_v; // nctrl x npoints
    Matrix basis_a; // nctrl x npoints
    u32 nctrl;
    u32 p;

    PvaBspline(u32 ncontrol, u32 npoints, u32 P, u32 dof);

    void compute_trajectory(Array&x, Matrix& task) override;
    u32 xlen() override;

    void compute_basis();
    void compute_control(Array&x, Matrix& task);
};




//------ FUNCTIONS ------------------------------------------------------------------------------------
inline Pva::Pva(u32 npoints, u32 njoints) :
    pos(njoints, npoints),
    vel(njoints, npoints),
    acc(njoints, npoints),
    t(npoints),
    joints(njoints),
    points(npoints)
{}

inline PvaBspline::PvaBspline(u32 ncontrol, u32 npoints, u32 P, u32 dof) :
    Pva(npoints, dof),
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
    Array knots(m+1);
    {
        for (u32 i = m; i > m-p-1; i--)
            knots[i] = 1.0f;
        const real du = 1.0f / (real)(m+1 - 2*(p+1) + 1);
        for (u32 i = p+1; i < m-p; i++)
            knots[i] = knots[i-1] + du;
    }

    Array N(m*(p+1)); // triangle basis function
    const real du = 1.0f / (points-1);
    real* basis_p_col = basis_p.data;
    real* basis_v_col = basis_v.data;
    real* basis_a_col = basis_a.data;
    for (u32 point=0; point < points; point++) {
        const real u = point*du;

        for (u32 i=0; i<m; i++)
            N[i] = u >= knots[i] && u < knots[i+1] ? 1.0f : 0.0f;
        if (point == points-1)
            N[nctrl-1] = 1.0f;
        for (u32 pi = 1; pi <= p; pi++) {
            for (u32 i = 0; i < m - pi; i++) {
                if (knots[i+pi] != knots[i])
                    N[m*pi + i] = N[m*(pi-1) + i] * (u - knots[i]) / (knots[i+pi] - knots[i]);
                else
                    N[m*pi + i] = 0.0f;
                if (knots[i+pi+1] != knots[i+1])
                    N[m*pi + i] += N[m*(pi-1) + i+1] * (knots[i+pi+1] - u) / (knots[i+pi+1] - knots[i+1]);
            }
        }

        // position basis functions
        for (u32 i = 0; i < nctrl; i++)
            basis_p_col[i] = N[m*p + i];

        // velocity basis functions
        for (u32 i = 0; i < nctrl-1; i++)
            basis_v_col[i] = -(real)p * N[m*(p-1) + i+1] / (knots[i+p+1] - knots[i+1]);
        for (u32 i=1; i < nctrl; i++)
            basis_v_col[i] += (real)p * N[m*(p-1) + i] / (knots[i+p] - knots[i]);

        // acceleration basis functions
        for (u32 i = 0; i < nctrl-2; i++)
            basis_a_col[i] = (p*(p-1)) * N[m*(p-2) + i+2] / ((knots[i+p+1] - knots[i+1]) * (knots[i+p+1] - knots[i+2]));
        for (u32 i = 1; i < nctrl-1; i++)
            basis_a_col[i] -= (p*(p-1)) * N[m*(p-2) + i+1] * (1.0f/(knots[i+p]-knots[i]) + 1.0f/(knots[i+p+1]-knots[i+1])) / (knots[i+p] - knots[i+1]);
        for (u32 i = 2; i < nctrl; i++)
            basis_a_col[i] += (p*(p-1)) * N[m*(p-2) + i] / ((knots[i+p-1] - knots[i]) * (knots[i+p] - knots[i]));

        // increment pointers
        basis_p_col += nctrl;
        basis_v_col += nctrl;
        basis_a_col += nctrl;
    }
}

inline void PvaBspline::compute_control(Array&x, Matrix& task) {
    const real T = x[x.size-1];
    const real du = 1.0f/(nctrl-p);
    const real T2 = T*T;
    const real one_over_T = 1/T;
    const real one_over_T2 = one_over_T*one_over_T;

    const auto p0 = &task(0, 0);
    const auto v0 = &task(0, 1);
    const auto a0 = &task(0, 2);
    const auto pf = &task(0, 3);
    const auto vf = &task(0, 4);
    const auto af = &task(0, 5);

    u32 ctr_i = 0;
    u32 x_i = 0;
    auto ctr = control.data;

    const real Kv = T*du/p;
    const real Ka = 2*T2*du*du / (p*(p-1));
    for (u32 joint = 0; joint < joints; joint++) {
        // Initial PVA
        ctr[ctr_i++] = p0[joint];
        ctr[ctr_i++] = Kv*v0[joint] + p0[joint];
        ctr[ctr_i++] = Ka*a0[joint] + 3*ctr[ctr_i - 1] - 2*p0[joint];

        // From optimization vector
        for (u32 i = 0; i < nctrl-6; i++)
            ctr[ctr_i++] = x[x_i++];

        // Final PVA
        const real Pn = pf[joint];
        const real Pn_minus_1 = Pn - Kv*vf[joint];
        const real Pn_minus_2 = Ka*af[joint] - 2*Pn + 3*Pn_minus_1;
        ctr[ctr_i++] = Pn_minus_2;
        ctr[ctr_i++] = Pn_minus_1;
        ctr[ctr_i++] = Pn;
    }
}

inline u32 PvaBspline::xlen() {
    return joints*(nctrl-6)+1;
}

inline void PvaBspline::compute_trajectory(Array& x, Matrix& task) {
    Assert(x.size == xlen());
    Assert(task.rows == joints);
    Assert(task.cols == 6);

    compute_control(x, task);

    const real one_over_T = 1/x[x.size-1];
    const real one_over_T2 = one_over_T*one_over_T;

    auto p = pos.data;
    auto v = vel.data;
    auto a = acc.data;
    for (u32 point = 0; point < points; point++) {
        const auto bp = &basis_p(0, point);
        const auto bv = &basis_v(0, point);
        const auto ba = &basis_a(0, point);
        for (u32 joint = 0; joint < joints; joint++) {

            const auto c = &control(0, joint);

            auto accum_p = _mm256_setzero_pd();
            auto accum_v = _mm256_setzero_pd();
            auto accum_a = _mm256_setzero_pd();
            for (u32 i = 0; i<nctrl; i += 4) {
                const auto c_v  = _mm256_loadu_pd(&c[i]);
                const auto bp_v  = _mm256_loadu_pd(&bp[i]);
                const auto bv_v  = _mm256_loadu_pd(&bv[i]);
                const auto ba_v  = _mm256_loadu_pd(&ba[i]);
                accum_p = _mm256_fmadd_pd(c_v, bp_v, accum_p);
                accum_v = _mm256_fmadd_pd(c_v, bv_v, accum_v);
                accum_a = _mm256_fmadd_pd(c_v, ba_v, accum_a);
            }

            *p = simd_hadd(accum_p);
            *v = simd_hadd(accum_v) * one_over_T;
            *a = simd_hadd(accum_a) * one_over_T2;
            p++;
            v++;
            a++;
        }
    }
}

} // namespace blast
