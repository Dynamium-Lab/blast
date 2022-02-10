#pragma once

#include "blast.hpp"

// todo: cleanup
#include <utility>
// SSE SIMD intrinsics
#include <xmmintrin.h>
// AVX SIMD intrinsics
#include <immintrin.h>


namespace blast {

// Constants that define the b-spline
struct BsplineDef {
    u32 nctrl;
    u32 npts;
    u32 p;
    u32 njoints;
};

// Container for the three basis function matrices (p, v, a).
// The basis functions do not change once the object is constructed so if you want different npts, njoints, nctrl or p, then make another object.
// The basis functions are stored in a `nctrl` x `npts` Matrices.
struct BsplineBasis {
    Matrix pos; // nctrl x npts
    Matrix vel;
    Matrix acc;

    BsplineBasis(BsplineDef def);
};

// Container for the position, velocity, and acceleration matrices and the time vector.
struct Pva {
    Matrix pos; // njoints x npts
    Matrix vel;
    Matrix acc;
    Array t;

    Pva(BsplineDef def);

    // compute the PVA using bsplines.
    Matrix bspline_control(BsplineDef def, Array& x, Matrix& task);
    void bspline(BsplineDef def, real T, Matrix& control, BsplineBasis& basis);

    u32 joints() const;
    u32 points() const;
};






//------ FUNCTIONS ------------------------------------------------------------------------------------
inline Pva::Pva(BsplineDef def) :
    pos(def.njoints, def.npts),
    vel(def.njoints, def.npts),
    acc(def.njoints, def.npts),
    t(def.npts)
{}

inline u32 Pva::joints() const {
    return pos.rows;
}

inline u32 Pva::points() const {
    return t.size;
}

inline BsplineBasis::BsplineBasis(BsplineDef def) :
    pos(def.nctrl, def.npts),
    vel(def.nctrl, def.npts),
    acc(def.nctrl, def.npts) {

    const auto npts = def.npts;
    const auto nctrl = def.nctrl;
    const auto p = def.p;

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
    const real du = 1.0f / (npts-1);
    real* basis_p_col = pos.data;
    real* basis_v_col = vel.data;
    real* basis_a_col = acc.data;
    for (u32 point=0; point < npts; point++) {
        const real u = point*du;

        for (u32 i=0; i<m; i++)
            N[i] = u >= knots[i] && u < knots[i+1] ? 1.0f : 0.0f;
        if (point == npts-1)
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

inline Matrix Pva::bspline_control(BsplineDef def, Array& x, Matrix& task) {
    const auto njoints = def.njoints;
    const auto nctrl = def.nctrl;
    const auto p = def.p;
    const auto npts = def.npts;

    Matrix control(nctrl, njoints);

    Assert(x.size == njoints*(nctrl-6) + 1 );
    Assert(task.rows == njoints);
    Assert(task.cols == 6);

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
    for (u32 joint = 0; joint < njoints; joint++) {
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

    return std::move(control);
}

// Horizontal sum of 4 lanes of the vector
inline float hadd_ps( __m128 r4 ) {
    // Add 4 values into 2
    const __m128 r2 = _mm_add_ps( r4, _mm_movehl_ps( r4, r4 ) );
    // Add 2 lower values into the final result
    const __m128 r1 = _mm_add_ss( r2, _mm_movehdup_ps( r2 ) );
    // Return the lowest lane of the result vector.
    // The intrinsic below compiles into noop, modern compilers return floats in the lowest lane of xmm0 register.
    return _mm_cvtss_f32( r1 );
}

// Horizontal sum of 8 lanes of the vector
inline float hadd_ps( __m256 r8 ) {
    const __m128 low = _mm256_castps256_ps128( r8 );
    const __m128 high = _mm256_extractf128_ps( r8, 1 );
    return hadd_ps( _mm_add_ps( low, high ) );
}

inline double hadd_pd(__m256d v) {
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
    vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128

    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
}

inline void Pva::bspline(BsplineDef def, real T, Matrix& control, BsplineBasis& basis) {
    // Note: performance critical.
    const auto njoints = def.njoints;
    const auto nctrl = def.nctrl;
    const auto npts = def.npts;
    const real one_over_T = 1/T;
    const real one_over_T2 = one_over_T*one_over_T;

    auto p = pos.data;
    auto v = vel.data;
    auto a = acc.data;
    for (u32 point = 0; point < npts; point++) {
        const auto bp = &basis.pos(0, point);
        const auto bv = &basis.vel(0, point);
        const auto ba = &basis.acc(0, point);
        for (u32 joint = 0; joint < njoints; joint++) {

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

            *p = hadd_pd(accum_p);
            *v = hadd_pd(accum_v) * one_over_T;
            *a = hadd_pd(accum_a) * one_over_T2;
            p++;
            v++;
            a++;
        }
    }
}

} // namespace blast
