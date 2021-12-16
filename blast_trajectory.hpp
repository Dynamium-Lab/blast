#pragma once

namespace blast {

// Constants that define the b-spline
struct BsplineDef {
    u32 nctrl;
    u32 npts;
    u32 p;
    u32 njoints;
};

// Container for the three basis function matrices (p, v, a).
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
};

// Compute the basis functions that will be multiplied by the control points to produce a
// trajectory. These basis functions only change if npts, njoints, nctrl, or p changes. The basis
// functions for position (`basis.pos`) are stored in a `nctrl` x `npts` Matrix. The basis functions
// for velocity and acceleration are the same sizes.x

// Note: The matrices in the output need to already be the right size.
// Note: This function is not explicitly optimized for performance because it should not be called
//       in a loop.
void bspline_basis_functions(BsplineDef def, BsplineBasis& basis);

// compute the position, velocity and acceleration for each point, for each joint directly from the
// optimization vector.
void bspline_pva(BsplineDef def, Array& x, Matrix& task, BsplineBasis& basis, Pva& pva);





//--------------------------------------------------------------------------------------------------------------------
inline Pva::Pva(BsplineDef def) :
    pos(def.njoints, def.npts),
    vel(def.njoints, def.npts),
    acc(def.njoints, def.npts),
    t(def.npts)
{}

inline BsplineBasis::BsplineBasis(BsplineDef def) :
    pos(def.nctrl, def.npts),
    vel(def.nctrl, def.npts),
    acc(def.nctrl, def.npts)
{}

inline void bspline_basis_functions(BsplineDef def, BsplineBasis& basis) {
    const auto npts = def.npts;
    const auto nctrl = def.nctrl;
    const auto p = def.p;
    Assert(basis.pos.rows == nctrl && basis.pos.cols == npts);
    Assert(basis.vel.rows == nctrl && basis.vel.cols == npts);
    Assert(basis.acc.rows == nctrl && basis.acc.cols == npts);

    u32 m = nctrl + p;
    Array knots(m+1);
    {
        for (int i = m; i > m-p-1; i--)
            knots[i] = 1.0f;
        const real du = 1.0f / (real)(m+1 - 2*(p+1) + 1);
        for (int i = p+1; i < m-p; i++)
            knots[i] = knots[i-1] + du;
    }

    Array N(m*(p+1)); // triangle basis function
    const real du = 1.0f / (npts-1);
    real* basis_p_col = basis.pos.data;
    real* basis_v_col = basis.vel.data;
    real* basis_a_col = basis.acc.data;
    for (int point=0; point < npts; point++) {
        const real u = point*du;

        for (int i=0; i<m; i++)
            N[i] = u >= knots[i] && u < knots[i+1] ? 1.0f : 0.0f;
        if (point == npts-1)
            N[nctrl-1] = 1.0f;
        for (int pi = 1; pi <= p; pi++) {
            for (int i = 0; i < m - pi; i++) {
                if (knots[i+pi] != knots[i])
                    N[m*pi + i] = N[m*(pi-1) + i] * (u - knots[i]) / (knots[i+pi] - knots[i]);
                else
                    N[m*pi + i] = 0.0f;
                if (knots[i+pi+1] != knots[i+1])
                    N[m*pi + i] += N[m*(pi-1) + i+1] * (knots[i+pi+1] - u) / (knots[i+pi+1] - knots[i+1]);
            }
        }

        // position basis functions
        for (int i = 0; i < nctrl; i++)
            basis_p_col[i] = N[m*p + i];

        // velocity basis functions
        for (int i = 0; i < nctrl-1; i++)
            basis_v_col[i] = -(real)p * N[m*(p-1) + i+1] / (knots[i+p+1] - knots[i+1]);
        for (int i=1; i < nctrl; i++)
            basis_v_col[i] += (real)p * N[m*(p-1) + i] / (knots[i+p] - knots[i]);

        // acceleration basis functions
        for (int i = 0; i < nctrl-2; i++)
            basis_a_col[i] = (p*(p-1)) * N[m*(p-2) + i+2] / ((knots[i+p+1] - knots[i+1]) * (knots[i+p+1] - knots[i+2]));
        for (int i = 1; i < nctrl-1; i++)
            basis_a_col[i] -= (p*(p-1)) * N[m*(p-2) + i+1] * (1.0f/(knots[i+p]-knots[i]) + 1.0f/(knots[i+p+1]-knots[i+1])) / (knots[i+p] - knots[i+1]);
        for (int i = 2; i < nctrl; i++)
            basis_a_col[i] += (p*(p-1)) * N[m*(p-2) + i] / ((knots[i+p-1] - knots[i]) * (knots[i+p] - knots[i]));

        // increment pointers
        basis_p_col += nctrl;
        basis_v_col += nctrl;
        basis_a_col += nctrl;
    }
}

inline void bspline_pva(BsplineDef def, Array& x, Matrix& task, BsplineBasis& basis, Pva& pva) {
    // Note: performance critical.
    const auto njoints = def.njoints;
    const auto nctrl = def.nctrl;
    const auto p = def.p;
    const auto npts = def.npts;

    static thread_local Matrix control(nctrl, njoints);

    Assert(x.size == njoints*(nctrl-6) + 1 );
    Assert(task.rows == njoints);
    Assert(task.cols == 6);

    const real T = x[x.size-1];
    const real du = 1.0/(nctrl-p);
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

    // Note: This function is performance critical.
    // Todo: SIMD
    zero(pva.pos);
    zero(pva.vel);
    zero(pva.acc);
    for (u32 point = 0; point < npts; point++) {
        for (u32 joint = 0; joint < njoints; joint++) {
            auto& p = pva.pos(joint, point);
            auto& v = pva.vel(joint, point);
            auto& a = pva.acc(joint, point);
            for (int i = 0; i < nctrl; i++) {
                const auto c = control(i, joint);
                p += c * basis.pos(i, point);
                v += c * basis.vel(i, point);
                a += c * basis.acc(i, point);
            }
            v *= one_over_T;
            a *= one_over_T2;
        }
    }
}

} // namespace blast
