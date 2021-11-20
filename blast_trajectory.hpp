#pragma once

namespace blast {

// Compute the basis functions that will be multiplied by the control points to produce a
// trajectory. These basis functions only change if the number of points, the number of control
// points or the order of the splines change. The basis functions for position (`basis_p`) are
// stored in a `nctrl` x `npts` Matrix.
// Note: This function is not explicitly optimized for performance because it should not be called in a loop.
// Note: The three output matrices need to already be the right size.
void bspline_basis_functions(u32 nctrl, u32 npts, u32 p, Matrix& basis_p, Matrix& basis_v, Matrix& basis_a);

// Compute the control points corresponding to the current optimization vector `x` that will ensure
// that the computed trajectory will respect the boundary conditions in `task`.
// - control is a nctrl x dof matrix
void bspline_control_points(u32 nctrl, u32 njoints, u32 p, Array& x, Matrix& task, Matrix& control);

// Compute the Position, velocity and acceleration for each point, for each joint, corresponding to
// the given control points and basis functions.
// - `control` is a nctrl x njoints matrix
// - `basis_p`, `basis_v`, `basis_a` are nctrl x npts matrices
// - `pos`, `vel`, `acc` are njoints x npts matrices
inline void bspline_pva(u32 nctrl, u32 npts, u32 njoints, real T, Matrix& control, Matrix& basis_p, Matrix& basis_v, Matrix& basis_a, Matrix& pos, Matrix& vel, Matrix& acc);












//--------------------------------------------------------------------------------------------------------------------
inline void bspline_control_points(u32 nctrl, u32 njoints, u32 p, Array& x, Matrix& task, Matrix& control) {
    // Note: performance critical.
    Assert(x.size == njoints*(nctrl-6) + 1 );
    Assert(task.rows == njoints);
    Assert(task.cols == 6);
    Assert(control.rows == nctrl);
    Assert(control.cols == njoints);

    float T = x[x.size-1];
    float du = 1.0f/(nctrl-p);
    float T2 = T*T;

    const auto p0 = &task(0, 0);
    const auto v0 = &task(0, 1);
    const auto a0 = &task(0, 2);
    const auto pf = &task(0, 3);
    const auto vf = &task(0, 4);
    const auto af = &task(0, 5);

    u32 ctr_i = 0;
    auto ctr = control.data;
    u32 x_idx = 0;

    const float Kv = T*du/p;
    const float Ka = 2*T2*du*du / (p*(p-1));
    for (u32 joint = 0; joint < njoints; joint++) {
        // Initial PVA
        ctr[ctr_i++] = p0[joint];
        ctr[ctr_i++] = Kv*v0[joint] + p0[joint];
        ctr[ctr_i++] = Ka*a0[joint] + 3*ctr[ctr_i - 1] - 2*p0[joint];

        // From optimization vector
        for (u32 i = 0; i < nctrl-6; i++)
            ctr[ctr_i++] = x[x_idx++];

        // Final PVA
        const float Pn = pf[joint];
        const float Pn_minus_1 = Pn - Kv*vf[joint];
        const float Pn_minus_2 = Ka*af[joint] - 2*Pn + 3*Pn_minus_1;
        ctr[ctr_i++] = Pn_minus_2;
        ctr[ctr_i++] = Pn_minus_1;
        ctr[ctr_i++] = Pn;
    }

//---- Old dof x nctrl shaped matrix ----
    // // copy first control point (p0)
    // for (u32 joint = 0; joint < njoints; joint++)
    //     ctr[ctr_i++] = p0[joint];
    // // compute control points for initial velocity
    // const float Kv = T*du/p;
    // for (u32 joint = 0; joint < njoints; joint++)
    //     ctr[ctr_i++] = Kv*v0[joint] + p0[joint];
    // // compute control points for initial acceleration
    // const float Ka = 2*T2*du*du / (p*(p-1));
    // for (u32 joint = 0; joint < njoints; joint++)
    //     ctr[ctr_i++] = Ka*a0[joint] + 3*ctr[ctr_i - njoints] - 2*p0[joint];
    // // insert all but the last element of the optimization vector
    // for (u32 i = 0; i < x.size-1; i++)
    //     ctr[ctr_i++] = x[i];
    // // compute control points for final pva
    // for (u32 joint = 0; joint < njoints; joint++) {
    //     const float Pn = pf[joint];
    //     const float Pn_minus_1 = Pn - Kv*vf[joint];
    //     const float Pn_minus_2 = Ka*af[joint] - 2*Pn + 3*Pn_minus_1;
    //     ctr[ctr_i] = Pn_minus_2;
    //     ctr[ctr_i+njoints] = Pn_minus_1;
    //     ctr[ctr_i+2*njoints] = Pn;
    //     ctr_i++;
    // }
}

inline void bspline_basis_functions(u32 nctrl, u32 npts, u32 p, Matrix& basis_p, Matrix& basis_v, Matrix& basis_a) {
    Assert(basis_p.rows == nctrl && basis_p.cols == npts);
    Assert(basis_v.rows == nctrl && basis_v.cols == npts);
    Assert(basis_a.rows == nctrl && basis_a.cols == npts);

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
    real* basis_p_col = basis_p.data;
    real* basis_v_col = basis_v.data;
    real* basis_a_col = basis_a.data;
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

inline void bspline_pva(u32 nctrl, u32 npts, u32 njoints, real T, Matrix& control, Matrix& basis_p, Matrix& basis_v, Matrix& basis_a, Matrix& pos, Matrix& vel, Matrix& acc) {
    // Note: This function is performance critical.
    // Todo: SIMD
    const auto one_over_T = 1/T;
    const auto one_over_T2 = one_over_T*one_over_T;
    zero(pos);
    zero(vel);
    zero(acc);
    for (u32 point = 0; point < npts; point++) {
        for (u32 joint = 0; joint < njoints; joint++) {
            auto& p = pos(joint, point);
            auto& v = vel(joint, point);
            auto& a = acc(joint, point);
            for (int i = 0; i < nctrl; i++) {
                const auto c = control(i, joint);
                p += c * basis_p(i, point);
                v += c * basis_v(i, point);
                a += c * basis_a(i, point);
            }
            v *= one_over_T;
            a *= one_over_T2;
        }
    }
}

} // namespace blast
