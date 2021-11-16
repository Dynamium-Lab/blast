#pragma once

namespace blast {


// Compute the basis functions that will be multiplied by the control points to produce a
// trajectory. These basis functions only change if the number of points, the number of control
// points or the order of the splines change. This function is not explicitly optimized for
// performance because it should not be called in a loop. NOTE: the three output matrices need to
// already be the right size.
void bspline_basis_functions(u32 nctrl, u32 npts, u32 p, Matrix& basis_p, Matrix& basis_v, Matrix& basis_a);





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

} // namespace blast