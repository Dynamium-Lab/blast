#include <stdio.h>
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "blast_optional_utilities.hpp"

using namespace blast;

void rosensuzuki_fun(Array &x, double &f, Array &g) {
    f = x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] - 5*x[0] - 5*x[1] - 21*x[2] + 7*x[3];
    g[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] - x[1] + x[2] - x[3] - 8;
    g[1] = x[0]*x[0] + 2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] - x[0] - x[3] - 10;
    g[2] = 2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] - x[1] - x[3] - 5;
}

// Activation function
//
//  m        = number of constraints
//  result   = constraints are inserted here
//  arr_act  = list of active constraints indexes
//  arr_nact = lits of non active constraints indexes
//  num_act    = number of active constraints
//  num_nact   = number of non active constraints
//  vmax     = value of the worst violated constraint
void active(Array &g, Array &arr_act, Array &arr_nact, double &vmax) {
    const real delta = 0.05;
    int num_act = 0;
    int num_nact = 0;

    vmax = array_max(g);
    const real eps = delta - vmax;
    for (u32 j = 0; j < g.size; j++) {
        if (g[j] >= -eps) {
            arr_act.data[num_act] = j;
            num_act++;
        }
        else {
            arr_nact.data[num_nact] = j;
            num_nact++;
        }
    }
}

// Active constraints gradients
//
// m = number of constraints
// n = number of variables
host_fn void cstr_manip_active(Array &g, const Array &x, Matrix &grad, double &f, double &vmax, Array &arr_act, Array &arr_nact, Matrix &grad_a, Matrix &grad_na, Array &cstr_a, Array &cstr_na) {
    active(g, arr_act, arr_nact, vmax);
    unsigned m_act = sizeof(arr_act) / sizeof(arr_act[0]);
    auto m_nact = g.size - m_act;

    // Optimisation *opt = (Optimisation *)f_data;
    // internal_cstr_manip_single(m, result, n, x, opt);

    const real eps = 1e-5;
    Array x_plus(x.size);
    Array g_plus;
    real f_plus;
    // todo: parallel?
    for (u32 j = 0; j < x.size; j++) {
        memcpy(x_plus.data, x.data, x.size * sizeof(real));
        x_plus[j] += eps;
        // internal_cstr_manip_single(m, g_plus.data, n, x_plus.data, opt);
        rosensuzuki_fun(x_plus, f_plus, g_plus);
        for (u32 i = 0; i < g.size; i++)
            grad(j, i) = (g_plus[i] - g[i]) / eps;
    }

    for (u32 j = 0; j < x.size; j++) {
        for (u32 i = 0; i < m_act; i++) {
            grad_a(j, i) = grad(j, arr_act[i]); //[arr_act[i] * n + j];
            cstr_a[i] = g[arr_act[i]];
        }
        for (u32 i = 0; i < m_nact; i++) {
            grad_na(j, i) = grad(j, arr_nact[i]);//[arr_nact[i] * n + j];
            cstr_na[i] = g[arr_nact[i]];
        }
    }
}

host_fn void SQP_nocedal(const unsigned iter_lim, double *g, double *x, double *f, double *grad) {
    for(u32 iter = 0; iter < iter_lim; iter++) {
        // todo: fval
        // todo: iterval

        // todo: stepsize(*a, *d, *vmax)
        // todo: dn = sqrt(d_transp*d) & test if dn == 0

        // todo: x0 = x & find new x = x + a*d;
        // todo: f, g, grad _old & f, g new

        // todo: xtol & test if change in opt vector

        // todo: cstr_manip_act
        // todo: H = hessien
        //todo: auto grad_a_t = transpose(grad_a) (needs to be Matrix)
        //todo: auto lambda = inv(grad_a_t)
    }
}
int main() {
    using namespace blast;

    // init
    int iter_lim = 1000;

    // Rosensuzuki Problem
    unsigned nv = 4;
    unsigned ncon = 3;
    Array x0(nv);
    for (u32 i = 0; i < nv; i++)
        x0[i] = get_random();
    Matrix H = eye(nv);
    double f0;
    Array g0(ncon);
    rosensuzuki_fun(x0, f0, g0);

    Matrix grad;
    Matrix grad_a;
    Matrix grad_na;
    Array cstr_a;
    Array cstr_na;
    Array arr_act;
    Array arr_nact;
    double vmax;

    cstr_manip_active(g0, x0, grad, f0, vmax, arr_act, arr_nact, grad_a, grad_na, cstr_a, cstr_na);

    auto grad_a_t = transpose(grad_a);
    //todo: auto lambda = inv(grad_a_t)

    //todo: SQP

    return 0;
}