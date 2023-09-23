#include <stdio.h>
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "blast_optional_utilities.hpp"

using namespace blast;

struct gradients {
    Array   grad_f;
    Matrix  grad_g;
    Matrix  grad_a;
    Matrix  grad_na;
    Array   arr_act;
    Array   arr_nact;
    Array   cstr_a;
    Array   cstr_na;
};

struct opt_result {
    Array       fval;
    unsigned    iter_val;
    Array       x_opt;
};


void rosensuzuki_fun(Array& x, real& f, Array& g) {
    f = x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] - 5*x[0] - 5*x[1] - 21*x[2] + 7*x[3];
    g[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] - x[1] + x[2] - x[3] - 8;
    g[1] = x[0]*x[0] + 2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] - x[0] - x[3] - 10;
    g[2] = 2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] - x[1] - x[3] - 5;
}

// Hessien approximation using BFGS
//
// H = hessien
// s(1, nv) = x - x0; (in main loop)
// q(1,nv) = (grad_f_t + lambda*grad_g_t) - (grad_f_t_old + lambda*grad_g_t_old)
// v(nv) = grad_g*g - grad_g_old*g_old
Matrix& BFGS_hess(Matrix& H, Matrix& s, Matrix& q, Matrix& v) {
    auto H_t = transpose(H);
    auto q_t = transpose(q); // q_t(nv,1)
    auto s_t = transpose(s);
    auto pos_def = q_t*s; // todo: Array/Matrix size 1 becomes a real ?

    if (pos_def(0, 0) >= 0) { // todo: proper fix
        while (pos_def(0, 0) < 0) { // todo: proper fix
            auto i_min = matrix_min_id(pw_mult(q, s));
            q(0, i_min) = q(0, i_min)/2;
            q_t = transpose(q);
            pos_def = q_t*s;

            if (pos_def(0, 0) >= -0.0001) { // todo: proper fix
                real w = 1;
                while (pos_def(0, 0) <= 0) { // todo: proper fix
                    for (u32 i = 0; i < q.size; i++) {
                        if (q(0, i)*w >= 0 && q(0, i)*s(0, i) >= 0)
                            v(0, i) = 0;
                    }
                    auto v_t = transpose(v);
                    q = q + v_t*w;
                    q_t = transpose(q);
                    pos_def = q_t*s;
                    w *= 2;
                }
            }
        }
    }
    H = H + (q*q_t)*pinv(q_t*s) - H*(s*s_t)*H_t*pinv(s_t*H*s);
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
void cstr_manip_active(const Array& x, real& f, Array& g, gradients& grad, real& vmax) {

    active(g, grad.arr_act, grad.arr_nact, vmax);
    unsigned m_act = sizeof(grad.arr_act) / sizeof(grad.arr_act[0]);
    auto m_nact = g.size - m_act;

    // Optimisation *opt = (Optimisation *)f_data;
    // internal_cstr_manip_single(m, result, n, x, opt);

    const real eps = 1e-5;
    Array x_plus(x.size);
    Array g_plus;
    real f_plus;
    // todo: parallel?
    for (u32 j = 0; j < x.size; j++) {
        x_plus = x;
        x_plus[j] += eps;
        // internal_cstr_manip_single(m, g_plus.data, n, x_plus.data, opt);
        rosensuzuki_fun(x_plus, f_plus, g_plus);
        grad.grad_f[j] = (f_plus - f) / eps;
        for (u32 i = 0; i < g.size; i++)
            grad.grad_g(j, i) = (g_plus[i] - g[i]) / eps;
    }

    for (u32 j = 0; j < x.size; j++) {
        for (u32 i = 0; i < m_act; i++) {
            grad.grad_a(j, i) = grad.grad_g(j, grad.arr_act[i]);
            grad.cstr_a[i] = g[grad.arr_act[i]];
        }
        for (u32 i = 0; i < m_nact; i++) {
            grad.grad_na(j, i) = grad.grad_g(j, grad.arr_nact[i]);
            grad.cstr_na[i] = g[grad.arr_nact[i]];
        }
    }
}

void SQP_nocedal(unsigned iter_lim, opt_result& opt, const Array& x, real& f, Array& g, gradients& grad, real& vmax, Array& lambda) {
    unsigned iter;
    Array s;
    Array q;
    for(iter = 0; iter < iter_lim; iter++) {

        opt.fval[iter] = f;
        print(f);
        // todo: opt.x_opt for optimal x

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
    opt.iter_val = iter;
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

    gradients  grad;
    opt_result opt;
    double vmax;

    cstr_manip_active(x0, f0, g0, grad, vmax);

    auto grad_a_t = transpose(grad.grad_a);
    auto lambda = pinv(grad_a_t)*(grad.grad_f + H*x0);

    //todo: SQP

    return 0;
}