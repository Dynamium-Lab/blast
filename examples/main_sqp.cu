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
    Array       x_init;
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
Matrix& BFGS_hess(Matrix& H, Array& s, Array& q, Array& v) {
    auto H_t = transpose(H);
    auto q_t = transpose(q);
    auto s_t = transpose(s);
    auto pos_def = dot(q, s);

    if (pos_def >= 0) {
        while (pos_def < 0) {
            real m_min;
            auto i_min = matrix_min_id(pw_mult(q, s), m_min); // todo: remove m_min when not used ?
            q[i_min] = q[i_min]/2;
            // q_t = transpose(q);
            pos_def = dot(q, s);

            if (pos_def >= -0.0001) {
                real w = 1;
                while (pos_def <= 0) {
                    for (u32 i = 0; i < q.size; i++) {
                        if (q[i]*w >= 0 && q[i]*s[i] >= 0)
                            v[i] = 0;
                    }
                    // auto v_t = transpose(v);
                    q = q + v*w;
                    // q_t = transpose(q);
                    pos_def = dot(q, s);
                    w *= 2;
                }
            }
        }
    }
    auto q_s = dot(q, s);
    auto s_Hs = dot(s, H*s);
    H = H + (q*q_t)/q_s - (H*(s*s_t)*H_t)/s_Hs;
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

// qp active set
Array qp_active_set(Matrix& H, Array& grad_f, Matrix& A, Array& b) {
    unsigned iter_lim = 1000;
    auto x0 = pinv(A)*b; // todo: Matrix of 1 col is Array ?
    auto x = x0;

    auto H_eye = eye(H.cols);
    if( H == H_eye) {
        x = -grad_f;
        return x;
    }

    // Active subset initialization
    auto Ax_b = A*x - b;
    auto vmax = array_max(Ax_b);
    auto eps = 0.95*vmax;

    Array arr_act;
    Array arr_nact;
    int num_act = 0;
    int num_nact = 0;
    for (int i = 0; i < b.size; i++) {
        if (Ax_b[i] >= -eps) {
            arr_act[num_act] = i;
            num_act++;
        }
        else {
            arr_nact[num_nact] = i;
            num_nact++;
        }
    }
    Matrix A_a;
    Matrix A_na;
    Array b_a;
    Array b_na;

    // Main loop qp
    for(u32 i = 0; i< iter_lim; i++) {
        // Initialize A_a, A_na, b_a, b_na active and non active sets
        for (u32 j = 0; j < x.size; j++) {
            for (u32 i = 0; i < arr_act.size; i++) {
                A_a(j, i) = A(j, arr_act[i]);
                b_a[i] = b[arr_act[i]];
            }
            for (u32 i = 0; i < arr_nact.size; i++) {
                A_na(j, i) = A(j, arr_nact[i]);
                b_na[i] = b[arr_nact[i]];
            }
        }

        // Schur-complement method to find direction  p_k
        auto A_a_Hinv   = A_a*pinv(H);
        auto g          =  grad_f + H*x;
        auto h          = A_a*x - b_a;
        auto A_a_t      =  transpose(A_a);
        auto lambda     = pinv(A_a_Hinv*A_a_t)*(A_a_Hinv*g - h);
        auto p_k        = pinv(H)*(A_a_t*lambda - g);

        // if p_k == 0
        if (norm(p_k) < 1e-6) {
            auto lambda_bar = pinv(A_a_t)*b; // todo: check is its ok
            real lambda_min;
            auto lambda_min_id = array_min_id(lambda_bar, lambda_min);
            if(lambda_min >= 0)
                return x;
            else {  // remove lambda_min_id from arr_act
                int num_act = 0;
                int num_nact = 0;
                for (u32  i = 0;  i < arr_act.size; i++) {
                    if (arr_act[i] != lambda_min_id) {
                        arr_act[num_act] = arr_act[i];
                        num_act++;
                    }
                }
                int j = 0;
                for (u32 i = 0; i < b.size; i++) {
                    if (arr_act[j] == i) {
                        j++;
                    }
                    else {
                        arr_nact[num_nact] = i;
                        num_nact++;
                    }
                }
            }
        }
        else {
            real alpha_k;
            real alpha_min;
            unsigned alpha_min_id;
            Array alpha_tmp;
            if (A_na.size == 0)
                alpha_k  = 1;
            else {
                for (u32 i = 0; i < A_na.size; i++) {
                    auto a_na = A_na.row(i);
                    alpha_tmp[i] = (b_na[i] - dot(a_na, x))/dot(a_na, p_k);
                }
                real alpha_min;
                alpha_min_id = array_min_id(alpha_tmp, alpha_min);
                alpha_k = alpha_min < 1 ? alpha_min : 1;
            }
            if (alpha_k < 0)
                alpha_k = 0;

            x = x + alpha_k*p_k;
            if (alpha_k < 1) {
                int num_act = 0;
                int num_nact = 0;
                for (u32  i = 0;  i < arr_nact.size; i++) {
                    if (arr_nact[i] != alpha_min_id) {
                        arr_nact[num_nact] = arr_nact[i];
                        num_nact++;
                    }
                }
                int j = 0;
                for (u32 i = 0; i < b.size; i++) {
                    if (arr_nact[j] == i) {
                        j++;
                    }
                    else {
                        arr_act[num_act] = i;
                        num_act++;
                    }
                }
            }
            else {
                return x;
            }
        }
    }
}

void SQP_nocedal(unsigned iter_lim, opt_result& opt, const Array& x, real& f, Array& g, gradients& grad, real& vmax, Array& lambda) {
    unsigned iter;
    Array s;
    Array q;
    opt.x_init = x;
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
        // todo: H = hessien (H, s, q, v)
        //todo: auto grad_a_t = transpose(grad_a) (needs to be Matrix)
        //todo: auto lambda = inv(grad_a_t)

    }
    opt.x_opt = x;
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