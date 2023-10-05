#include <stdio.h>
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "blast_optional_utilities.hpp"

using namespace blast;

struct gradients {
    Array           grad_f;
    Matrix          grad_g;
    Matrix          grad_a;
    Matrix          grad_na;
    vector<real>    arr_act;
    vector<real>    arr_nact;
    Array           cstr_a;
    Array           cstr_na;
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
    return H;
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
void active(Array& g, vector<real>& arr_act, vector<real>& arr_nact, double& vmax) {
    const real delta = 0.05;

    vmax = array_max(g);
    const real eps = delta - vmax;
    for (u32 j = 0; j < g.size; j++) {
        if (g[j] >= -eps)
            arr_act.push_back(j);
        else
            arr_nact.push_back(j);
    }
}

// Active constraints gradients
//
// m = number of constraints
// n = number of variables
void cstr_manip_active(Array& x, real& f, Array& g, gradients& grad, real& vmax) {
    vector<real> arr_act;
    vector<real> arr_nact;
    active(g, arr_act, arr_nact, vmax);
    unsigned m_act = arr_act.size();
    auto m_nact = g.size - m_act;

    // Optimisation *opt = (Optimisation *)f_data;
    // internal_cstr_manip_single(m, result, n, x, opt);

    const real eps = 1e-5;
    Array x_plus(x.size);
    Array g_plus(g.size);
    real f_plus;

    Array grad_f(x.size);
    Matrix grad_g(x.size, g.size);
    Matrix grad_a(x.size, m_act);
    Matrix grad_na(x.size, m_nact);
    Array cstr_a(m_act);
    Array cstr_na(m_nact);
    // todo: parallel?
    for (u32 j = 0; j < x.size; j++) {
        x_plus = x;
        x_plus[j] += eps;
        // internal_cstr_manip_single(m, g_plus.data, n, x_plus.data, opt);
        rosensuzuki_fun(x_plus, f_plus, g_plus);
        grad_f[j] = (f_plus - f) / eps;
        for (u32 i = 0; i < g.size; i++)
            grad_g(j, i) = (g_plus[i] - g[i]) / eps;
    }

    for (u32 j = 0; j < x.size; j++) {
        for (u32 i = 0; i < m_act; i++) {
            grad_a(j, i) = grad_g(j, arr_act[i]);
            cstr_a[i] = g[arr_act[i]];
        }
        for (u32 i = 0; i < m_nact; i++) {
            grad_na(j, i) = grad_g(j, arr_nact[i]);
            cstr_na[i] = g[arr_nact[i]];
        }
    }

    grad.grad_f = grad_f;
    grad.grad_g = grad_g;
    grad.grad_a = grad_a;
    grad.grad_na = grad_na;
    grad.arr_act = arr_act;
    grad.arr_nact = arr_nact;
    grad.cstr_a = cstr_a;
    grad.cstr_na = cstr_na;
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

    vector<real> arr_act;
    vector<real> arr_nact;
    for (int i = 0; i < b.size; i++) {
        if (Ax_b[i] >= -eps) {
            arr_act.push_back(i);
        }
        else {
            arr_nact.push_back(i);
        }
    }
    unsigned m_act = sizeof(arr_act) / sizeof(arr_act[0]);
    auto m_nact = b.size - m_act;

    Matrix A_a(x.size, m_act);
    Matrix A_na(x.size, m_nact);
    Array b_a(m_act);
    Array b_na(m_nact);

    // Main loop qp
    for(u32 i = 0; i< iter_lim; i++) {
        // Initialize A_a, A_na, b_a, b_na active and non active sets
        for (u32 j = 0; j < x.size; j++) {
            for (u32 i = 0; i < m_act; i++) {
                A_a(j, i) = A(j, arr_act[i]);
                b_a[i] = b[arr_act[i]];
            }
            for (u32 i = 0; i < m_nact; i++) {
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
                arr_act.clear();
                arr_nact.clear();
                for (u32  i = 0;  i < m_act; i++) {
                    if (arr_act[i] != lambda_min_id)
                        arr_act.push_back(arr_act[i]);
                }
                int j = 0;
                for (u32 i = 0; i < b.size; i++) {
                    if (arr_act[j] == i)
                        j++;
                    else
                        arr_nact.push_back(i);
                }
            }
        }
        else {
            real alpha_k;
            real alpha_min;
            unsigned alpha_min_id;
            Array alpha_tmp(A_na.rows);
            if (A_na.size == 0)
                alpha_k  = 1;
            else {
                for (u32 i = 0; i < A_na.rows; i++) {
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
                arr_act.clear();
                arr_nact.clear();
                for (u32  i = 0;  i < m_nact; i++) {
                    if (arr_nact[i] != alpha_min_id)
                        arr_nact.push_back(arr_nact[i]);
                }
                int j = 0;
                for (u32 i = 0; i < b.size; i++) {
                    if (arr_nact[j] == i) {
                        j++;
                    }
                    else {
                        arr_act.push_back(i);
                    }
                }
            }
            else {
                return x;
            }
        }
    }
    return x;
}

real step_size(Matrix& H, gradients& grad, Array& d, Array& x, real& f, Array& g, real& vmax) {
    d = qp_active_set(H, grad.grad_f, -grad.grad_g, g);

    auto x_new = x + d;
    real f_new;
    Array g_new(g.size);
    gradients grad_new;
    rosensuzuki_fun(x_new, f_new, g_new);
    cstr_manip_active(x_new, f_new, g_new, grad_new, vmax);
    auto g_new_max = array_max(g_new);
    auto g_max = array_max(g);

    auto grad_na = grad_new.grad_na;
    auto g_na = grad_new.cstr_na;

    real a_bar;
    if(g_new_max <= g_max)
        a_bar = 1;
    else {
        auto a_d = array_min(grad_na*d);
        vector<real> a_tmp;
        if (a_d < 0) {
            for (u32 i = 0; i < grad_na.size; i++) {
                auto a_na = grad_na.row(i);
                a_tmp.push_back((g_na[i] - dot(a_na, x))/dot(a_na, d));
            }
        }
        else {
            a_tmp.push_back(1);
        }
        unsigned m_tmp_size = sizeof(a_tmp) / sizeof(a_tmp[0]);
        for (u32 i = 0; i < m_tmp_size; i++) {
            if (a_tmp[i] < 0)
                a_tmp[i] = INF_REAL;
        }
        a_bar = array_min(a_tmp);
    }
    a_bar = a_bar < 1 ? 1: a_bar;

    Array r(g.size);
    Array r_new(g.size);
    for(u32 i = 0; i < g.size; i++) {
        r[i] = norm(grad.grad_f)/norm(grad.grad_g.row(i));
        r_new[i] = norm(grad_new.grad_f)/norm(grad_new.grad_g.row(i));
    }

    auto grad_a_t = transpose(grad.grad_a);
    auto lambda = pinv(grad_a_t)*(grad.grad_f + H*x);
    real w = 0.5;
    unsigned iter_lim = 100;
    real a = a_bar;
    for (u32 iter = 0; iter < iter_lim; iter++) {
        x_new = x + a*d;
        rosensuzuki_fun(x_new, f_new, g_new);
        cstr_manip_active(x_new, f_new, g_new, grad_new, vmax);
        auto grad_a_t_new = transpose(grad_new.grad_a);
        auto lambda_new = pinv(grad_a_t_new)*(grad_new.grad_f + H*x_new);

        Array g_array_max(g.size);
        Array g_array_max_new(g.size);
        Array r_g(g.size);
        Array r_g_new(g.size);
        for (u32 i = 0; i < g.size; i++) {
            g_array_max[i] = 0.0 <= g[i] ? g[i]: 0.0;
            g_array_max_new[i] = 0.0 <= g_new[i] ? g_new[i]: 0.0;

            r_g[i] = r[i]*g_array_max[i];
            r_g_new[i] = r_new[i]*g_array_max_new[i];

            r[i] = lambda[i] < (r[i] + lambda[i])/2 ? (r[i] + lambda[i])/2 : lambda[i];
            r_new[i] = lambda_new[i] < (r_new[i] + lambda_new[i])/2 ? (r_new[i] + lambda_new[i])/2 : lambda_new[i];
        }

        auto Psi = f + sum(r_g);
        auto Psi_new = f_new + sum(r_g_new);
        auto grad_f_t = transpose(grad.grad_f);
        real coeff = 0.5;
        auto Psi_comp = Psi + coeff*(a*grad_f_t)*d; // todo: Psi_comp 1d Array ?
        if (Psi_new <= Psi_comp[0])
            return a;
        a *= w;
    }
    return a;
}

void SQP_nocedal(unsigned iter_lim, opt_result& opt, Array& x, real& f, Array& g, Matrix& H, gradients& grad, real& vmax, Array& lambda) {
    unsigned iter;
    Array s(x.size);
    Matrix q(x.size, 1);
    Array v(x.size);
    Array d(x.size);
    opt.x_init = x;
    for(iter = 0; iter < iter_lim; iter++) {
        print(f);

        real a = step_size(H, grad, d, x, f, g, vmax);
        auto dn = norm(d);
        if (dn == 0) {
            break; // todo: will this work ?
        }

        auto x0 = x;
        x = x + a*d;
        auto f_old = f;
        auto g_old = g;
        auto grad_g_old = grad.grad_g;
        auto grad_f_old = grad.grad_f;
        rosensuzuki_fun(x, f, g);

        auto xtol = norm(x)*1e-6 + 1e-12;
        if (norm(x - x0) < xtol) {
            break; // todo: will this work ?
        }

        auto ftol = abs(f)*1e-6 + 1e-12;
        auto tolG = 1e-6;
        if (abs(f - f_old) < ftol && vmax < tolG) {
            break; // todo: will this work ?
        }

        cstr_manip_active(x, f, g, grad, vmax);
        auto grad_g = grad.grad_g;
        auto grad_f = grad.grad_f;

        s = x - x0;
        q = (transpose(grad_f) + lambda*transpose(grad_g)) - (transpose(grad_f_old) + lambda*transpose(grad_g_old));
        v = grad_g*g - grad_g_old*g_old;
        H = BFGS_hess(H, s, q.row(1), v); // todo: row or col for q ?

        auto grad_a_t = transpose(grad.grad_a);
        auto lambda = pinv(grad_a_t)*(grad.grad_f + H*x);
    }
    opt.fval = f;
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
    x0[0] = 1;
    x0[1] = 1;
    x0[2] = 1;
    x0[3] = 1;
    // for (u32 i = 0; i < nv; i++)
    //     x0[i] = get_random();
    Matrix H = eye(nv);
    double f0;
    Array g0(ncon);
    rosensuzuki_fun(x0, f0, g0);

    gradients  grad;
    opt_result opt;
    double vmax;

    cstr_manip_active(x0, f0, g0, grad, vmax);

    // todo: fix lambda pinv
    Array lambda(g0.size);
    auto grad_a_t = transpose(grad.grad_a);
    auto grad_a = grad.grad_a;
    auto lu_a = LU_decomp(grad_a*grad_a_t);
    // auto p_inv_a = pinv(grad_a_t);
    // print(grad_a);
    print(lu_a);

    // auto p_inv_a_array = p_inv_a.col(1);S
    // auto lambda_tmp = dot(p_inv_a_array, (grad.grad_f + H*x0));
    // for(u32 i = 0; i < grad.arr_act.size(); i++) {
    //     lambda[grad.arr_act[i]] = lambda_tmp; //pinv(grad.grad_a)*(grad.grad_f + H*x0);
    // }

    // SQP_nocedal(iter_lim, opt, x0, f0, g0, H, grad, vmax, lambda);

    return 0;
}