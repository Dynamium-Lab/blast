#include <stdio.h>
#include "blast.hpp"
#include "blast_optimization.hpp"
#include "blast_optional_utilities.hpp"
#include "nlopt.h"

using namespace blast;

struct gradients {

    Array           grad_f;
    Matrix          grad_g;
    Matrix          grad_a;
    Matrix          grad_na;
    vector<int>    arr_act;
    vector<int>    arr_nact;
    Array           cstr_a;
    Array           cstr_na;
};

struct opt_result {
    real       fval;
    int    iter_val = 0;
    Array       x_opt;
    Array       x_init;
};


void rosensuzuki_fun(Array& x, real& f, Array& g) {
    f = x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] - 5*x[0] - 5*x[1] - 21*x[2] + 7*x[3];
    if(g.size != 0) {

        g[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] - x[1] + x[2] - x[3] - 8;
        g[1] = x[0]*x[0] + 2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] - x[0] - x[3] - 10;
        g[2] = 2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] - x[1] - x[3] - 5;
    }
}

// Hessien approximation using BFGS
//
// H = hessien
// s(1, nv) = x - x0; (in main loop)
// q(1,nv) = (grad_f_t + lambda*grad_g_t) - (grad_f_t_old + lambda*grad_g_t_old)
// v(nv) = grad_g*g - grad_g_old*g_old
Matrix BFGS_hess(Matrix& H, Array& s, Array& q, Array& v) {
    Matrix H_new(H.rows, H.cols);
    auto H_t = transpose(H);
    auto q_t = transpose(q);
    auto s_t = transpose(s);
    auto pos_def = dot(q, s);

    if (pos_def <= 0) {
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
    H_new = H + (q*q_t)/q_s - (H*(s*s_t)*H_t)/s_Hs;
    return H_new;
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
void active(Array& g, vector<int>& arr_act, vector<int>& arr_nact, double& vmax) {
    arr_act.clear();
    arr_nact.clear();
    vmax = array_max(g);

    const real eps = vmax - 0.05*abs(vmax);
    for (u32 j = 0; j < g.size; j++) {
        if (g[j] >= eps)
            arr_act.push_back(j);
        else
            arr_nact.push_back(j);
    }
}

// Active constraints gradients
//
// m = number of constraints
// n = number of variables
gradients cstr_manip_active(Array& x, real& f, Array& g, real& vmax) {
    gradients grad_r;
    vector<int> arr_act;
    vector<int> arr_nact;
    active(g, arr_act, arr_nact, vmax);
    unsigned m_act = arr_act.size();
    unsigned m_nact = g.size - m_act;

    // Optimisation *opt = (Optimisation *)f_data;
    // internal_cstr_manip_single(m, result, n, x, opt);
    Array cst(x.size);
    constant(cst, 1e-6);
    auto eps = 0.001*x + cst;
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
        x_plus[j] += eps[j];
        // internal_cstr_manip_single(m, g_plus.data, n, x_plus.data, opt);
        rosensuzuki_fun(x_plus, f_plus, g_plus);
        grad_f[j] = (f_plus - f) / eps[j];
        for (u32 i = 0; i < g.size; i++)
            grad_g(j, i) = (g_plus[i] - g[i]) / eps[j];
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

    grad_r.grad_f = grad_f;
    grad_r.grad_g = grad_g;
    grad_r.grad_a = grad_a;
    grad_r.grad_na = grad_na;
    grad_r.arr_act = arr_act;
    grad_r.arr_nact = arr_nact;
    grad_r.cstr_a = cstr_a;
    grad_r.cstr_na = cstr_na;

    return grad_r;
}

host_fn double obj_rs(unsigned n, const double* x, double* grad, void*) {
    using namespace blast;
    double result;

    Array g;
    Array g_plus;
    Array x_alias ((real*)x, n);
    Assert(x_alias.is_alias);

    rosensuzuki_fun(x_alias, result, g);

    Array cst(x_alias.size);
    constant(cst, 1e-6);
    auto eps = 0.001*x_alias + cst;

    Array x_plus(x_alias.size);
    real f_plus;

    if (grad) {
        for (u32 j = 0; j < x_alias.size; j++) {

            x_plus = x_alias;
            x_plus[j] += eps[j];
            // internal_cstr_manip_single(m, g_plus.data, n, x_plus.data, opt);
            rosensuzuki_fun(x_plus, f_plus, g_plus);
            grad[j] = (f_plus - result) / eps[j];
        }
    }
    return result;
}

// Simple manipulator defined constraints function.
//  The given manipulator's constraints function is called and no additionnal constraints are added.
//
//  m      = number of constraints
//  result = constraints are inserted here
//  n      = number of elements in optimization vector
//  x      = optimization vector
//  grad   = if not NULL, gradient is inserted here
//  data   = is cast to Optimization struct
//
// Returns nothing
host_fn void cstr_manip_rs(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    // Optimisation* opt = (Optimisation*)f_data;

    // internal_cstr_manip_single(m, result, n, x, opt);
    // f_data = x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] - 5*x[0] - 5*x[1] - 21*x[2] + 7*x[3];
    result[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[0] - x[1] + x[2] - x[3] - 8;
    result[1] = x[0]*x[0] + 2*x[1]*x[1] + x[2]*x[2] + 2*x[3]*x[3] - x[0] - x[3] - 10;
    result[2] = 2*x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + 2*x[0] - x[1] - x[3] - 5;

    if (grad) {
        const real eps = 1e-5;
        Array x_plus(n);
        Array r_plus(m);
        // todo: parallel?
        for (u32 j = 0; j < n; j++) {
            memcpy(x_plus.data, x, n * sizeof(real));
            x_plus[j] += eps;
            // internal_cstr_manip_single(m, r_plus.data, n, x_plus.data, opt);
            r_plus[0] = x_plus[0]*x_plus[0] + x_plus[1]*x_plus[1] + x_plus[2]*x_plus[2] + x_plus[3]*x_plus[3] + x_plus[0] - x_plus[1] + x_plus[2] - x_plus[3] - 8;
            r_plus[1] = x_plus[0]*x_plus[0] + 2*x_plus[1]*x_plus[1] + x_plus[2]*x_plus[2] + 2*x_plus[3]*x_plus[3] - x_plus[0] - x_plus[3] - 10;
            r_plus[2] = 2*x_plus[0]*x_plus[0] + x_plus[1]*x_plus[1] + x_plus[2]*x_plus[2] + 2*x_plus[0] - x_plus[1] - x_plus[3] - 5;
            for (u32 i = 0; i < m; i++)
                grad[i*n + j] = (r_plus[i]-result[i])/eps;
        }
    }
}

// qp active set
Array qp_active_set(Matrix& H, Array& grad_f, Matrix& A, Array& b) {
    unsigned iter_lim = 1000;
    auto A_t = transpose(A);
    auto x0 = pinv_svd(A_t)*b;
    auto x = x0;

    // auto H_eye = eye(H.cols);
    // if( H == H_eye) {
    //     x = -grad_f;
    //     return x;
    // }

    // Active subset initialization
    Array Ax_b(b.size);
    Ax_b = A_t*x - b;
    real vmax;

    vector<int> arr_act;
    vector<int> arr_nact;
    active(Ax_b, arr_act, arr_nact, vmax);
    unsigned m_act = arr_act.size();
    unsigned m_nact = b.size - m_act;
    u32 i = 0;
    Matrix A_a(x.size, b.size);
    Matrix A_na(x.size, b.size);
    Array b_a(b.size);
    Array b_na(b.size);

    // Main loop qp
    for(i = 0; i< iter_lim; i++) {
        // Initialize A_a, A_na, b_a, b_na active and non active sets
        A_a.resize(x.size, m_act);
        A_na.resize(x.size, m_nact);
        b_a.resize(m_act);
        b_na.resize(m_nact);
        for (u32 j = 0; j < x.size; j++) {
            for (u32 i = 0; i < m_act; i++) {
                A_a(j, i) = -A(j, arr_act[i]);
                b_a[i] = -b[arr_act[i]];
            }
            for (u32 i = 0; i < m_nact; i++) {
                A_na(j, i) = -A(j, arr_nact[i]);
                b_na[i] = -b[arr_nact[i]];
            }
        }

        // Schur-complement method to find direction  p_k
        auto A_a_t = transpose(A_a);
        auto A_a_Hinv   = A_a_t*pinv_svd(H);
        auto g          =  grad_f + H*x;
        auto h          = A_a_t*x - b_a;
        auto lambda     = pinv_svd(A_a_Hinv*A_a)*(A_a_Hinv*g - h);
        // auto test = A_a*transpose(lambda);
        auto p_k        = pinv_svd(H)*(A_a*lambda - g);

        // if p_k == 0
        if (norm(p_k) < 1e-6) {
            auto lambda_bar = pinv_svd(A_a_t)*b_a;
            real lambda_min;
            auto lambda_min_id = array_min_id(lambda_bar, lambda_min);
            if(lambda_min >= 0) {

                arr_act.clear();
                arr_nact.clear();
                return x;
            }
            else {  // remove lambda_min_id from arr_act
                arr_nact.push_back(arr_act[lambda_min_id]);
                arr_act.erase(arr_act.begin() + lambda_min_id);
                if (arr_act.empty() == false)
                    sort(arr_act.begin(), arr_act.end());
                if (arr_nact.empty() == false)
                    sort(arr_nact.begin(), arr_nact.end());
            }
        }
        else {
            real alpha_k;
            // real alpha_min;
            unsigned alpha_min_id;
            if (A_na.size == 0)
                alpha_k  = 1;
            else {
                Array alpha_tmp(A_na.cols);
                real alpha_min;
                for (u32 i = 0; i < A_na.cols; i++) {
                    auto a_na = A_na.col(i);
                    alpha_tmp[i] = (b_na[i] - dot(a_na, x))/dot(a_na, p_k);
                }
                alpha_min_id = array_min_id(alpha_tmp, alpha_min);
                alpha_k = clamp(alpha_min, 0, 1);
            }
            x = x + alpha_k*p_k;
            if (alpha_k < 1) {
                arr_act.push_back(arr_nact[alpha_min_id]);
                arr_nact.erase(arr_nact.begin() + alpha_min_id);
                if (arr_act.empty() == false)
                    sort(arr_act.begin(), arr_act.end());
                if (arr_nact.empty() == false)
                    sort(arr_nact.begin(), arr_nact.end());
            }
            else {
                arr_act.clear();
                arr_nact.clear();
                return x;
            }
        }
        m_act = arr_act.size();
        m_nact = b.size - m_act;
    }
    arr_act.clear();
    arr_nact.clear();
    return x;
}

real step_size(Matrix& H, gradients& grad, Array& d, Array& x, real& f, Array& g, real& vmax) {
    auto negative_grad_g = - grad.grad_g;
    auto T1 = get_tick_us();
    d = qp_active_set(H, grad.grad_f, negative_grad_g, g);
    auto T2 = get_tick_us();
    auto time = (T2 - T1)/1000.0;
    // printf("%f\n", time);
    auto d_k = d;

    auto x_new = x + d_k;
    real f_new;
    Array g_new(g.size);
    gradients grad_new;
    rosensuzuki_fun(x_new, f_new, g_new);
    grad_new = cstr_manip_active(x_new, f_new, g_new, vmax);
    auto g_new_max = array_max(g_new);
    auto g_max = array_max(g);

    auto grad_na = grad_new.grad_na;
    auto g_na = grad_new.cstr_na;

    real a_bar = 1;
    if(g_new_max <= g_max)
        a_bar = 1;
    else {
        auto a_d = array_min(transpose(grad_na)*d_k);
        if (a_d < 0) {
            Array a_tmp(grad_na.cols);
            unsigned m_tmp_size;
            for (u32 i = 0; i < grad_na.cols; i++) {
                auto a_na = grad_na.col(i);
                a_tmp[i] = (g_na[i] - dot(a_na, x))/dot(a_na, d);
            }
            m_tmp_size = sizeof(a_tmp) / sizeof(a_tmp[0]);
            for (u32 i = 0; i < m_tmp_size; i++) {
                if (a_tmp[i] < 0)
                    a_tmp[i] = INF_REAL;
            }
            a_bar = array_min(a_tmp);
        }
        else {
            a_bar = 1;
        }


    }
    a_bar = clamp(a_bar, 0, 1);

    Array r(g.size);
    Array r_new(g.size);
    auto grad_g_t = transpose(grad.grad_g);
    auto grad_new_g_t = transpose(grad_new.grad_g);
    for(u32 i = 0; i < g.size; i++) {
        r[i] = norm(grad.grad_f)/norm(grad.grad_g.col(i));
        r_new[i] = norm(grad_new.grad_f)/norm(grad_new.grad_g.col(i));
    }

    Array lambda(g.size);
    Array lambda_new(g.size);
    auto grad_a_inv = pinv_svd(grad.grad_a);

    vector<int> arr_act;
    vector<int> arr_nact;
    active(g, arr_act, arr_nact, vmax);
    unsigned m_act = arr_act.size();
    // unsigned m_nact = g.size - m_act;

    for (u32 i =0; i < m_act; i++) {
        lambda[arr_act[i]] = dot(grad_a_inv.row(i), (grad.grad_f + H*x));
    }

    const real w = 0.5;
    const real coeff = 0.5;
    unsigned iter_lim = 100;
    real a = a_bar;
    for (u32 iter = 0; iter < iter_lim; iter++) {
        x_new = x + a*d_k;
        rosensuzuki_fun(x_new, f_new, g_new);
        grad_new = cstr_manip_active(x_new, f_new, g_new, vmax);
        // auto grad_a_t_new = transpose(grad_new.grad_a);
        auto grad_a_inv_new = pinv_svd(transpose(grad_new.grad_a));

        active(g_new, arr_act, arr_nact, vmax);
        m_act = arr_act.size();
        // unsigned m_nact = g.size - m_act;
        auto tst = H*x_new;

        zero(lambda_new);
        for (u32 i =0; i < m_act; i++) {
            auto g_row = grad_a_inv_new.col(i);
            lambda_new[arr_act[i]] = dot(grad_a_inv_new.col(i), (grad_new.grad_f + H*x_new));
        }
        Array g_array_max(g.size);
        Array g_array_max_new(g.size);
        Array r_g(g.size);
        Array r_g_new(g.size);
        for (u32 i = 0; i < g.size; i++) {
            // g_array_max[i] = 0.0 <= g[i] ? g[i]: 0.0;
            // g_array_max_new[i] = 0.0 <= g_new[i] ? g_new[i]: 0.0;
            g_array_max[i] = max(0.0, g[i]);
            g_array_max_new[i] = max(0.0, g_new[i]);

            r_g[i] = r[i]*g_array_max[i];
            r_g_new[i] = r_new[i]*g_array_max_new[i];

            r[i] = max(lambda[i], (r[i] + lambda[i])/2);
            r_new[i] = max(lambda_new[i], (r_new[i] + lambda_new[i])/2);
        }

        auto Psi = f + sum(r_g);
        auto Psi_new = f_new + sum(r_g_new);

        auto Psi_comp = Psi + coeff*a*dot(grad.grad_f, d_k);
        if (Psi_new <= Psi_comp) {
            break;
        }
        a *= w;
    }
    return a;
}

opt_result SQP_nocedal(unsigned iter_lim, Array& x, real& f, Array& g, Matrix& H, gradients& grad, real& vmax) {
    opt_result opt_r;
    unsigned iter;
    Array lambda(g.size);
    Matrix s(x.size, 1);
    Matrix q(x.size, 1);
    Array v(x.size);
    Array d(x.size);

    auto grad_a_t = transpose(grad.grad_a);
    auto grad_a = grad.grad_a;
    auto grad_a_inv = pinv_svd(grad_a);

    auto arr_act = grad.arr_act;
    auto arr_nact = grad.arr_nact;
    unsigned m_act = arr_act.size();

    for (u32 i =0; i < m_act; i++) {
        lambda[arr_act[i]] = dot(grad_a_inv.row(i), (grad.grad_f + H*x));
    }
    auto lambda_t = transpose(lambda);

    opt_r.x_init = x;
    for(iter = 0; iter < iter_lim; iter++) {
        // auto T1_step = get_tick_us();
        real a = step_size(H, grad, d, x, f, g, vmax);
        // auto T2_step = get_tick_us();
        // auto time_step = (T2_step - T1_step)/1000.0;
        // printf("time step %f\n", time_step);
        auto dn = norm(d);
        if (dn == 0) {
            // std::cout<<"dn == 0"<<std::endl;
            // std::cout<<dn<<std::endl;
            break;
        }

        auto x0 = x;
        x = x + a*d;
        auto f_old = f;
        auto g_old = g;
        auto grad_g_old = grad.grad_g;
        auto grad_f_old = grad.grad_f;
        rosensuzuki_fun(x, f, g);

        auto xtol = norm(x)*1e-12 + 1e-15;
        if (norm(x - x0) < xtol) {
            // std::cout<<"norm(x - x0)"<<std::endl;
            // std::cout<<norm(x - x0)<<std::endl;
            // std::cout<<xtol<<std::endl;
            // std::cout<<vmax<<std::endl;
            x = x0;
            f = f_old;
            break;
        }

        auto ftol = abs(f)*1e-6 + 1e-12;
        auto tolG = 1e-6;
        if (abs(f - f_old) < ftol && vmax < tolG) {
            // std::cout<<"abs(f - f_old)"<<std::endl;
            // std::cout<<abs(f - f_old)<<std::endl;
            break;
        }

        gradients grad_new;
        grad_new = cstr_manip_active(x, f, g, vmax);
        auto grad_g = grad_new.grad_g;
        auto grad_f = grad_new.grad_f;

        s = x - x0;
        s = transpose(s);
        q =  (transpose(grad_f) + lambda_t*transpose(grad_g)) - (transpose(grad_f_old) + lambda_t*transpose(grad_g_old));
        v = grad_g*g - grad_g_old*g_old;
        auto s_row = s.row(0);
        auto q_row = q.row(0);
        H = BFGS_hess(H, s_row, q_row, v);

        grad_a_t = transpose(grad_new.grad_a);
        zero(lambda);
        arr_act = grad_new.arr_act;
        m_act = arr_act.size();
        for (u32 i =0; i < m_act; i++) {
            lambda[arr_act[i]] = dot(pinv_svd(grad_new.grad_a).row(i), (grad_new.grad_f + H*x));
        }
        grad = grad_new;
    }
    opt_r.fval = f;
    opt_r.x_opt = x;
    opt_r.iter_val = iter;

    return opt_r;
}

int main() {
    using namespace blast;
    // init
    int iter_lim = 1000;
    // Rosensuzuki Problem
    unsigned nv = 4;
    unsigned ncon = 3;
    Array x0 = random_array(nv, 10);
    // Array x0(nv);
    // x0[0] = 0;
    // x0[1] = 1;
    // x0[2] = 2;
    // x0[3] = -1;
    // for (u32 i = 0; i < nv; i++)
    //     x0[i] = get_random();
    print(x0);
    Matrix H = eye(nv);
    double f0;
    Array g0(ncon);
    rosensuzuki_fun(x0, f0, g0);

    gradients  grad;
    double vmax;

    grad = cstr_manip_active(x0, f0, g0, vmax);

    auto T1_own = get_tick_us();
    auto opt = SQP_nocedal(iter_lim, x0, f0, g0, H, grad, vmax);
    auto T2_own = get_tick_us();
    auto time_own = (T2_own - T1_own)/1000.0;

    auto f_val = opt.fval;
    auto x_val = opt.x_opt;
    print(x_val);
    std::cout<<f_val<<std::endl;

    opt_result optim;
    nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, nv);
    nlopt_result result;
    Array con_tol(ncon, 0.001);
    Array xtol(nv, 0.000001);

    result = nlopt_add_inequality_mconstraint(o, ncon, cstr_manip_rs, &optim, con_tol.data);
    Assert(result == NLOPT_SUCCESS);
    result = nlopt_set_min_objective(o, obj_rs, &optim.fval);
    Assert(result == NLOPT_SUCCESS);
    result = nlopt_set_lower_bound(o, ncon, -100); //todo: ncon or nv ?
    Assert(result == NLOPT_SUCCESS);
    result = nlopt_set_upper_bound(o, ncon, 100); //todo: ncon or nv ?
    Assert(result == NLOPT_SUCCESS);
    result = nlopt_set_ftol_abs(o, 0.001);
    Assert(result == NLOPT_SUCCESS);
    result = nlopt_set_xtol_abs(o, xtol.data);
    Assert(result == NLOPT_SUCCESS);
    result = nlopt_set_maxtime(o, 5.0);
    Assert(result == NLOPT_SUCCESS);
    result = nlopt_set_maxeval(o, 10000);
    Assert(result == NLOPT_SUCCESS);

    double f;
    auto T1_nlopt = get_tick_us();
    result = nlopt_optimize(o, x0.data, &f);
    auto T2_nlopt = get_tick_us();
    auto time_nlopt = (T2_nlopt - T1_nlopt)/1000.0;

    auto frac_time = time_own / time_nlopt;

    printf("%f\n", frac_time);
    printf("%f\n", f);
    print(x0);

    return 0;
}