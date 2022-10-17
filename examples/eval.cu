
#include <stdio.h>
#include "blast.hpp"
#include "blast_optional_utilities.hpp"
#include "nlopt.h"





// todo: move out of GUI
struct Optimisation {
    blast::Manipulator*    manip;
    blast::Matrix*         task;
    blast::Pva*            pva;
};

double objective(unsigned n, const double* x, double* grad, void*) {
    using namespace blast;
    double result = x[n-1];
    if (grad) {
        for(u32 i = 0; i < n-1; i++)
            grad[i] = 0;
        grad[n-1] = 1;
    }
    return result;
}

// Constraints function
// m = number of constraints
// n = number of optim var (xlen)
void constraints(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data) {
    using namespace blast;
    Optimisation* opt = (Optimisation*)f_data;
    Array xv(n);
    std::copy_n(x, n, xv.data);

    // compute
    opt->pva->compute_trajectory(xv, *opt->task);
    opt->manip->constraints(*opt->pva, result);

    if (grad) {
        const real eps = 1e-5;
        Array x_plus(n);
        Array r_plus(m);
        for (u32 j = 0; j < n; j++) {
            memcpy(x_plus.data, x, n * sizeof(real));
            x_plus[j] += eps;

            opt->pva->compute_trajectory(x_plus, *opt->task);
            opt->manip->constraints(*opt->pva, r_plus.data);

            for (u32 i = 0; i < m; i++)
                grad[i*n + j] = (r_plus[i]-result[i])/eps;
        }
    }
}


int main() {
    using namespace blast;
    const u32 npts = 100;
    const u32 nctrl = 8;
    const u32 p = 5;
    const u32 noptim = 50;

    // setup manipulator
    Gen3_7DOF manip;

    // setup trajectory
    PvaBspline pva(nctrl, npts, p, manip.joints);

    //--- setup task
    Matrix task(manip.joints, 6);
    {
        // init position
        auto tmp = task.col(0);
        Assert(tmp.is_alias);
        tmp = {1.0, -2.13274122871835, 3.07081139897158, -1.32741228718346, 0.0, -0.159265358979322, 2.03540569948580};
        // final position
        tmp = task.col(3);
        Assert(tmp.is_alias);
        tmp = {-2.28318530717959, -0.849555921538759, 1.23895832717571, 1.10621709845737, 0.0354056994857963, -0.309709437440564, 2.03540569948580};
        // validate
        if (!manip.validate_task(task)) {
            printf("The required task is NOT valid\n");
            return -1;
        }
    }

    // prep optimization data
    Optimisation optim;
    {
        optim.manip = &manip;
        optim.pva = &pva;
        optim.task = &task;
    }

    // prep optimization
    nlopt_result result;
    // nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LN_COBYLA, pva.xlen());
    nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, pva.xlen());
    {
        // objective function
        result = nlopt_set_min_objective(o, objective, &optim);
        Assert(result ==NLOPT_SUCCESS);
        // nonlinear constraints
        const u32 ncon = (5 + 2 + 7*2)*npts;
        Array con_tol(ncon);
        constant(con_tol, 0.001);
        result = nlopt_add_inequality_mconstraint(o, ncon, constraints, &optim, con_tol.data);
        Assert(result ==NLOPT_SUCCESS);
        // lower bounds
        result = nlopt_set_lower_bound(o, pva.xlen()-1, 0.1);
        Assert(result ==NLOPT_SUCCESS);
        // upper bounds
        result = nlopt_set_upper_bound(o, pva.xlen()-1, 10.0);
        Assert(result ==NLOPT_SUCCESS);
        // ftol
        result = nlopt_set_ftol_abs(o, 0.001);
        Assert(result ==NLOPT_SUCCESS);
        // xtol
        Array xtol(pva.xlen(), 0.000001);
        result = nlopt_set_xtol_abs(o, xtol.data);
        Assert(result ==NLOPT_SUCCESS);
        // max time
        result = nlopt_set_maxtime(o, 5.0);
        Assert(result ==NLOPT_SUCCESS);
        // max eval
        result = nlopt_set_maxeval(o, 10000);
        Assert(result ==NLOPT_SUCCESS);
    }

    for (int iter = 0; iter < noptim; iter++) {
        auto T1 = get_tick_us();

        // random optimization vector
        Array x(pva.xlen());
        {
            fill_random(x, 1);
            x.back() = abs(x.back()) * 5 + 0.1;
        }

        // launch optimization
        double f;
        result = nlopt_optimize(o, x.data, &f);
        auto T2 = get_tick_us();
        double time = (T2-T1) / 1000.0;

        // check solution
        pva.compute_trajectory(x, task);
        auto max_con = array_max(manip.constraints(pva));
        bool is_valid = max_con < 0.01;
        printf("code: %d, result: %f, time: %fms, %s.\n", result, x.back(), time, is_valid? "valid" : "not valid");
    }

    return 0;
}
