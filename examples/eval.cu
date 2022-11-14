
#include <stdio.h>
#include "blast.hpp"
#include "blast_optional_utilities.hpp"
#include "nlopt.h"
#include <thread>




int main() {
    using namespace blast;
    const u32 npts = 100;
    const u32 nctrl = 8;
    const u32 p = 5;
    const u32 noptim = 50;

    Gen3_7DOF manip;
    PvaBspline pva(nctrl, npts, p, manip.joints);
    Matrix task(manip.joints, 6);
    {
        // init position
        auto tmp = task.col(0); Assert(tmp.is_alias);
        tmp = {1.0, -2.13274122, 3.07081139, -1.32741228, 0.0, -0.159265358, 2.03540569};
        // final position
        tmp = task.col(3); Assert(tmp.is_alias);
        tmp = {-2.28318530, -0.8495559, 1.238958, 1.106217098, 0.0354056994, -0.3097094374, 2.03540569};
        Assert(manip.validate_task(task));
    }

    // prep optimization
    Optimisation optim {&manip, &task, &pva};
    nlopt_result result;
    nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, pva.xlen());
    {
        const u32 ncon = (5 + 2 + 7*2)*npts;
        Array con_tol(ncon, 0.001);
        Array xtol(pva.xlen(), 0.000001);
        result = nlopt_add_inequality_mconstraint(o, ncon, cstr_manip, &optim, con_tol.data); Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_min_objective(o, obj_time, &optim); Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_lower_bound(o, pva.xlen()-1, 0.1);  Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_upper_bound(o, pva.xlen()-1, 10.0); Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_ftol_abs(o, 0.001);                 Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_xtol_abs(o, xtol.data);             Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_maxtime(o, 5.0);                    Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_maxeval(o, 10000);                  Assert(result ==NLOPT_SUCCESS);
    }

    auto start_time = get_tick_us();
    for (int iter = 0; iter < noptim; iter++) {
        auto T1 = get_tick_us();

        // random optimization vector
        // auto x = guess_shot_max(manip, pva, task, 50);
        // auto x = guess_shot_sum(manip, pva, task, 50);
        auto x = guess_random(manip, pva);

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
    auto total_time = (get_tick_us() - start_time) / 1000.0;
    printf("total time: %f\n", total_time);

    return 0;
}
