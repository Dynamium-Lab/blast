
#include <stdio.h>
#include "blast.hpp"
#include "blast_optional_utilities.hpp"
#include "nlopt.h"
#include "json.hpp"
#include <fstream>
#include <thread>
#include <mutex>

const int thread_count = 4;


struct OptimConfig {
    blast::u32     npts;
    blast::u32     nctrl;
    blast::u32     p;
    blast::u32     noptim;
    blast::u32     nshot;
    blast::real    m;
    blast::Matrix  task;
};

struct Result {
    blast::real     compute_time;
    OptimConfig     config;
    blast::Array    x;
};

std::vector<OptimConfig> config_list;
int config_index = 0;
std::mutex mut;

void eval_function() {
    using namespace blast;

    for(;;) {
        mut.lock();
        printf("hello from thread %d\n", std::this_thread::get_id());
        if (config_index >= config_list.size()) {
            mut.unlock();
            break;
        }
        auto config = config_list[config_index];
        config_index++;
        mut.unlock();

        const u32 npts      = config.npts;
        const u32 nctrl     = config.nctrl;
        const u32 p         = config.p;
        const u32 noptim    = config.noptim;

        Gen3_7DOF manip;
        manip.set_payload(config.m);
        PvaBspline pva(nctrl, npts, p, manip.joints);

        // prep optimization
        Optimisation optim {&manip, &config.task, &pva};
        u32     ncon = manip.ncon(npts);
        Array   con_tol(ncon, 0.001);
        Array   xtol(pva.xlen(), 0.000001);

        nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, pva.xlen());
        nlopt_result result;
        result = nlopt_add_inequality_mconstraint(o, ncon, cstr_manip, &optim, con_tol.data); Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_min_objective(o, obj_time, &optim); Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_lower_bound(o, pva.xlen()-1, 0.1);  Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_upper_bound(o, pva.xlen()-1, 10.0); Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_ftol_abs(o, 0.001);                 Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_xtol_abs(o, xtol.data);             Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_maxtime(o, 5.0);                    Assert(result ==NLOPT_SUCCESS);
        result = nlopt_set_maxeval(o, 10000);                  Assert(result ==NLOPT_SUCCESS);

        auto start_time = get_tick_us();
        for (int iter = 0; iter < noptim; iter++) {
            auto T1 = get_tick_us();

            // random optimization vector
            auto x = config.nshot == 1 ?
                     guess_random(manip, pva) :
                     guess_shot_mean(manip, pva, config.task, config.nshot);

            // launch optimization
            double f;
            result = nlopt_optimize(o, x.data, &f);
            auto T2 = get_tick_us();
            double time = (T2-T1) / 1000.0;

            // check solution
            pva.compute_trajectory(x, config.task);
            auto max_con = array_max(manip.constraints(pva));
            bool is_valid = max_con < 0.01;


            // printf("code: %d, result: %f, time: %fms, %s.\n", result, x.back(), time, is_valid? "valid" : "not valid");
        }
        auto total_time = (get_tick_us() - start_time) / 1000.0;
        // printf("total time: %f\n", total_time);
    }
}

int main() {
    using namespace blast;

    const u32   nctrl_list[] {8, 10, 12, 14, 16};
    const u32   npts_list[]  {100, 128, 256, 512};
    const real  m_list[]     {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const u32   nshot_list[] {25, 50, 100, 200};
    Matrix      task_list[9];

    auto start_time = get_tick_us();

    nlohmann::json j_file;
    std::ifstream file_handle("task_binpick.json");
    file_handle >> j_file;
    u32 i = 0;
    for (auto& task : task_list) {
        task.resize(7, 6);
        auto task_tmp = j_file["task"][i++];
        auto& pi = task_tmp["pi"];
        auto& vi = task_tmp["vi"];
        auto& ai = task_tmp["ai"];
        auto& pf = task_tmp["pf"];
        auto& vf = task_tmp["vf"];
        auto& af = task_tmp["af"];
        for (u32 j = 0; j < 7; j++) {
            task(j, 0) = pi[j];
            task(j, 1) = vi[j];
            task(j, 2) = ai[j];
            task(j, 3) = pf[j];
            task(j, 4) = vf[j];
            task(j, 5) = af[j];
        }
    }

    for (auto npts : npts_list) {
        for (auto nctrl : nctrl_list) {
            for (auto m : m_list) {
                for (auto nshot : nshot_list) {
                    for (auto task : task_list) {
                        OptimConfig c;
                        c.npts = npts;
                        c.nctrl = nctrl;
                        c.p = 5;
                        c.m = m;
                        c.nshot = nshot;
                        c.noptim = 50;
                        c.task = task;
                        config_list.push_back(c);
                    }
                }
            }
        }
    }


    std::vector<std::thread> workers;

    for (int i =0; i < thread_count; i++) {
        workers.push_back(std::thread(eval_function));
    }


    for (auto& t : workers)
        t.join();

    auto end_time = get_tick_us();
    auto total_time = (end_time - start_time) / 1000.0;

    printf("Total computation time: %fms\n", total_time);


    return 0;
}
