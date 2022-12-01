
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
    blast::u32 npts;
    blast::u32 nctrl;
    blast::u32 p;
    blast::u32 noptim;
    blast::u32 nshot;
    blast::real m;
    blast::Matrix task;
};
std::vector<OptimConfig> config_list;
int config_index = 0;
std::mutex mut_config;

struct Result {
    blast::real compute_time;
    OptimConfig config;
    blast::Array x;
    blast::Array x0;
};
std::vector<Result> result_list;
int result_index = 0;
std::mutex mut_result;

void eval_function() {
    using namespace blast;

    for (;;) {
        mut_config.lock();
        printf("hello from thread %d\n", std::this_thread::get_id());
        if (config_index >= config_list.size()) {
            mut_config.unlock();
            break;
        }
        auto config = config_list[config_index];
        config_index++;
        mut_config.unlock();

        const u32 npts = config.npts;
        const u32 nctrl = config.nctrl;
        const u32 p = config.p;
        const u32 noptim = config.noptim;
        std::vector<Result> tmp_result_list(noptim);

        Gen3_7DOF manip;
        manip.set_payload_without_gripper(config.m);
        PvaBspline pva(nctrl, npts, p, manip.joints);

        // prep optimization
        Optimisation optim{&manip, &config.task, &pva};
        u32 ncon = manip.ncon(npts);
        Array con_tol(ncon, 0.001);
        Array xtol(pva.xlen(config.task), 0.000001);

        nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, pva.xlen(config.task));
        nlopt_result result;
        result = nlopt_add_inequality_mconstraint(o, ncon, cstr_manip, &optim, con_tol.data);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_min_objective(o, obj_time, &optim);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_lower_bound(o, pva.xlen(config.task) - 1, 0.1);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_upper_bound(o, pva.xlen(config.task) - 1, 10.0);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_ftol_abs(o, 0.001);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_xtol_abs(o, xtol.data);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_maxtime(o, 5.0);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_maxeval(o, 10000);
        Assert(result == NLOPT_SUCCESS);

        auto start_time = get_tick_us();
        for (int iter = 0; iter < noptim; iter++) {
            auto T1 = get_tick_us();

            // random optimization vector
            auto x = config.nshot == 1 ? guess_random(manip, pva, config.task) : guess_shot_mean(manip, pva, config.task, config.nshot);
            tmp_result_list[iter].x0 = x;

            // launch optimization
            double f;
            result = nlopt_optimize(o, x.data, &f);
            auto T2 = get_tick_us();
            double time = (T2 - T1) / 1000.0;

            // check solution
            pva.compute_trajectory(x, config.task);
            auto max_con = array_max(manip.constraints(pva));
            bool is_valid = max_con < 0.01;

            // printf("code: %d, result: %f, time: %fms, %s.\n", result, x.back(), time, is_valid? "valid" : "not valid");
            tmp_result_list[iter].compute_time = time;
            tmp_result_list[iter].config = config;
            tmp_result_list[iter].x = x;

        }
        mut_result.lock();
        result_list.insert(result_list.end(), tmp_result_list.begin(), tmp_result_list.end());
        mut_result.unlock();

        auto total_time = (get_tick_us() - start_time) / 1000.0;
        // printf("total time: %f\n", total_time);
    }
}

int main() {
    using namespace blast;

    // const u32 nctrl_list[] {8, 10, 12, 14, 16};
    // const u32 npts_list[] {100, 128, 256, 512};
    // const real m_list[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // const u32 nshot_list[] {25, 50, 100, 200};
    // Matrix task_list[9];
    const u32 nctrl_list[] {8};
    const u32 npts_list[] {100};
    const real m_list[] {2};
    const u32 nshot_list[] {50};
    Matrix task_list[3];

    auto start_time = get_tick_us();

    nlohmann::json j_file;
    std::ifstream file_handle("task_binpick.json");
    file_handle >> j_file;
    u32 i = 0;
    for (auto &task : task_list) {
        task.resize(7, 6);
        auto task_tmp = j_file["task"][i++];
        auto &pi = task_tmp["pi"];
        auto &vi = task_tmp["vi"];
        auto &ai = task_tmp["ai"];
        auto &pf = task_tmp["pf"];
        auto &vf = task_tmp["vf"];
        // auto &af = task_tmp["af"];
        for (u32 j = 0; j < 7; j++) {
            task(j, 0) = pi[j];
            task(j, 1) = vi[j];
            task(j, 2) = ai[i];
            task(j, 3) = pf[j];
            task(j, 4) = vf[j];
            task(j, 5) = NAN_REAL;
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

    result_list.reserve(config_list.size() * config_list[0].noptim);

    std::vector<std::thread> workers;

    for (int i = 0; i < thread_count; i++) {
        workers.push_back(std::thread(eval_function));
    }

    for (auto &t : workers)
        t.join();

    auto end_time = get_tick_us();
    auto total_time = (end_time - start_time) / 1000.0;

    printf("Total computation time: %fms\n", total_time);

    nlohmann::json j_result;
    for(u32 i = 0; i < result_list.size(); i++) {
        auto& j_current = j_result["results"][i];
        j_current["compute time"] = result_list[i].compute_time;
        auto& j_config = j_current["config"];
        auto& result_config = result_list[i].config;
        j_config["npts"] = result_config.npts;
        j_config["nctrl"] = result_config.nctrl;
        j_config["p"] = result_config.p;
        j_config["noptim"] = result_config.noptim;
        j_config["nshot"] = result_config.nshot;
        j_config["m"] = result_config.m;
        j_config["task"]["rows"] = result_config.task.rows;
        j_config["task"]["cols"] = result_config.task.cols;
        j_config["task"]["size"] = result_config.task.size;
        for(u32 j = 0; j < result_config.task.size; j++)
            j_config["task"]["data"][j] = result_config.task.data[j];
        for(u32 j = 0; j < result_list[i].x.size; j++)
            j_current["x"][j] = result_list[i].x[j];
        for(u32 j = 0; j < result_list[i].x0.size; j++)
            j_current["x0"][j] = result_list[i].x0[j];
    }

    std::ofstream file_handle_result("result_list.json");
    file_handle_result << j_result;


    return 0;
}
