#include <stdio.h>
#include "blast.hpp"
#include "utilities/blast_utilities.hpp"
#include "nlopt.h"
#include "json.hpp"
#include <fstream>

struct OptimConfig {
    blast::u32 npts;
    blast::u32 nctrl;
    blast::u32 p;
    blast::u32 noptim;
    blast::u32 nshot;
    blast::real m;
    blast::u32 task_id;
    blast::Matrix task;
};
std::vector<OptimConfig> config_list;
int config_index = 0;

struct Result {
    bool success;
    bool success_false;
    blast::real compute_time;
    OptimConfig config;
    blast::Array x;
    blast::Array x0;
};
std::vector<Result> result_list;
int result_index = 0;

int main() {
    using namespace blast;

    const real m_list[] {0};
    const u32 nctrl_list[] {10};
    const u32 npts_list[] {50};
    const u32 nshot_list[] {50};
    Matrix task_list[6];

    auto start_time = get_tick_us();

    // read tasks
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
        auto &af = task_tmp["af"];
        for (u32 j = 0; j < 7; j++) {
            task(j, 0) = pi[j].is_null() ? NAN_REAL : pi[j];
            task(j, 1) = vi[j].is_null() ? NAN_REAL : vi[j];
            task(j, 2) = ai[i].is_null() ? NAN_REAL : ai[j];
            task(j, 3) = pf[j].is_null() ? NAN_REAL : pf[j];
            task(j, 4) = vf[j].is_null() ? NAN_REAL : vf[j];
            task(j, 5) = af[i].is_null() ? NAN_REAL : af[j];
        }
    }

    for (auto npts : npts_list)
        for (auto nctrl : nctrl_list)
            for (auto m : m_list)
                for (auto nshot : nshot_list)
                    for (int i = 0; i < sizeof(task_list)/sizeof(Matrix); i++) {
                        OptimConfig c;
                        c.npts = npts;
                        c.nctrl = nctrl;
                        c.p = 5;
                        c.m = m;
                        c.nshot = nshot;
                        c.noptim = 50;
                        c.task = task_list[i];
                        c.task_id = i;
                        config_list.push_back(c);
                    }

    result_list.reserve(config_list.size() * config_list[0].noptim);

    objlist world;
    add_OBB({0.25, 0.2, 0.05}, {0.05, 0.05, 0.05}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world);

    for (;;) {
        if (config_index >= config_list.size()) {
            break;
        }
        auto config = config_list[config_index];
        config_index++;

        const u32 npts = config.npts;
        const u32 nctrl = config.nctrl;
        const u32 p = config.p;
        const u32 noptim = config.noptim;
        std::vector<Result> tmp_result_list(noptim);

        Gen3_7DOF manip;
        manip.set_payload_without_gripper(config.m);
        Bspline bspline(nctrl, npts, p, manip.joints);

        Gen3_7DOF manip_more_points;
        manip_more_points.set_payload_without_gripper(config.m);
        Bspline bspline_more_points(nctrl, 10000, p, manip_more_points.joints);

        // prep optimization
        Optimisation optim{&manip, &config.task, &bspline, &world};
        u32 ncon = manip.ncon(npts) + optim.n_collision_constraints; // todo: good ?
        Array con_tol(ncon, 0.001);
        const int xlen = bspline.xlen(config.task);
        Array xtol(xlen, 0.000001);

        nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, xlen);
        nlopt_result result;
        result = nlopt_add_inequality_mconstraint(o, ncon, cstr_world_gen3, &optim, con_tol.data);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_min_objective(o, obj_time, &optim);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_lower_bound(o, xlen - 1, 0.1);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_upper_bound(o, xlen - 1, 60.0);
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
            auto x = config.nshot == 1 ? guess_random(manip, bspline, config.task) : guess_shot_mean(manip, bspline, config.task, config.nshot);
            tmp_result_list[iter].x0 = x;

            // launch optimization
            double f;
            result = nlopt_optimize(o, x.data, &f);
            auto T2 = get_tick_us();
            double time = (T2 - T1) / 1000.0;

            // check solution
            bspline.compute_trajectory(x, config.task);
            Array const_result(ncon);
            cstr_world_gen3(ncon, const_result.data, xlen, x.data, NULL, &optim);
            auto max_con = array_max(const_result);
            bool is_valid = max_con < 0.01;

            bspline_more_points.compute_trajectory(x, config.task);
            auto max_con_more_points = array_max(manip_more_points.constraints(bspline_more_points.traj)); // todo: check obstacle avoidance with more points ?
            bool is_valid_more_points = max_con_more_points < 0.05;

            tmp_result_list[iter].success = is_valid && is_valid_more_points;
            tmp_result_list[iter].success_false = is_valid && !is_valid_more_points;
            tmp_result_list[iter].compute_time = time;
            tmp_result_list[iter].config = config;
            tmp_result_list[iter].x = x;

        }
        result_list.insert(result_list.end(), tmp_result_list.begin(), tmp_result_list.end());

        auto total_time = (get_tick_us() - start_time) / 1000.0;
        // printf("total time: %f\n", total_time);
    }

    auto end_time = get_tick_us();
    auto total_time = (end_time - start_time) / 1000.0;

    printf("Total computation time: %fms\n", total_time);

    nlohmann::json j_result;
    for(u32 i = 0; i < result_list.size(); i++) {
        auto& j_current = j_result["results"][i];
        j_current["success"] = result_list[i].success;
        j_current["success_false"] = result_list[i].success_false;
        j_current["compute time"] = result_list[i].compute_time;
        auto& j_config = j_current["config"];
        auto& result_config = result_list[i].config;
        j_config["npts"] = result_config.npts;
        j_config["nctrl"] = result_config.nctrl;
        j_config["p"] = result_config.p;
        j_config["noptim"] = result_config.noptim;
        j_config["nshot"] = result_config.nshot;
        j_config["m"] = result_config.m;
        j_config["task_id"] = result_config.task_id;
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