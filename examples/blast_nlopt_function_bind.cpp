#include <blast>

#include "nlopt.h"
#include "json.hpp"

#include <stdio.h>
#include <functional>
#include <fstream>
#include <iostream>
#include <thread>
#include <mutex>


const int thread_count = 8;

struct OptimConfig {
    blast::u32 npts;
    blast::u32 nctrl;
    blast::u32 p;
    blast::u32 noptim;
    blast::u32 nshot;
    blast::real m;
    blast::Matrix task;
    blast::u32 task_id;
    blast::World world;
};
std::vector<OptimConfig> config_list;
int config_index = 0;
std::mutex mut_config;

struct Result {
    bool success;
    bool success_false;
    blast::real compute_time;
    OptimConfig config;
    blast::Array x;
    blast::Array x0;
};
std::vector<Result> result_list;
std::mutex mut_result;

std::ofstream file_handle_result("./result_list.json");

std::ofstream file_handle_traj("./trajectory.csv");
blast::Trajectory traj_csv;
blast::real time_csv = 0;
std::mutex mut_csv;

void eval_function() {
    using namespace blast;

    for (;;) {
        mut_config.lock();
        if (config_index >= config_list.size()) {
            mut_config.unlock();
            break;
        }
        auto config = config_list[config_index];
        std::cout << "Worker " << std::this_thread::get_id() << ", working on unit " << config_index << " of " << config_list.size() << std::endl;
        config_index++;
        mut_config.unlock();

        // init
        const u32 npts = config.npts;
        const u32 nctrl = config.nctrl;
        const u32 p = config.p;
        const u32 noptim = config.noptim;
        std::vector<Result> tmp_result_list(noptim);

        Link6 manip;
        manip.set_payload(config.m);
        Bspline bspline(nctrl, npts, p, manip.joints);

        Link6 manip_more_points;
        manip_more_points.set_payload(config.m);


        // prep optimization
        Optimization<Link6> opt{&manip, &config.task, &bspline, &config.world};
        u32 ncon = manip.ncon(npts);
        const Array con_tol(ncon, 0.001);
        Array xtol(bspline.xlen(config.task), 0.000001);
        using namespace std::placeholders;
        // auto cstr_world = std::bind(Optimization<Link6>::constraints, _1, _2, _3, _4, _5, _6, &opt);

        nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, bspline.xlen(config.task));
        nlopt_result result;
        result = nlopt_add_inequality_mconstraint(o, ncon, cstr_world<Link6>, &opt, con_tol.data);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_min_objective(o, obj_time, &opt);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_lower_bound(o, bspline.xlen(config.task) - 1, 0.1);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_upper_bound(o, bspline.xlen(config.task) - 1, 60.0);
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
        for (u32 iter = 0; iter < noptim; iter++) {
            auto T1 = get_tick_us();
            // random optimization vector
            auto x = config.nshot == 1 ? guess_random(bspline, config.task) : guess_shot_mean(opt, config.nshot);
            tmp_result_list[iter].x0 = x;

            // launch optimization
            double f;
            result = nlopt_optimize(o, x.data, &f);
            auto T2 = get_tick_us();
            double time = (T2 - T1) / 1000.0;

            // check solution
            bspline.compute_trajectory(x, config.task);
            Array constraints_points(manip.ncon(bspline.traj.t.size));
            manip.internal_constraints(bspline.traj, constraints_points.data);
            auto max_con = max(constraints_points);
            bool is_valid = max_con < 0.01;
            // traj_csv = is_valid ? bspline.traj : traj_csv;

            // Calculated from traj time (dt = 1ms)
            int points_1ms = (u32)ceil(x.data[x.size-1] * 1000 + 1); // todo: is this good ?
            Bspline bspline_more_points(nctrl, points_1ms, p, manip_more_points.joints);

            bspline_more_points.compute_trajectory(x, config.task);
            Array constraints_more_points(manip.ncon(bspline_more_points.traj.t.size));
            manip_more_points.internal_constraints(bspline_more_points.traj, constraints_more_points.data);
            auto max_con_more_points = max(constraints_more_points);
            bool is_valid_more_points = max_con_more_points < 0.05;
            traj_csv = is_valid_more_points ? bspline_more_points.traj : traj_csv;

            tmp_result_list[iter].success = is_valid && is_valid_more_points;
            tmp_result_list[iter].success_false = is_valid && !is_valid_more_points;
            tmp_result_list[iter].compute_time = time;
            tmp_result_list[iter].config = config;
            tmp_result_list[iter].x = x;
        }
        mut_result.lock();
        result_list.insert(result_list.end(), tmp_result_list.begin(), tmp_result_list.end());
        mut_result.unlock();
        manip_more_points.dynamics(traj_csv);
        auto trajectory_torque = manip._efforts;

        // mut_csv.lock();
        // for (u32 i = 0; i < traj_csv.t.size; i++)
        //     file_handle_traj << traj_csv.t[i] + time_csv
        //                      << "," << traj_csv.pos(0, i) << "," << traj_csv.pos(1, i) << "," << traj_csv.pos(2, i) << "," << traj_csv.pos(3, i) << "," << traj_csv.pos(4, i) << "," << traj_csv.pos(5, i)
        //                      << "," << traj_csv.vel(0, i) << "," << traj_csv.vel(1, i) << "," << traj_csv.vel(2, i) << "," << traj_csv.vel(3, i) << "," << traj_csv.vel(4, i) << "," << traj_csv.vel(5, i)
        //                      << "," << traj_csv.acc(0, i) << "," << traj_csv.acc(1, i) << "," << traj_csv.acc(2, i) << "," << traj_csv.acc(3, i) << "," << traj_csv.acc(4, i) << "," << traj_csv.acc(5, i)
        //                      << "," << trajectory_torque(0, i) << "," << trajectory_torque(1, i) << "," << trajectory_torque(2, i) << "," << trajectory_torque(3, i) << "," << trajectory_torque(4, i) << "," << trajectory_torque(5, i)
        //                      << std::endl;
        // mut_csv.unlock();
        // time_csv += traj_csv.t[traj_csv.t.size - 1];

        auto total_time = (get_tick_us() - start_time) / 1000.0;
    }
}


int main() {
    using namespace blast;
    Link6 manip;

    // const real m_list[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    const real m_list[] {0};
    // const u32 nctrl_list[] {10};
    const u32 nctrl_list[] {8};
    // const u32 npts_list[] {50};
    const u32 npts_list[] {100};
    const u32 nshot_list[] {50};
    Matrix task_list[1];
    bool robots_is_OBB  = true; // false -> capsules
    bool demo_is_simple = true; // false -> complex

    auto start_time = get_tick_us();

    // read tasks
    /*
    // nlohmann::json j_file;
    // std::ifstream file_handle("./examples/task_binpick.json");

    // file_handle.clear();
    // file_handle.seekg(0);
    // file_handle >> j_file;
    // u32 i = 0;
    // for (auto &task : task_list) {
    //     task.resize(6, 6);
    //     auto task_tmp = j_file["task"][i++];
    //     auto &pi = task_tmp["pi"];
    //     auto &vi = task_tmp["vi"];
    //     auto &ai = task_tmp["ai"];
    //     auto &pf = task_tmp["pf"];
    //     auto &vf = task_tmp["vf"];
    //     auto &af = task_tmp["af"];
    //     for (u32 j = 0; j < manip.joints; j++) {
    //         task(j, 0) = pi[j].is_null() ? NAN_REAL : (real)pi[j];
    //         task(j, 1) = vi[j].is_null() ? NAN_REAL : (real)vi[j];
    //         task(j, 2) = ai[j].is_null() ? NAN_REAL : (real)ai[j];
    //         task(j, 3) = pf[j].is_null() ? NAN_REAL : (real)pf[j];
    //         task(j, 4) = vf[j].is_null() ? NAN_REAL : (real)vf[j];
    //         task(j, 5) = af[j].is_null() ? NAN_REAL : (real)af[j];
    //     }
    // }*/

    // Read tasks Demo 1
    Array pos_zero = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    Array pos_black_1 = deg2rad({-51.44, -56.88, 33.31, 186.09, 172.32, 126.22});
    Array pos_black_1_10cm = deg2rad({-59.08, -39.48, 66.98, 176.31, 158.49, 126.22});

    Array pos_white_4 = deg2rad({-60.01, -52.99, 63.76, 176.33, 151.74, 126.21});
    Array pos_white_4_10cm = deg2rad({-40.26, -52.16, 67.02, 182.93, 144.11, 142.73});

    Array pos_white_2 = deg2rad({-55.32, -50.63, 45.20, 179.97, 170.65, 126.22});
    Array pos_white_2_10cm = deg2rad({-40.53, -37.18, 66.71, 183.99, 155.30, 142.73});

    Array pos_black_5 = deg2rad({-43.47, -46.36, 78.22, 175.84, 144.41, 138.33});
    Array pos_black_5_10cm = deg2rad({-55.44, -62.33, 44.39, 179.98, 162.66, 126.22});

    Array pos_white_3 = deg2rad({-50.81, -68.03, 33.30, 186.10, 161.69, 126.21});
    Array pos_white_3_10cm = deg2rad({-52.59, -38.51, 97.11, 180.56, 134.22, 130.07});

    Array pos_black_6 = deg2rad({30.50, -12.64, 102.90, 172.29, 145.76, 195.22});
    Array pos_black_6_10cm = deg2rad({-52.59, -17.11, 97.12, 181.98, 148.05, 130.07});

    Array pos_white_dbox = deg2rad({13.80, -62.21, 13.09, 178.88, 181.38, 195.22});
    Array pos_black_dbox = deg2rad({-45.71, -29.32, 78.21, 180.71, 155.91, 138.34});

    // For every bloc, the sequence is bloc -> +10cm -> dbox -> next bloc's +10cm

    std::vector<Matrix> tasks;
    Matrix task_zero_start(6, 6);
    Matrix task_zero_finish(6, 6);
    Matrix task_b1_0(6, 6);
    Matrix task_b1_1(6, 6);
    Matrix task_b1_2(6, 6);
    Matrix task_b1_w4(6, 6);
    Matrix task_w4_0(6, 6);
    Matrix task_w4_1(6, 6);
    Matrix task_w4_2(6, 6);
    Matrix task_w4_w2(6, 6);
    Matrix task_w2_0(6, 6);
    Matrix task_w2_1(6, 6);
    Matrix task_w2_2(6, 6);
    Matrix task_w2_b5(6, 6);
    Matrix task_b5_0(6, 6);
    Matrix task_b5_1(6, 6);
    Matrix task_b5_2(6, 6);
    Matrix task_b5_w3(6, 6);
    Matrix task_w3_0(6, 6);
    Matrix task_w3_1(6, 6);
    Matrix task_w3_2(6, 6);
    Matrix task_w3_b6(6, 6);
    Matrix task_b6_0(6, 6);
    Matrix task_b6_1(6, 6);
    Matrix task_b6_2(6, 6);
    for (u32 i = 0; i < 6; i++) {
        // Zero to start
        task_zero_start(i, 0) = pos_zero[i];
        task_zero_start(i, 3) = pos_black_1_10cm[i];
        // Black 1
        task_b1_0(i, 0) = pos_black_1_10cm[i];
        task_b1_0(i, 3) = pos_black_1[i];
        task_b1_1(i, 0) = pos_black_1[i];
        task_b1_1(i, 3) = pos_black_1_10cm[i];
        task_b1_2(i, 0) = pos_black_1_10cm[i];
        task_b1_2(i, 3) = pos_black_dbox[i];
        // Move to next
        task_b1_w4(i, 0) = pos_black_dbox[i];
        task_b1_w4(i, 3) = pos_white_4_10cm[i];
        // White 4
        task_w4_0(i, 0) = pos_white_4_10cm[i];
        task_w4_0(i, 3) = pos_white_4[i];
        task_w4_1(i, 0) = pos_white_4[i];
        task_w4_1(i, 3) = pos_white_4_10cm[i];
        task_w4_2(i, 0) = pos_white_4_10cm[i];
        task_w4_2(i, 3) = pos_white_dbox[i];
        // Move to next
        task_w4_w2(i, 0) = pos_white_dbox[i];
        task_w4_w2(i, 3) = pos_white_2_10cm[i];
        // White 2
        task_w2_0(i, 0) = pos_white_2_10cm[i];
        task_w2_0(i, 3) = pos_white_2[i];
        task_w2_1(i, 0) = pos_white_2[i];
        task_w2_1(i, 3) = pos_white_2_10cm[i];
        task_w2_2(i, 0) = pos_white_2_10cm[i];
        task_w2_2(i, 3) = pos_white_dbox[i];
        // Move to next
        task_w2_b5(i, 0) = pos_white_dbox[i];
        task_w2_b5(i, 3) = pos_black_5_10cm[i];
        // Black 5
        task_b5_0(i, 0) = pos_black_5_10cm[i];
        task_b5_0(i, 3) = pos_black_5[i];
        task_b5_1(i, 0) = pos_black_5[i];
        task_b5_1(i, 3) = pos_black_5_10cm[i];
        task_b5_2(i, 0) = pos_black_5_10cm[i];
        task_b5_2(i, 3) = pos_black_dbox[i];
        // Move to next
        task_b5_w3(i, 0) = pos_black_dbox[i];
        task_b5_w3(i, 3) = pos_white_3_10cm[i];
        // White 3
        task_w3_0(i, 0) = pos_white_3_10cm[i];
        task_w3_0(i, 3) = pos_white_3[i];
        task_w3_1(i, 0) = pos_white_3[i];
        task_w3_1(i, 3) = pos_white_3_10cm[i];
        task_w3_2(i, 0) = pos_white_3_10cm[i];
        task_w3_2(i, 3) = pos_white_dbox[i];
        // Move to next
        task_w3_b6(i, 0) = pos_white_dbox[i];
        task_w3_b6(i, 3) = pos_black_6_10cm[i];
        // Black 6
        task_b6_0(i, 0) = pos_black_6_10cm[i];
        task_b6_0(i, 3) = pos_black_6[i];
        task_b6_1(i, 0) = pos_black_6[i];
        task_b6_1(i, 3) = pos_black_6_10cm[i];
        task_b6_2(i, 0) = pos_black_6_10cm[i];
        task_b6_2(i, 3) = pos_black_dbox[i];
        // Zero to finish
        task_zero_finish(i, 0) = pos_black_dbox[i];
        task_zero_finish(i, 3) = pos_zero[i];
    }
    // Sequence tasks
    tasks.push_back(task_zero_start);
    tasks.push_back(task_b1_0);
    tasks.push_back(task_b1_1);
    tasks.push_back(task_b1_2);
    tasks.push_back(task_b1_w4);
    tasks.push_back(task_w4_0);
    tasks.push_back(task_w4_1);
    tasks.push_back(task_w4_2);
    tasks.push_back(task_w4_w2);
    tasks.push_back(task_w2_0);
    tasks.push_back(task_w2_1);
    tasks.push_back(task_w2_2);
    tasks.push_back(task_w2_b5);
    tasks.push_back(task_b5_0);
    tasks.push_back(task_b5_1);
    tasks.push_back(task_b5_2);
    tasks.push_back(task_b5_w3);
    tasks.push_back(task_w3_0);
    tasks.push_back(task_w3_1);
    tasks.push_back(task_w3_2);
    tasks.push_back(task_w3_b6);
    tasks.push_back(task_b6_0);
    tasks.push_back(task_b6_1);
    tasks.push_back(task_b6_2);
    tasks.push_back(task_zero_finish);

    World world;
    // todo: populate world

    // Tables
    add_box({0.4826, 0.0237, -0.4372}, {2.0, 2.0, 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world);
    add_box({0.0, 0.0, -0.4091}, {0.1397, 0.1905, 0.4091}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); // link6 base

    // Computer screens & controlers
    add_box({-0.5654, -0.8145, 0.3248}, {0.381, 0.635, 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); // next to link6
    add_box({0.4506, -1.3479, 0.3248}, {0.381, 0.635, 0.381}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); // next to gen3lite

    // Add Obstacle
    add_box({0.6477, -0.1956, -0.0562}, {0.375/2, 0.01, 0.375/2}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); // vertical plate

    add_box({0, -0.75, 1.0}, {0.10, 0.15, 2.0}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); // cables

    if (robots_is_OBB) {
        // Robots as OBB
        add_box({0.6665, 1.1286, 0.0}, {0.508, 0.3175, 1.0457}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); // UR5
        add_box({1.3011, 0.0, -0.0435}, {0.1774, 0.8636, 1.1903}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); // Gen3 exagerated on the whole linear axis
    }
    else {
        // Robots as caps
        add_capsule({0.6665, 1.1286, 0.0}, {0.6665, 1.1286, 1.0457}, 0.3175, &world); // UR5
        add_capsule({1.3011, 0.0, -0.0435}, {1.3011, 0.0, 1.2301}, 0.1774, &world); // Gen3 exactly accross from link6
    }

    if (demo_is_simple) {
        // Bins for demo 1 (simple)
        add_box({0.4114, -0.3784, 0.01785}, {0.2125, 0.155, 0.075}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world);
        add_box({0.8716, -0.3784, 0.01785}, {0.2125, 0.155, 0.075}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world);
    }
    else {
        // Bins for demo 1 (complex)
        // bin 1
        add_box({0.4114, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //bottom
        add_box({0.2114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //front
        add_box({0.6114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //back
        add_box({0.4114, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //left
        add_box({0.4114, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //right

        // bin 2
        add_box({0.8716, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //bottom
        add_box({0.6716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //front
        add_box({1.0716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //back
        add_box({0.8716, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //left
        add_box({0.8716, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //right
    }

    for (auto npts : npts_list)
        for (auto nctrl : nctrl_list)
            for (auto m : m_list)
                for (auto nshot : nshot_list)
                    for (int i = 0; i < tasks.size(); i++) {
                        OptimConfig c;
                        c.npts = npts;
                        c.nctrl = nctrl;
                        c.p = 5;
                        c.m = m;
                        c.nshot = nshot;
                        c.noptim = 50;
                        c.task = tasks[i];
                        c.task_id = i;
                        c.world = world;
                        config_list.push_back(c);
                    }

    result_list.reserve(config_list.size() * config_list[0].noptim);

    std::vector<std::thread> workers;

    for (int i = 0; i < thread_count; i++)
        workers.push_back(std::thread(eval_function));

    /*for (;;) {

        if (config_index >= config_list.size()) {
            break;
        }
        auto config = config_list[config_index];
        config_index++;

        // init
        const u32 npts = config.npts;
        const u32 nctrl = config.nctrl;
        const u32 p = config.p;
        const u32 noptim = config.noptim;
        std::vector<Result> tmp_result_list(noptim);

        Link6 manip;
        manip.set_payload(config.m);
        Bspline bspline(nctrl, npts, p, manip.joints);

        Link6 manip_more_points;
        manip_more_points.set_payload(config.m);
        Bspline bspline_more_points(nctrl, 10000, p, manip_more_points.joints);

        // prep optimization
        Optimization<Link6> opt{&manip, &config.task, &bspline, &world};
        u32 ncon = manip.ncon(npts);
        Array con_tol(ncon, 0.001);
        Array xtol(bspline.xlen(config.task), 0.000001);

        std::cout << "Setting up Optimization" << std::endl;
        nlopt_opt o = nlopt_create(nlopt_algorithm::NLOPT_LD_SLSQP, bspline.xlen(config.task));
        nlopt_result result;
        result = nlopt_add_inequality_mconstraint(o, ncon, cstr_world<Link6>, &opt, con_tol.data);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_min_objective(o, obj_time, &opt);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_lower_bound(o, bspline.xlen(config.task) - 1, 0.1);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_upper_bound(o, bspline.xlen(config.task) - 1, 60.0);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_ftol_abs(o, 0.001);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_xtol_abs(o, xtol.data);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_maxtime(o, 5.0);
        Assert(result == NLOPT_SUCCESS);
        result = nlopt_set_maxeval(o, 10000);
        Assert(result == NLOPT_SUCCESS);

        std::cout << "Start optimization iterations" << std::endl;
        auto start_time = get_tick_us();
        for (u32 iter = 0; iter < noptim; iter++) {
            auto T1 = get_tick_us();
            // random optimization vector
            auto x = config.nshot == 1 ? guess_random(bspline, config.task) : guess_shot_mean(opt, config.nshot);
            tmp_result_list[iter].x0 = x;

            // launch optimization
            double f;
            result = nlopt_optimize(o, x.data, &f);
            auto T2 = get_tick_us();
            double time = (T2 - T1) / 1000.0;

            // check solution
            bspline.compute_trajectory(x, config.task);
            Array constraints_points(manip.ncon(bspline.traj.t.size));
            manip.internal_constraints(bspline.traj, constraints_points.data);
            auto max_con = max(constraints_points);
            bool is_valid = max_con < 0.01;

            bspline_more_points.compute_trajectory(x, config.task);
            Array constraints_more_points(manip.ncon(bspline_more_points.traj.t.size));
            manip_more_points.internal_constraints(bspline_more_points.traj, constraints_more_points.data);
            auto max_con_more_points = max(constraints_more_points);
            bool is_valid_more_points = max_con_more_points < 0.05;

            tmp_result_list[iter].success = is_valid && is_valid_more_points;
            tmp_result_list[iter].success_false = is_valid && !is_valid_more_points;
            tmp_result_list[iter].compute_time = time;
            tmp_result_list[iter].config = config;
            tmp_result_list[iter].x = x;
        }
        result_list.insert(result_list.end(), tmp_result_list.begin(), tmp_result_list.end());

        auto total_time = (get_tick_us() - start_time) / 1000.0;
        std::cout << "Total time: " << total_time << std::endl;
    }*/

    for (auto &t : workers)
        t.join();

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
        j_config["world"]["capslist"]["size"] = result_config.world.capsules.size();
        for (u32 j = 0; j < result_config.world.capsules.size(); j++) {
            auto& j_caps = j_config["world"]["obblist"][j];
            j_caps["c"]["x"] = result_config.world.boxes[j].c.x;
            j_caps["c"]["y"] = result_config.world.boxes[j].c.y;
            j_caps["c"]["z"] = result_config.world.boxes[j].c.z;
            j_caps["e"]["y"] = result_config.world.boxes[j].e.y;
            j_caps["e"]["z"] = result_config.world.boxes[j].e.z;
            j_caps["e"]["x"] = result_config.world.boxes[j].e.x;
            for (u32 k = 0; k < 9; k++)
                j_caps["R"]["data"] = result_config.world.boxes[j].R.data[k];
        }
        for(u32 j = 0; j < result_list[i].x.size; j++)
            j_current["x"][j] = result_list[i].x[j];
        for(u32 j = 0; j < result_list[i].x0.size; j++)
            j_current["x0"][j] = result_list[i].x0[j];
    }

    file_handle_result << j_result;


    return 0;
}