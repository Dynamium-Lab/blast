#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include "catch2/catch.hpp"
#include "blast"
#include "test_helper/test_helper.hpp"
#include "test_helper/test_functions.hpp"

using namespace blast;

struct BenchCase {
    u32 nctrl;
    u32 npoints;
    u32 p;
};

TEST_CASE("benchmark phase 1 gradient acceleration  (P,V,A,Tau)", "[Acceleration]") {
    using namespace blast;

    std::string input_file_initial_guess = "../../../examples/acceleration/files/link6_initial_guess.csv";

    // blast::begin_profile();
    // --- Define Manipulator ---
    auto link6 = get_generic_Link6();

    // --- Define Tasks ---
    u32 task_size = 7;
    std::vector<Array> positions(task_size + 1);
    positions[0] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // zero
    positions[1] = deg2rad({-40.4, -26.9, 83.6, 1.5, 20.0, -42.2});     // w1 + 10
    positions[2] = deg2rad({51.9, -13.6, 107.9, 3.6, 33.1, 51.2});      // wb1
    positions[3] = deg2rad({-33.9, -35.6, 70.9, -0.7, 16.9, -31.8});    // w3 + 10
    positions[4] = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});       // wb3
    positions[5] = deg2rad({-20.0, -19.3, 95.9, -3.0, 23.5, -17.6});    // w4 + 10
    positions[6] = deg2rad({53.9, -6.2, 121.9, 5.3, 34.8, 45.6});       // wb4
    positions[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // zero

    std::vector<Matrix> tasks;
    for (u32 i = 0; i < positions.size() - 1; i++) {
        Matrix tmp(6, 6);

        for (u32 j = 0; j < link6.joints; j++) {
            tmp(j, 0) = positions[i][j];
            tmp(j, 3) = positions[i+1][j];
        }
        tasks.push_back(tmp);
    }

    // --- Define Guess ---
    Guess guess;
    // std::vector<Guess> guess_predefined;

    // auto xlen = bspline.xlen(tasks[0]);
    // guess_predefined = read_array_from_csv(input_file_initial_guess, xlen);
    // Assert(!guess_predefined.empty());


    // -- Define Constraints ---
    Constraints<GenericManipulator> constraints;
    constraints.show_info           = false;

    constraints.position            = true;
    constraints.velocity            = true;
    constraints.acceleration        = true;
    constraints.torque              = true;

    constraints.tcp_speed           = false;
    constraints.self_collisions     = false;
    constraints.external_collisions = false;
    constraints.n_collision_constraints = 100;
    constraints.n_collision_skip = 2;

    // --- Define Objective ---
    Objective<GenericManipulator> objective;
    objective.K_time = 1;

    u32 bench_counter = 0;
    std::vector<real> traj_time(task_size);
    std::vector<real> traj_time_acc(task_size);
    std::vector<real> comp_time(task_size);
    std::vector<real> comp_time_acc(task_size);

    std::vector<BenchCase> cases;
    const u32 nctrl[] = {12, 16, 20};
    const u32 npoints[] = {100, 200, 256};
    const u32 p[] = {3, 5};
    for (const auto& p_val : p) {
        for (const auto& n_ctrl : nctrl) {
            for (const auto& n_points : npoints) {
                cases.push_back({n_ctrl, n_points, p_val});
            }
        }
    }

    for (u32 c = 0; c < cases.size(); c++) {

        // --- Define Bspline ---
        Bspline bspline(cases[c].nctrl, cases[c].npoints, cases[c].p, 6);

        for (u32 t = 0; t < tasks.size(); t++) {
            // get underlying buffer
            std::streambuf* orig_buf = cout.rdbuf();

            // Base Line Optimization
            Optimization<GenericManipulator> opt(link6, tasks[t]);

            opt.set_bspline(bspline);
            opt.set_guess(guess);
            opt.set_constraints(constraints);
            opt.set_objective(objective);

            bench_counter = 0;
            // set null
            cout.rdbuf(NULL);
            BENCHMARK("Optimize") {
                auto results = optimize(&opt);
                traj_time[t] += results.x[results.x.size - 1];
                comp_time[t] += results.compute_time;
                bench_counter++;
            };
            traj_time[t] = traj_time[t] / bench_counter;
            comp_time[t] = comp_time[t] / bench_counter;
            // restore buffer
            cout.rdbuf(orig_buf);

            // Acceleration 1 Optimization
            Optimization<GenericManipulator> opt_acc(link6, tasks[t]);
            opt_acc.set_bspline(bspline);
            opt_acc.set_guess(guess);
            opt_acc.set_constraints(constraints);
            opt_acc.set_objective(objective);

            bench_counter = 0;
            // set null
            cout.rdbuf(NULL);
            BENCHMARK("Optimize Acceleration") {
                auto result_acc = optimize_acc(&opt_acc);
                traj_time_acc[t] += result_acc.x[result_acc.x.size - 1];
                comp_time_acc[t] += result_acc.compute_time;
                bench_counter++;
            };
            traj_time_acc[t] = traj_time_acc[t] / bench_counter;
            comp_time_acc[t] = comp_time_acc[t] / bench_counter;
            // restore buffer
            cout.rdbuf(orig_buf);

        }

        std::cout << std::endl;
        std::cout << "Benchmark Results for n_ctrl = " << bspline.nctrl << ", n_points = " << bspline.points << " & p = " << bspline.p << std::endl;
        for (u32 i = 0; i < traj_time.size(); i++) {
            auto acceleration = (comp_time[i] - comp_time_acc[i]) / comp_time[i] * 100;
            auto quality = 100 - (traj_time[i] - traj_time_acc[i]) / traj_time[i] * 100;

            std::cout << "Traj t:        " << i << std::endl;
            std::cout << "Acceleration:  " << acceleration << "%" << std::endl;
            std::cout << "Quality:       " << quality << "%" << std::endl;
            std::cout << "Traj time:     " << traj_time[i] <<     " Comp time:     " << comp_time[i] << std::endl;
            std::cout << "Traj time acc: " << traj_time_acc[i] << " Comp time acc: " << comp_time_acc[i] << std::endl;
        }

    }


}
