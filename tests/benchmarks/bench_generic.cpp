#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING

#include "catch2/catch.hpp"
#include <blast>
#include "test_helper/test_helper.hpp"

using namespace blast;

/*inline void add_task(Matrix& task, Array pos1, Array pos2, Array vel1 = {}, Array vel2 = {}, Array acc1 = {}, Array acc2 = {}) {
    Assert(pos1.size == pos2.size);
    task.resize(pos1.size, 6);
    for (u32 j = 0; j < task.rows; j++) {
        task(j, 0) = pos1[j];
        task(j, 3) = pos2[j];
        if (vel1.size != 0 && vel2.size != 0) {
            Assert(vel1.size == vel2.size == task.rows);
            task(j, 1) = vel1[j];
            task(j, 4) = vel2[j];
        }
        if (acc1.size != 0 && acc2.size != 0) {
            Assert(acc1.size == acc2.size == task.rows);
            task(j, 2) = acc1[j];
            task(j, 5) = acc2[j];
        }
    }
}

TEST_CASE("Link6 Optimization", "[Generic]") {
    using namespace blast;

    // --- Test characteristics ---
    bool position = true;
    bool velocity = true;
    bool acceleration = true;
    bool torque = true;
    bool tcp_speed = true;
    bool self_collisions = true;
    bool external_collisions = true;


    // manipulator data (todo: could be imported and exported to file - json ?)
    GenericManipulator link6 = get_generic_Link6();

    // Task
    Matrix task(link6.joints, 6);
    // Array pos_zero = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // Array pos_home = {0, -0.349, 1.92, 0, 0.698, 0};

    // for (u32 j = 0; j < link6.joints; j++) {
    //     task(j, 0) = pos_zero[j];
    //     task(j, 3) = pos_home[j];
    // }

    Array pos_wb1 = deg2rad({51.851436614990234,
                             -13.578636169433594,
                             107.87167358398438,
                             3.6194305419921875,
                             33.133209228515625,
                             51.21833801269531});
    Array pos_w1_10cm = deg2rad({-40.445762634277344,
                                 -26.876392364501953,
                                 83.60868835449219,
                                 1.49383544921875,
                                 19.951095581054688,
                                 -42.22943115234375});
    for (int i = 0; i < link6.joints; i++) {
        task(i, 0) = pos_wb1[i];
        task(i, 3) = pos_w1_10cm[i];
    }

    // B-spline
    Bspline bspline(link6.joints);

    // Guess
    Guess guess; // using default guess

    // Constraints
    Constraints<GenericManipulator> constraints;
    constraints.position = position;
    constraints.velocity = velocity;
    constraints.acceleration = acceleration;
    constraints.torque = torque;
    constraints.tcp_speed = tcp_speed;
    constraints.self_collisions = self_collisions;
    constraints.external_collisions = external_collisions;

    // Objective
    Objective<GenericManipulator> objective;
    objective.K_time = 1;

    // World
    World world;
    add_static_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // vertical plate (no coll)
    add_static_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); //table

    // --- Hard-coded Link6 ---

    // B-spline
    Bspline bspline_hc(link6.joints);

    // Constraints
    Link6 hc_link6;
    Constraints<Link6> constraints_hc;
    constraints_hc.position = position;
    constraints_hc.velocity = velocity;
    constraints_hc.acceleration = acceleration;
    constraints_hc.torque = torque;
    constraints_hc.tcp_speed = tcp_speed;
    constraints_hc.self_collisions = self_collisions;
    constraints_hc.external_collisions = external_collisions;

    // Objective
    Objective<Link6> objective_hc;
    objective_hc.K_time = 1;

    real comp_time = 0.0;
    real traj_time = 0.0;
    u32 counter = 0;
    u32 sucess_counter = 0;

    // Generic Manipulator Link6
    Optimization opt_gen(link6, task);

    opt_gen.set_bspline(bspline);
    opt_gen.set_guess(guess);
    opt_gen.set_world(world);
    opt_gen.set_objective(objective);
    opt_gen.set_constraints(constraints);

    std::vector<Matrix> tasks = get_Link6_demo1_tasks();

    BENCHMARK("Generic Link6 Optimization") {
        for (int i = 0; i < tasks.size(); i++) {
            opt_gen.set_task(tasks[i]);
            auto result_gen = blast::optimize(&opt_gen);
            comp_time += result_gen.compute_time;
            traj_time += result_gen.x[result_gen.x.size - 1];
            counter++;
            if(result_gen.success && !result_gen.success_false)
                sucess_counter++;
        }
    };

    real comp_time_hc = 0.0;
    real traj_time_hc = 0.0;
    u32 counter_hc = 0;
    u32 sucess_counter_hc = 0;
    // HC Link6
    Optimization opt_hc(hc_link6, task);
    opt_hc.set_bspline(bspline_hc);
    opt_hc.set_guess(guess);
    opt_hc.set_world(world);
    opt_hc.set_constraints(constraints_hc);
    opt_hc.set_objective(objective_hc);
    BENCHMARK("Hard Coded Link6 Optimization") {
        for (int i = 0; i < tasks.size(); i++) {
            opt_hc.set_task(tasks[i]);
            auto result_hc = blast::optimize(&opt_hc);
            comp_time_hc += result_hc.compute_time;
            traj_time_hc += result_hc.x[result_hc.x.size - 1];
            counter_hc++;
            if(result_hc.success && !result_hc.success_false)
                sucess_counter_hc++;
        }
    };


    std::cout << " " << std::endl;
    std::cout << "Compute time Generic Link6: " << comp_time/counter << std::endl;
    std::cout << "Traj time Generic Link6: " << traj_time/counter << std::endl;
    std::cout << "Success Rate Generic Link6: " << sucess_counter/counter << std::endl;
    std::cout << "Compute time Hard Coded Link6: " << comp_time_hc/counter_hc << std::endl;
    std::cout << "Traj time Hard Coded Link6: " << traj_time_hc/counter_hc << std::endl;
    std::cout << "Success Rate Hard Coded Link6: " << sucess_counter_hc/counter << std::endl;

}
*/

TEST_CASE("Link6 Generic vs Hard-coded", "[]") {
    Link6 manip_hc;
    int n_tests = 1;

    auto opt_gen = get_generic_link6_opt();
    Optimization<Link6> opt_hc = get_hardcoded_link6_opt();

    World world;
    add_static_box({0.67, -0.1475, -0.0562}, {0.35, 0.025, 0.4}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); // vertical plate (no coll)
    add_static_box({0.6415, 0.0237, -0.53815}, {2.0, 2.0, 0.381}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world); //table

    opt_gen.set_world(world);
    opt_hc.set_world(world);

    std::vector<Matrix> tasks = get_Link6_demo1_tasks();

    std::vector<real> constraints_time_gen;
    std::vector<real> constraints_time_hc;
    std::vector<real> objectives_time_gen;
    std::vector<real> objectives_time_hc;

    int max_loop = n_tests > tasks.size() ? tasks.size() : n_tests;
    for (int i = 0; i < tasks.size(); i++) {
        int n = opt_gen.bspline.xlen(tasks[i]);

        opt_gen.set_task(tasks[i]);
        opt_hc.set_task(tasks[i]);

        ncon(&opt_gen);
        ncon(&opt_hc);

        Array x(n);
        fill_random(x, PI);

        Array result_gen(opt_gen.constraints.n_constraints);
        Array result_hc(opt_hc.constraints.n_constraints);
        Array grad_gen(opt_gen.constraints.n_constraints*n);
        Array grad_hc(opt_hc.constraints.n_constraints*n);
        BENCHMARK_ADVANCED("Constraints function generic") (Catch::Benchmark::Chronometer meter) {
            meter.measure([&] {
                auto T1 = get_tick_us();
                nlopt_constraints<GenericManipulator>(opt_gen.constraints.n_constraints, result_gen.data, n, x.data, grad_gen.data, &opt_gen);
                auto T2 = get_tick_us();
                constraints_time_gen.push_back((T2-T1)/1000);
                return result_gen;
            });
        };
        BENCHMARK_ADVANCED("Constraints function hard-coded") (Catch::Benchmark::Chronometer meter) {
            meter.measure([&] {
                auto T1 = get_tick_us();
                nlopt_constraints<Link6>(opt_hc.constraints.n_constraints, result_hc.data, n, x.data, grad_hc.data, &opt_hc);
                auto T2 = get_tick_us();
                constraints_time_hc.push_back((T2-T1)/1000);
                return result_hc;
            });
        };
        CHECK(is_close(result_gen, result_hc));
        CHECK(is_close(grad_gen, grad_hc));

        real obj_result_gen;
        real obj_result_hc;
        Array obj_grad_gen(n);
        Array obj_grad_hc(n);
        BENCHMARK_ADVANCED("Objective function generic") (Catch::Benchmark::Chronometer meter) {
            meter.measure([&] {
                auto T1 = get_tick_us();
                obj_result_gen = objective_function<GenericManipulator>(n, x.data, obj_grad_gen.data, &opt_gen);
                auto T2 = get_tick_us();
                objectives_time_gen.push_back((T2-T1)/1000);
                return obj_result_gen;
            });
        };
        BENCHMARK_ADVANCED("Objective function hard-coded") (Catch::Benchmark::Chronometer meter) {
            meter.measure([&] {
                auto T1 = get_tick_us();
                obj_result_hc = objective_function<Link6>(n, x.data, obj_grad_hc.data, &opt_hc);
                auto T2 = get_tick_us();
                objectives_time_hc.push_back((T2-T1)/1000);
                return obj_result_hc;
            });
        };
        CHECK(obj_result_gen == Approx(obj_result_hc));
        CHECK(is_close(obj_grad_gen, obj_grad_hc));
    }

    std::cout << "Constraints time generic : " << mean(constraints_time_gen) << " ms" << std::endl;
    std::cout << "Constraints time hard-coded : " << mean(constraints_time_hc) << " ms" << std::endl;
    std::cout << "Objectives time generic : " << mean(objectives_time_gen) << " ms" << std::endl;
    std::cout << "Objectives time hard-coded : " << mean(objectives_time_gen) << " ms" << std::endl;
}
