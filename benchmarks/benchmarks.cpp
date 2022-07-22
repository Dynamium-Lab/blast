#include "blast.hpp"
#include "blast_optional_utilities.hpp"

#include <benchmark/benchmark.h>

static void BM_Splines(benchmark::State& state) {
    // init code here

    using namespace blast;

    const auto npts = 256;
    const auto njoints = 7;
    const auto nctrl = 8*3;
    const auto p = 5;

    // random task
    real amp = 10;
    Matrix task(njoints, 6);
    for (u32 i = 0; i < task.rows; i++) {
        for (u32 j = 0; j < task.cols; j++) {
            task(i, j) = amp * get_random();
        }
    }

    // random optimization vector
    Array x(njoints*(nctrl-6) + 1);
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x[x.size-1] = std::abs(x[x.size-1]);
    // printf("T = %f\n", x[x.size-1]);

    // Compute basis functions
    PvaBspline pva(nctrl, npts, p, njoints);

    for (auto _ : state) {
        pva.compute_trajectory(x, task);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_Splines)->Unit(benchmark::kMicrosecond);

static void BM_Dynamics_Gen3_7dof(benchmark::State& state) {
    using namespace blast;

    Gen3_7DOF manip;

    const auto npts = 256;
    const auto njoints = manip.joints;
    const auto nctrl = 8*3;
    const auto p = 5;
    PvaBspline pva(nctrl, npts, p, njoints);

    // random task
    real amp = 10;
    Matrix task(njoints, 6);
    for (u32 i = 0; i < task.rows; i++) {
        for (u32 j = 0; j < task.cols; j++) {
            task(i, j) = amp * get_random();
        }
    }

    // random optimization vector
    Array x(pva.xlen());
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x[x.size-1] = std::abs(x[x.size-1]);

    pva.compute_trajectory(x, task);
    Matrix efforts(njoints, npts);

    for(auto _ : state) {
        manip.dynamics(pva, efforts);
    }
}
BENCHMARK(BM_Dynamics_Gen3_7dof)->Unit(benchmark::kMicrosecond);
