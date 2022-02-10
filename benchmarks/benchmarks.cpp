#include "blast.hpp"
#include "blast_optional_utilities.hpp"
#include <random>

#include <benchmark/benchmark.h>

// static void BM_StringCreation(benchmark::State& state) {
//     for (auto _ : state)
//         std::string empty_string;
// }
// BENCHMARK(BM_StringCreation);

// static void BM_StringCopy(benchmark::State& state) {
//     std::string x = "hello";
//     for (auto _ : state)
//         std::string copy(x);
// }
// BENCHMARK(BM_StringCopy);


static blast::real get_random() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 e2(rd());
    static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
    return dis(e2);
}

static void BM_Splines(benchmark::State& state) {
    // init code here

    using namespace blast;

    BsplineDef def;

    def.nctrl = 8*3;
    def.npts = 256;
    def.p = 5;
    def.njoints = 7;

    const auto nctrl = def.nctrl;
    const auto npts = def.npts;
    const auto njoints = def.njoints;
    const auto p = def.p;

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
    Pva pva(def);
    BsplineBasis basis(def);

    for (auto _ : state) {
        auto control = pva.bspline_control(def, x, task);
        pva.bspline(def, x[x.size-1], control, basis);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_Splines)->Unit(benchmark::kMicrosecond);
