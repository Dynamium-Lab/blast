


#define BLAST_USE_DOUBLES 0
#include "blast.hpp"
#include "blast_optional_utilities.hpp"

#include <benchmark/benchmark.h>

#ifndef __NVCC__
#error "NVCC is required to run these benchmarks"
#endif

const blast::u32 npoints = 1024*80;
const blast::u32 nblocks = 1024*4;
static_assert(npoints % nblocks == 0);

// static void BM_Mat4(benchmark::State& state) {
//     using namespace blast;
//     const auto npts = npoints;
//     const auto njoints = 7;
//     const auto nctrl = 8*3;
//     const auto p = 5;
//     // random task
//     real amp = 10;
//     Matrix task(njoints, 6);
//     for (u32 i = 0; i < task.rows; i++)
//         for (u32 j = 0; j < task.cols; j++)
//             task(i, j) = amp * get_random();
//     // random optimization vector
//     Array x(njoints*(nctrl-6) + 1);
//     for (u32 i = 0; i < x.size; i++)
//         x[i] = amp * get_random();
//     x[x.size-1] = std::abs(x[x.size-1]);
//     PvaBspline pva(nctrl, npts, p, njoints);
//     pva.compute_trajectory(x, task);
//     Gen3_7DOF manip;
//     for (auto _ : state) {
//         // for (int i = 0; i < npts; i++) {
//         //     manip.forward_kinematics(pva.pos.col(i));
//         // }
//         // note: this is more than 2x faster than calling the above function 'npts' times
//         manip.forward_kinematics(pva.pos);
//         benchmark::ClobberMemory();
//     }
// }
// BENCHMARK(BM_Mat4)->Unit(benchmark::kMicrosecond);

static void BM_CPU_Constraints_PVA(benchmark::State& state) {
    using namespace blast;
    const auto npts = npoints;
    const auto njoints = 7;
    const auto nctrl = 12;
    const auto p = 5;
    // Compute basis functions
    Gen3_7DOF manip;
    PvaBspline pva(nctrl, npts, p, njoints);
    for (auto _ : state) {
        // random task
        real amp = 10;
        Matrix task(njoints, 6);
        for (u32 i = 0; i < task.rows; i++)
            for (u32 j = 0; j < task.cols; j++)
                task(i, j) = amp * get_random();
        // random optimization vector
        Array x(njoints*(nctrl-6) + 1);
        for (u32 i = 0; i < x.size; i++)
            x[i] = amp * get_random();
        x[x.size-1] = std::abs(x[x.size-1]);
        // compute trajectory
        pva.compute_trajectory(x, task);
        manip.validate(pva);
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_CPU_Constraints_PVA)->Unit(benchmark::kMicrosecond);

static void BM_Cuda_Constraints_PVA(benchmark::State& state) {
    using namespace blast;
    const auto npts = npoints;
    const auto njoints = 7;
    const auto nctrl = 12;
    const auto p = 5;

    cuPvaBspline pva;
    cuGen3_7DOF manip;
    pva.init(npts, njoints, p, nctrl);
    manip.init(0, npts);

    cuda_check( cudaDeviceSynchronize() );

    for (auto _ : state) {
        // random task
        const real amp = 10;
        Matrix task(njoints, 6);
        for (u32 i = 0; i < task.rows; i++)
            for (u32 j = 0; j < task.cols; j++)
                task(i, j) = amp * get_random();

        // random optimization vector
        auto x = blast::random_array(pva.host->xlen(), amp);
        x.back() = std::abs(x.back());

        // compute trajectory
        pva.compute_control_and_send(x, task);
        pva_constraints_kernel<<< nblocks, npts/nblocks >>>(pva);
        cuda_check_kernel;
        pva.fetch_pva();
        manip.fetch_constraints(npts);
    }
    cudaDeviceSynchronize();
}
BENCHMARK(BM_Cuda_Constraints_PVA)->Unit(benchmark::kMicrosecond);

static void BM_Cuda_Constraints(benchmark::State& state) {
    using namespace blast;
    const auto npts = npoints;
    const auto njoints = 7;
    const auto nctrl = 12;
    const auto p = 5;

    cuPvaBspline pva;
    cuGen3_7DOF manip;
    pva.init(npts, njoints, p, nctrl);
    manip.init(0, npts);

    cuda_check( cudaDeviceSynchronize() );

    for (auto _ : state) {
        // random task
        real amp = 10;
        Matrix task(njoints, 6);
        for (u32 i = 0; i < task.rows; i++)
            for (u32 j = 0; j < task.cols; j++)
                task(i, j) = amp * get_random();

        // random optimization vector
        auto x = blast::random_array(pva.host->xlen(), amp);
        x.back() = std::abs(x.back());

        // compute trajectory
        pva.compute_control_and_send(x, task);
        constraints_no_pva_kernel<<< nblocks, npts/nblocks >>>(pva);
        cuda_check_kernel;
        manip.fetch_constraints(npts);
    }
    cudaDeviceSynchronize();
}
BENCHMARK(BM_Cuda_Constraints)->Unit(benchmark::kMicrosecond);

// static void BM_Cuda_Splines_NoObjects(benchmark::State& state) {
//     using namespace blast;
//     const auto npts = npoints;
//     const auto njoints = 7;
//     const auto nctrl = 8*3;
//     const auto p = 5;
//     // Compute basis functions
//     cuPvaBspline pva;
//     pva.init(npts, njoints, p, nctrl);
//     cudaDeviceSynchronize();
//     for (auto _ : state) {
//         // random task
//         real amp = 10;
//         Matrix task(njoints, 6);
//         for (u32 i = 0; i < task.rows; i++)
//             for (u32 j = 0; j < task.cols; j++)
//                 task(i, j) = amp * get_random();
//         // random optimization vector
//         Array x(njoints*(nctrl-6) + 1);
//         for (u32 i = 0; i < x.size; i++)
//             x[i] = amp * get_random();
//         x[x.size-1] = std::abs(x[x.size-1]);
//         // compute trajectory
//         pva.compute_control_and_send(x, task);
//         no_object_kernel<<< 1, npts >>>(pva.dt, pva.one_over_T, pva.one_over_T2, pva.joints, pva.ncontrol, pva.device_basis_p, pva.device_basis_v, pva.device_basis_a, pva.device_control, pva.device_pos, pva.device_vel, pva.device_acc, pva.device_t);
//         benchmark::ClobberMemory();
//     }
//     cudaDeviceSynchronize();
// }
// BENCHMARK(BM_Cuda_Splines_NoObjects)->Unit(benchmark::kMicrosecond);

// static void BM_Dynamics_Gen3_7dof(benchmark::State& state) {
//     using namespace blast;
//     Gen3_7DOF manip;
//     const auto npts = npoints;
//     const auto njoints = manip.joints;
//     const auto nctrl = 8*3;
//     const auto p = 5;
//     PvaBspline pva(nctrl, npts, p, njoints);
//     // random task
//     real amp = 10;
//     Matrix task(njoints, 6);
//     for (u32 i = 0; i < task.rows; i++)
//         for (u32 j = 0; j < task.cols; j++)
//             task(i, j) = amp * get_random();
//     // random optimization vector
//     Array x(pva.xlen());
//     for (u32 i = 0; i < x.size; i++)
//         x[i] = amp * get_random();
//     x[x.size-1] = std::abs(x[x.size-1]);
//     pva.compute_trajectory(x, task);
//     Matrix efforts(njoints, npts);
//     for(auto _ : state) {
//         manip.dynamics(pva, efforts);
//     }
// }
// BENCHMARK(BM_Dynamics_Gen3_7dof)->Unit(benchmark::kMicrosecond);
