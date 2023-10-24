#define BLAST_ENABLE_TESTS
#include "blast.hpp"
#include "gpu/gpu_manipulator.hpp"
#include "gpu/gpu_trajectory.hpp"




// Benchmarks
// TEST_CASE("GPU trajectory computation speed", "[Trajectory]") {
//     using namespace blast;
//     const u32 points = 256;
//     const u32 joints = 7;
//     const u32 p = 5;
//     const u32 ncontrol = 24;
//     const u32 ntrajectories = 80;
//     cuMultiBsplines gpu_bsplines(points, joints, p, ncontrol, ntrajectories);
//     Bspline host_bspline(ncontrol, points, p, joints);
//     // random task
//     Matrix task(joints, 6);
//     real amp = 10;
//     for (u32 i = 0; i < task.rows; i++)
//         for (u32 j = 0; j < task.cols; j++)
//             task(i, j) = amp * get_random();
//     // random optimization vector
//     Array x(joints*(ncontrol-6) + 1);
//     for (u32 i = 0; i < x.size; i++)
//         x[i] = amp * get_random();
//     x.back() = abs(x.back());
//     // copy the 'x' Array 'ntrajectories' times for the gpu version
//     const u32 xlen = host_bspline.xlen(task);
//     Array x_multi(ntrajectories * host_bspline.xlen(task));
//     for (int i = 0; i < ntrajectories; i++) {
//         for (int j = 0; j < xlen; j++) {
//             x_multi[i*xlen + j] = x[j];
//         }
//     }
//     {
//         std::string msg = "Computing " + std::to_string(ntrajectories) + " trajectories on the CPU";
//         BENCHMARK(msg.c_str()) {
//             for (int i = 0; i < ntrajectories; i++)
//                 host_bspline.compute_trajectory(x, task);
//             return host_bspline.traj.pos(0, 0);
//         };
//     }
//     {
//         std::string msg = "Computing " + std::to_string(ntrajectories) + " trajectories on the GPU";
//         BENCHMARK(msg.c_str()) {
//             gpu_bsplines.compute_trajectories(x_multi, task);
//             gpu_bsplines.fetch_trajectories();
//             return gpu_bsplines.pos[0];
//         };
//     }
// }

// TEST_CASE("Collision detection speed", "[World]") {
//     using namespace blast;
//     BENCHMARK("Capsule - Sphere (4 objects)") {
//         real dist = 0;
//         for (auto t : test) {
//             dist = distmin(t.caps, t.sph);
//         }
//         return dist;
//     };
//     BENCHMARK("Capsule - Capsule (5 objects)") {
//         real dist = 0;
//         for (auto t : test_caps) {
//             dist += distmin(t.caps1, t.caps2);
//         }
//         return dist;
//     };
//     BENCHMARK("Capsule - Cylinder (5 objects)") {
//         for (auto t : test_cyl) {
//             real dist = distmin(t.caps, t.cyl);
//         }
//     };
//     BENCHMARK("Capsule - OBB without GJK (30 objects)") {
//         real dist = 0;
//         for (auto t : test_obb) {
//             dist = distmin(t.box, t.caps);
//         }
//         return dist;
//     };
//     BENCHMARK("Capsule - OBB with GJK (30 objects)") {
//         gjkresult res;
//         for (auto t : test_obb) {
//             res = GJK_solve_gjk_simple(t.caps, t.box);
//         }
//         return res;
//     };
//     BENCHMARK("Capsule - OBB with boolean GJK (30 objects)") {
//         bool res;
//         for (auto t : test_obb) {
//             res = GJK_bool(t.caps, t.box);
//         }
//         return res;
//     };
// }
