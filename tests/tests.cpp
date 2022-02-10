#include "blast.hpp"
#include "blast_optional_utilities.hpp"
#include <cmath>
#include <random>

#include <gtest/gtest.h>



TEST(SplineTest, TrajectoryCorrectness) {

    using namespace blast;

    const u32 nctrl = 8*3;
    const u32 npts = 4000;
    const u32 p = 5;
    const u32 njoints = 7;

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
    x[x.size-1] = 3.f;
    printf("T = %f\n", x[x.size-1]);

    // Compute basis functions

    // Compute trajectory
    PvaBspline pva(nctrl, npts, p, njoints);
    pva.compute_trajectory(x, task);

    double init_max_pos_error = 0;
    double init_max_vel_error = 0;
    double init_max_acc_error = 0;
    double final_max_pos_error = 0;
    double final_max_vel_error = 0;
    double final_max_acc_error = 0;
    double max_vel_error = 0;
    double max_acc_error = 0;
    double dt = x[x.size-1] / (double)(npts-1);
    for (u32 i = 0; i < njoints; i++) {
        // boundary conditions
        init_max_pos_error  = std::max(init_max_pos_error,  std::abs((double)pva.pos(i, 0) - (double)task(i, 0)));
        init_max_vel_error  = std::max(init_max_vel_error,  std::abs((double)pva.vel(i, 0) - (double)task(i, 1)));
        init_max_acc_error  = std::max(init_max_acc_error,  std::abs((double)pva.acc(i, 0) - (double)task(i, 2)));
        final_max_pos_error = std::max(final_max_pos_error, std::abs((double)pva.pos(i, npts-1) - (double)task(i, 3)));
        final_max_vel_error = std::max(final_max_vel_error, std::abs((double)pva.vel(i, npts-1) - (double)task(i, 4)));
        final_max_acc_error = std::max(final_max_acc_error, std::abs((double)pva.acc(i, npts-1) - (double)task(i, 5)));

        // derivatives
        for (u32 j = 1; j < npts-1; j++) {
            double diff_p = ((double)pva.pos(i, j+1) - (double)pva.pos(i, j-1)) / (2.0*dt);
            max_vel_error = std::max(max_vel_error, std::abs(diff_p - (double)pva.vel(i, j)));
            double diff_v = ((double)pva.vel(i, j+1) - (double)pva.vel(i, j-1)) / (2.0*dt);
            max_acc_error = std::max(max_acc_error, std::abs(diff_v - (double)pva.acc(i, j)));
        }
    }

    EXPECT_TRUE(init_max_pos_error < 0.2);
    EXPECT_TRUE(init_max_vel_error < 0.2);
    EXPECT_TRUE(init_max_acc_error < 0.2);
    EXPECT_TRUE(final_max_pos_error < 0.2);
    EXPECT_TRUE(final_max_vel_error < 0.2);
    EXPECT_TRUE(final_max_acc_error < 0.2);
    EXPECT_TRUE(max_vel_error < 0.2);
    EXPECT_TRUE(max_acc_error < 0.9);
}
