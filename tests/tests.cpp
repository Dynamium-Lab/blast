#include "blast.hpp"
#include <cmath>
#include <random>

#include <gtest/gtest.h>


static blast::real get_random() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 e2(rd());
    static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
    return dis(e2);
}

TEST(SimdTests, test) {
    auto a = _mm256_set_pd(1, 2, 3, 4);
    auto b = _mm256_set_pd(4, 3, 2, 1);

    auto r = blast::hadd_pd(_mm256_mul_pd(a, b));

    EXPECT_DOUBLE_EQ(r, 20.0);
}

TEST(SplineTest, TrajectoryCorrectness) {

    using namespace blast;

    BsplineDef def;

    def.nctrl = 8*3;
    def.npts = 4000;
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
    x[x.size-1] = 3.f;
    printf("T = %f\n", x[x.size-1]);

    // Compute basis functions
    BsplineBasis basis(def);

    // Compute trajectory
    Pva pva(def);
    pva.bspline(def, x[x.size-1], pva.bspline_control(def, x, task), basis);

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
