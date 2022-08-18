#include "blast.hpp"
#include "blast_optional_utilities.hpp"

#include <gtest/gtest.h>

TEST(Math, ArrayOperations) {
    using namespace blast;

    Array a(8);
    a = { 1, 1, 1, 1, 1, 1, 1, 1 };
    Array b(8);
    b = { 1, 1, 1, 1, 1, 1, 2, 2 };

    auto r = dot(a, b);
    EXPECT_EQ(r, (real)10.0);
}

TEST(Math, Mat3Operations) {
    using namespace blast;
    Mat3 m1 = {1, 0, 0,     0, 1, 0,    0, 0, 1};
    Mat3 m2 = {0, 1, 0,     1, 0, 0,    0, 0, 1};

    m1*m2;
}

TEST(Math, Mat4Operations) {
    using namespace blast;

    Mat4 m1;
    m1 = {1, 2, 3, 4,  1, 2, 3, 4,    1, 2, 2, 4,    1, 2, 3, 4};
    Mat4 m2;
    m2 = {1, 2, 3, 4,  1, 2, 2, 4,    1, 2, 3, 4,    1, 2, 3, 100};
    for (u32 i = 0; i < 16; i++) {
        m1[i] = 10*get_random();
        m2[i] = 5*get_random();
    }

    Mat4 m3;
    for (u32 i = 0; i < 4; i++) {
        for (u32 j = 0; j < 4; j++) {
            m3(i, j) = m1(i, 0)*m2(0, j) + m1(i, 1)*m2(1, j) + m1(i, 2)*m2(2, j) + m1(i, 3)*m2(3, j);
        }
    }

    Mat4 m4 = m1 * m2;

    m1 *= m2;


    for (u32 i = 0; i < 16; i++) {
        EXPECT_FLOAT_EQ((float)m1[i], (float)m3[i]);
        EXPECT_FLOAT_EQ((float)m1[i], (float)m4[i]);
    }
}

TEST(Math, TwoSegmentDist) {
    auto dist_sqr_test1 = blast::two_segment_distance_sqr({1, 1, 1}, {2, 2, 2}, {1, 0, 0}, {2, 0, 0});
    auto dist_sqr_test2 = blast::two_segment_distance_sqr({1, 1, 1}, {2, 2, 2}, {1, 0, 0}, {2, 3, 2});
    auto dist_sqr_test3 = blast::two_segment_distance_sqr({1.5, 3, 2}, {7, 0, 4}, {8.2, 0, 5}, {2, 2, 0});
    auto dist_sqr_test4 = blast::two_segment_distance_sqr({1, 1, 1}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4});

    EXPECT_FLOAT_EQ((float)dist_sqr_test1, (float)(1.414213562373095*1.414213562373095));
    EXPECT_FLOAT_EQ((float)dist_sqr_test2, (float)(0.408248290463863*0.408248290463863));
    EXPECT_FLOAT_EQ((float)dist_sqr_test3, (float)(0.277660159805987*0.277660159805987));
    EXPECT_FLOAT_EQ((float)dist_sqr_test4, (float)(1.732050807568877*1.732050807568877));
}

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


TEST(BlastManip, SelfCollision) {
    using namespace blast;
    blast::Gen3_7DOF manip;
    blast::Array theta1(7);
    theta1 = {0, 15, 180, 230, 360, 55, 90};
    theta1 = deg2rad(theta1);
    blast::Array theta2(7);
    theta2 = {0, 0, 0, 0, 0, 0, 0};
    theta2 = deg2rad(theta2);
    blast::Array theta3(7);
    theta3 = {0, 16, 180, 221, 358, 284, 88};
    theta3 = deg2rad(theta3);
    blast::Array theta4(7);
    theta4 = {347, 47, 158, 212, 341, 300, 8};
    theta4 = deg2rad(theta4);

    auto dist_sqr_min_1 = manip.collision_dist_sqr(theta1);
    auto dist_sqr_min_2 = manip.collision_dist_sqr(theta2);
    auto dist_sqr_min_3 = manip.collision_dist_sqr(theta3);
    auto dist_sqr_min_4 = manip.collision_dist_sqr(theta4);

    EXPECT_TRUE(dist_sqr_min_1[0] > 0 && dist_sqr_min_1[1] > 0);
    EXPECT_TRUE(dist_sqr_min_2[0] > 0 && dist_sqr_min_2[1] > 0);
    EXPECT_TRUE(dist_sqr_min_3[0] < 0 && dist_sqr_min_3[1] < 0);
    EXPECT_TRUE(dist_sqr_min_4[0] < 0 && dist_sqr_min_4[1] > 0);
}
