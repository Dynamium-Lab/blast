#include <blast>
#include "test_helper/test_helper.hpp"
#include "test_helper/test_functions.hpp"

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch2/catch.hpp"

using namespace blast;

// does not pass  todo: (andre) works for me
TEST_CASE("Basis function test", "[bspline]") {
    int n_ctrl = 12;
    int p = 5;
    int n_joints = 3;
    int n_pts = 100;

    Bspline spline_test_base(n_ctrl, n_pts, p, n_joints);
    spline_test_base.compute_basis();

    real total;
    for (int i = 0; i < n_pts; i++) {
        total = 0;
        for (int j = 0; j < n_ctrl; j++) {
            total += spline_test_base.basis_p(j, i);
        }
        CHECK(is_close(total, 1.0));
    }

    SECTION("CHANGING BSPLINE ORDER p") {
        p = 3;
        Bspline spline_test_p_less(n_ctrl, n_pts, p, n_joints);
        spline_test_p_less.compute_basis();
        
        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_p_less.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
    
        p = 7;
        Bspline spline_test_p_more(n_ctrl, n_pts, p, n_joints);
        spline_test_p_more.compute_basis();
    
        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_p_more.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
    }

    SECTION("CHANGING NUMBER OF CONTROL POINTS n_ctrl") {
        p = 5;
        n_ctrl = 10;
        Bspline spline_test_ctrl_less(n_ctrl, n_pts, p, n_joints);
        spline_test_ctrl_less.compute_basis();

        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_ctrl_less.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
        
        n_ctrl = 16;
        Bspline spline_test_ctrl_more(n_ctrl, n_pts, p, n_joints);
        spline_test_ctrl_more.compute_basis();
        
        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_ctrl_more.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
    }

    SECTION("CHANGING NUMBER OF POINTS n_pts") {
        n_ctrl = 12;
        n_pts = 50;
        Bspline spline_test_pts_less(n_ctrl, n_pts, p, n_joints);
        spline_test_pts_less.compute_basis();

        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_pts_less.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
        
        n_pts = 250;
        Bspline spline_test_pts_more(n_ctrl, n_pts, p, n_joints);
        spline_test_pts_more.compute_basis();
        
        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_pts_more.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
    }

    SECTION("CHANGING NUMBER OF JOINTS n_joints") {
        n_pts = 100;
        n_joints = 1;
        Bspline spline_test_joints_less(n_ctrl, n_pts, p, n_joints);
        spline_test_joints_less.compute_basis();
        
        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_joints_less.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
        
        n_joints = 6;
        Bspline spline_test_joints_more(n_ctrl, n_pts, p, n_joints);
        spline_test_joints_more.compute_basis();
        
        for (int i = 0; i < n_pts; i++) {
            total = 0;
            for (int j = 0; j < n_ctrl; j++) {
                total += spline_test_joints_more.basis_p(j, i);
            }
            CHECK(is_close(total, 1.0));
        }
    }

    SECTION("compute_basis_derivative() function test") {
        n_pts = 1000;
        n_joints = 1;
        Bspline spline_test(n_ctrl, n_pts, p, n_joints);

        spline_test.compute_basis();
        spline_test.compute_basis_derivative(2);

        CHECK(is_close(spline_test.basis_p, spline_test.basis[0]));
        CHECK(is_close(spline_test.basis_v, spline_test.basis[1]));
        CHECK(is_close(spline_test.basis_a, spline_test.basis[2]));
    }
}

// todo: Very large error, confirm this test is ok
TEST_CASE("Bspline velocity and acceleration test", "[bspline]") {
    int n_ctrl = 12;
    int p = 5;
    int n_joints = 1;
    int n_pts = 5;
    // int n_pts = 2001;

    Bspline spline_test(n_ctrl, n_pts, p, n_joints);

    Matrix task(n_joints, 6);
    fill_random(task, PI);

    Array x_test(spline_test.x_len(task));
    fill_random(x_test, 2*PI);
    x_test.back() = 2.0; // time must be positive

    spline_test.compute_trajectory(x_test, task);
    Matrix vel_approx(n_joints, n_pts);
    Matrix acc_approx(n_joints, n_pts);
    real dt = x_test.back() / (n_pts-1); // dt = 0.001 with n_pts = 2001 and T = 2.0
    for (int i = 1; i < n_pts; i++) {
        for (int j = 0; j < n_joints; j++) {
            vel_approx(j, i) = (spline_test.traj.pos(j, i) - spline_test.traj.pos(j, i-1)) / dt;
            acc_approx(j, i) = (spline_test.traj.vel(j, i) - spline_test.traj.vel(j, i-1)) / dt;
        }
    }
    print(vel_approx);
    print(spline_test.traj.vel);
    CHECK(is_close(vel_approx, spline_test.traj.vel, 1e-1));
    CHECK(is_close(acc_approx, spline_test.traj.acc, 1e-1));
}