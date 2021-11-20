#include "blast.hpp"
#include <cmath>
#include <random>



static blast::real get_random() {
    static thread_local std::random_device rd;
    static thread_local std::mt19937 e2(rd());
    static thread_local std::uniform_real_distribution<blast::real> dis(-1, 1);
    return dis(e2);
}

int main() {
    using namespace blast;

    const u32 nctrl = 21;
    const u32 npts = 4000;
    const u32 p = 5;
    const u32 njoints = 1000;

    // random task
    real amp = 10;
    Matrix task(njoints, 6);
    for (u32 i = 0; i < task.rows; i++) {
        for (u32 j = 0; j < task.cols; j++) {
            task(i, j) = amp * get_random();
        }
    }

    Matrix basis_p(nctrl, npts);
    Matrix basis_v(nctrl, npts);
    Matrix basis_a(nctrl, npts);
    bspline_basis_functions(nctrl, npts, p, basis_p, basis_v, basis_a);

    Array x(njoints*(nctrl-6) + 1);
    real T = amp * std::abs(get_random());
    x[x.size-1] = T;
    real dt = T / (npts-1);
    Matrix ctrl(nctrl, njoints);
    bspline_control_points(nctrl, njoints, p, x, task, ctrl);

    Matrix pos(njoints, npts);
    Matrix vel(njoints, npts);
    Matrix acc(njoints, npts);
    bspline_pva(nctrl, npts, njoints, T, ctrl, basis_p, basis_v, basis_a, pos, vel, acc);

    real init_max_pos_error = 0;
    real init_max_vel_error = 0;
    real init_max_acc_error = 0;
    real final_max_pos_error = 0;
    real final_max_vel_error = 0;
    real final_max_acc_error = 0;
    real max_vel_error = 0;
    real max_acc_error = 0;
    for (u32 i = 0; i < njoints; i++) {
        // boundary conditions
        init_max_pos_error = std::max(init_max_pos_error, std::abs(pos(i, 0) - task(i, 0)));
        init_max_vel_error = std::max(init_max_vel_error, std::abs(vel(i, 0) - task(i, 1)));
        init_max_acc_error = std::max(init_max_acc_error, std::abs(acc(i, 0) - task(i, 2)));
        final_max_pos_error = std::max(final_max_pos_error, std::abs(pos(i, npts-1) - task(i, 3)));
        final_max_vel_error = std::max(final_max_vel_error, std::abs(vel(i, npts-1) - task(i, 4)));
        final_max_acc_error = std::max(final_max_acc_error, std::abs(acc(i, npts-1) - task(i, 5)));

        // derivatives
        for (u32 j = 1; j < npts-1; j++) {
            real diff_p = (pos(i, j+1) - pos(i, j-1)) / (2*dt);
            max_vel_error = std::max(max_vel_error, std::abs(diff_p - vel(i, j)));
            real diff_v = (vel(i, j+1) - vel(i, j-1)) / (2*dt);
            max_acc_error = std::max(max_acc_error, std::abs(diff_v - acc(i, j)));
        }
    }
    printf("Boundary conditions:\n");
    printf("Initial position max error = %f\n", init_max_pos_error);
    printf("Initial velocity max error = %f\n", init_max_vel_error);
    printf("Initial acceleration max error = %f\n", init_max_acc_error);
    printf("Final position max error = %f\n", final_max_pos_error);
    printf("Final velocity max error = %f\n", final_max_vel_error);
    printf("Final acceleration max error = %f\n", final_max_acc_error);
    printf("\nDerivative errors:\n");
    printf("Max velocity error: %f\n", max_vel_error);
    printf("Max acceleration error: %f\n", max_acc_error);
    return 0;
}
