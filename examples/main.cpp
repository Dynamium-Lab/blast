
#include <stdio.h>
#include "blast.hpp"
#include "blast_optional_utilities.hpp"

int main() {
    using namespace blast;
    ManipulatorGeneric manip(7);

    // init
    const u32 njoints = 7;
    const u32 npts = 4000;
    const u32 nctrl = 8*3;
    const u32 p = 5;
    PvaBspline pva(nctrl, npts, p, njoints);

    // random task
    real amp = 2;
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
    x[x.size-1] = 0.9f;

    // Compute trajectory
    pva.compute_trajectory(x, task);

    double init_max_pos_error = 0;
    double init_max_vel_error = 0;
    double init_max_acc_error = 0;
    double final_max_pos_error = 0;
    double final_max_vel_error = 0;
    double final_max_acc_error = 0;
    double max_vel_error = 0;
    double max_acc_error = 0;
    double dt = (double)x[x.size-1] / (double)(npts-1);
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

    printf("Normal:\n");
    printf("\tMax initial position error: %f\n", init_max_pos_error);
    printf("\tMax initial velocity error: %f\n", init_max_vel_error);
    printf("\tMax initial acceleration error: %f\n", init_max_acc_error);
    printf("\tMax final position error: %f\n", final_max_pos_error);
    printf("\tMax final velocity error: %f\n", final_max_vel_error);
    printf("\tMax final acceleration error: %f\n", final_max_acc_error);
    printf("\tMax derived velocity error: %f\n", max_vel_error);
    printf("\tMax derived acceleration error: %f\n", max_acc_error);

    return 0;
}
