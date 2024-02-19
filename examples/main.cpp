
#include <stdio.h>
#include "blast.hpp"

int main() {
    using namespace blast;
    Gen3_7DOF manip;

    // init
    const u32 njoints = 7;
    const u32 npts = 101;
    const u32 nctrl = 8;
    const u32 p = 5;
    Bspline bspline(nctrl, npts, p, njoints);

    // random task
    real amp = 2;
    Matrix task(njoints, 6);
    for (u32 i = 0; i < task.rows; i++)
        for (u32 j = 0; j < task.cols; j++)
            task(i, j) = amp * get_random();

    Optimisation opt{&manip, &task, &bspline};

    // random optimization vector
    Array x(bspline.xlen(task));
    for (u32 i = 0; i < x.size; i++)
        x[i] = amp * get_random();
    x[x.size-1] = 5;

    {
        Array cstr(manip.ncon(npts));
        auto f = obj_time(x.size, x.data, nullptr, nullptr);
        cstr_manip(manip.ncon(npts), cstr.data, x.size, x.data, nullptr, &opt);
        // add penalty for values in cstr that are positive
    }

    // Compute trajectory
    bspline.compute_trajectory(x, task);
    auto &traj = bspline.traj;

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
        init_max_pos_error  = std::max(init_max_pos_error,  std::abs((double)traj.pos(i, 0) - (double)task(i, 0)));
        init_max_vel_error  = std::max(init_max_vel_error,  std::abs((double)traj.vel(i, 0) - (double)task(i, 1)));
        init_max_acc_error  = std::max(init_max_acc_error,  std::abs((double)traj.acc(i, 0) - (double)task(i, 2)));
        final_max_pos_error = std::max(final_max_pos_error, std::abs((double)traj.pos(i, npts-1) - (double)task(i, 3)));
        final_max_vel_error = std::max(final_max_vel_error, std::abs((double)traj.vel(i, npts-1) - (double)task(i, 4)));
        final_max_acc_error = std::max(final_max_acc_error, std::abs((double)traj.acc(i, npts-1) - (double)task(i, 5)));

        // derivatives
        for (u32 j = 1; j < npts-1; j++) {
            double diff_p = ((double)traj.pos(i, j+1) - (double)traj.pos(i, j-1)) / (2.0*dt);
            max_vel_error = std::max(max_vel_error, std::abs(diff_p - (double)traj.vel(i, j)));
            double diff_v = ((double)traj.vel(i, j+1) - (double)traj.vel(i, j-1)) / (2.0*dt);
            max_acc_error = std::max(max_acc_error, std::abs(diff_v - (double)traj.acc(i, j)));
        }
    }

    printf("\tMax initial position error: %f\n",     init_max_pos_error);
    printf("\tMax initial velocity error: %f\n",     init_max_vel_error);
    printf("\tMax initial acceleration error: %f\n", init_max_acc_error);
    printf("\tMax final position error: %f\n",       final_max_pos_error);
    printf("\tMax final velocity error: %f\n",       final_max_vel_error);
    printf("\tMax final acceleration error: %f\n",   final_max_acc_error);
    printf("\tMax derived velocity error: %f\n",     max_vel_error);
    printf("\tMax derived acceleration error: %f\n", max_acc_error);

    return 0;
}
