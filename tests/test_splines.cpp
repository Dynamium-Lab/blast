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

    BsplineDef def;

    def.nctrl = 21;
    def.npts = 4000;
    def.p = 5;
    def.njoints = 100;

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
    printf("T = %f\n", x[x.size-1]);

    // Compute basis functions
    BsplineBasis basis(def);

    // Compute trajectory
    Pva pva(def);
    pva.bspline(def, x, task, basis);

    real init_max_pos_error = 0;
    real init_max_vel_error = 0;
    real init_max_acc_error = 0;
    real final_max_pos_error = 0;
    real final_max_vel_error = 0;
    real final_max_acc_error = 0;
    real max_vel_error = 0;
    real max_acc_error = 0;
    real dt = x[x.size-1] / (npts-1);
    for (u32 i = 0; i < njoints; i++) {
        // boundary conditions
        init_max_pos_error = std::max(init_max_pos_error, std::abs(pva.pos(i, 0) - task(i, 0)));
        init_max_vel_error = std::max(init_max_vel_error, std::abs(pva.vel(i, 0) - task(i, 1)));
        init_max_acc_error = std::max(init_max_acc_error, std::abs(pva.acc(i, 0) - task(i, 2)));
        final_max_pos_error = std::max(final_max_pos_error, std::abs(pva.pos(i, npts-1) - task(i, 3)));
        final_max_vel_error = std::max(final_max_vel_error, std::abs(pva.vel(i, npts-1) - task(i, 4)));
        final_max_acc_error = std::max(final_max_acc_error, std::abs(pva.acc(i, npts-1) - task(i, 5)));

        // derivatives
        for (u32 j = 1; j < npts-1; j++) {
            real diff_p = (pva.pos(i, j+1) - pva.pos(i, j-1)) / (2*dt);
            max_vel_error = std::max(max_vel_error, std::abs(diff_p - pva.vel(i, j)));
            real diff_v = (pva.vel(i, j+1) - pva.vel(i, j-1)) / (2*dt);
            max_acc_error = std::max(max_acc_error, std::abs(diff_v - pva.acc(i, j)));
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
