
#include <stdio.h>
#include "blast.hpp"
#include "blast_optional_utilities.hpp"

int main() {
    using namespace blast;

    // init
    const u32 njoints = 7;
    const u32 npts = 256;
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
    printf("T = %f\n", x[x.size-1]);

    // Compute trajectory
    pva.compute_trajectory(x, task);


    return 0;
}
