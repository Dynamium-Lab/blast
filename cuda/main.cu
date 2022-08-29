
#include "blast.hpp"
#include "optional/blast_optional_utilities.hpp"

int main() {
    using namespace blast;
    print_device_properties();

    const u32 points = 256;
    const u32 joints = 7;
    const u32 p = 5;
    const u32 ncontrol = 12;

    // prep manip
    cuGen3_7DOF manip;
    manip.init(0, points);

    //-- prep trajectory
    cuPvaBspline pva;
    pva.init(points, joints, p, ncontrol);

    // optim vector
    Array x(pva.host->xlen());
    for (int i = 0; i < x.size-1; i++)
        x[i] = get_random();
    x.back() = 5;

    // task
    Matrix task(7, 6);
    task.col(0) = {1, 1, 1, 1, 1, 1, 1};
    task.col(3) = {2, 2, 2, 2, 2, 2, 2};

    // compute control points
    pva.compute_control_and_send(x, task);

    // launch kernel
    pva_constraints_kernel<<< 1, pva.points >>>(pva);
    cuda_check_kernel;

    // collect results
    pva.fetch_pva();
    // printf("GPU computed trajectory:\n");
    // print(pva.host->vel);

    // pva.host->compute_trajectory(x, task);
    // printf("CPU computed trajectory:\n");
    // print(pva.host->vel);
    pva.clear();
}
