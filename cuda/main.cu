
#include "cuda/blast_cuda.cuh"

#include "blast.hpp"
#include "optional/blast_optional_utilities.hpp"

int main() {
    using namespace blast;
    print_device_properties();

    // prep manip
    cuGen3_7DOF manip;
    manip.init(0);
    manip.map_to_gpu();

    //-- prep trajectory
    cuPvaBspline pva;
    pva.init(11, 7, 5, 8);

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
    test_kernel<<< 1, pva.points >>>(pva);

    // collect results
    pva.fetch_pva();
    printf("GPU computed trajectory:\n");
    print(pva.host->vel);

    pva.host->compute_trajectory(x, task);
    printf("CPU computed trajectory:\n");
    print(pva.host->vel);
    pva.clear();
}
