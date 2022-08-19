
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
    pva.init(100, 7, 5, 8);

    // optim vector
    Array x(pva.host_pva->xlen());
    for (int i = 0; i < x.size-1; i++)
        x[i] = get_random();
    x.back() = 5;

    // task
    Matrix task(7, 6);
    task.col(0) = {1, 1, 1, 1, 1, 1, 1};
    task.col(3) = {2, 2, 2, 2, 2, 2, 2};
    print(task);

    // compute control points
    pva.compute_control_and_send(x, task);

    printf("host control points:\n");
    print(pva.host_pva->control);
    fflush(stdout);

    // launch kernel
    test_kernal<<< 1, pva.points >>>(pva);

    // collect results
    pva.fetch_pva();
    // print(pva.host_pva->pos);
    // pva.clear();
}
