
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
    auto c = task.col(0);
    Assert(c.is_alias);

    c = {1, 2, 3, 4, 5, 6, 7};

    print(task);
    // compute control points
    pva.compute_control_and_send(x, task);

    // launch kernel
    test_kernal<<< pva.points, 1 >>>(pva);

    // collect results
    pva.fetch_pva();
    print(pva.host_pva->pos);
    pva.clear();
}
