
#include "cuda/blast_cuda.cuh"

#include "blast.hpp"

int main() {
    using namespace blast;
    print_device_properties();

    cuGen3_7DOF manip;
    manip.init(0);
    manip.map_to_gpu(manip_broadcast_arena);
    test_kernal<<< 1, 1 >>>();
}
