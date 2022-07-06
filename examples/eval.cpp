
#include <stdio.h>
#include "blast.hpp"
#include "blast_optional_utilities.hpp"

int main() {
    using namespace blast;

    // setup manipulator
    ManipulatorUR5 manip;
    manip.init_dynamics(2); // 2 kg point mass payload

    // setup trajectory
    const u32 njoints = manip.joints;
    const u32 npts = 256;
    const u32 nctrl = 8*3;
    const u32 p = 5;
    PvaBspline pva(nctrl, npts, p, njoints);

    // setup task
    Matrix task(njoints, 6);
    // init pos
    task(0, 0) = 0  * pi / 180;
    task(1, 0) = 15 * pi / 180;
    task(2, 0) = 15 * pi / 180;
    task(3, 0) = 15 * pi / 180;
    task(4, 0) = 15 * pi / 180;
    task(5, 0) = 0  * pi / 180;
    // final pos
    task(0, 3) = 30 * pi / 180;
    task(1, 3) = 30 * pi / 180;
    task(2, 3) = 30 * pi / 180;
    task(3, 3) = 45 * pi / 180;
    task(4, 3) = 30 * pi / 180;
    task(5, 3) = 15 * pi / 180;

    // setup optimization vector
    Array x(pva.xlen());
    x.back() = 2;

    auto objective = [](Array& x, Matrix& task, Pva& pva, Manipulator& manip) {
        return x.back();
    };
    auto constraints = [](Array& x, Matrix& task, Pva& pva, Manipulator& manip) {
        Matrix efforts(pva.joints, pva.points); // todo: construct every iteration?

        pva.compute_trajectory(x, task);
        manip.dynamics(pva, efforts);

        Array result(pva.points * pva.joints * 2);
        u32 r_index = 0;
        for (u32 point = 0; point < pva.points; point++) {
            for (u32 j = 0; j < pva.joints; j++) {
                result[r_index++] = (manip.tau_min[j] - efforts(j, point)) / manip.tau_max[j];
                result[r_index++] = (efforts(j, point) - manip.tau_max[j]) / manip.tau_max[j];
            }
        }
        return result;
    };

    printf("obj: %f\n", objective(x, task, pva, manip));
    printf("constraints: \n");
    print(constraints(x, task, pva, manip));
    return 0;
}
