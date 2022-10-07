
#include <stdio.h>
#include "blast.hpp"
#include "blast_optional_utilities.hpp"

int main() {
    using namespace blast;

    // setup manipulator
    Gen3_7DOF manip;

    // setup trajectory
    const u32 npts = 50;
    const u32 nctrl = 8;
    const u32 p = 5;
    PvaBspline pva(nctrl, npts, p, manip.joints);

    // setup task
    Matrix task(manip.joints, 6);
    // init pos
    task(0, 0) = 0  * pi / 180;
    task(1, 0) = 15 * pi / 180;
    task(2, 0) = 15 * pi / 180;
    task(3, 0) = 15 * pi / 180;
    task(4, 0) = 15 * pi / 180;
    task(5, 0) = 0  * pi / 180;
    task(6, 0) = 0  * pi / 180;
    // final pos
    task(0, 3) = 30 * pi / 180;
    task(1, 3) = 30 * pi / 180;
    task(2, 3) = 30 * pi / 180;
    task(3, 3) = 45 * pi / 180;
    task(4, 3) = 30 * pi / 180;
    task(5, 3) = 15 * pi / 180;
    task(6, 3) = 15 * pi / 180;

    // setup optimization vector
    Array x(pva.xlen());
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x[3] = 4;
    x.back() = 2;

    Optimisation optim;
    optim.manip = &manip;
    optim.pva = &pva;
    optim.task = &task;

    const u32 nconstraints = (5 + 2 + 7*2)*npts;
    Array c(nconstraints);
    Array g(nconstraints*x.size);
    constraints(nconstraints, c.data, x.size, x.data, g.data, &optim);

    return 0;
}
