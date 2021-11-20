#include "blast.hpp"

int main() {
    using namespace blast;

    const u32 nctrl = 9;
    const u32 npts = 30;
    const u32 p = 5;
    const u32 njoints = 6;

    Matrix basis_p(nctrl, npts);
    Matrix basis_v(nctrl, npts);
    Matrix basis_a(nctrl, npts);
    bspline_basis_functions(nctrl, npts, p, basis_p, basis_v, basis_a);

    Array x(njoints*(nctrl-6) + 1);
    real T = 2;
    x[x.size-1] = T;
    Matrix task(njoints, 6);
    task(0, 3) = 1;
    task(1, 3) = 2;
    task(2, 3) = 3;
    task(3, 3) = 4;
    task(4, 3) = 5;
    task(5, 3) = 6;
    Matrix ctrl(nctrl, njoints);
    bspline_control_points(nctrl, njoints, p, x, task, ctrl);

    Matrix pos(njoints, npts);
    Matrix vel(njoints, npts);
    Matrix acc(njoints, npts);
    bspline_pva(nctrl, npts, njoints, T, ctrl, basis_p, basis_v, basis_a, pos, vel, acc);

    return 0;
}
