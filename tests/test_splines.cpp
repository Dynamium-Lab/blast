#include "blast.hpp"

int main() {
    using namespace blast;

    u32 nctrl = 21;
    u32 npts = 256;
    u32 p = 5;

    Matrix basis_p(nctrl, npts);
    Matrix basis_v(nctrl, npts);
    Matrix basis_a(nctrl, npts);

    bspline_basis_functions(nctrl, npts, p, basis_p, basis_v, basis_a);

    return 0;
}