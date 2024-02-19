
#include "blast.hpp"
#include <iostream>

using std::cout;
using std::endl;

int main() {
    using namespace blast;

    // init
    const u32 njoints = 6;
    const u32 npts = 256;
    const u32 nctrl = 7;
    const u32 p = 3;
    Bspline bspline(nctrl, npts, p, njoints);

    print_to_csv(bspline.basis_p, "basis_p.csv");
    print_to_csv(bspline.basis_v, "basis_v.csv");
    print_to_csv(bspline.basis_a, "basis_a.csv");
}