#include "blast.h"
#include "blast_error.h"

namespace blast {

double obj_time(unsigned n, const double* x, double* grad, void*) {
    using namespace blast;
    double result = x[n-1];
    if (grad) {
        for(u32 i = 0; i < n-1; i++)
            grad[i] = 0;
        grad[n-1] = 1;
    }
    return result;
}

Array guess_random(Bspline& bspline, Matrix& task) {
    Array x(bspline.xlen(task));
    fill_random(x, 1);
    x.back() = std::abs(x.back()) * 5 + 0.1;
    return x;
}
}