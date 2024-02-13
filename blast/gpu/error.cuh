#pragma once


// check the output of a cuda function call and abort if there was an error
#define cuda_check(expr) \
{ \
    auto code = (expr); \
    if (code != cudaSuccess){ \
        fprintf(stderr, "Cuda error: %s. In file: %s (%d)\n", cudaGetErrorString(code), __FILE__,__LINE__);\
        abort();\
    } \
}

#define cuda_check_kernel cuda_check( cudaPeekAtLastError() )

