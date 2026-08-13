#pragma once


// check the output of a cuda function call and abort if there was an error
#define CUDA_CHECK(expr)                                                                                   \
  {                                                                                                        \
    auto code = (expr);                                                                                    \
    if (code != cudaSuccess) {                                                                             \
      fprintf(stderr, "Cuda error: %s. In file: %s (%d)\n", cudaGetErrorString(code), __FILE__, __LINE__); \
      abort();                                                                                             \
    }                                                                                                      \
  }

#define CUDA_CHECK_KERNEL CUDA_CHECK(cudaPeekAtLastError())
