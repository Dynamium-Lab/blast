#pragma once


#if defined(__CUDA_ARCH__) && (defined(BLAST_DEBUG) || defined(_DEBUG))
#define Assert(expr) \
    if (!(expr)){\
    printf("Assertion failed in function: %s. File: %s(%d).\n", __FUNCTION__, __FILE__, __LINE__); \
    __trap(); \
}\

#elif (defined(BLAST_DEBUG) || defined(_DEBUG))
#define Assert(expr) \
    if (!(expr)){\
        fprintf(stderr, "Assertion failed in function: %s. File: %s(%d).\n", __FUNCTION__, __FILE__, __LINE__); \
        abort(); \
}
#else
#define Assert(expr) do{} while(0)
#endif

template <typename ToCheck, std::size_t buffer_size, std::size_t real_size = sizeof(ToCheck)>
void assert_buffer_size() {
    static_assert(buffer_size >= real_size, "Buffer is not large enough!");
}

template <typename type1, typename type2, std::size_t size1 = sizeof(type1), std::size_t size2 = sizeof(type2)>
void assert_size() {
    static_assert(size1 == size2, "Types to not have the same size!");
}








//-- functions that are only visible to nvcc GPU compiler -----------------------
#ifdef __NVCC__



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



#endif
