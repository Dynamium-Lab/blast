#pragma once


#if defined(__CUDA_ARCH__) && (defined(BLAST_DEBUG) || defined(_DEBUG))
#define Assert(expr) \
    if (!(expr)){\
    printf("Assertion failed in function: %s. File: %s(%d).\n", __FUNCTION__, __FILE__, __LINE__); \
    __trap(); \
}
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
