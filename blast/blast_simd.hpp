#pragma once

#include <xmmintrin.h>
#include <immintrin.h>



// Horizontal sum of 4 lanes of the vector
inline float simd_hadd( __m128 r4 ) {
    // Add 4 values into 2
    const __m128 r2 = _mm_add_ps( r4, _mm_movehl_ps( r4, r4 ) );
    // Add 2 lower values into the final result
    const __m128 r1 = _mm_add_ss( r2, _mm_movehdup_ps( r2 ) );
    // Return the lowest lane of the result vector.
    // The intrinsic below compiles into noop, modern compilers return floats in the lowest lane of xmm0 register.
    return _mm_cvtss_f32( r1 );
}

// Horizontal sum of 8 lanes of the vector
inline float simd_hadd( __m256 r8 ) {
    const __m128 low = _mm256_castps256_ps128( r8 );
    const __m128 high = _mm256_extractf128_ps( r8, 1 );
    return simd_hadd( _mm_add_ps( low, high ) );
}

inline double simd_hadd(__m256d v) {
#if defined(__MINGW64__) || defined(__MINGW64__) || defined(__MINGW32__)
    return v[0] + v[1] + v[2] + v[3];
#else
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
    vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128

    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));  // reduce to scalar
#endif
}