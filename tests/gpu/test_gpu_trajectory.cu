// GPU trajectory correctness test
// Compiles with nvcc and links against Catch2 + blast.
//
// The actual TEST_CASE is defined inside trajectory.cuh behind
// #ifdef BLAST_ENABLE_TESTS — we just need to define that macro,
// provide the Catch2 main, and include the header.

#define BLAST_ENABLE_TESTS

// Catch2 v2 single-translation-unit main
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

// Pull in the GPU trajectory implementation and its embedded test.
// trajectory.cuh is guarded by #ifdef __NVCC__ so it only compiles
// when nvcc is the active compiler — exactly what we want here.
#include "blast/gpu/trajectory.cuh"
