// GPU manipulator constraint correctness test
// Compiles with nvcc and links against Catch2 + blast.
//
// The actual TEST_CASE is defined inside manipulator.cuh behind
// #ifdef BLAST_ENABLE_TESTS — we just need to define that macro,
// provide the Catch2 main, and include the header.

#define BLAST_ENABLE_TESTS

// Catch2 v2 single-translation-unit main
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

// Pull in the GPU manipulator implementation and its embedded test.
#include "blast/gpu/manipulator.cuh"
