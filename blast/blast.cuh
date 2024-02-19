#pragma once

#ifdef BLAST_ENABLE_TESTS
#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch/catch.hpp"
#endif

#ifndef BLAST_USE_DOUBLES
#define BLAST_USE_DOUBLES 1
#endif

#define blast_fn __host__ __device__ inline
#define host_fn __host__ inline
#define dev_fn __device__ inline

#include "utilities/utilities.hpp"

#include "blast_error.hpp"
#include "blast_math.hpp"
#include "blast_trajectory.hpp"
#include "blast_optimization.hpp"
#include "blast_world.hpp"

#include "gpu/error.cuh"
#include "gpu/manipulator.hpp"
#include "gpu/trajectory.hpp"
