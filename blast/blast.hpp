#pragma once
#include <stdint.h>
#include <iostream>
#include <math.h>

#ifdef BLAST_ENABLE_TESTS
#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch/catch.hpp"
#endif

#include "Eigen/Dense"

// use doubles by default
#ifndef BLAST_USE_DOUBLES
#define BLAST_USE_DOUBLES 1
#endif

// determine if the compilation enables building for nvidia GPUs
#ifdef __NVCC__
#define blast_fn __host__ __device__ inline
#define host_fn __host__ inline
#else
#define blast_fn inline
#define host_fn inline
#endif


#include "utilities/utilities.hpp"

#include "blast_error.hpp"
#include "blast_math.hpp"
#include "blast_trajectory.hpp"
#include "blast_manipulator.hpp"
#include "blast_optimization.hpp"
#include "manipulators/kinova_gen3.h"
#include "collisions/world.h"
#include "optimization/gwo.hpp"
#include "optimization/pso.hpp"
#include "collisions/Collision_pso.hpp"
#include "collisions/Collision_gwo.hpp"
#include "collisions/Collision_ga.hpp"

