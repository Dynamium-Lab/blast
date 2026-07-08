#pragma once

#include <blast>
#include "blast_manipulator.hpp"
#include "blast_trajectory.hpp"

#if defined(_MSC_VER)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#undef NOMINMAX
#undef WIN32_LEAN_AND_MEAN
#else
#include <time.h>
#endif
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include "json.hpp"

namespace blast {
struct ConstraintSelection;
struct Objective;
struct Guess;
struct Optimization;
struct Result;

inline host_fn int64_t get_tick_us() {
#if defined(_MSC_VER)
  LARGE_INTEGER start, frequency;
  QueryPerformanceFrequency(&frequency);
  QueryPerformanceCounter(&start);
  return (start.QuadPart * 1000000) / frequency.QuadPart;
#else
  struct timespec start;
  clock_gettime(CLOCK_MONOTONIC, &start);
  return (start.tv_sec * 1000000LLU) + (start.tv_nsec / 1000);
#endif
}

inline blast_fn void print(Vec3 v);

inline blast_fn void print(Mat3 m);

inline blast_fn void print(const Array& a);

inline blast_fn void print(const Matrix& m);

inline blast_fn void print(const double* data, unsigned size);

inline blast_fn void print(const float* data, unsigned size);

inline blast_fn void print(const double* data, unsigned rows, unsigned cols);

inline blast_fn void print(const float* data, unsigned rows, unsigned cols);

inline host_fn void print_to_csv(const Matrix& m, const std::string& filename);

inline host_fn void print_to_csv(const Trajectory& traj, const std::string& filename);

inline host_fn void print_to_csv(const std::vector<Trajectory>& traj, const std::string& filename);

inline host_fn blast::Matrix read_csv_matrix(const std::string& filename, const char* csv_sep);

inline host_fn blast::Matrix read_csv_matrix_no_header(const std::string& filename, const char* csv_sep);

inline host_fn Trajectory read_csv_trajectory_no_header(const std::string& filename, const char* csv_sep);

template<typename T, std::size_t N>
host_fn bool is_close(const std::array<T, N>& a1, const std::array<T, N>& a2, real eps = BLAST_EPSILON);

template<typename T>
host_fn bool is_close(const ObjMatrix<T>& a1, const ObjMatrix<T>& a2, real eps = BLAST_EPSILON);

template<typename T>
host_fn bool is_close(const std::vector<T>& a1, const std::vector<T>& a2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(real r1, real r2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const Box& box1, const Box& box2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const DynamicBox& box1, const DynamicBox& box2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const Sphere& sph1, const Sphere& sph2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const DynamicSphere& sph1, const DynamicSphere& sph2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const Capsule& capsule1, const Capsule& capsule2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const DynamicCapsule& capsule1, const DynamicCapsule& capsule2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const World& world1, const World& world2, real eps = BLAST_EPSILON);

inline host_fn bool is_close(const ManipulatorTempData& manip_data1, const ManipulatorTempData& manip_data2, const u32 n_joints, const u32 n_caps, real eps = BLAST_EPSILON);

inline host_fn bool operator==(const CollisionModelCapsule& a, const CollisionModelCapsule& b);
inline host_fn bool operator==(const Manipulator& a, const Manipulator& b);
inline host_fn bool operator==(const Bspline& a, const Bspline& b);
inline host_fn bool operator==(const ConstraintSelection& a, const ConstraintSelection& b);
inline host_fn bool operator==(const Objective& a, const Objective& b);
inline host_fn bool operator==(const Guess& a, const Guess& b);
inline host_fn bool operator==(const Optimization& a, const Optimization& b);
inline host_fn bool operator==(const Result& a, const Result& b);

} // namespace blast
