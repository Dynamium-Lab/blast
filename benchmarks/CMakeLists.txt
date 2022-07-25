cmake_minimum_required(VERSION 3.5)

project(blast_tests)

set(CMAKE_INCLUDE_CURRENT_DIR   ON)
set(CMAKE_VERBOSE_MAKEFILE      OFF)
set(CMAKE_COLOR_MAKEFILE        ON)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


include_directories(../blast/)
include_directories(../blast/optional)

# todo: generalize this
find_package(benchmark REQUIRED)

add_executable(bench benchmarks.cpp)
target_link_libraries(bench benchmark::benchmark benchmark::benchmark_main)

if (MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /Qpar")
endif (MSVC)