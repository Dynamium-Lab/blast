cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(blast_benchmarks LANGUAGES CXX CUDA)

set(CMAKE_INCLUDE_CURRENT_DIR   ON)
set(CMAKE_VERBOSE_MAKEFILE      OFF)
set(CMAKE_COLOR_MAKEFILE        ON)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


include_directories(../blast/)
include_directories(../blast/optional)

# todo: generalize this
find_package(benchmark REQUIRED)

add_executable(bench benchmarks.cu)
target_link_libraries(bench benchmark::benchmark benchmark::benchmark_main)
set_target_properties(bench PROPERTIES CUDA_ARCHITECTURES "75;86")
target_compile_options(bench PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xptxas -v>)
# target_compile_options(bench PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -maxrregcount 64>)

# target_link_options(bench "/NODEFAULTLIB:libcmt.lib /NODEFAULTLIB:libcmtd.lib /NODEFAULTLIB:msvcrtd.lib")
# set_target_properties(bench PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

if (MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /Qpar /MT")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:libcmt.lib")
endif (MSVC)
