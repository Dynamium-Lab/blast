include(FetchContent)
FetchContent_Declare(
    catch
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v2.13.10
)
FetchContent_MakeAvailable(catch)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tests/CMake/Modules/")
include(ParseAndAddCatchTests)

# Usage: add_blast_target(target_name source_file [TEST | BENCH])
function(add_blast_target target_name source_file type)
  add_executable(${target_name} ${source_file})
  target_compile_features(${target_name} PRIVATE cxx_std_17)

  # PRIVATE: machine-specific tuning for dev targets only, never leaks to
  # downstream consumers of the blast library.
  if(MSVC)
    target_compile_options(${target_name} PRIVATE /arch:AVX2)
  else()
    target_compile_options(${target_name} PRIVATE -march=native)
  endif()

  # tests/ is the include root for test_helper/, etc.
  target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  # Base dependencies
  target_link_libraries(${target_name} PRIVATE Catch2::Catch2 blast)

  # Conditional Tracy logic
  if (ENABLE_TRACY)
    target_link_libraries(${target_name} PRIVATE TracyClient)
    target_compile_definitions(${target_name} PRIVATE TRACY_ENABLE TRACY_NO_EXIT)
  endif ()

  # If marked as a TEST, register it with CTest
  if ("${type}" STREQUAL "TEST")
    add_test(NAME ${target_name}_test COMMAND ${target_name})
  endif ()
endfunction()

# note: Add new benchmark files as follows :
# NAME_OF_BENCHMARK BENCHMARK_TARGET_LOCATION
# List of benchmarks with their source files
set(BENCHMARKS
    # Acceleration
    #  bench_gradients_acceleration acceleration/bench_gradients_acceleration.cpp
    # bench_generic benchmarks/bench_generic.cpp
    # bench_test_manip benchmarks/bench_test_manip.cpp
    #    bench_manipulator benchmarks/bench_manipulator.cpp
    #    bench_accelerations acceleration/bench_accelerations.cpp
    #    bench_inplace acceleration/bench_inplace.cpp
    #    bench_paper_gradients benchmarks/bench_paper_gradients.cpp
    #    bench_constraints benchmarks/bench_constraints.cpp
    #    bench_functions benchmarks/bench_functions.cpp

    # bench_optimization_methods benchmarks/bench_optimization_methods.cpp
)

# note: Add new tests files as follows :
# NAME_OF_TEST TEST_TARGET_LOCATION
# List of tests with their source files
set(TESTS
    # Acceleration
    # test_gradients_acceleration acceleration/test_gradients_acceleration.cpp
    # test_constraints_acceleration acceleration/test_constraints_acceleration.cpp
    #    test_accelerations acceleration/test_accelerations.cpp
    #    test_inplace acceleration/test_inplace.cpp

    # Container
    test_objmatrix container/test_objmatrix.cpp

    # Manipulator
    test_functions manipulators/test_functions.cpp # note: empty test
    test_generic_gen3 manipulators/test_generic_gen3.cpp
    test_generic_ur5e manipulators/test_generic_ur5e.cpp

    # Optimization
    # test_optimization optimization/test_optimization.cpp
    # test_optimization_guess optimization/test_initial_guess.cpp
    # test_optimization_objective optimization/test_objective.cpp
    # test_optimization_constraints optimization/test_constraints.cpp
    # test_optimization_tasks optimization/test_task.cpp

    # Math
    # test_math_matrix math/test_matrix.cpp
    test_misc math/test_misc.cpp

    test_math_vec3 math/test_vec3.cpp
    # World
    test_primitives world/test_primitives.cpp

    # Trajectory
    test_bspline trajectory/test_bsplines.cpp
)

# Loop through and add benchmarks
while (BENCHMARKS)
  list(POP_FRONT BENCHMARKS target_name)
  list(POP_FRONT BENCHMARKS source_file)
  add_blast_target(${target_name} ${source_file} "BENCH")
endwhile ()

# Loop through and add tests
while (TESTS)
  list(POP_FRONT TESTS target_name)
  list(POP_FRONT TESTS source_file)
  add_blast_target(${target_name} ${source_file} "TEST")
endwhile ()