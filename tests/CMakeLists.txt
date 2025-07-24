include(FetchContent)
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.10)
FetchContent_MakeAvailable(catch)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/tests/CMake/Modules/)
include(ParseAndAddCatchTests)

# Define a function to automate adding benchamrks
function(automate_add_benchmarks target_name source_file)
  add_executable(${target_name} ${source_file})
  target_compile_features(${target_name} PUBLIC cxx_std_17)
  target_link_libraries(${target_name} PRIVATE Catch2::Catch2 blast)

  # Set runtime output directories for each configuration
  set_target_properties(${target_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG       "${CMAKE_CURRENT_BINARY_DIR}/Debug"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE     "${CMAKE_CURRENT_BINARY_DIR}/Release"
    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_CURRENT_BINARY_DIR}/RelWithDebInfo"
  )
endfunction()

# Define a function to automate adding tests
function(automate_add_tests target_name source_file)
  add_executable(${target_name} ${source_file})
  target_compile_features(${target_name} PUBLIC cxx_std_17)
  target_link_libraries(${target_name} PRIVATE Catch2::Catch2 blast)
  add_test(NAME ${target_name}_test COMMAND ${target_name})

  # Set runtime output directories for each configuration
  set_target_properties(${target_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG       "${CMAKE_CURRENT_BINARY_DIR}/Debug"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE     "${CMAKE_CURRENT_BINARY_DIR}/Release"
    RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_CURRENT_BINARY_DIR}/RelWithDebInfo"
  )
endfunction()

# note: Add new benchmark files as follows :
# NAME_OF_BENCHMARK BENCHMARK_TARGET_LOCATION
# List of benchmarks with their source files
set(BENCHMARKS
  # Acceleration
  bench_gradients_acceleration acceleration/bench_gradients_acceleration.cpp
  # bench_generic benchmarks/bench_generic.cpp
  # bench_test_manip benchmarks/bench_test_manip.cpp
  # bench_constraints benchmarks/bench_constraints.cpp
  # bench_functions benchmarks/bench_functions.cpp
)

# note: Add new tests files as follows :
# NAME_OF_TEST TEST_TARGET_LOCATION
# List of tests with their source files
set(TESTS
  # Acceleration
  # test_gradients_acceleration acceleration/test_gradients_acceleration.cpp
  # test_constraints_acceleration acceleration/test_constraints_acceleration.cpp

  # Container
  # test_container_objmatrix container/test_objmatrix.cpp
  # test_container_objdblarray container/test_objdblarray.cpp
  # Manipulator
  # test_manipulator_functions manipulators/test_functions.cpp
  # test_generic_gen3 manipulators/test_generic_gen3.cpp
  # test_generic_link6 manipulators/test_generic_link6.cpp
  # test_generic_crx25 manipulators/test_generic_crx25.cpp
  # test_generic_ur5e manipulators/test_generic_ur5e.cpp
  # Optimization
  # test_optimization optimization/test_optimization.cpp
  # test_optimization_guess optimization/test_initial_guess.cpp
  # test_optimization_objective optimization/test_objective.cpp
  # test_optimization_constraints optimization/test_constraints.cpp
  # test_optimization_tasks optimization/test_task.cpp
  # test_optimization_hierarchical optimization/test_hierarchical_optimization.cpp
  # Math
  # test_math_matrix math/test_matrix.cpp
  # World
  # test_primitives world/test_primitives.cpp
  # Trajectory
  # test_bspline trajectory/test_bsplines.cpp
)

# Loop through and add benchmarks
while(BENCHMARKS)
  list(POP_FRONT BENCHMARKS target_name)
  list(POP_FRONT BENCHMARKS source_file)
  automate_add_benchmarks(${target_name} ${source_file})
endwhile()

# Loop through and add tests
while(TESTS)
  list(POP_FRONT TESTS target_name)
  list(POP_FRONT TESTS source_file)
  automate_add_tests(${target_name} ${source_file})
endwhile()