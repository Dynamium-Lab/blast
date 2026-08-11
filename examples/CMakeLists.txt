# Define a function to automate adding benchmarks
function(automate_add_examples target_name source_file)
    add_executable(${target_name} ${source_file})
    target_compile_features(${target_name} PRIVATE cxx_std_17)
    target_link_libraries(${target_name} PRIVATE blast)
    # PRIVATE: machine-specific tuning for dev targets only, never leaks to
    # downstream consumers of the blast library.
    if (MSVC)
        target_compile_options(${target_name} PRIVATE /arch:AVX2)
    else ()
        target_compile_options(${target_name} PRIVATE -march=native)
    endif ()
endfunction()

set(EXAMPLES
        example_01_forward_kinematics example_01_forward_kinematics.cpp
        example_02_trajectory_optimization example_02_trajectory_optimization.cpp
        example_03_collision_avoidance example_03_collision_avoidance.cpp

        print_bsplines_and_basis_functions extras/print_bsplines_and_basis_functions.cpp

        ICRA2026_benchmark_speedups ICRA2026/ICRA2026_benchmark_speedups.cpp
)

# Loop through and add examples
while (EXAMPLES)
    list(POP_FRONT EXAMPLES target_name)
    list(POP_FRONT EXAMPLES source_file)
    automate_add_examples(${target_name} ${source_file})
endwhile ()
