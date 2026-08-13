include(FetchContent)
FetchContent_Declare(
        catch
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v2.13.10
)
FetchContent_MakeAvailable(catch)

# nanobench powers the (optional) benchmarks; only fetched when building them.
if (BLAST_BUILD_BENCHMARKS)
    FetchContent_Declare(
            nanobench
            GIT_REPOSITORY https://github.com/martinus/nanobench.git
            GIT_TAG v4.3.11
    )
    FetchContent_MakeAvailable(nanobench)
endif ()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tests/CMake/Modules/")
include(ParseAndAddCatchTests)

# Usage: add_blast_target(target_name source_file [TEST | BENCH])
function(add_blast_target target_name source_file type)
    add_executable(${target_name} ${source_file})
    target_compile_features(${target_name} PRIVATE cxx_std_17)

    # PRIVATE: machine-specific tuning for dev targets only, never leaks to
    # downstream consumers of the blast library.
    if (MSVC)
        target_compile_options(${target_name} PRIVATE /arch:AVX2)
    else ()
        target_compile_options(${target_name} PRIVATE -march=native)
    endif ()

    # tests/ is the include root for test_helper/, etc.
    target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

    # Base dependencies. Tests use Catch2; benchmarks use nanobench.
    target_link_libraries(${target_name} PRIVATE blast)
    if ("${type}" STREQUAL "BENCH")
        target_link_libraries(${target_name} PRIVATE nanobench)
    else ()
        target_link_libraries(${target_name} PRIVATE Catch2::Catch2)
    endif ()

    # Conditional Tracy logic
    if (ENABLE_TRACY)
        target_link_libraries(${target_name} PRIVATE TracyClient)
        target_compile_definitions(${target_name} PRIVATE TRACY_ENABLE TRACY_NO_EXIT)
    endif ()

    # The external-NLopt backend builds nlopt as a shared library; copy it beside
    # the executable so it is found at runtime. No-op for the in-house SQP backend
    # (no nlopt target exists there).
    if (NOT BLAST_USE_NATIVE_SQP)
        add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:nlopt> $<TARGET_FILE_DIR:${target_name}>)
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
        bench_core benchmarks/bench_core.cpp
)

# note: Add new tests files as follows :
# NAME_OF_TEST TEST_TARGET_LOCATION
# List of tests with their source files
set(TESTS

        # Container
        test_objmatrix container/test_objmatrix.cpp

        # Manipulator
        test_generic_gen3 manipulators/test_generic_gen3.cpp
        test_generic_ur5e manipulators/test_generic_ur5e.cpp
        test_manipulator_constructor manipulators/test_manipulator_constructor.cpp

        # Optimization
        test_constraints optimization/test_constraints.cpp
        test_task optimization/test_task.cpp
        test_optimization_basic optimization/test_optimization_basic.cpp
        test_optimization_task_violation optimization/test_optimization_task_violation.cpp

        # Math
        test_math_matrix math/test_matrix.cpp
        test_math_mat3 math/test_mat3.cpp
        test_math_array math/test_array.cpp
        test_misc math/test_misc.cpp
        test_math_vec3 math/test_vec3.cpp

        # World
        test_primitives world/test_primitives.cpp

        # Trajectory
        test_bspline trajectory/test_bsplines.cpp
        test_polynomial trajectory/test_polynomial.cpp
)

# Loop through and add benchmarks (opt-in via BLAST_BUILD_BENCHMARKS)
if (BLAST_BUILD_BENCHMARKS)
    while (BENCHMARKS)
        list(POP_FRONT BENCHMARKS target_name)
        list(POP_FRONT BENCHMARKS source_file)
        add_blast_target(${target_name} ${source_file} "BENCH")
    endwhile ()
endif ()

# Loop through and add tests
while (TESTS)
    list(POP_FRONT TESTS target_name)
    list(POP_FRONT TESTS source_file)
    add_blast_target(${target_name} ${source_file} "TEST")
endwhile ()

# --------------------------------------------------------------------------
# GPU tests (only built when BLAST_ENABLE_CUDA=ON)
# --------------------------------------------------------------------------
if(BLAST_ENABLE_CUDA)
    # add_blast_gpu_target(target_name source_file)
    # Compiles a .cu file with nvcc, links blast + Catch2 + CUDA::cudart,
    # and registers it as a CTest test.
    function(add_blast_gpu_target target_name source_file)
        add_executable(${target_name} ${source_file})

        # Require C++17 for both host and device code
        target_compile_features(${target_name} PRIVATE cxx_std_17)
        set_target_properties(${target_name} PROPERTIES
            CUDA_STANDARD 17
            CUDA_STANDARD_REQUIRED ON
            # Separate compilation lets .cu files call device functions
            # defined in other translation units (not needed today, but
            # avoids a painful retrofit later).
            # NOTE: CUDA_SEPARABLE_COMPILATION is intentionally disabled: with
            # GCC 11/12 on Ubuntu 22.04 + CUDA 12.3, enabling it triggers a bug
            # where the cudafe-generated stub references CUDA RNG globals
            # (precalc_xorwow_matrix, mrg32k3aM1, etc.) that GCC fails to find.
            # Since all device code lives in a single .cu file, we don't need it.
            CUDA_SEPARABLE_COMPILATION OFF
        )

        # nvcc flags:
        #   --expt-relaxed-constexpr  — allows constexpr functions in device code
        #                               (blast headers mark helpers constexpr)
        #   -lineinfo                 — embeds source line info for cuda-memcheck /
        #                               Nsight without the full -G overhead
        #   -allow-unsupported-compiler — suppress warnings about GCC 12 with CUDA 12.3
        #   -DBLAST_NO_CURAND         — prevents curand_kernel.h from being included via
        #                               Array.hpp; without this guard CUDA 12.3 generates
        #                               an invalid host stub referencing ::precalc_xorwow_matrix
        #                               etc. that GCC 11/12 cannot resolve (CUDA bug).
        target_compile_options(${target_name} PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:
                --expt-relaxed-constexpr
                -lineinfo
                -allow-unsupported-compiler
                -DBLAST_NO_CURAND
            >
        )

        # tests/ is the include root for test_helper.hpp etc.
        target_include_directories(${target_name} PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            # blast/gpu/*.cuh uses #include <blast> which resolves through
            # the blast INTERFACE include dirs, but nvcc also needs the
            # project root so that #include "blast/gpu/..." paths resolve
            # correctly — the "blast" segment in the path is the source
            # subdirectory, so the root must be one level above it.
            ${CMAKE_SOURCE_DIR}
        )

        target_link_libraries(${target_name} PRIVATE
            blast
            Catch2::Catch2
            CUDA::cudart
        )

        # Tracy (if enabled globally)
        if(ENABLE_TRACY)
            target_link_libraries(${target_name} PRIVATE TracyClient)
            target_compile_definitions(${target_name} PRIVATE TRACY_ENABLE TRACY_NO_EXIT)
        endif()

        # Copy nlopt DLL next to the executable (external-NLopt backend only)
        if(NOT BLAST_USE_NATIVE_SQP)
            add_custom_command(TARGET ${target_name} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:nlopt> $<TARGET_FILE_DIR:${target_name}>)
        endif()

        add_test(NAME ${target_name}_test COMMAND ${target_name})
    endfunction()

    # ----- GPU test targets -----
    # Add new GPU tests here in the same NAME  SOURCE_FILE pattern.
    add_blast_gpu_target(test_gpu_trajectory gpu/test_gpu_trajectory.cu)
    add_blast_gpu_target(test_gpu_manipulator gpu/test_gpu_manipulator.cu)
    add_blast_gpu_target(test_gpu_manipulator_ur5e gpu/test_gpu_manipulator_ur5e.cu)

    # ----- GPU benchmark (not a CTest test, just run directly) -----
    add_executable(bench_gpu_trajectory gpu/bench_gpu_trajectory.cu)
    target_compile_features(bench_gpu_trajectory PRIVATE cxx_std_17)
    set_target_properties(bench_gpu_trajectory PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF
    )
    target_compile_options(bench_gpu_trajectory PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:
            --expt-relaxed-constexpr
            -lineinfo
            -allow-unsupported-compiler
            -DBLAST_NO_CURAND
        >
    )
    target_include_directories(bench_gpu_trajectory PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}
    )
    target_link_libraries(bench_gpu_trajectory PRIVATE
        blast
        CUDA::cudart
    )
endif() # BLAST_ENABLE_CUDA