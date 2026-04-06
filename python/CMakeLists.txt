find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
    _blast_core
    STABLE_ABI
    NB_STATIC
    src/module.cpp
    src/bind_math.cpp
    src/bind_world.cpp
    src/bind_manipulator.cpp
    src/bind_task.cpp
    src/bind_trajectory.cpp
    src/bind_optimization.cpp
)

# Blast is header-only with non-inline function definitions in its headers.
# Compiling multiple TUs that each include <blast> produces duplicate symbols.
# Unity build merges all source files into one TU, resolving the ODR issue.
set_target_properties(_blast_core PROPERTIES UNITY_BUILD ON)

target_link_libraries(_blast_core PRIVATE Blast::blast)

# nanobind's NB_STATIC adds -Os (size-optimize) which overrides CMake's -O3.
# Re-assert -O3 here so the heavy optimization kernel runs at full speed.
if(NOT MSVC)
    target_compile_options(_blast_core PRIVATE -O3)
    # Architecture: fall back to minimum AVX2+FMA if the caller hasn't set -march=.
    if(NOT CMAKE_CXX_FLAGS MATCHES "-march=")
        target_compile_options(_blast_core PRIVATE -mavx2 -mfma)
    endif()
else()
    target_compile_options(_blast_core PRIVATE /O2 /arch:AVX2)
endif()

target_compile_definitions(_blast_core PRIVATE BLAST_USE_DOUBLES=1 BLAST_USE_NATIVE_SQP=1)

install(TARGETS _blast_core LIBRARY DESTINATION blast)
