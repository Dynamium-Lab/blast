cmake_minimum_required(VERSION 3.21)

project(
    Blast
    VERSION 0.1.0
    DESCRIPTION "Blast: A fast trajectory optimization library for robot manipulators"
    LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------
option(ENABLE_TRACY "Enable Tracy Profiler support" OFF)

# Scalar precision of blast::real. ON = double (default), OFF = float.
option(BLAST_USE_DOUBLES "Use double precision for blast::real (OFF = float)" ON)

# Solver backend. OFF = external NLopt (double-only), ON = in-house SQP.
option(BLAST_USE_NATIVE_SQP "Use the in-house SQP solver instead of external NLopt" OFF)

# Developer benchmarks (only built when Blast is the top-level project).
option(BLAST_BUILD_BENCHMARKS "Build benchmarks" OFF)

# A float build cannot use external NLopt: its C API is double-only. Float
# therefore requires the in-house SQP solver. Fail early with a clear message
# rather than a deep template/linker error.
if(NOT BLAST_USE_DOUBLES AND NOT BLAST_USE_NATIVE_SQP)
    message(FATAL_ERROR
        "BLAST_USE_DOUBLES=OFF (float) requires BLAST_USE_NATIVE_SQP=ON: "
        "external NLopt is double-only.")
endif()

# --------------------------------------------------------------------------
# Tracy (optional profiling dependency, fetched only when requested)
# --------------------------------------------------------------------------
if(ENABLE_TRACY)
    message(STATUS "Tracy Profiler is enabled — fetching dependency...")
    include(FetchContent)
    FetchContent_Declare(
        tracy
        GIT_REPOSITORY https://github.com/wolfpld/tracy.git
        GIT_TAG        v0.12.2
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   TRUE
    )
    FetchContent_MakeAvailable(tracy)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    message(STATUS "Building with Tracy enabled — make sure the server is launched when executing.")
endif()

# --------------------------------------------------------------------------
# Core library
# --------------------------------------------------------------------------
add_subdirectory(blast)

# --------------------------------------------------------------------------
# Python bindings (optional)
# --------------------------------------------------------------------------
option(BLAST_BUILD_PYTHON "Build Python bindings via nanobind" OFF)
if(BLAST_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# --------------------------------------------------------------------------
# Developer targets — tests, benchmarks, examples
# Only built when Blast is the top-level project; silently absent when
# consumed via add_subdirectory() or find_package().
# --------------------------------------------------------------------------
if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    add_subdirectory(examples)

    enable_testing()
    add_subdirectory(tests)
endif()

# --------------------------------------------------------------------------
# Install + Export
# --------------------------------------------------------------------------
install(EXPORT BlastTargets
    FILE        BlastTargets.cmake
    NAMESPACE   Blast::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Blast
)

configure_package_config_file(
    cmake/BlastConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/BlastConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Blast
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/BlastConfigVersion.cmake
    VERSION     ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/BlastConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/BlastConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Blast
)