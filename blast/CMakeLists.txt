# Guard against duplicate registration when Blast is consumed as a submodule
# that is itself a dependency of a larger project that also uses Blast.
if(TARGET blast)
    return()
endif()

include(GNUInstallDirs)

# --------------------------------------------------------------------------
# Vendored dependency: NLopt
# Disable language bindings and tests — we only need the C library.
# EXCLUDE_FROM_ALL keeps nlopt out of the default install; we export it
# explicitly below so downstream find_package(Blast) sees it.
# --------------------------------------------------------------------------
set(NLOPT_PYTHON  OFF CACHE BOOL "NLopt Python bindings"  FORCE)
set(NLOPT_OCTAVE  OFF CACHE BOOL "NLopt Octave bindings"  FORCE)
set(NLOPT_MATLAB  OFF CACHE BOOL "NLopt MATLAB bindings"  FORCE)
set(NLOPT_GUILE   OFF CACHE BOOL "NLopt Guile bindings"   FORCE)
set(NLOPT_SWIG    OFF CACHE BOOL "NLopt SWIG wrappers"    FORCE)
set(NLOPT_TESTS   OFF CACHE BOOL "NLopt unit tests"       FORCE)

if(NOT TARGET nlopt)
    # Pre-initialize variables that nlopt's CMakeLists reads before setting,
    # to silence --warn-uninitialized in the nlopt subdirectory (we cannot
    # modify extern/ since it is a vendored dependency).
    set(NLOPT_DEFINITIONS      "" CACHE STRING "NLopt compile definitions" FORCE)
    set(CPACK_SOURCE_IGNORE_FILES "" CACHE STRING "" FORCE)
    # nlopt line 325 uses Python_VERSION_MAJOR/MINOR unconditionally, even
    # when NLOPT_PYTHON=OFF (find_package(Python) is skipped so they are never
    # set by CMake).  Provide harmless fallbacks so the variable is initialized.
    if(NOT DEFINED Python_VERSION_MAJOR)
        set(Python_VERSION_MAJOR "3")
    endif()
    if(NOT DEFINED Python_VERSION_MINOR)
        set(Python_VERSION_MINOR "0")
    endif()
    
    # ======================================================================
    # FIX: Force NLopt to build as a static library to prevent missing DLLs
    # ======================================================================
    set(NLOPT_CXX ON CACHE BOOL "Build NLopt with C++ compiler" FORCE)
    
    set(OLD_BUILD_SHARED ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libraries" FORCE)

    add_subdirectory(extern/nlopt EXCLUDE_FROM_ALL)
    
    # Restore the original global setting so you don't force the outer 
    # project's other dependencies to be static unexpectedly
    set(BUILD_SHARED_LIBS ${OLD_BUILD_SHARED} CACHE BOOL "Build static libraries" FORCE)

    # Keep the interface include fix from earlier
    target_include_directories(nlopt INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extern/nlopt/src/api>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/extern/nlopt>
    )
endif()

# --------------------------------------------------------------------------
# Blast — header-only interface library
# --------------------------------------------------------------------------
add_library(blast INTERFACE)
add_library(Blast::blast ALIAS blast)

# Build-tree consumers see the source directories directly.
# Install-tree consumers see the installed prefix layout.
#
# Install layout (relative to CMAKE_INSTALL_PREFIX):
#   include/blast/          — blast headers (blast, blast_*.hpp, subdirs)
#   include/blast-vendor/   — bundled vendor headers (VCL, json.hpp)
#   include/blast-vendor/vcl/
target_include_directories(blast
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extern/vcl>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extern>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/blast>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/blast-vendor/vcl>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/blast-vendor>
)

target_compile_features(blast INTERFACE cxx_std_17)

target_link_libraries(blast INTERFACE nlopt)

# Tracy: propagate the correct macro to all consumers of blast.
# TRACY_ENABLE is Tracy's own compile-time switch (checked in blast/blast).
# Do NOT use add_compile_definitions — that would leak to every target.
if(ENABLE_TRACY)
    target_compile_definitions(blast INTERFACE TRACY_ENABLE)
    target_link_libraries(blast INTERFACE TracyClient)
endif()

# -march=native is machine-specific *tuning*. It must NOT be INTERFACE — that
# would force it on every downstream consumer, breaking cross-compilation and
# reproducible builds. The dev targets (examples/, tests/) apply it per-target
# as PRIVATE, choosing /arch:AVX2 for MSVC and -march=native otherwise.
#
# AVX2+FMA, by contrast, is a hard *minimum-ISA requirement*: Blast's headers
# call AVX2/FMA intrinsics directly (e.g. _mm256_fmadd_pd), so consumers cannot
# compile without them. Propagate the feature flags so find_package(Blast)
# consumers just work. Feature flags are additive — they never downgrade a
# consumer's -march=native. MSVC is intentionally excluded: its /arch: is a
# single mutually-exclusive level, so forcing /arch:AVX2 would downgrade an
# /arch:AVX512 consumer (warning D9025); MSVC is covered by the #error guard in
# blast/blast instead.
target_compile_options(blast INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-mavx2;-mfma>
)

# --------------------------------------------------------------------------
# Install rules
# --------------------------------------------------------------------------

# Library targets (nlopt is a compiled dependency; blast is INTERFACE)
install(TARGETS blast nlopt
    EXPORT      BlastTargets
    ARCHIVE     DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY     DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME     DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Main header file (no extension — included as #include <blast>)
install(FILES blast
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/blast
)

# blast_*.hpp + subdirectory headers, preserving directory structure.
# Excludes extern/ (handled separately below) and CMakeLists.txt.
install(
    DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/blast
    FILES_MATCHING PATTERN "*.hpp"
    PATTERN "extern"        EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
)

# VCL vendor headers (needed by consumers — blast/blast includes vectorclass.h)
install(
    DIRECTORY   extern/vcl/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/blast-vendor/vcl
    FILES_MATCHING PATTERN "*.h"
)

# nlopt public headers.
# The exported Blast::nlopt target declares INTERFACE_INCLUDE_DIRECTORIES
# pointing to ${CMAKE_INSTALL_INCLUDEDIR}, so we must install nlopt.h and
# the generated nlopt.hpp there.  (nlopt's own install rules are suppressed
# by EXCLUDE_FROM_ALL, so we handle them explicitly.)
install(FILES extern/nlopt/src/api/nlopt.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/extern/nlopt/nlopt.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# JSON single-header (blast_utilities.hpp includes "json.hpp")
install(FILES extern/json.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/blast-vendor
)