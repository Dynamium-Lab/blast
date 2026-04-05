cmake_minimum_required(VERSION 3.16)

project(
    Blast
    VERSION 1.0
    DESCRIPTION "Blast: A fast trajectory optimization library for robot manipulators"
    LANGUAGES CXX
)
include(FetchContent)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

option(ENABLE_TRACY "Enable Tracy Profiler support" OFF)

if(ENABLE_TRACY)
  message(STATUS "Tracy Profiler is enabled - fetching dependency...")

  include(FetchContent)
  FetchContent_Declare(
          tracy
          GIT_REPOSITORY https://github.com/wolfpld/tracy.git
          GIT_TAG v0.12.2
          GIT_SHALLOW TRUE
          GIT_PROGRESS TRUE
  )
  FetchContent_MakeAvailable(tracy)

  # Enable LTO/IPO as requested
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

  # Optional: Add a preprocessor macro to your project so your code
  # knows whether to include Tracy headers or not.
  add_compile_definitions(BLAST_TRACY_ENABLE)
  message(STATUS "Building with tracy enabled - make sure the server is launched when executing.")
endif()


add_subdirectory(blast)

# Only build examples, docs, and tests if main projects, not if included as a dependency to another project.
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # build examples
  add_subdirectory(examples)

  # add tests
  enable_testing()
  add_subdirectory(tests)
endif ()

