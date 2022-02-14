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

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(
    tests 
    tests.cpp
)
target_link_libraries(
  tests
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(tests)
