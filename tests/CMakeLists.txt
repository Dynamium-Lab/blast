include(FetchContent)
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.10)
FetchContent_MakeAvailable(catch)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/tests/CMake/Modules/)
include(ParseAndAddCatchTests)

add_executable(profiler test_profiler.cpp)
target_link_libraries(profiler PRIVATE blast)
target_compile_definitions(profiler PRIVATE PROFILER=1)

add_executable(link6_test test_manipulator/test_Link6.cpp)
target_link_libraries(link6_test PRIVATE blast)