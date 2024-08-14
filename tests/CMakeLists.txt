include(FetchContent)
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.10)
FetchContent_MakeAvailable(catch)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/tests/CMake/Modules/)
include(ParseAndAddCatchTests)

set(AVX_FLAGS)
if (MSVC)
    set(AVX_FLAGS "/arch:AVX2")
else ()
    set(AVX_FLAGS " -mavx2")
endif ()
message(avxflags="${AVX_FLAGS}")
