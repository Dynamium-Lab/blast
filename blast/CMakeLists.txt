set(BUILD_SHARED_LIBS OFF)
add_subdirectory(extern/nlopt)

add_library(blast INTERFACE
    blast)

target_include_directories(blast INTERFACE
    ./
    extern
    extern/vcl
    Manipulator
    Math
    Optimization
    Trajectory
    World
)

target_compile_features(blast INTERFACE cxx_std_17)
target_link_libraries(blast INTERFACE nlopt TracyClient)


if(MSVC)
target_compile_options(blast INTERFACE /arch:AVX2)
else()
target_compile_options(blast INTERFACE -march=native)
endif()
