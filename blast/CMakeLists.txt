add_subdirectory(extern/nlopt)


add_library(blast INTERFACE
        world/CoDO.hpp)

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
target_link_libraries(blast INTERFACE nlopt)


if(MSVC)
target_compile_options(blast INTERFACE /arch:AVX2)
else()
target_compile_options(blast INTERFACE -march=native)
endif()
