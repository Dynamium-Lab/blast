add_subdirectory(extern/nlopt)

add_library(blast
    manipulator/Gen3.cpp
    manipulator/Link6.cpp
    manipulator/R2.cpp
    optimization/optimization.cpp
    trajectory/bspline.cpp
    trajectory/polynomial_5order.cpp
    world/world.cpp
    world/collision_pso.cpp
    world/collision_gwo.cpp
    world/collision_ga.cpp
)

target_include_directories(blast PUBLIC
    ../include
    ./
    Manipulator
    Math
    Optimization
    Trajectory
    World
)


target_compile_features(blast PUBLIC cxx_std_17)
target_link_libraries(blast PUBLIC nlopt)
target_include_directories(blast PUBLIC extern/taskflow extern/vcl)

if(MSVC)
target_compile_options(blast PUBLIC /arch:AVX2)
else()
target_compile_options(blast PUBLIC -march=native)
endif()
