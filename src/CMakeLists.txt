add_library(blast
    
    Manipulator/Gen3.cpp
    Manipulator/Link6.cpp
    Manipulator/R2.cpp

    Math/math.cpp
    Math/Vec3.cpp
    Math/Mat3.cpp
    Math/Array.cpp
    Math/Matrix.cpp

    Optimization/Optimization.cpp

    Trajectory/Bspline.cpp
    Trajectory/polynomial_5order.cpp

    World/world.cpp

    World/collision_pso.cpp
    World/collision_gwo.cpp
    World/collision_ga.cpp

    blast_utilities.cpp
)

target_include_directories(blast PUBLIC
    ../include
    ./
    Manipulator
    Math
    Optimization
    Trajectory
    World
    ../extern/mipp
)

target_compile_features(blast PUBLIC
    cxx_std_17
)


set(AVX_FLAGS)
if (MSVC)
    set(AVX_FLAGS "/arch:AVX2")
else ()
    set(AVX_FLAGS " -mavx2")
endif ()
message(avxflags="${AVX_FLAGS}")
target_compile_options(blast PUBLIC ${AVX_FLAGS})