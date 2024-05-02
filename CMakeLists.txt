cmake_minimum_required(VERSION 3.10)

project(
    Blast
    VERSION 1.0
    DESCRIPTION "Blast: A trajectory optimization library for robot manipulators"
    LANGUAGES CXX
)


add_executable(main 
    main.cpp
    src/Math/math.cpp
    src/Math/Vec3.cpp
    src/Math/Mat3.cpp
    src/Math/Array.cpp
    src/Math/Matrix.cpp

    src/World/world.cpp
)

target_include_directories(main PUBLIC
    include
    src
    src/Math
    src/World
    extern/mipp
)

target_compile_features(main PUBLIC
    cxx_std_17
)


# Only build examples, docs, and tests if main projects, not if included as a dependancy to another project.
# if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
#     # build examples
#     add_subdirectory(examples)

#     # Let's ensure -std=c++xx instead of -std=g++xx
#     set(CMAKE_CXX_EXTENSIONS OFF)

#     # Let's nicely support folders in IDEs
#     set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#     # add tests
#     include(CTest)

#     # Generate documentation
#     find_package(Doxygen)
#     if(Doxygen_FOUND)
#     add_subdirectory(docs)
#     else()
#     message(STATUS "Doxygen not found, not building docs")
#     endif()
# endif()

# add dependancies
# add_subdirectory(extern/eigen)
# add_subdirectory(extern/nlopt)

# add source code
# add_subdirectory(src)

# if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING)
#    AND BUILD_TESTING)
#   add_subdirectory(test)
# endif()
