cmake_minimum_required(VERSION 3.5)

project(blast_tests)

set(CMAKE_INCLUDE_CURRENT_DIR   ON)
set(CMAKE_VERBOSE_MAKEFILE      OFF)
set(CMAKE_COLOR_MAKEFILE        ON)

set(CMAKE_CXX_STANDARD 17)

include_directories(../blast/)
include_directories(../blast/optional)

set(NLOPT_DIR "${PROJECT_SOURCE_DIR}/nlopt")

include_directories(${NLOPT_DIR}/include)
include_directories(${NLOPT_DIR}/include/api)
include_directories(${NLOPT_DIR}/include/algs)
link_directories(${NLOPT_DIR}/lib/)

add_executable(
    main 
    main.cpp test.cpp
)
add_executable(
    eval 
    eval.cpp
)
target_link_libraries(eval PRIVATE nlopt.lib)

