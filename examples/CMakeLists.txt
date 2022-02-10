cmake_minimum_required(VERSION 3.5)

project(blast_tests)

set(CMAKE_INCLUDE_CURRENT_DIR   ON)
set(CMAKE_VERBOSE_MAKEFILE      OFF)
set(CMAKE_COLOR_MAKEFILE        ON)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 11)

include_directories(../)
include_directories(../optional)

add_executable(
    main 
    main.cpp test.cpp
)
