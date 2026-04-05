cmake_minimum_required(VERSION 3.21)
project(blast_test_package CXX)

find_package(Blast REQUIRED)

add_executable(test_consumer main.cpp)
target_link_libraries(test_consumer PRIVATE Blast::blast)