cmake_minimum_required(VERSION 3.16)

project(
    Blast
    VERSION 1.0
    DESCRIPTION "Blast: A fast trajectory optimization library for robot manipulators"
    LANGUAGES CXX
)

set(CMAKE_INCLUDE_CURRENT_DIR   ON)

add_subdirectory(blast)


# Only build examples, docs, and tests if main projects, not if included as a dependancy to another project.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # build examples
    # add_subdirectory(examples)


    # add tests
    enable_testing()
    add_subdirectory(tests)

    # # Generate documentation
    # find_package(Doxygen)
    # if(Doxygen_FOUND)
    # add_subdirectory(docs)
    # else()
    # message(STATUS "Doxygen not found, not building docs")
    # endif()
endif()

