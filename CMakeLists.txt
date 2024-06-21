cmake_minimum_required(VERSION 3.10)

project(
    Blast
    VERSION 1.0
    DESCRIPTION "Blast: A fast trajectory optimization library for robot manipulators"
    LANGUAGES CXX
)

# add_compile_options("/MT")
set(CMAKE_INCLUDE_CURRENT_DIR   ON)

add_subdirectory(src)


# Only build examples, docs, and tests if main projects, not if included as a dependancy to another project.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # build examples
    add_subdirectory(examples)

    add_executable(main main.cpp)  
    target_compile_features(main PUBLIC cxx_std_17)
    target_link_libraries(main PRIVATE blast)


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

# add dependancies
# add_subdirectory(extern/eigen)
# add_subdirectory(extern/nlopt)

# if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR MODERN_CMAKE_BUILD_TESTING)
#    AND BUILD_TESTING)
#   add_subdirectory(test)
# endif()
