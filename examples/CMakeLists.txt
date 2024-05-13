add_executable(examples
    blast_nlopt_function_bind.cpp
)

target_include_directories(examples PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../extern/nlopt/include/api/
    ${CMAKE_CURRENT_SOURCE_DIR}/../extern/nlopt/include/util/
    ${CMAKE_CURRENT_SOURCE_DIR}/../extern/nlopt/include/algs/
)

target_compile_features(examples PUBLIC
    cxx_std_17
)

target_link_libraries(examples PRIVATE blast ${CMAKE_CURRENT_SOURCE_DIR}/../extern/nlopt/lib/nlopt.lib)