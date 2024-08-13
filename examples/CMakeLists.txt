add_executable(examples
    blast_nlopt_function_bind.cpp
)

target_compile_features(examples PUBLIC
    cxx_std_17
)

target_link_libraries(examples PRIVATE blast nlopt)