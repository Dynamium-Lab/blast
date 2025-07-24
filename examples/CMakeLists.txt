add_executable(examples
    main.cpp
)

target_compile_features(examples PUBLIC
    cxx_std_17
)

target_link_libraries(examples PRIVATE blast)