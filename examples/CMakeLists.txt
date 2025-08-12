# Define a function to automate adding benchamrks
function(automate_add_examples target_name source_file)
    add_executable(${target_name} ${source_file})
    target_compile_features(${target_name} PUBLIC cxx_std_17)
    target_link_libraries(${target_name} PRIVATE blast nlopt)

    # Set runtime output directories for each configuration
    set_target_properties(${target_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_BINARY_DIR}/Debug"
            RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_BINARY_DIR}/Release"
            RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_CURRENT_BINARY_DIR}/RelWithDebInfo"
    )
endfunction()

set(EXAMPLES
        example_manip example_manipulator.cpp
        #        example_demo1 demo1_link6.cpp
        #        example_demo2 demo2.cpp
        #        example_demo2_link6 demo2_link6.cpp
        #        example_demo2_gen3 demo2_gen3.cpp
        #        example_generic_urdf_parser generic_urdf_parser.cpp
        #        example_generic_xacro_processor generic_xacro_processor.cpp
        #        example_generic_gen3_nodep generic_gen3_nodep.cpp
        #        example_demo_ptech demo_ptech.cpp
        #        example_generate_trajectory example_generate_trajectory.cpp

        #        example_acceleration_rrt_connect acceleration/example_rrt_connect.cpp
        #        example_acceleration_gradients acceleration/gradients.cpp
        #        eval_acceleration_gradients acceleration/eval_gradients.cpp
        #        profile_acceleration_gradients acceleration/profile_gradients.cpp
)

# Loop through and add examples
while (EXAMPLES)
    list(POP_FRONT EXAMPLES target_name)
    list(POP_FRONT EXAMPLES source_file)
    automate_add_examples(${target_name} ${source_file})
endwhile ()