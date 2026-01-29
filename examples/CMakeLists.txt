# Define a function to automate adding benchmarks
function(automate_add_examples target_name source_file)
  add_executable(${target_name} ${source_file})
  target_compile_features(${target_name} PUBLIC cxx_std_17)
  target_link_libraries(${target_name} PRIVATE blast nlopt)
endfunction()

set(EXAMPLES
    example_manipulator example_manipulator.cpp
    eval_manipulator ICRA2026/eval_manipulator.cpp
    get_bspline_and_basis_functions get_bspline_and_basis_functions.cpp
    get_start_and_end_pos get_start_and_end_pos.cpp

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
