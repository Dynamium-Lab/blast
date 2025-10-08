#include <stdio.h>
#include <blast>
#include "blast_math.h"
#include "blast_manipulator.h"
#include "collisions/world.h"
#include "nlopt.h"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <sstream>
#include <vector>
#include <string>

using namespace blast;
// This is given to the test function and contains all the informations neeeded to perform the tests.
struct OptimTest {
    std::string csvInput;
    std::string csvOutput;
    std::string csv_sep = ",";
    Array N_individuals;
    Array N_iterations;
    int N_tests = 0;
    Array input_dimensions;
    bool PSO = FALSE;
    bool GWO = FALSE;
    bool GA = FALSE;
    bool caps1obj = FALSE;
    bool capsallobj = FALSE;
    bool robotallobj = FALSE;
    objlist* world;
};

// Reads a csv file from the string filename
Matrix read_csv_matrix (const std::string& filename, Array dimensions) {
    std::ifstream file(filename);
    Assert(file.is_open());
    Matrix pos(dimensions[0], dimensions[1]);
    std::string line;
    std::getline(file, line); // discard header
    for (int i = 0; i < pos.cols; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        for (int k = 0; k < pos.rows; k++) {
            std::string value;
            if (std::getline(iss, value, ',')) {
                pos(k,i)= std::stod(value);
            }
        }
    }
    return pos;
}

// Tests algorithms as specified by Optim_information and writes the answers in a csv file
void test_algorithms(OptimTest& Optim_information) {
    // ------------------------------
    // ---     Initialization     ---
    // ------------------------------
    int n_collision_results = 1;
    Gen3_7DOF manip;
    real r = 0.0;
    std::string csv_sep = Optim_information.csv_sep;
    const int n_tests = Optim_information.N_tests;
    objlist* world = Optim_information.world;

    // Reading trajectory
    Matrix joint_pos;
    joint_pos = read_csv_matrix(Optim_information.csvInput, Optim_information.input_dimensions);
    int njoint = joint_pos.rows;
    // ---  Initialization ends   ---

    // Get cartesian positions
    Matrix cart_pos = manip.robot_capsules(joint_pos, 1);
    int npoint = cart_pos.cols;

    // Create capsules
    capslist capsules;
    capsules.caps.reserve(npoint*(njoint-1));
    for (int i = 0; i < npoint; i++) {
        for (int j = 0; j < cart_pos.rows/3 - 1; j++) {
            capsule temp;
            temp.p1 = { cart_pos(3*j, i),   cart_pos(3*j+1, i), cart_pos(3*j+2, i) };
            temp.p2 = { cart_pos(3*j+3, i), cart_pos(3*j+4, i), cart_pos(3*j+5, i) };
            temp.r = r;
            capsules.caps.push_back(temp);
        }
    }

    // Open output csv file
    std::ofstream csvFile(Optim_information.csvOutput);
    csvFile.clear();
    Assert(csvFile.is_open());

    // Write headers to the CSV file
    csvFile  << "Number of Particles" << csv_sep << "Number of iterations" << csv_sep << "Primitive Distance" << csv_sep << "Primitive test time (ms)" << csv_sep;
    if (Optim_information.PSO == true) {
        if (Optim_information.caps1obj == true) {
            csvFile << "PSO 1 capsule, 1 object test time (ms)"  <<  csv_sep << "PSO 1 capsule, 1 object rel. error" << csv_sep;
        }
        if (Optim_information.capsallobj == true) {
            csvFile << "PSO 1 capsule, all objects test time (ms)"  <<  csv_sep << "PSO 1 capsule, all objects rel. error" << csv_sep;
        }
        if (Optim_information.robotallobj == true) {
            csvFile << "PSO full robot, all objects test time (ms)"  <<  csv_sep << "PSO full robot, all objects rel. error" << csv_sep;
        }
    }
    if (Optim_information.GWO == true) {
        if (Optim_information.caps1obj == true) {
            csvFile << "GWO 1 capsule, 1 object test time (ms)"  <<  csv_sep << "GWO 1 capsule, 1 object rel. error" << csv_sep;
        }
        if (Optim_information.capsallobj == true) {
            csvFile << "GWO 1 capsule, all objects test time (ms)"  <<  csv_sep << "GWO 1 capsule, all objects rel. error" << csv_sep;
        }
        if (Optim_information.robotallobj == true) {
            csvFile << "GWO full robot, all objects test time (ms)"  <<  csv_sep << "GWO full robot, all objects rel. error" << csv_sep;
        }
    }
    if (Optim_information.GA == true) {
        if (Optim_information.caps1obj == true) {
            csvFile << "GA 1 capsule, 1 object test time (ms)"  <<  csv_sep << "GA 1 capsule, 1 object rel. error" << csv_sep;
        }
        if (Optim_information.capsallobj == true) {
            csvFile << "GA 1 capsule, all objects test time (ms)"  <<  csv_sep << "GA 1 capsule, all objects rel. error" << csv_sep;
        }
        if (Optim_information.robotallobj == true) {
            csvFile << "GA full robot, all objects test time (ms)"  <<  csv_sep << "GA full robot, all objects rel. error" << csv_sep;
        }
    }
    csvFile << std::endl;

    // Declaring all possible arrays
    Array times_PSO_caps1obj(n_tests);
    Array result_PSO_caps1obj(n_tests);
    Array error_PSO_caps1obj(n_tests);
    Array times_PSO_capsallobj(n_tests);
    Array result_PSO_capsallobj(n_tests);
    Array error_PSO_capsallobj(n_tests);
    Array times_PSO_robotallobj(n_tests);
    Array result_PSO_robotallobj(n_tests);
    Array error_PSO_robotallobj(n_tests);
    Array times_GWO_caps1obj(n_tests);
    Array result_GWO_caps1obj(n_tests);
    Array error_GWO_caps1obj(n_tests);
    Array times_GWO_capsallobj(n_tests);
    Array result_GWO_capsallobj(n_tests);
    Array error_GWO_capsallobj(n_tests);
    Array times_GWO_robotallobj(n_tests);
    Array result_GWO_robotallobj(n_tests);
    Array error_GWO_robotallobj(n_tests);
    Array times_GA_caps1obj(n_tests);
    Array result_GA_caps1obj(n_tests);
    Array error_GA_caps1obj(n_tests);
    Array times_GA_capsallobj(n_tests);
    Array result_GA_capsallobj(n_tests);
    Array error_GA_capsallobj(n_tests);
    Array times_GA_robotallobj(n_tests);
    Array result_GA_robotallobj(n_tests);
    Array error_GA_robotallobj(n_tests);
    
    // Theoretical minimum
    // For now, just using PSO with a lot of iterations and particles. 
    // todo: Add gradient descent function to find true minimum
    Optim_solution pso_min = test_collision_pso_OBB_data(cart_pos, world, 1000, 1000);
    // real theo_min = grad_desc()

    csvFile << Optim_information.N_tests << csv_sep << "tests" << csv_sep << "Optimal solution : " << csv_sep <<pso_min.best_f << std::endl;

    std::vector<real> result_primitives(n_collision_results);
    for (int i = 0; i < Optim_information.N_individuals.size; i++) {
        for (int j = 0; j < Optim_information.N_iterations.size; j++) {
            // ------------------------------
            // ---  Test with primitives  ---
            // ------------------------------
            auto start_prim = get_tick_us();
            result_primitives = test_collision(&capsules, world, n_collision_results);
            auto stop_prim = get_tick_us();
            real times_prim = (stop_prim - start_prim) / 1000.0;
            for (int test = 0; test < n_tests; test ++) {
                // // note: This will be able to be tested once GJK is up and running
                // // ------------------------------
                // // ---     Test with GJK      ---
                // // ------------------------------
                // auto start_GJK = get_tick_us();
                // // std::vector<real> result_GJK = test_collision_GJK(&capsules, &world, n_collision_results);
                // auto stop_GJK = get_tick_us();
                // times_GJK[test] = (stop_GJK - start_GJK) / 1000.0;

                // ------------------------------
                // ---     Test with PSO      ---
                // ------------------------------
                if (Optim_information.PSO == true) {
                    if (Optim_information.caps1obj == true) {
                        auto start_PSO = get_tick_us();
                        result_PSO_caps1obj[test] = test_collision_pso_OBB(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_PSO = get_tick_us();
                        times_PSO_caps1obj[test] = (stop_PSO - start_PSO) / 1000.0;
                        error_PSO_caps1obj[test] = (result_primitives[0] - result_PSO_caps1obj[test])/result_primitives[0];
                    }
                    if (Optim_information.capsallobj == true) {
                        auto start_PSO = get_tick_us();
                        result_PSO_capsallobj[test] = test_collision_pso_world_1caps(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_PSO = get_tick_us();
                        times_PSO_capsallobj[test] = (stop_PSO - start_PSO) / 1000.0;
                        error_PSO_capsallobj[test] = (result_primitives[0] - result_PSO_capsallobj[test])/result_primitives[0];
                    }
                    if (Optim_information.robotallobj == true) {
                        auto start_PSO = get_tick_us();
                        result_PSO_robotallobj[test] = test_collision_pso_world_full_robot(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_PSO = get_tick_us();
                        times_PSO_robotallobj[test] = (stop_PSO - start_PSO) / 1000.0;
                        error_PSO_robotallobj[test] = (result_primitives[0] - result_PSO_robotallobj[test])/result_primitives[0];
                    }
                }
                // ------------------------------
                // ---     Test with GWO      ---
                // ------------------------------
                if (Optim_information.GWO == true) {
                    if (Optim_information.caps1obj == true) {
                        auto start_GWO = get_tick_us();
                        result_GWO_caps1obj[test] = test_collision_gwo_OBB(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_GWO = get_tick_us();
                        times_GWO_caps1obj[test] = (stop_GWO - start_GWO) / 1000.0;
                        error_GWO_caps1obj[test] = (result_primitives[0] - result_GWO_caps1obj[test])/result_primitives[0];
                    }
                    if (Optim_information.capsallobj == true) {
                        auto start_GWO = get_tick_us();
                        result_GWO_capsallobj[test] = test_collision_gwo_world_1caps(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_GWO = get_tick_us();
                        times_GWO_capsallobj[test] = (stop_GWO - start_GWO) / 1000.0;
                        error_GWO_capsallobj[test] = (result_primitives[0] - result_GWO_capsallobj[test])/result_primitives[0];
                    }
                    if (Optim_information.robotallobj == true) {
                        auto start_GWO = get_tick_us();
                        result_GWO_robotallobj[test] = test_collision_ga_world_full_robot(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_GWO = get_tick_us();
                        times_GWO_robotallobj[test] = (stop_GWO - start_GWO) / 1000.0;
                        error_GWO_robotallobj[test] = (result_primitives[0] - result_GWO_robotallobj[test])/result_primitives[0];
                    }
                }
                // ------------------------------
                // ---     Test with GA       ---
                // ------------------------------
                if (Optim_information.GA == true) {
                    if (Optim_information.caps1obj == true) {
                        auto start_GA = get_tick_us();
                        result_GA_caps1obj[test] = test_collision_ga_OBB(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_GA = get_tick_us();
                        times_GA_caps1obj[test] = (stop_GA - start_GA) / 1000.0;
                        error_GA_caps1obj[test] = (result_primitives[0] - result_GA_caps1obj[test])/result_primitives[0];
                    }
                    if (Optim_information.capsallobj == true) {
                        auto start_GA = get_tick_us();
                        result_GA_capsallobj[test] = test_collision_ga_world_1caps(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_GA = get_tick_us();
                        times_GA_capsallobj[test] = (stop_GA - start_GA) / 1000.0;
                        error_GA_capsallobj[test] = (result_primitives[0] - result_GA_capsallobj[test])/result_primitives[0];
                    }
                    if (Optim_information.robotallobj == true) {
                        auto start_GA = get_tick_us();
                        result_GA_robotallobj[test] = test_collision_ga_world_full_robot(cart_pos, world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_GA = get_tick_us();
                        times_GA_robotallobj[test] = (stop_GA - start_GA) / 1000.0;
                        error_GA_robotallobj[test] = (result_primitives[0] - result_GA_robotallobj[test])/result_primitives[0];
                    }
                }
            }

            // Output information
            csvFile << Optim_information.N_individuals[i] << csv_sep << Optim_information.N_iterations[j] << csv_sep << result_primitives[0] << csv_sep << times_prim << csv_sep;
            if (Optim_information.PSO == true) {
                if (Optim_information.caps1obj == true) {
                    csvFile << sum(times_PSO_caps1obj)/n_tests << csv_sep << sum(error_PSO_caps1obj)/n_tests << csv_sep;
                }
                if (Optim_information.capsallobj == true) {
                    csvFile << sum(times_PSO_capsallobj)/n_tests << csv_sep << sum(error_PSO_capsallobj)/n_tests << csv_sep;
                }
                if (Optim_information.robotallobj == true) {
                    csvFile << sum(times_PSO_robotallobj)/n_tests << csv_sep << sum(error_PSO_robotallobj)/n_tests << csv_sep;
                }
            }
            if (Optim_information.GWO == true) {
                if (Optim_information.caps1obj == true) {
                    csvFile << sum(times_GWO_caps1obj)/n_tests << csv_sep << sum(error_GWO_caps1obj)/n_tests << csv_sep;
                }
                if (Optim_information.capsallobj == true) {
                    csvFile << sum(times_GWO_capsallobj)/n_tests << csv_sep << sum(error_GWO_capsallobj)/n_tests << csv_sep;
                }
                if (Optim_information.robotallobj == true) {
                    csvFile << sum(times_GWO_robotallobj)/n_tests << csv_sep << sum(error_GWO_robotallobj)/n_tests << csv_sep;
                }
            }
            if (Optim_information.GA == true) {
                if (Optim_information.caps1obj == true) {
                    csvFile << sum(times_GA_caps1obj)/n_tests << csv_sep << sum(error_GA_caps1obj)/n_tests << csv_sep;
                }
                if (Optim_information.capsallobj == true) {
                    csvFile << sum(times_GA_capsallobj)/n_tests << csv_sep << sum(error_GA_capsallobj)/n_tests << csv_sep;
                }
                if (Optim_information.robotallobj == true) {
                    csvFile << sum(times_GA_robotallobj)/n_tests << csv_sep << sum(error_GA_robotallobj)/n_tests << csv_sep;
                }
            }

            csvFile << std::endl;
        }
    }
    // note : This is not used for now but might be in the future for statistics
    // // -----------------------------------------
    // // ---     Quartiles Classification      ---
    // // -----------------------------------------
    // // Open csv file again to overwrite
    // // std::ofstream csvFile(Optim_information.csvOutput);
    // csvFile.clear();
    // // Write headers to the CSV file
    // csvFile << "Test Number" << csv_sep << "PSO rel. error" << csv_sep << "GWO rel. error" << csv_sep << "Quartile PSO" << csv_sep << "Quartile GWO" << csv_sep << std::endl;
    // // Sort error vectors
    // std::sort(error_PSO_values.begin(), error_PSO_values.end());
    // std::sort(error_GWO_values.begin(), error_GWO_values.end());
    // // Calculate quartiles
    // real Q0_GWO = 0;
    // real Q1_GWO = error_GWO_values[n_tests / 4];
    // real Q2_GWO = error_GWO_values[n_tests / 2];
    // real Q3_GWO = error_GWO_values[(3 * n_tests) / 4];
    // real Q0_PSO = 0;
    // real Q1_PSO = error_PSO_values[n_tests / 4];
    // real Q2_PSO = error_PSO_values[n_tests / 2];
    // real Q3_PSO = error_PSO_values[(3 * n_tests) / 4];
    // // Assign quartile labels
    // for (int test = 0; test < n_tests; test ++) {
    //     real error_GWO = error_GWO_values[test];
    //     real error_PSO = error_PSO_values[test];
    //     std::string quartile_GWO, quartile_PSO;
    //     if (error_GWO <= Q0_GWO) {
    //         quartile_GWO = "Q0";
    //     } else if (error_GWO <= Q1_GWO) {
    //         quartile_GWO = "Q1";
    //     } else if (error_GWO <= Q2_GWO) {
    //         quartile_GWO = "Q2";
    //     } else if (error_GWO <= Q3_GWO) {
    //         quartile_GWO = "Q3";
    //     } else {
    //         quartile_GWO = "Q4";
    //     }
    //     if (error_PSO <= Q0_PSO) {
    //         quartile_PSO = "Q0";
    //     } else if (error_PSO <= Q1_PSO) {
    //         quartile_PSO = "Q1";
    //     } else if (error_PSO <= Q2_PSO) {
    //         quartile_PSO = "Q2";
    //     } else if (error_PSO <= Q3_PSO) {
    //         quartile_PSO = "Q3";
    //     } else {
    //         quartile_PSO = "Q4";
    //     }
    //     csvFile << test << csv_sep << error_PSO << csv_sep << error_GWO << csv_sep << quartile_PSO << csv_sep << quartile_GWO << csv_sep << std::endl;
    // }
}

int main() {
    using namespace blast;
    OptimTest Optim_information;
    objlist world;

    //  --- CREATING WORLD ---

    // Figure box (from the project)
    // add_box({0.9, -0.7, 0.1}, {0.1, 0.1, 0.1}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world);

    // DEMO 1 - COMMENT ONE BOX TYPE (simple or complex)
    // // Bins for demo 1 (simple)
    // add_box({0.4114, -0.3784, 0.01785}, {0.2125, 0.155, 0.075}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world);
    // add_box({0.8716, -0.3784, 0.01785}, {0.2125, 0.155, 0.075}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world);

    // Bins for demo 1 (complex)
    // bin 1
    add_box({0.4114, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //bottom
    add_box({0.2114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //front
    add_box({0.6114, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //back
    add_box({0.4114, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //left
    add_box({0.4114, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //right

    // // bin 2
    add_box({0.8716, -0.3784, -0.04465}, {0.2125, 0.155, 0.0125}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //bottom
    add_box({0.6716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //front
    add_box({1.0716, -0.3784, 0.03035}, {0.0125, 0.155, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //back
    add_box({0.8716, -0.2359, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //left
    add_box({0.8716, -0.5209, 0.03035}, {0.1875, 0.0125, 0.0625}, Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1), &world); //right

    Optim_information.world = &world;

    //  --- SETTING OPTIMIZATION PARAMETERS ---
    Array individuals = {10, 20, 50, 100, 150, 250};
    Optim_information.N_individuals = individuals;
    Array iterations = {10, 100, 500};
    Optim_information.N_iterations = iterations;

    //  --- DECLARING ALGORITHMS THAT WILL BE TESTED ---
    Optim_information.PSO = TRUE;
    // Optim_information.GWO = TRUE;
    // Optim_information.GA = TRUE;
    
    //  --- DECLARING TYPE OF OPTIMIZATION ---
    Optim_information.caps1obj = TRUE;
    Optim_information.capsallobj = TRUE;
    // Optim_information.robotallobj = TRUE;


    //  --- SETTING INPUT / OUTPUT PARAMETERS ---
    Optim_information.csv_sep = ";";
    Optim_information.N_tests = 100;
    Optim_information.csvInput = "examples/build/trajectory_100_pts.csv";
    Optim_information.csvOutput = "C:/Users/thoma/Desktop/data_collision_check_traj1.csv";
    Array dimensions = {6, 100};    // SIZE OF THE CSV INPUT 
    Optim_information.input_dimensions = dimensions;

    //  --- TESTING ALGORITHMS ---
    test_algorithms(Optim_information);

    return 0;
}

