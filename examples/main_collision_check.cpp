#include <stdio.h>
#include "blast.hpp"
#include "blast_math.hpp"
#include "blast_manipulator.hpp"
#include "collisions/Collision_pso.hpp"
#include "collisions/Collision_gwo.hpp"
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
struct OptimTest {
    std::string csvInput;
    std::string csvOutput;
    std::string csv_sep = ",";
    Array N_individuals;
    Array N_iterations;
    int N_tests = 0;
    bool PSO = FALSE;
    bool GWO = FALSE;
    objlist world;
};

Matrix read_csv_matrix (const std::string& filename) { 
    std::ifstream file(filename);
    Assert(file.is_open());
    Matrix pos(7, 15995);
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    for (int i = 0; i < pos.rows; i++) {
        for (int k = 0; k < pos.cols; k++) {
            std::string value;
            if (std::getline(iss, value, ',')) {
                pos(i,k)= std::stod(value);
            } 
        }
    }
    return pos;
}

void test_algorithms(OptimTest& Optim_information) { 
    // ------------------------------
    // ---     Initialization     ---
    // ------------------------------
    int n_collision_results = 1;
    Gen3_7DOF manip;
    real r = 0.00875;
    // std::ofstream* csvFile = &Optim_information.csvOutput;
    std::string csv_sep = Optim_information.csv_sep;
    const int n_tests = Optim_information.N_tests;
    objlist world = Optim_information.world;

    // Reading trajectory
    Matrix joint_pos;
    joint_pos = read_csv_matrix(Optim_information.csvInput);
    int njoint = joint_pos.rows;
    int npoint = joint_pos.cols;
    // ---  Initialization ends   ---
    
    // Get cartesian positions
    Matrix cart_pos = manip.robot_capsules(joint_pos, 1);

    // Create capsules
    capslist capsules;
    capsules.caps.reserve(npoint*(njoint-1));
    for (int i = 0; i < npoint; i++) {
        for (int j = 0; j < cart_pos.rows/3 - 1; j++) {
            capsule temp;
            temp.p1 = { cart_pos(3*j, i),   cart_pos(3*j+1, i), cart_pos(3*j+2, i) };
            temp.p1 = { cart_pos(3*j+3, i), cart_pos(3*j+4, i), cart_pos(3*j+5, i) };
            temp.r = r;
            capsules.caps.push_back(temp);
        }
    }

    // Open csv file
    std::ofstream csvFile(Optim_information.csvOutput); 
    csvFile.clear();
    Assert(csvFile.is_open());

    // Write headers to the CSV file
    if (Optim_information.GWO == true && Optim_information.PSO == true) {
        csvFile << "Test Number" << csv_sep << "Number of Particles "<< csv_sep << "Number of iterations" << csv_sep << "Primitive Distance" << csv_sep << "Primitive test time (ms)" << csv_sep << "PSO test time (ms)" << csv_sep << "GWO test time (ms)"  <<  csv_sep << "PSO rel. error" << csv_sep << "GWO rel. error" << csv_sep << std::endl;
    }
    else {
        if (Optim_information.GWO == true)
            csvFile << "Test Number" << csv_sep << "Number of Particles "<< csv_sep << "Number of iterations" << csv_sep << "Primitive Distance" << csv_sep << "Primitive test time (ms)" << csv_sep << "PSO test time (ms)"  <<  csv_sep << "PSO rel. error" << csv_sep << std::endl;
        if (Optim_information.GWO == true)
            csvFile << "Test Number" << csv_sep << "Number of Particles "<< csv_sep << "Number of iterations" << csv_sep << "Primitive Distance" << csv_sep << "Primitive test time (ms)" << csv_sep << "GWO test time (ms)"  <<  csv_sep << "GWO rel. error" << csv_sep << std::endl;
    }
    // Do multiple tests for statistics
    Array times_prim(n_tests);
    Array times_GJK(n_tests);
    Array times_PSO(n_tests);
    Array times_GWO(n_tests);
    for (int i = 0; i < Optim_information.N_individuals.size; i++) {
        for (int j = 0; j < Optim_information.N_iterations.size; j++) {
            for (int test = 0; test < n_tests; test ++) {
                // ------------------------------
                // ---  Test with primitives  ---
                // ------------------------------
                auto start_prim = get_tick_us();
                std::vector<real> result_primitives = test_collision(&capsules, &world, n_collision_results);
                auto stop_prim = get_tick_us();
                times_prim[test] = (stop_prim - start_prim) / 1000.0;

                // // todo: fix GJK
                // // ------------------------------
                // // ---     Test with GJK      ---
                // // ------------------------------
                // auto start_GJK = get_tick_us();
                // // std::vector<real> result_GJK = test_collision_GJK(&capsules, &world, n_collision_results);
                // auto stop_GJK = get_tick_us();
                // times_GJK[test] = (stop_GJK - start_GJK) / 1000.0;
                
                if (Optim_information.PSO == true && Optim_information.GWO == true) {
                    // ------------------------------
                    // ---     Test with PSO      ---
                    // ------------------------------
                    auto start_PSO = get_tick_us();
                    real result_PSO = test_collision_pso(cart_pos, &world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                    auto stop_PSO = get_tick_us();
                    times_PSO[test] = (stop_PSO - start_PSO) / 1000.0;
                    // ------------------------------
                    // ---     Test with GWO      ---
                    // ------------------------------
                    auto start_GWO = get_tick_us();
                    real result_GWO = test_collision_gwo(cart_pos, &world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                    auto stop_GWO = get_tick_us();
                    times_GWO[test] = (stop_GWO - start_GWO) / 1000.0;
                    csvFile << test << csv_sep << Optim_information.N_individuals[i] << csv_sep << Optim_information.N_iterations[j] << csv_sep << result_primitives[0] << csv_sep << times_prim[test] << csv_sep << times_PSO[test] << csv_sep <<times_GWO[test] << csv_sep << (result_primitives[0] - result_PSO)/result_primitives[0] << csv_sep << (result_primitives[0] - result_GWO)/result_primitives[0] << csv_sep << std::endl;
                }
                else {
                    if (Optim_information.PSO == true) {
                        // ------------------------------
                        // ---     Test with PSO      ---
                        // ------------------------------
                        auto start_PSO = get_tick_us();
                        real result_PSO = test_collision_pso(cart_pos, &world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_PSO = get_tick_us();
                        times_PSO[test] = (stop_PSO - start_PSO) / 1000.0;
                        csvFile << test << csv_sep << Optim_information.N_individuals[i] << csv_sep << Optim_information.N_iterations[j] << csv_sep << result_primitives[0] << csv_sep << times_prim[test] << csv_sep << times_PSO[test] << csv_sep << (result_primitives[0] - result_PSO)/result_primitives[0]  << csv_sep << std::endl;
                    }
                    if (Optim_information.GWO == true) {
                        // ------------------------------
                        // ---     Test with GWO      ---
                        // ------------------------------
                        auto start_GWO = get_tick_us();
                        real result_GWO = test_collision_gwo(cart_pos, &world, Optim_information.N_individuals[i], Optim_information.N_iterations[j]) - r;
                        auto stop_GWO = get_tick_us();
                        times_GWO[test] = (stop_GWO - start_GWO) / 1000.0;
                        csvFile << test << csv_sep << Optim_information.N_individuals[i] << csv_sep << Optim_information.N_iterations[j] << csv_sep << result_primitives[0] << csv_sep << times_prim[test] << csv_sep << times_GWO[test] << csv_sep << (result_primitives[0] - result_GWO)/result_primitives[0]  << csv_sep << std::endl;
                    }
                }
                // ------------------------------
                // ---     Test with GA       ---
                // ------------------------------
                // ------------------------------
                // ---     Test with SQP      ---
                // ------------------------------
            }
        }
    }

    csvFile.close();
    std::cout << "CSV file updated successfully." << std::endl;
}

int main() {
    using namespace blast;
    OptimTest Optim_information;
    objlist world;
    add_OBB({0.9, -0.7, 0.1}, {0.1, 0.1, 0.1}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world);
    Array individuals = {1, 2, 5, 10, 20, 50, 100, 150, 250};
    Array iterations = {1, 2, 3, 4, 5, 10, 15, 25, 50};
    Optim_information.N_individuals = individuals;
    Optim_information.N_iterations = iterations;
    Optim_information.PSO = TRUE;
    Optim_information.GWO = TRUE;
    Optim_information.world = world;
    Optim_information.csv_sep = ";";
    Optim_information.N_tests = 1;
    Optim_information.csvInput = "trajectory1.csv";
    Optim_information.csvOutput = "C:/Users/thoma/Desktop/data_collision_check.csv";
 
    test_algorithms(Optim_information);

    return 0;
}

