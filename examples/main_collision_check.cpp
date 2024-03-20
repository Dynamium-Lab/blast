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
Matrix read_csv_matrix (const std::string& filename) { 
    std::ifstream file(filename);
    Assert(!file.is_open());
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

int main() {
    using namespace blast;
    // ------------------------------
    // ---     Initialization     ---
    // ------------------------------
    int n_collision_results = 1;
    const int n_tests = 100;
    bool write_csv = TRUE;
    std::string csv_sep = ";";
    Gen3_7DOF manip;

    // Reading trajectory
    Matrix joint_pos;
    joint_pos = read_csv_matrix("trajectory1.csv");
    int njoint = joint_pos.rows;
    int npoint = joint_pos.cols;

    // Creating world
    objlist world;
    add_OBB({0.9, -0.7, 0.1}, {0.1, 0.1, 0.1}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, &world);
    real r = 0.00875;
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

    // Do multiple tests for statistics
    long times_prim[n_tests];
    long times_GJK[n_tests];
    long times_PSO[n_tests];
    long times_GWO[n_tests];
    for (int test = 0; test < n_tests; test ++) {
        // ------------------------------
        // ---  Test with primitives  ---
        // ------------------------------
        auto start_prim = get_tick_us();
        std::vector<real> result_primitves = test_collision(&capsules, &world, n_collision_results);
        auto stop_prim = get_tick_us();
        times_prim[test] = (stop_prim - start_prim) / 1000.0;

        // todo: fix GJK
        // ------------------------------
        // ---     Test with GJK      ---
        // ------------------------------
        auto start_GJK = get_tick_us();
        // std::vector<real> result_GJK = test_collision_GJK(&capsules, &world, n_collision_results);
        auto stop_GJK = get_tick_us();
        times_GJK[test] = (stop_GJK - start_GJK) / 1000.0;
        
        // ------------------------------
        // ---     Test with PSO      ---
        // ------------------------------
        auto start_PSO = get_tick_us();
        real result_PSO = test_collision_pso(cart_pos, &world) - r;
        auto stop_PSO = get_tick_us();
        times_PSO[test] = (stop_PSO - start_PSO) / 1000.0;

        // ------------------------------
        // ---     Test with GWO      ---
        // ------------------------------
        auto start_GWO = get_tick_us();
        real result_GWO = test_collision_gwo(cart_pos, &world) - r;
        auto stop_GWO = get_tick_us();
        times_GWO[test] = (stop_GWO - start_GWO) / 1000.0;

        // ------------------------------
        // ---     Test with GA       ---
        // ------------------------------
        // ------------------------------
        // ---     Test with SQP      ---
        // ------------------------------
    }

    if (write_csv) {
        #include <fstream>
        // Open a CSV file for writing
        std::ofstream csvFile("C:/Users/thoma/Desktop/data_collision_check.csv");
        // Check if the file is open
        if (!csvFile.is_open()) {
            std::cerr << "Error opening file." << std::endl;
        }
        // Write headers to the CSV file
        csvFile << "Test" << csv_sep << "Primitive test time (ms)" << csv_sep << "GJK test time (ms)" << csv_sep << "PSO test time (ms)" <<  csv_sep << "GWO test time (ms)" << std::endl;
        // Write data
        for (int i = 0; i < n_tests; i++) {
            csvFile << i+1 << csv_sep << times_prim[i] << csv_sep << times_GJK[i] << csv_sep << times_PSO[i] << csv_sep << times_GWO[i] << std::endl;
        }
        // Close the file
        csvFile.close();
        std::cout << "CSV file created successfully." << std::endl;
    }

    return 0;
}

