#include <stdio.h>
#include "blast.hpp"
#include "blast_math.hpp"
#include "blast_manipulator.hpp"
#include "optimization/Collision_pso.hpp"
#include "collisions/world.h"
#include "nlopt.h"
#include "json.hpp"
#include <chrono>
#include <iostream>
#include <fstream>
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
    const int n_tests = 1;
    bool write_csv = TRUE;
    std::string csv_sep = ";";
    Gen3_7DOF manip;

    // Reading trajectory
    Matrix joint_pos;
    joint_pos = read_csv_matrix("trajectories.csv");
    int njoint = joint_pos.rows;
    int npoint = joint_pos.cols;

    // Creating world
    objlist* world;
    add_OBB({1, 1, 0}, {0.2, 0.2, 0.2}, {1, 0, 0, 0, 1, 0, 0, 0, 1}, world);
    int r = 0.0875;
    // ---  Initialization ends   ---
    
    // Get cartesian positions
    Matrix cart_pos = manip.robot_capsules(joint_pos, 1);

    // Create capsules
    capslist* capsules;
    (*capsules).caps.reserve(npoint*(njoint-1));
    for (int i = 0; i < npoint; i++) {
        for (int j = 0; j < njoint - 1; j++) {
            capsule temp;
            temp.p1 = { cart_pos(3*j, i),   cart_pos(3*j+1, i), cart_pos(3*j+2, i) };
            temp.p1 = { cart_pos(3*j+3, i), cart_pos(3*j+4, i), cart_pos(3*j+5, i) };
            temp.r = r;
            (*capsules).caps.push_back(temp);
        }
    }

    // Do multiple tests for statistics
    long times_prim[n_tests];
    long times_GJK[n_tests];
    long times_PSO[n_tests];
    std::chrono::duration<long, std::micro> int_usec;
    for (int test = 0; test < n_tests; test ++) {
        // ------------------------------
        // ---  Test with primitives  ---
        // ------------------------------
        
        auto start_prim = std::chrono::high_resolution_clock::now();
        std::vector<real> result_primitves = test_collision(capsules, world, n_collision_results);     
        auto stop_prim = std::chrono::high_resolution_clock::now();
        auto duration_prim = std::chrono::duration_cast<std::chrono::microseconds>(stop_prim - start_prim);
        int_usec = duration_prim;
        times_prim[test] = int_usec.count();

        // ------------------------------
        // ---     Test with GJK      ---
        // ------------------------------
        auto start_GJK = std::chrono::high_resolution_clock::now();
        std::vector<real> result_GJK = test_collision_GJK(capsules, world, n_collision_results);
        auto stop_GJK = std::chrono::high_resolution_clock::now();
        auto duration_GJK = std::chrono::duration_cast<std::chrono::microseconds>(stop_GJK - start_GJK);
        int_usec = duration_GJK;
        times_GJK[test] = int_usec.count();
        
        // ------------------------------
        // ---     Test with PSO      ---
        // ------------------------------
        auto start_PSO = std::chrono::high_resolution_clock::now();
        real result_pso = INF_REAL;
        for (int i = 0; i < world->OBBlist.size(); i++) {
            real result_tmp = collision_pso(cart_pos, world->OBBlist[i]);
            result_pso = result_tmp < result_pso ? result_tmp : result_pso;
        }
        auto stop_PSO = std::chrono::high_resolution_clock::now();
        auto duration_PSO = std::chrono::duration_cast<std::chrono::microseconds>(stop_PSO - start_PSO);
        int_usec = duration_PSO;
        times_PSO[test] = int_usec.count();

        // ------------------------------
        // ---     Test with GWO      ---
        // ------------------------------
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
            std::cerr << "Error opening file!" << std::endl;
        }
        // Write headers to the CSV file
        csvFile << "Test, Primitive test time (ms), GJK test time (ms), PSO test time (ms)" << std::endl;
        // Write data
        for (int i = 0; i < n_tests; i++) {
            // csvFile << i+1 << csv_sep << times_prim(i) << csv_sep << times_GJK(i) << csv_sep << times_PSO(i)<< std::endl;
            csvFile << i+1 << csv_sep << times_prim[i] << csv_sep << times_GJK[i]  << csv_sep << times_PSO[i] << std::endl;
        }
        // Close the file
        csvFile.close();
        std::cout << "CSV file created successfully." << std::endl;
    }

    return 0;
}
