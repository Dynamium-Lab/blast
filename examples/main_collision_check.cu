#include <stdio.h>
#include "blast.hpp"
#include "blast_math.hpp"
#include "nlopt.h"
#include "json.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace blast;
Matrix read_csv_matrix (const std::string& filename) { 
    std::ifstream file(filename);
    Assert(!file.is_open()) 
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

// Array read_csv_array (const std::string& filename) { 
//     std::ifstream file(filename);
//     Assert(file.is_open()) 
//     Array pos(7*15996);
//     std::string line;
//    if (std::getline(file, line)) {
//         std::istringstream iss(line);
//         for (int i = 0; i < pos.size; i++) {
//             std::string value;
//             if (std::getline(iss, value, ',')) {
//                 pos[i] = std::stod(value);
//             }
//             else {
//                 break;
//             };
//         }
//     }
//     return pos;
// }

int main() {
    using namespace blast;
    // Array posi(7*15995);
    Matrix matrice;
    matrice = read_csv_matrix("trajectories.csv");
    //posi = read_csv_array("trajectories.csv");
    return 0;
}




    // Assert(!file.is_open()) 
    // std::vector<std::vector<std::string>> data; // 2D vector to store CSV data
    // std::string line;
    // while (std::getline(file, line)) {
    //     std::vector<std::string> row;
    //     std::stringstream ss(line);
    //     std::string cell;

    //     while (std::getline(ss, cell, ',')) {
    //         row.push_back(cell);
    //     }

    //     data.push_back(row);
    // }
    // // for (const auto& row : data) {
    // //     for (const auto& cell : row) {
    // //         std::cout << cell << "\t";
    // //     }
    // //     std::cout << std::endl;
    // // }
    // file.close();
    // return data;