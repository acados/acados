/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "test/test_utils/read_matrix.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <locale>
#include <string>

#define MAXBUFSIZE  ((int) 1e6)

// Trim from start (in place)
static inline void ltrim(std::string *s) {
    s->erase(s->begin(), std::find_if(s->begin(), s->end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// Trim from end (in place)
static inline void rtrim(std::string *s) {
    s->erase(std::find_if(s->rbegin(), s->rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s->end());
}

// Trim from both ends (in place)
static inline void trim(std::string *s) {
    ltrim(s);
    rtrim(s);
}

// Read space delimited file into Eigen matrix
Eigen::MatrixXd readMatrix(const std::string filename) {
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename.c_str());
    if (!infile.is_open()) {
        std::string err_message;
        err_message + "Unable to open file " + filename;
        std::cout << err_message << std::endl;
        throw new std::runtime_error(err_message);
    }
    while (infile.is_open() && !infile.eof()) {
        std::string line;
        getline(infile, line);
        trim(&line);
        int temp_cols = 0;
        std::stringstream stream(line);
        while (!stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i, j) = buff[ cols*i+j ];

    return result;
}

Eigen::MatrixXd readMatrixFromFile(std::string filename, int_t rows, int_t cols) {
    Eigen::MatrixXd M = readMatrix(filename);
    return M;
}

Eigen::VectorXd readVectorFromFile(std::string filename, int_t length) {
    Eigen::VectorXd v = readMatrix(filename);
    return v;
}
