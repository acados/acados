/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
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

#include "catch/include/catch.hpp"

#define MAXBUFSIZE ((int) 1e6)

// Trim from start (in place)
static inline void ltrim(std::string *s)
{
    s->erase(s->begin(),
             std::find_if(s->begin(), s->end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// Trim from end (in place)
static inline void rtrim(std::string *s)
{
    s->erase(std::find_if(s->rbegin(), s->rend(), std::not1(std::ptr_fun<int, int>(std::isspace)))
                 .base(),
             s->end());
}

// Trim from both ends (in place)
static inline void trim(std::string *s)
{
    ltrim(s);
    rtrim(s);
}

// Read space delimited file into Eigen matrix
Eigen::MatrixXd readMatrix(const std::string &filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename.c_str());
    if (!infile.is_open())
    {
        std::string err_message;
        err_message + "Unable to open file " + filename;
        std::cout << err_message << std::endl;
        throw new std::runtime_error(err_message);
    }
    while (infile.is_open() && !infile.eof())
    {
        std::string line;
        getline(infile, line);
        trim(&line);
        int temp_cols = 0;
        std::stringstream stream(line);
        while (!stream.eof()) stream >> buff[cols * rows + temp_cols++];

        if (temp_cols == 0) continue;

        if (cols == 0) cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) result(i, j) = buff[cols * i + j];

    return result;
}

Eigen::MatrixXd readMatrixFromFile(const std::string &filename, int_t rows, int_t cols)
{
    Eigen::MatrixXd M = readMatrix(filename);
    REQUIRE(M.rows() == rows);
    REQUIRE(M.cols() == cols);
    return M;
}

Eigen::VectorXd readVectorFromFile(const std::string &filename, int_t length)
{
    Eigen::VectorXd v = readMatrix(filename);
    REQUIRE(v.rows() == length);
    return v;
}
