#include "read_matrix.hpp"
#include <cctype>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <functional>
#include <locale>
#include <exception>

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
Eigen::MatrixXd readMatrix(const char *filename) {
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename);
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
