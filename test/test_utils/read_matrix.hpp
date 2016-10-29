#ifndef TEST_TEST_UTILS_READ_MATRIX_HPP_
#define TEST_TEST_UTILS_READ_MATRIX_HPP_

#include <string>
#include "acados/utils/types.h"

/* Ignore compiler warnings from Eigen */
#if defined(__clang__)
#include <Eigen/Dense>
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#endif

Eigen::MatrixXd readMatrix(const std::string filename);

Eigen::MatrixXd readMatrixFromFile(std::string filename, int_t rows, int_t cols);

Eigen::VectorXd readVectorFromFile(std::string filename, int_t length);

#endif  // TEST_TEST_UTILS_READ_MATRIX_HPP_
