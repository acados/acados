#ifndef TEST_TEST_UTILS_READ_MATRIX_HPP_
#define TEST_TEST_UTILS_READ_MATRIX_HPP_

/* Ignore compiler warnings from Eigen */
#if defined(__clang__)
#include <Eigen/Dense>
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
#endif


Eigen::MatrixXd readMatrix(const char *filename);

#endif  // TEST_TEST_UTILS_READ_MATRIX_HPP_
