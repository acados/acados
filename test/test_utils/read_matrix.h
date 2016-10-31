#ifndef TEST_TEST_UTILS_READ_MATRIX_H_
#define TEST_TEST_UTILS_READ_MATRIX_H_

#include <string>
#include "test/test_utils/eigen.h"
#include "acados/utils/types.h"

Eigen::MatrixXd readMatrix(const std::string filename);

Eigen::MatrixXd readMatrixFromFile(std::string filename, int_t rows, int_t cols);

Eigen::VectorXd readVectorFromFile(std::string filename, int_t length);

#endif  // TEST_TEST_UTILS_READ_MATRIX_H_
