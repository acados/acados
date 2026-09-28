/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef TEST_TEST_UTILS_READ_MATRIX_H_
#define TEST_TEST_UTILS_READ_MATRIX_H_

#include <string>

#include "acados/utils/types.h"
#include "test/test_utils/eigen.h"

Eigen::MatrixXd readMatrix(const std::string &filename);

Eigen::MatrixXd readMatrixFromFile(const std::string &filename, int_t rows, int_t cols);

Eigen::VectorXd readVectorFromFile(const std::string &filename, int_t length);

#endif  // TEST_TEST_UTILS_READ_MATRIX_H_
