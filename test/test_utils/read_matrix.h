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

#ifndef TEST_TEST_UTILS_READ_MATRIX_H_
#define TEST_TEST_UTILS_READ_MATRIX_H_

#include <string>

#include "acados/utils/types.h"
#include "test/test_utils/eigen.h"

Eigen::MatrixXd readMatrix(const std::string filename);

Eigen::MatrixXd readMatrixFromFile(std::string filename, int_t rows, int_t cols);

Eigen::VectorXd readVectorFromFile(std::string filename, int_t length);

#endif  // TEST_TEST_UTILS_READ_MATRIX_H_
