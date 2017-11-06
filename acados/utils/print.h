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

#ifndef ACADOS_UTILS_PRINT_H_
#define ACADOS_UTILS_PRINT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

void print_matrix(char *file_name, const real_t *matrix, const int_t nrows,
                  const int_t ncols);

void print_matrix_name(char *file_name, char *name, const real_t *matrix,
                       const int_t nrows, const int_t ncols);

void print_int_matrix(char *file_name, const int_t *matrix, const int_t nrows,
                      const int_t ncols);

void print_array(char *file_name, real_t *array, int_t size);

void print_int_array(char *file_name, const int_t *array, int_t size);

void read_matrix(const char *file_name, real_t *array, const int_t nrows,
                 const int_t ncols);

void print_ocp_qp(ocp_qp_in *qp);

void print_ocp_qp_to_file(ocp_qp_in *qp);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_PRINT_H_
