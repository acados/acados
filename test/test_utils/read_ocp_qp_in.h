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

#ifndef TEST_TEST_UTILS_READ_OCP_QP_IN_H_
#define TEST_TEST_UTILS_READ_OCP_QP_IN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

int_t read_int_vector_from_txt(int_t *vec, int_t n, const char *filename);
int_t read_double_vector_from_txt(real_t *vec, int_t n, const char *filename);
int_t read_double_matrix_from_txt(real_t *mat, int_t m, int_t n, const char *filename);
int_t write_double_vector_to_txt(real_t *vec, int_t n, const char *fname);
int_t write_int_vector_to_txt(int_t *vec, int_t n, const char *fname);

void print_ocp_qp_in(ocp_qp_in const in);

ocp_qp_in *read_ocp_qp_in(const char *fpath_, int_t BOUNDS, int_t INEQUALITIES, int_t MPC,
                          int_t QUIET);

void write_ocp_qp_in_to_txt(ocp_qp_in *const in, const char *dir);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* TEST_TEST_UTILS_READ_OCP_QP_IN_H_ */
