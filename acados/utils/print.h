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

#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common_frontend.h"
#include "acados/utils/types.h"

void print_matrix(char *file_name, const real_t *matrix, const int_t nrows, const int_t ncols);

void print_matrix_name(char *file_name, char *name, const real_t *matrix, const int_t nrows,
                       const int_t ncols);

void print_int_matrix(char *file_name, const int_t *matrix, const int_t nrows, const int_t ncols);

void print_array(char *file_name, real_t *array, int_t size);

void print_int_array(char *file_name, const int_t *array, int_t size);

void read_matrix(const char *file_name, real_t *array, const int_t nrows, const int_t ncols);

void write_double_vector_to_txt(real_t *vec, int_t n, const char *fname);

// ocp nlp
// TODO(andrea): inconsistent naming
void ocp_nlp_dims_print(ocp_nlp_dims *dims);
// TODO(andrea): inconsistent naming
void ocp_nlp_out_print(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out);
// TODO(andrea): inconsistent naming
void ocp_nlp_res_print(ocp_nlp_dims *dims, ocp_nlp_res *nlp_res);

// ocp qp
void print_ocp_qp_dims(ocp_qp_dims *dims);

void print_dense_qp_dims(dense_qp_dims *dims);

void print_ocp_qp_in(ocp_qp_in *qp_in);

void print_ocp_qp_out(ocp_qp_out *qp_out);

void print_ocp_qp_in_to_string(char string_out[], ocp_qp_in *qp_in);

void print_ocp_qp_out_to_string(char string_out[], ocp_qp_out *qp_out);

void print_ocp_qp_res(ocp_qp_res *qp_res);

void print_colmaj_ocp_qp_in(colmaj_ocp_qp_in *qp);

void print_colmaj_ocp_qp_in_to_file(colmaj_ocp_qp_in *qp);

void print_colmaj_ocp_qp_out(char *filename, colmaj_ocp_qp_in *qp, colmaj_ocp_qp_out *out);

void print_dense_qp_in(dense_qp_in *qp_in);

void print_ocp_qp_info(ocp_qp_info *info);

// blasfeo
void print_blasfeo_target();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_PRINT_H_
