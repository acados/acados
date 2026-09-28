/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
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

// void print_matrix(char *file_name, const real_t *matrix, const int_t nrows, const int_t ncols);

// void print_matrix_name(char *file_name, char *name, const real_t *matrix, const int_t nrows,
//                        const int_t ncols);

// void print_int_matrix(char *file_name, const int_t *matrix, const int_t nrows, const int_t ncols);

// void print_array(char *file_name, real_t *array, int_t size);

// void print_int_array(char *file_name, const int_t *array, int_t size);

void read_matrix(const char *file_name, real_t *array, const int_t nrows, const int_t ncols);

void write_double_vector_to_txt(real_t *vec, int_t n, const char *fname);

// ocp nlp
void print_ocp_nlp_dims(ocp_nlp_dims *dims);

void print_ocp_nlp_out(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out);

void print_ocp_nlp_res(ocp_nlp_dims *dims, ocp_nlp_res *nlp_res);

// ocp qp
void print_ocp_qp_dims(ocp_qp_dims *dims);

// void print_dense_qp_dims(dense_qp_dims *dims);

void print_ocp_qp_in(ocp_qp_in *qp_in);

void print_ocp_qp_in_to_file(FILE *file, ocp_qp_in *qp_in);

void print_ocp_qp_out(ocp_qp_out *qp_out);

void print_ocp_qp_out_to_file(FILE *file, ocp_qp_out *qp_out);

void print_ocp_qp_res(ocp_qp_res *qp_res);

void print_dense_qp_in(dense_qp_in *qp_in);
// void print_ocp_qp_in_to_string(char string_out[], ocp_qp_in *qp_in);

// void print_ocp_qp_out_to_string(char string_out[], ocp_qp_out *qp_out);

// void print_colmaj_ocp_qp_in(colmaj_ocp_qp_in *qp);

// void print_colmaj_ocp_qp_in_to_file(colmaj_ocp_qp_in *qp);

// void print_colmaj_ocp_qp_out(char *filename, colmaj_ocp_qp_in *qp, colmaj_ocp_qp_out *out);

void print_qp_info(qp_info *info);

// void acados_warning(char warning_string[]);

// void acados_error(char error_string[]);

// void acados_not_implemented(char feature_string[]);

// blasfeo
// void print_blasfeo_target();

void print_debug_output(char* message, int print_level, int required_print_level);
//
void print_debug_output_double(char* message, double value, int print_level, int required_print_level);

const char* status_to_string(return_values_t status);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_PRINT_H_
