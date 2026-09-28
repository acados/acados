/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
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
