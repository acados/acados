/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
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
