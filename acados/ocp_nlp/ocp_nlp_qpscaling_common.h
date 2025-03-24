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


/// \ingroup ocp_nlp
/// @{

/// \defgroup ocp_nlp_reg ocp_nlp_reg
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_qpscaling_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_qpscaling_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"



/* dims */

// same as qp_dims
typedef struct
{
    int *nx;
    int *nu;
    int N;
} ocp_nlp_qpscaling_dims;

//
acados_size_t ocp_nlp_qpscaling_dims_calculate_size(int N);
//
ocp_nlp_qpscaling_dims *ocp_nlp_qpscaling_dims_assign(int N, void *raw_memory);
//
void ocp_nlp_qpscaling_dims_set(void *config_, ocp_nlp_qpscaling_dims *dims, int stage, char *field, int* value);



/* config */

typedef struct
{
    /* dims */
    acados_size_t (*dims_calculate_size)(int N);
    ocp_nlp_qpscaling_dims *(*dims_assign)(int N, void *raw_memory);
    void (*dims_set)(void *config, ocp_nlp_qpscaling_dims *dims, int stage, char *field, int *value);
    /* opts */
    acados_size_t (*opts_calculate_size)(void);
    void *(*opts_assign)(void *raw_memory);
    void (*opts_initialize_default)(void *config, ocp_nlp_qpscaling_dims *dims, void *opts);
    void (*opts_set)(void *config, void *opts, const char *field, void* value);
    /* memory */
    acados_size_t (*memory_calculate_size)(void *config, ocp_nlp_qpscaling_dims *dims, void *opts);
    void *(*memory_assign)(void *config, ocp_nlp_qpscaling_dims *dims, void *opts, void *raw_memory);
    /* functions */
    void (*scale_qp)(void *config, ocp_nlp_qpscaling_dims *dims, void *opts, void *memory, ocp_qp_in *qp_in);
    void (*rescale_solution)(void *config, ocp_nlp_qpscaling_dims *dims, void *opts, void *memory, ocp_qp_in *qp_in, ocp_qp_out *qp_out);
} ocp_nlp_qpscaling_config;

//
acados_size_t ocp_nlp_qpscaling_config_calculate_size(void);
//
void *ocp_nlp_qpscaling_config_assign(void *raw_memory);

void ocp_qp_scale_objective(ocp_qp_in *qp_in, double factor);
//
void compute_gershgorin_max_abs_eig_estimate(int n, struct blasfeo_dmat *A, double *out);
//
void ocp_qp_out_scale_duals(ocp_qp_out *qp_out, double factor);

#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_qpscaling_COMMON_H_
/// @}
/// @}
