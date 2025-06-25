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


#ifndef ACADOS_OCP_NLP_OCP_NLP_QPSCALING_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_QPSCALING_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"


/* dims */

// same as qp_dims
typedef struct
{
    ocp_qp_dims *qp_dim;
} ocp_nlp_qpscaling_dims;

//
acados_size_t ocp_nlp_qpscaling_dims_calculate_size(int N);
//
ocp_nlp_qpscaling_dims *ocp_nlp_qpscaling_dims_assign(int N, void *raw_memory);

/************************************************
 * options
 ************************************************/
typedef struct
{
    double ub_max_abs_eig;
    double lb_norm_inf_grad_obj;
    int print_level;

    qpscaling_scale_objective_type scale_qp_objective;
    ocp_nlp_qpscaling_constraint_type scale_qp_constraints;
} ocp_nlp_qpscaling_opts;
// use all functions just through config pointers


typedef struct {
    int status;
    double obj_factor;
    struct blasfeo_dvec *constraints_scaling_vec;
    ocp_qp_in *scaled_qp_in;
    ocp_qp_out *scaled_qp_out;
} ocp_nlp_qpscaling_memory;

acados_size_t ocp_nlp_qpscaling_opts_calculate_size(void);
void ocp_nlp_qpscaling_opts_initialize_default(ocp_nlp_qpscaling_dims *dims, void *opts_);
void *ocp_nlp_qpscaling_opts_assign(void *raw_memory);
void ocp_nlp_qpscaling_opts_set(void *opts_, const char *field, void* value);

//

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_qpscaling_memory_calculate_size(ocp_nlp_qpscaling_dims *dims, void *opts_, ocp_qp_dims *orig_qp_dim);
void *ocp_nlp_qpscaling_memory_assign(ocp_nlp_qpscaling_dims *dims, void *opts_, ocp_qp_dims *orig_qp_dim, void *raw_memory);
void *ocp_nlp_qpscaling_get_constraints_scaling_ptr(void *memory_, void* opts_);
void ocp_nlp_qpscaling_memory_get(ocp_nlp_qpscaling_dims *dims, void *mem_, const char *field, int stage, void* value);


/************************************************
 * functionality
 ************************************************/
void ocp_nlp_qpscaling_precompute(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out);
//
void ocp_nlp_qpscaling_scale_qp(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in);

void ocp_nlp_qpscaling_rescale_solution(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out);

#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_QPSCALING_COMMON_H_
