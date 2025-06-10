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


#ifndef ACADOS_OCP_NLP_OCP_NLP_qpscaling_OBJECTIVE_GERSHGORIN_H_
#define ACADOS_OCP_NLP_OCP_NLP_qpscaling_OBJECTIVE_GERSHGORIN_H_

#ifdef __cplusplus
extern "C" {
#endif


// blasfeo
#include "blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_qpscaling_common.h"



/************************************************
 * dims
 ************************************************/

// use the functions in ocp_nlp_qpscaling_common


/************************************************
 * options
 ************************************************/
typedef struct
{
    double ub_max_abs_eig; // upper bound
    double ub_norm_inf_grad_obj; // upper bound
    double lb_norm_inf_grad_obj; // lower bound

    bool scale_qp_objective;
    bool scale_qp_constraints;
} ocp_nlp_qpscaling_obj_gershgorin_opts;
// use all functions just through config pointers


typedef struct {
    double obj_factor;
    struct blasfeo_dvec *constraints_scaling_vec;
} ocp_nlp_qpscaling_obj_gershgorin_memory;

//
void ocp_nlp_qpscaling_obj_gershgorin_config_initialize_default(ocp_nlp_qpscaling_config *config);
//
void *ocp_nlp_qpscaling_obj_gershgorin_get_constraints_scaling_ptr(void *memory_, void* opts_);


#ifdef __cplusplus
}
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_qpscaling_OBJECTIVE_GERSHGORIN_H_
