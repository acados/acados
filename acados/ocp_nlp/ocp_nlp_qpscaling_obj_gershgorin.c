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


#include "acados/ocp_nlp/ocp_nlp_qpscaling_obj_gershgorin.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "acados/ocp_nlp/ocp_nlp_qpscaling_common.h"
#include "acados/utils/math.h"

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"



/************************************************
 * opts
 ************************************************/
typedef struct
{
    double ub_max_abs_eig; // upper bound
    double ub_norm_inf_grad_obj; // upper bound
    double lb_norm_inf_grad_obj; // lower bound
} ocp_nlp_qpscaling_obj_gershgorin_opts;


acados_size_t ocp_nlp_qpscaling_obj_gershgorin_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_qpscaling_obj_gershgorin_opts);
}



void *ocp_nlp_qpscaling_obj_gershgorin_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_qpscaling_obj_gershgorin_opts_initialize_default(void *config_, ocp_nlp_qpscaling_dims *dims, void *opts_)
{
    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;

    opts->lb_norm_inf_grad_obj = 1e-4;
    opts->ub_norm_inf_grad_obj = 1e2;
    opts->ub_max_abs_eig = 1e5;

    return;
}

void ocp_nlp_qpscaling_obj_gershgorin_opts_set(void *config_, void *opts_, const char *field, void* value)
{

    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;

    if (!strcmp(field, "ub_max_abs_eig"))
    {
        double *d_ptr = value;
        opts->ub_max_abs_eig = *d_ptr;
        printf("ub_max_eig_abs: %2.e\n", opts->ub_max_abs_eig);
    }
    else if (!strcmp(field, "ub_norm_inf_grad_obj"))
    {
        double *d_ptr = value;
        opts->ub_norm_inf_grad_obj = *d_ptr;
        printf("ub_norm_inf_grad_obj : %2.e\n", opts->ub_norm_inf_grad_obj);
    }
    else if (!strcmp(field, "lb_norm_inf_grad_obj"))
    {
        double *d_ptr = value;
        opts->lb_norm_inf_grad_obj = *d_ptr;
        printf("lb_norm_inf_grad_obj : %2.e\n", opts->lb_norm_inf_grad_obj);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_qpscaling_obj_gershgorin_opts_set\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_qpscaling_obj_gershgorin_memory_calculate_size(void *config_, ocp_nlp_qpscaling_dims *dims, void *opts_)
{

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_qpscaling_obj_gershgorin_memory);

    size += (2*(N+1)+N)*sizeof(struct blasfeo_dmat *); // RSQrq BAbt DCt
    size += (3*(N+1)+2*N)*sizeof(struct blasfeo_dvec *); // rq b ux pi lam
    size += (N+1)*sizeof(int *); // idxb

    return size;
}



void *ocp_nlp_qpscaling_obj_gershgorin_memory_assign(void *config_, ocp_nlp_qpscaling_dims *dims, void *opts_, void *raw_memory)
{
    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_qpscaling_obj_gershgorin_memory *mem = (ocp_nlp_qpscaling_obj_gershgorin_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_qpscaling_obj_gershgorin_memory);

    assert((char *)mem + ocp_nlp_qpscaling_obj_gershgorin_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_qpscaling_obj_gershgorin_scale_qp(void *config, ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    double max_abs_eig = 0.0;
    double tmp;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *ns = qp_in->dim->ns;
    ocp_nlp_qpscaling_obj_gershgorin_memory *memory = mem_;
    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;

    struct blasfeo_dmat *RSQrq = qp_in->RSQrq;
    for (int stage = 0; stage <= dims->N; stage++)
    {
        compute_gershgorin_max_abs_eig_estimate(nx[stage]+nu[stage], RSQrq+stage, &tmp);
        max_abs_eig = fmax(max_abs_eig, tmp);
        // take Z into account
        blasfeo_dvecnrm_inf(2*ns[stage], qp_in->Z+stage, 0, &tmp);
        max_abs_eig = fmax(max_abs_eig, tmp);
    }

    if (max_abs_eig < opts->ub_max_abs_eig)
    {
        max_abs_eig = 1.0;
    }
    memory->obj_factor = 1.0 / max_abs_eig;


    // scale QP cost
    // print_ocp_qp_in(qp_in);
    if (memory->obj_factor != 1.0)
    {
        ocp_qp_scale_objective(qp_in, memory->obj_factor);
    }
    // printf("AFTER SCALING\n");
    // print_ocp_qp_in(qp_in);
    // printf("ocp_nlp_qpscaling_obj_gershgorin_scale_qp: computed obj_factor = %e\n", memory->obj_factor);
    // exit(1);
}


void ocp_nlp_qpscaling_obj_gershgorin_rescale_solution(void *config, ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    ocp_nlp_qpscaling_obj_gershgorin_memory *memory = mem_;

    // printf("ocp_nlp_qpscaling_obj_gershgorin_rescale_solution: qp_out before rescaling\n");
    // print_ocp_qp_out(qp_out);
    if (memory->obj_factor != 1.0)
    {
        ocp_qp_out_scale_duals(qp_out, 1/memory->obj_factor);
        ocp_qp_scale_objective(qp_in, 1/memory->obj_factor);
    }
    // printf("ocp_nlp_qpscaling_obj_gershgorin_rescale_solution: qp_out after rescaling\n");
    // print_ocp_qp_out(qp_out);
    // printf("ocp_nlp_qpscaling_obj_gershgorin_rescale_solution: using obj_factor = %e\n", memory->obj_factor);

    return;
}


void ocp_nlp_qpscaling_obj_gershgorin_config_initialize_default(ocp_nlp_qpscaling_config *config)
{
    // dims
    config->dims_calculate_size = &ocp_nlp_qpscaling_dims_calculate_size;
    config->dims_assign = &ocp_nlp_qpscaling_dims_assign;
    config->dims_set = &ocp_nlp_qpscaling_dims_set;
    // opts
    config->opts_calculate_size = &ocp_nlp_qpscaling_obj_gershgorin_opts_calculate_size;
    config->opts_assign = &ocp_nlp_qpscaling_obj_gershgorin_opts_assign;
    config->opts_initialize_default = &ocp_nlp_qpscaling_obj_gershgorin_opts_initialize_default;
    config->opts_set = &ocp_nlp_qpscaling_obj_gershgorin_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_qpscaling_obj_gershgorin_memory_calculate_size;
    config->memory_assign = &ocp_nlp_qpscaling_obj_gershgorin_memory_assign;
    // functions
    config->scale_qp = &ocp_nlp_qpscaling_obj_gershgorin_scale_qp;
    config->rescale_solution = &ocp_nlp_qpscaling_obj_gershgorin_rescale_solution;
}

