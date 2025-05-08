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
#include "acados/utils/mem.h"
#include "acados/utils/math.h"

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"


/************************************************
 * opts
 ************************************************/
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

    opts->scale_qp_objective = false;
    opts->scale_qp_dynamics = false;
    opts->scale_qp_constraints = false;

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
    else if (!strcmp(field, "scale_qp_objective"))
    {
        bool *d_ptr = value;
        opts->scale_qp_objective = *d_ptr;
    }
    else if (!strcmp(field, "scale_qp_dynamics"))
    {
        bool *d_ptr = value;
        opts->scale_qp_dynamics = *d_ptr;
    }
    else if (!strcmp(field, "scale_qp_constraints"))
    {
        bool *d_ptr = value;
        opts->scale_qp_constraints = *d_ptr;
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
    int i;
    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;
    // int *nx = dims->nx;
    // int *nu = dims->nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_qpscaling_obj_gershgorin_memory);

    // What is this used for?
    size += (2*(N+1)+N)*sizeof(struct blasfeo_dmat *); // RSQrq BAbt DCt
    size += (3*(N+1)+2*N)*sizeof(struct blasfeo_dvec *); // rq b ux pi lam
    size += (N+1)*sizeof(int *); // idxb

    // scaling vectors
    for (i = 0; i <= N; i++)
    {
        if (opts->scale_qp_constraints)
            size += 1 * blasfeo_memsize_dvec(2 * dims->ni[i]);
        if (i < N)
        {
            if (opts->scale_qp_dynamics)
                size += 1 * blasfeo_memsize_dvec(dims->nx[i + 1]);
        }
    }
    size += 1 * (N + 1) * sizeof(struct blasfeo_dvec);  // constraints_scaling_vec
    size += 1 * N * sizeof(struct blasfeo_dvec);  // dynamics_scaling_vec

    return size;
}



void *ocp_nlp_qpscaling_obj_gershgorin_memory_assign(void *config_, ocp_nlp_qpscaling_dims *dims, void *opts_, void *raw_memory)
{
    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;
    char *c_ptr = (char *) raw_memory;
    int N = dims->N;

    ocp_nlp_qpscaling_obj_gershgorin_memory *mem = (ocp_nlp_qpscaling_obj_gershgorin_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_qpscaling_obj_gershgorin_memory);

    assert((char *)mem + ocp_nlp_qpscaling_obj_gershgorin_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    if (opts->scale_qp_dynamics)
    {
        assign_and_advance_blasfeo_dvec_structs(dims->N, &mem->dynamics_scaling_vec, &c_ptr);
    }
    if (opts->scale_qp_constraints)
    {
        assign_and_advance_blasfeo_dvec_structs(dims->N + 1, &mem->constraints_scaling_vec, &c_ptr);
    }

    if (opts->scale_qp_dynamics)
    {
        for (int i = 0; i < dims->N; ++i)
        {
            assign_and_advance_blasfeo_dvec_mem(dims->nx[i + 1], mem->dynamics_scaling_vec + i, &c_ptr);
        }
    }
    if (opts->scale_qp_constraints)
    {
        for (int i = 0; i <= dims->N; ++i)
        {
            assign_and_advance_blasfeo_dvec_mem(2 * dims->ni[i], mem->constraints_scaling_vec + i, &c_ptr);
        }
    }

    // initialize with one
    if (opts->scale_qp_dynamics)
    {
        for(int i=0; i<dims->N; i++)
        {
            blasfeo_dvecse(dims->nx[i+1], 1.0, mem->dynamics_scaling_vec+i, 0);
        }

    }

    if (opts->scale_qp_constraints)
    {
        for(int i=0; i<dims->N; i++)
        {
            blasfeo_dvecse(dims->ni[i], 1.0, mem->constraints_scaling_vec+i, 0);
        }
        blasfeo_dvecse(dims->ni[N], 1.0, mem->constraints_scaling_vec+N, 0);
    }

    return mem;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_qpscaling_scale_qp_objective(void *config, ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    double max_abs_eig = 0.0;
    double tmp;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *ns = qp_in->dim->ns;
    ocp_nlp_qpscaling_obj_gershgorin_memory *memory = mem_;
    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;

    struct blasfeo_dmat *RSQrq = qp_in->RSQrq;
    double nrm_inf_grad_obj = 0.0;
    for (int stage = 0; stage <= dims->N; stage++)
    {
        compute_gershgorin_max_abs_eig_estimate(nx[stage]+nu[stage], RSQrq+stage, &tmp);
        max_abs_eig = fmax(max_abs_eig, tmp);
        // take Z into account
        blasfeo_dvecnrm_inf(2*ns[stage], qp_in->Z+stage, 0, &tmp);
        max_abs_eig = fmax(max_abs_eig, tmp);

        // norm gradient
        blasfeo_dvecnrm_inf(nx[stage]+nu[stage], qp_in->rqz+stage, 0, &tmp);
        nrm_inf_grad_obj = fmax(nrm_inf_grad_obj, fabs(tmp));
    }

    if (max_abs_eig < opts->ub_max_abs_eig)
    {
        max_abs_eig = 1.0;
    }
    else
    {
        // dividing by max_value helps that gradient does not get too small
        max_abs_eig = max_abs_eig/opts->ub_max_abs_eig;
    }
    memory->obj_factor = 1.0 / (max_abs_eig);
    printf("Scaling factor objective: %.2e\n",memory->obj_factor);

    // scale QP cost
    // print_ocp_qp_in(qp_in);
    if (memory->obj_factor != 1.0)
    {
        ocp_qp_scale_objective(qp_in, memory->obj_factor);
    }

    if (memory->obj_factor*nrm_inf_grad_obj <= opts->lb_norm_inf_grad_obj)
    {
        printf("lb_norm_inf_grad_obj violated! %.2e\n", opts->lb_norm_inf_grad_obj);
        printf("Gradient is very small! %.2e\n", memory->obj_factor*nrm_inf_grad_obj);
    }
    // printf("AFTER SCALING\n");
    // print_ocp_qp_in(qp_in);
    // printf("ocp_nlp_qpscaling_obj_gershgorin_scale_qp: computed obj_factor = %e\n", memory->obj_factor);
}

void ocp_nlp_qpscaling_scale_qp_dynamics(void *config, ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    printf("TBD");
}

void ocp_nlp_qpscaling_scale_qp_constraints(void *config, ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    printf("TBD");
}

void ocp_nlp_qpscaling_obj_gershgorin_scale_qp(void *config, ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;
    if (opts->scale_qp_objective)
    {
        ocp_nlp_qpscaling_scale_qp_objective(config, dims, opts_, mem_, qp_in);
    }
    if (opts->scale_qp_dynamics)
    {
        ocp_nlp_qpscaling_scale_qp_dynamics(config, dims, opts_, mem_, qp_in);
    }
    if (opts->scale_qp_constraints)
    {
        ocp_nlp_qpscaling_scale_qp_constraints(config, dims, opts_, mem_, qp_in);
    }
}


void ocp_nlp_qpscaling_obj_gershgorin_rescale_solution(void *config, ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    ocp_nlp_qpscaling_obj_gershgorin_memory *memory = mem_;
    ocp_nlp_qpscaling_obj_gershgorin_opts *opts = opts_;

    // printf("ocp_nlp_qpscaling_obj_gershgorin_rescale_solution: qp_out before rescaling\n");
    // print_ocp_qp_out(qp_out);
    if (memory->obj_factor != 1.0 && opts->scale_qp_objective)
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

