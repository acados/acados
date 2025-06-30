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


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "acados/utils/math.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"

#include "acados/ocp_nlp/ocp_nlp_qpscaling.h"

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
#include "blasfeo_d_aux_ext_dep.h"



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_nlp_qpscaling_dims_calculate_size(int N)
{
    acados_size_t size = sizeof(ocp_nlp_qpscaling_dims);
    // qp_dim
    size += sizeof(ocp_qp_dims);
    size += d_ocp_qp_dim_memsize(N);
    size += 8; // align
    make_int_multiple_of(8, &size);

    return size;
}



ocp_nlp_qpscaling_dims *ocp_nlp_qpscaling_dims_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // dims
    ocp_nlp_qpscaling_dims *dims = (ocp_nlp_qpscaling_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_qpscaling_dims);

    // qp_dim
    dims->qp_dim = (ocp_qp_dims *) c_ptr;
    c_ptr += sizeof(ocp_qp_dims);

    align_char_to(8, &c_ptr);

    // qp_dim
    d_ocp_qp_dim_create(N, dims->qp_dim, c_ptr);
    c_ptr += d_ocp_qp_dim_memsize(N);

    assert((char *) raw_memory + ocp_nlp_qpscaling_dims_calculate_size(N) >= c_ptr);

    return dims;
}



/************************************************
 * opts
 ************************************************/
acados_size_t ocp_nlp_qpscaling_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_qpscaling_opts);
}



void *ocp_nlp_qpscaling_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_qpscaling_opts_initialize_default(ocp_nlp_qpscaling_dims *dims, void *opts_)
{
    ocp_nlp_qpscaling_opts *opts = opts_;

    opts->lb_norm_inf_grad_obj = 1e-4;
    opts->ub_max_abs_eig = 1e5;
    opts->print_level = 0;

    opts->scale_qp_objective = NO_OBJECTIVE_SCALING;
    opts->scale_qp_constraints = NO_CONSTRAINT_SCALING;

    return;
}

void ocp_nlp_qpscaling_opts_set(void *opts_, const char *field, void* value)
{

    ocp_nlp_qpscaling_opts *opts = opts_;

    if (!strcmp(field, "ub_max_abs_eig"))
    {
        double *d_ptr = value;
        opts->ub_max_abs_eig = *d_ptr;
    }
    else if (!strcmp(field, "lb_norm_inf_grad_obj"))
    {
        double *d_ptr = value;
        opts->lb_norm_inf_grad_obj = *d_ptr;
    }
    else if (!strcmp(field, "scale_objective"))
    {
        qpscaling_scale_objective_type *d_ptr = value;
        opts->scale_qp_objective = *d_ptr;
    }
    else if (!strcmp(field, "scale_constraints"))
    {
        ocp_nlp_qpscaling_constraint_type *d_ptr = value;
        opts->scale_qp_constraints = *d_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_qpscaling_opts_set\n", field);
        exit(1);
    }

    return;

}


/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_qpscaling_memory_calculate_size(ocp_nlp_qpscaling_dims *dims, void *opts_, ocp_qp_dims *orig_qp_dim)
{
    int N = orig_qp_dim->N;
    int i;
    ocp_nlp_qpscaling_opts *opts = opts_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_qpscaling_memory);

    if (opts->scale_qp_objective != NO_OBJECTIVE_SCALING ||
        opts->scale_qp_constraints != NO_CONSTRAINT_SCALING)
    {
        size += ocp_qp_in_calculate_size(orig_qp_dim);
        size += ocp_qp_out_calculate_size(orig_qp_dim);
    }

    // constraints_scaling_vec
    if (opts->scale_qp_constraints)
    {
        size += (N + 1) * sizeof(struct blasfeo_dvec);  // constraints_scaling_vec
        for (i = 0; i <= N; i++)
        {
            size += blasfeo_memsize_dvec(orig_qp_dim->ng[i]);
        }
    }

    return size;
}



void *ocp_nlp_qpscaling_memory_assign(ocp_nlp_qpscaling_dims *dims, void *opts_, ocp_qp_dims *orig_qp_dim, void *raw_memory)
{
    ocp_nlp_qpscaling_opts *opts = opts_;
    char *c_ptr = (char *) raw_memory;
    int N = orig_qp_dim->N;

    // setup qp dimensions
    d_ocp_qp_dim_copy_all(orig_qp_dim, dims->qp_dim);

    ocp_nlp_qpscaling_memory *mem = (ocp_nlp_qpscaling_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_qpscaling_memory);

    mem->status = ACADOS_SUCCESS;

    if (opts->scale_qp_objective != NO_OBJECTIVE_SCALING ||
        opts->scale_qp_constraints != NO_CONSTRAINT_SCALING)
    {
        mem->scaled_qp_in = ocp_qp_in_assign(orig_qp_dim, c_ptr);
        c_ptr += ocp_qp_in_calculate_size(orig_qp_dim);
        mem->scaled_qp_out = ocp_qp_out_assign(orig_qp_dim, c_ptr);
        c_ptr += ocp_qp_out_calculate_size(orig_qp_dim);
    }

    if (opts->scale_qp_constraints)
    {
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->constraints_scaling_vec, &c_ptr);
        for (int i = 0; i <= N; ++i)
        {
            assign_and_advance_blasfeo_dvec_mem(orig_qp_dim->ng[i], mem->constraints_scaling_vec + i, &c_ptr);
            blasfeo_dvecse(orig_qp_dim->ng[i], 1.0, mem->constraints_scaling_vec+i, 0);
        }
    }
    assert((char *)mem + ocp_nlp_qpscaling_memory_calculate_size(dims, opts_, orig_qp_dim) >= c_ptr);

    return mem;
}

/************************************************
 * getter functions
 ************************************************/
void *ocp_nlp_qpscaling_get_constraints_scaling_ptr(void *memory_, void* opts_)
{
    ocp_nlp_qpscaling_memory *mem = memory_;
    ocp_nlp_qpscaling_opts *opts = opts_;
    if (opts->scale_qp_constraints)
        return mem->constraints_scaling_vec;
    else
        return NULL;
}


void ocp_nlp_qpscaling_memory_get(ocp_nlp_qpscaling_dims *dims, void *mem_, const char *field, int stage, void* value)
{
    ocp_nlp_qpscaling_memory *mem = mem_;

    if (!strcmp(field, "constr"))
    {
        double *ptr = value;
        blasfeo_unpack_dvec(dims->qp_dim->ng[stage], mem->constraints_scaling_vec + stage, 0, ptr, 1);
    }
    else if (!strcmp(field, "obj"))
    {
        double *ptr = value;
        *ptr = mem->obj_factor;
    }
    else if (!strcmp(field, "scaled_qp_in"))
    {
        ocp_qp_in **ptr = value;
        *ptr = mem->scaled_qp_in;
    }
    else if (!strcmp(field, "scaled_qp_out"))
    {
        ocp_qp_out **ptr = value;
        *ptr = mem->scaled_qp_out;
    }
    else if (!strcmp(field, "constraints_scaling_vec"))
    {
        struct blasfeo_dvec **ptr = value;
        *ptr = mem->constraints_scaling_vec + stage;
    }
    else if (!strcmp(field, "status"))
    {
        int *ptr = value;
        *ptr = mem->status;
    }
    else
    {
        printf("\nerror: ocp_nlp_qpscaling_memory_get: field %s not available\n", field);
        exit(1);
    }
}


/************************************************
 * helper functions
 ************************************************/
static double norm_inf_matrix_col(int col_idx, int col_length,  struct blasfeo_dmat *At)
{
    double norm = 0.0;
    for (int j = 0; j < col_length; ++j)
    {
        double tmp = BLASFEO_DMATEL(At, j, col_idx);
        norm = MAX(norm, fabs(tmp));
    }
    return norm;
}


static void rescale_solution_constraint_scaling(ocp_nlp_qpscaling_opts *opts, ocp_nlp_qpscaling_memory *mem, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    int *nb = qp_out->dim->nb;
    int *ng = qp_out->dim->ng;
    int *nx = qp_out->dim->nx;
    int *nu = qp_out->dim->nu;
    int *ns = qp_out->dim->ns;
    int N = qp_out->dim->N;
    double scaling_factor;
    int s_idx;

    for (int i = 0; i <= N; i++)
    {
        // copy ux;
        blasfeo_dveccp(nx[i]+nu[i]+2*ns[i], mem->scaled_qp_out->ux+i, 0, qp_out->ux+i, 0);

        if (opts->scale_qp_objective == NO_OBJECTIVE_SCALING)
        {
            // setup multipliers
            blasfeo_dveccp(2*nb[i]+2*ng[i]+2*ns[i], mem->scaled_qp_out->lam+i, 0, qp_out->lam+i, 0);
            if (i < N)
            {
                blasfeo_dveccp(nx[i+1], mem->scaled_qp_out->pi+i, 0, qp_out->pi+i, 0);
            }
        }

        // scale constraint multipliers
        for (int j = 0; j < ng[i]; ++j)
        {
            scaling_factor = BLASFEO_DVECEL(mem->constraints_scaling_vec+i, j);

            // scale lam of lower bound
            BLASFEO_DVECEL(qp_out->lam+i, nb[i]+j) *= scaling_factor;

            // scale lam of upper bound
            BLASFEO_DVECEL(qp_out->lam+i, 2*nb[i]+ng[i]+j) *= scaling_factor;

            s_idx = qp_in->idxs_rev[i][nb[i] + j];  // index of slack corresponding to this constraint
            if (s_idx >= 0)
            {
                // scale slack bound multipliers
                BLASFEO_DVECEL(qp_out->lam+i, 2*nb[i]+2*ng[i]+s_idx) *= scaling_factor;
                BLASFEO_DVECEL(qp_out->lam+i, 2*nb[i]+2*ng[i]+ns[i]+s_idx) *= scaling_factor;
                // scale slack variables
                BLASFEO_DVECEL(qp_out->ux+i, nx[i]+nu[i]+s_idx) /= scaling_factor;
                BLASFEO_DVECEL(qp_out->ux+i, nx[i]+nu[i]+ns[i]+s_idx) /= scaling_factor;
            }
        }
    }
}


static void ocp_qp_scale_objective(ocp_nlp_qpscaling_memory *mem, ocp_qp_in *qp_in, double factor)
{
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *ns = qp_in->dim->ns;

    for (int stage = 0; stage <= qp_in->dim->N; stage++)
    {
        // scale cost
        blasfeo_dveccpsc(nx[stage]+nu[stage]+2*ns[stage], factor, qp_in->rqz+stage, 0, mem->scaled_qp_in->rqz+stage, 0);
        blasfeo_dveccpsc(2*ns[stage], factor, qp_in->Z+stage, 0, mem->scaled_qp_in->Z+stage, 0);
        blasfeo_dgecpsc(nx[stage]+nu[stage], nx[stage]+nu[stage], factor, qp_in->RSQrq+stage, 0, 0, mem->scaled_qp_in->RSQrq+stage, 0, 0);
    }
}



static void ocp_qp_out_scale_duals(ocp_nlp_qpscaling_memory *mem, ocp_qp_out *qp_out, double factor)
{
    struct d_ocp_qp_dim *qp_dim = qp_out->dim;
    for (int i = 0; i <= qp_dim->N; i++)
    {
        blasfeo_dveccpsc(2*(qp_dim->nb[i]+qp_dim->ng[i]+qp_dim->ns[i]), factor, mem->scaled_qp_out->lam+i, 0, qp_out->lam+i, 0);
        if (i < qp_dim->N)
        {
            blasfeo_dveccpsc(qp_dim->nx[i+1], factor, mem->scaled_qp_out->pi+i, 0, qp_out->pi+i, 0);
        }
    }
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_qpscaling_precompute(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    ocp_nlp_qpscaling_opts *opts = opts_;
    ocp_nlp_qpscaling_memory *mem = mem_;

    // alias stuff that is the same between scaled and unscaled versions
    if (opts->scale_qp_constraints == NO_CONSTRAINT_SCALING && opts->scale_qp_objective == NO_OBJECTIVE_SCALING)
    {
        mem->scaled_qp_in = qp_in;
        mem->scaled_qp_out = qp_out;
    }
    else if (opts->scale_qp_constraints == NO_CONSTRAINT_SCALING)
    {
        /* qp_in */
        // only cost scaling -> alias everything not touched
        mem->scaled_qp_in->b = qp_in->b;
        mem->scaled_qp_in->BAbt = qp_in->BAbt;
        mem->scaled_qp_in->d = qp_in->d;
        mem->scaled_qp_in->d_mask = qp_in->d_mask;
        mem->scaled_qp_in->diag_H_flag = qp_in->diag_H_flag;
        mem->scaled_qp_in->DCt = qp_in->DCt;
        mem->scaled_qp_in->dim = qp_in->dim;
        mem->scaled_qp_in->idxb = qp_in->idxb;
        mem->scaled_qp_in->idxe = qp_in->idxe;
        mem->scaled_qp_in->idxs_rev = qp_in->idxs_rev;
        mem->scaled_qp_in->m = qp_in->m;
        // NOT aliased: rqz, RSQrq, Z
        /* qp_out */
        mem->scaled_qp_out->ux = qp_out->ux;
        mem->scaled_qp_out->misc = qp_out->misc;
        mem->scaled_qp_out->dim = qp_out->dim;
        mem->scaled_qp_out->t = qp_out->t;
        // NOT aliased: lam, pi
    }
    else
    {
        /* qp_in */
        // constraint scaling (& maybe cost scaling) -> alias everything not touched
        mem->scaled_qp_in->b = qp_in->b;
        mem->scaled_qp_in->BAbt = qp_in->BAbt;
        mem->scaled_qp_in->d_mask = qp_in->d_mask;
        mem->scaled_qp_in->diag_H_flag = qp_in->diag_H_flag;
        mem->scaled_qp_in->dim = qp_in->dim;
        mem->scaled_qp_in->idxb = qp_in->idxb;
        mem->scaled_qp_in->idxe = qp_in->idxe;
        mem->scaled_qp_in->idxs_rev = qp_in->idxs_rev;
        mem->scaled_qp_in->m = qp_in->m;
        // NOT aliased: rqz, Z, d, DCt
        if (opts->scale_qp_objective == NO_OBJECTIVE_SCALING)
        {
            // alias RSQrq, otherwise it is set up in objective scaling
            mem->scaled_qp_in->RSQrq = qp_in->RSQrq;
        }

        /* qp_out */
        mem->scaled_qp_out->misc = qp_out->misc;
        mem->scaled_qp_out->dim = qp_out->dim;
        mem->scaled_qp_out->t = qp_out->t;
        // NOT aliased: lam, pi, ux
    }
}

/**
 * @brief Scales the objective function of an OCP QP using Gershgorin eigenvalue estimates.
 *
 * - estimate max. abs. eigenvalue using gershgorin circles as max_abs_eig
 * - obj_factor = min(1.0, ub_max_abs_eig/max_abs_eig)
 */
void ocp_nlp_qpscaling_compute_obj_scaling_factor(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    double max_abs_eig = 0.0;
    double tmp, max_upscale_factor, lb_grad_norm_factor;

    ocp_qp_dims *dim = dims->qp_dim;
    ocp_nlp_qpscaling_memory *mem = mem_;
    ocp_nlp_qpscaling_opts *opts = opts_;

    int *nx = dim->nx;
    int *nu = dim->nu;
    int *ns = dim->ns;

    struct blasfeo_dmat *RSQrq = qp_in->RSQrq;
    double nrm_inf_grad_obj = 0.0;
    for (int stage = 0; stage <= dim->N; stage++)
    {
        compute_gershgorin_max_abs_eig_estimate(nx[stage]+nu[stage], RSQrq+stage, &tmp);
        max_abs_eig = MAX(max_abs_eig, tmp);
        // take Z into account
        blasfeo_dvecnrm_inf(2*ns[stage], qp_in->Z+stage, 0, &tmp);
        max_abs_eig = MAX(max_abs_eig, tmp);

        // norm gradient
        blasfeo_dvecnrm_inf(nx[stage]+nu[stage]+2*ns[stage], qp_in->rqz+stage, 0, &tmp);
        nrm_inf_grad_obj = MAX(nrm_inf_grad_obj, fabs(tmp));
    }

    if (max_abs_eig < opts->ub_max_abs_eig)
    {
        mem->obj_factor = 1.0;
        max_upscale_factor = opts->ub_max_abs_eig / max_abs_eig;
    }
    else
    {
        // scale objective down
        mem->obj_factor = opts->ub_max_abs_eig / max_abs_eig;
        max_upscale_factor = mem->obj_factor;
    }

    if (mem->obj_factor*nrm_inf_grad_obj <= opts->lb_norm_inf_grad_obj)
    {
        // grad norm would become too small -> scale cost up
        if (opts->print_level > 0)
        {
            printf("lb_norm_inf_grad_obj violated! %.2e\n", opts->lb_norm_inf_grad_obj);
            printf("Gradient is very small! %.2e\n", mem->obj_factor*nrm_inf_grad_obj);
        }
        lb_grad_norm_factor = opts->lb_norm_inf_grad_obj / nrm_inf_grad_obj;
        tmp = MIN(max_upscale_factor, lb_grad_norm_factor);
        mem->obj_factor = MAX(mem->obj_factor, tmp);
        mem->status = ACADOS_QPSCALING_BOUNDS_NOT_SATISFIED;
    }
    if (opts->print_level > 0)
    {
        printf("Scaling factor objective: %.2e\n", mem->obj_factor);
    }
}


// calculate scaling factors for all inequality constraints (except bounds) of the QP.
// The scaling factor is calculated as the maximum of the linear coefficients of the constraint and the maximum of the bounds.
void ocp_nlp_qpscaling_scale_constraints(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *ns = qp_in->dim->ns;
    int N = qp_in->dim->N;
    int i, j, s_idx;
    double mask_value_lower, mask_value_upper;
    ocp_nlp_qpscaling_memory *mem = mem_;
    ocp_nlp_qpscaling_opts *opts = opts_;
    double coeff_norm, scaling_factor;
    ocp_qp_in *scaled_qp_in = mem->scaled_qp_in;

    for (i = 0; i <= N; i++)
    {
        // setup all non-aliased stuff
        blasfeo_dveccp(2*(nb[i]+ng[i]+ns[i]), qp_in->d+i, 0, scaled_qp_in->d+i, 0);
        // copy cost, if not already done via cost scaling
        if (opts->scale_qp_objective == NO_OBJECTIVE_SCALING)
        {
            blasfeo_dveccp(nx[i]+nu[i]+2*ns[i], qp_in->rqz+i, 0, scaled_qp_in->rqz+i, 0);
            blasfeo_dveccp(2*ns[i], qp_in->Z+i, 0, scaled_qp_in->Z+i, 0);
        }
        // setup DCt, modify d in place
        for (j = 0; j < ng[i]; j++)
        {
            coeff_norm = norm_inf_matrix_col(j, nu[i]+nx[i], qp_in->DCt+i);
            mask_value_lower = BLASFEO_DVECEL(qp_in->d_mask+i, nb[i]+j);
            mask_value_upper = BLASFEO_DVECEL(qp_in->d_mask+i, 2*nb[i]+ng[i]+j);

            // calculate scaling factor from row norm
            double bound_max = MAX(fabs(mask_value_lower * BLASFEO_DVECEL(qp_in->d+i, nb[i]+j)),
                                    fabs(mask_value_upper * BLASFEO_DVECEL(qp_in->d+i, 2*nb[i]+ng[i]+j)));
            // only scale down.
            scaling_factor = 1.0 / MAX(1.0, MAX(bound_max, coeff_norm));

            // store scaling factor
            BLASFEO_DVECEL(mem->constraints_scaling_vec+i, j) = scaling_factor;

            // scale the constraint
            blasfeo_dgecpsc(nu[i]+nx[i], 1, scaling_factor, qp_in->DCt+i, 0, j, mem->scaled_qp_in->DCt+i, 0, j);

            s_idx = qp_in->idxs_rev[i][nb[i] + j];  // index of slack corresponding to this constraint
            if (s_idx != -1)
            {
                // printf("Scaling slack %d for constraint %d at stage %d with factor %.2e\n", s_idx, j, i, scaling_factor);
                // scale associated slack cost
                // lower
                BLASFEO_DVECEL(mem->scaled_qp_in->rqz+i, nu[i]+nx[i]+s_idx) /= scaling_factor;
                BLASFEO_DVECEL(mem->scaled_qp_in->Z+i, s_idx) /= (scaling_factor*scaling_factor);
                // upper
                BLASFEO_DVECEL(mem->scaled_qp_in->rqz+i, nu[i]+nx[i]+ns[i]+s_idx) /= scaling_factor;
                BLASFEO_DVECEL(mem->scaled_qp_in->Z+i, ns[i]+s_idx) /= (scaling_factor*scaling_factor);
                // scale slack bounds
                BLASFEO_DVECEL(mem->scaled_qp_in->d+i, 2*(nb[i]+ng[i])+s_idx) *= scaling_factor;
                BLASFEO_DVECEL(mem->scaled_qp_in->d+i, 2*(nb[i]+ng[i])+ns[i]+s_idx) *= scaling_factor;
            }

            // scale lower bound
            if (mask_value_lower == 1.0)
            {
                BLASFEO_DVECEL(mem->scaled_qp_in->d+i, nb[i]+j) *= scaling_factor;
            }

            // scale upper bound
            if (mask_value_upper == 1.0)
            {
                BLASFEO_DVECEL(mem->scaled_qp_in->d+i, 2*nb[i]+ng[i]+j) *= scaling_factor;
            }
        }
    }
}


static void print_qp_scaling_factors_constr(ocp_nlp_qpscaling_dims *dims, ocp_nlp_qpscaling_opts *opts, ocp_nlp_qpscaling_memory *mem)
{
    if (opts->scale_qp_constraints)
    {
        printf("Scaling factors for constraints:\n");
        for (int i = 0; i <= dims->qp_dim->N; i++)
        {
            printf("Stage %d: ", i);
            for (int j = 0; j < dims->qp_dim->ng[i]; j++)
            {
                printf("%.2e ", BLASFEO_DVECEL(mem->constraints_scaling_vec+i, j));
            }
            printf("\n");
        }
    }
}

void ocp_nlp_qpscaling_scale_qp(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    ocp_nlp_qpscaling_opts *opts = opts_;
    ocp_nlp_qpscaling_memory *mem = mem_;

    mem->status = ACADOS_SUCCESS;
    if (opts->scale_qp_objective)
    {
        ocp_nlp_qpscaling_compute_obj_scaling_factor(dims, opts_, mem_, qp_in);
        ocp_qp_scale_objective(mem, qp_in, mem->obj_factor);
    }
    else
    {
        // set the obj_factor to 1.0, for consinstency
        mem->obj_factor = 1.0;
    }

    if (opts->scale_qp_constraints)
    {
        ocp_nlp_qpscaling_scale_constraints(dims, opts_, mem_, qp_in);
    }
    if (opts->print_level > 0)
    {
        print_qp_scaling_factors_constr(dims, opts, mem);
    }
    // printf("qp_in AFTER SCALING\n");
    // print_ocp_qp_in(qp_in);
}


void ocp_nlp_qpscaling_rescale_solution(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    ocp_nlp_qpscaling_memory *mem = mem_;
    ocp_nlp_qpscaling_opts *opts = opts_;

    if (opts->scale_qp_objective)
    {
        ocp_qp_out_scale_duals(mem, qp_out, 1.0/mem->obj_factor);
    }
    if (opts->scale_qp_constraints)
    {
        rescale_solution_constraint_scaling(opts, mem, qp_in, qp_out);
        qp_info *info = (qp_info *) qp_out->misc;
        info->t_computed = 0;  // t needs to be recomputed if needed.
    }
    // printf("qp_out AFTER RESCALING\n");
    // print_ocp_qp_out(qp_out);
    return;
}
