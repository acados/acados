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

    size += 3*(N+1)*sizeof(int); // nx nu ng

    return size;
}



ocp_nlp_qpscaling_dims *ocp_nlp_qpscaling_dims_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // dims
    ocp_nlp_qpscaling_dims *dims = (ocp_nlp_qpscaling_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_qpscaling_dims);
    // nx
    dims->nx = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // nu
    dims->nu = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // ng
    dims->ng = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);

    dims->N = N;

    // initialize to zero by default
    int ii;
    // nx
    for(ii=0; ii<=N; ii++)
        dims->nx[ii] = 0;
    // nu
    for(ii=0; ii<=N; ii++)
        dims->nu[ii] = 0;

    assert((char *) raw_memory + ocp_nlp_qpscaling_dims_calculate_size(N) >= c_ptr);

    return dims;
}



void ocp_nlp_qpscaling_dims_set(ocp_nlp_qpscaling_dims *dims, int stage, char *field, int* value)
{
    if (!strcmp(field, "nx"))
    {
        dims->nx[stage] = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu[stage] = *value;
    }
    else if (!strcmp(field, "ng"))
    {
        dims->ng[stage] = *value;
    }
    else
    {
        printf("\nerror: field %s not available in module ocp_nlp_qpscaling_dims_set\n", field);
        exit(1);
    }

    return;
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

    opts->scale_qp_objective = NO_COST_SCALING;
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
        // printf("ub_max_eig_abs: %2.e\n", opts->ub_max_abs_eig);
    }
    else if (!strcmp(field, "lb_norm_inf_grad_obj"))
    {
        double *d_ptr = value;
        opts->lb_norm_inf_grad_obj = *d_ptr;
        // printf("lb_norm_inf_grad_obj : %2.e\n", opts->lb_norm_inf_grad_obj);
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

acados_size_t ocp_nlp_qpscaling_memory_calculate_size(ocp_nlp_qpscaling_dims *dims, void *opts_)
{

    int N = dims->N;
    int i;
    ocp_nlp_qpscaling_opts *opts = opts_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_qpscaling_memory);

    // constraints_scaling_vec
    if (opts->scale_qp_constraints)
    {
        size += (N + 1) * sizeof(struct blasfeo_dvec);  // constraints_scaling_vec
        for (i = 0; i <= N; i++)
        {
            size += blasfeo_memsize_dvec(dims->ng[i]);
        }
    }

    return size;
}



void *ocp_nlp_qpscaling_memory_assign(ocp_nlp_qpscaling_dims *dims, void *opts_, void *raw_memory)
{
    ocp_nlp_qpscaling_opts *opts = opts_;
    char *c_ptr = (char *) raw_memory;
    int N = dims->N;

    ocp_nlp_qpscaling_memory *mem = (ocp_nlp_qpscaling_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_qpscaling_memory);

    assert((char *)mem + ocp_nlp_qpscaling_memory_calculate_size(dims, opts_) >= c_ptr);

    if (opts->scale_qp_constraints)
    {
        assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->constraints_scaling_vec, &c_ptr);
        for (int i = 0; i <= N; ++i)
        {
            assign_and_advance_blasfeo_dvec_mem(dims->ng[i], mem->constraints_scaling_vec + i, &c_ptr);
            blasfeo_dvecse(dims->ng[i], 1.0, mem->constraints_scaling_vec+i, 0);
        }
    }

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
        blasfeo_unpack_dvec(dims->ng[stage], mem->constraints_scaling_vec + stage, 0, ptr, 1);
    }
    else if (!strcmp(field, "obj"))
    {
        double *ptr = value;
        *ptr = mem->obj_factor;
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
/*
The interesting matrices are stored transposed
*/
static double norm_inf_matrix_row(int row, int n_col,  struct blasfeo_dmat *At)
{
    double norm = 0.0;
    for (int j = 0; j < n_col; ++j)
    {
        double tmp = BLASFEO_DMATEL(At, j, row);
        norm = fmax(norm, fabs(tmp));
    }
    return norm;
}

/*
The interesting matrices are stored transposed
*/
// static double norm_2_matrix_row(int row, int n_col, struct blasfeo_dmat *At)
// {
//     double norm = 0.0;
//     for (int j = 0; j < n_col; ++j)
//     {
//         double tmp = BLASFEO_DMATEL(At, j, row);
//         norm += tmp * tmp;
//     }
//     return sqrt(norm);
// }

/*
The interesting matrices are stored transposed
*/
static void scale_matrix_row(int row, int n_col, struct blasfeo_dmat *At, double scaling_factor)
{
    for (int j = 0; j < n_col; ++j)
    {
        BLASFEO_DMATEL(At, j, row) = BLASFEO_DMATEL(At, j, row) * scaling_factor;
    }
}

static void rescale_solution_constraint_scaling(ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_nlp_qpscaling_memory *mem)
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


static void ocp_qp_scale_objective(ocp_qp_in *qp_in, double factor)
{
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *ns = qp_in->dim->ns;

    struct blasfeo_dvec *rqz = qp_in->rqz;
    struct blasfeo_dmat *RSQrq = qp_in->RSQrq;
    struct blasfeo_dvec *Z = qp_in->Z;

    for (int stage = 0; stage <= qp_in->dim->N; stage++)
    {
        // scale cost
        blasfeo_dvecsc(nx[stage]+nu[stage]+2*ns[stage], factor, rqz+stage, 0);
        blasfeo_dvecsc(2*ns[stage], factor, Z+stage, 0);
        blasfeo_dgecpsc(nx[stage]+nu[stage], nx[stage]+nu[stage], factor, RSQrq+stage, 0, 0, RSQrq+stage, 0, 0);
    }
}



static void ocp_qp_out_scale_duals(ocp_qp_out *qp_out, double factor)
{
    struct d_ocp_qp_dim *qp_dim = qp_out->dim;
    struct blasfeo_dvec *pi = qp_out->pi;
    struct blasfeo_dvec *lam = qp_out->lam;
    // print_ocp_qp_dims(qp_dim);
    for (int stage = 0; stage <= qp_dim->N; stage++)
    {
        blasfeo_dvecsc(2*qp_dim->nb[stage]+2*qp_dim->ng[stage]+2*qp_dim->ns[stage], factor, lam+stage, 0);
        if (stage < qp_dim->N)
            blasfeo_dvecsc(qp_dim->nx[stage+1], factor, pi+stage, 0);
    }
}






/************************************************
 * functions
 ************************************************/
/**
 * @brief Scales the objective function of an OCP QP using Gershgorin eigenvalue estimates.
 *
 * - estimate max. abs. eigenvalue using gershgorin circles as max_abs_eig
 * - obj_factor = min(1.0, ub_max_abs_eig/max_abs_eig)
 */
void ocp_nlp_qpscaling_scale_objective(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    double max_abs_eig = 0.0;
    double tmp;
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *ns = qp_in->dim->ns;
    ocp_nlp_qpscaling_memory *memory = mem_;
    ocp_nlp_qpscaling_opts *opts = opts_;

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
        blasfeo_dvecnrm_inf(nx[stage]+nu[stage]+2*ns[stage], qp_in->rqz+stage, 0, &tmp);
        nrm_inf_grad_obj = fmax(nrm_inf_grad_obj, fabs(tmp));
    }

    double max_upscale_factor, lb_grad_norm_factor;
    if (max_abs_eig < opts->ub_max_abs_eig)
    {
        memory->obj_factor = 1.0;
        max_upscale_factor = opts->ub_max_abs_eig / max_abs_eig;
    }
    else
    {
        // scale objective down
        memory->obj_factor = opts->ub_max_abs_eig / max_abs_eig;
        max_upscale_factor = memory->obj_factor;
    }
    // printf("Scaling factor objective: %.2e based on ub_max_abs_eig\n", memory->obj_factor);

    if (memory->obj_factor*nrm_inf_grad_obj <= opts->lb_norm_inf_grad_obj)
    {
        // grad norm would become too small -> upscale cost
        // printf("lb_norm_inf_grad_obj violated! %.2e\n", opts->lb_norm_inf_grad_obj);
        // printf("Gradient is very small! %.2e\n", memory->obj_factor*nrm_inf_grad_obj);
        lb_grad_norm_factor = opts->lb_norm_inf_grad_obj / nrm_inf_grad_obj;
        tmp = fmin(max_upscale_factor, lb_grad_norm_factor);
        memory->obj_factor = fmax(memory->obj_factor, tmp);
        // TODO: return some status code here?
    }
    // printf("Scaling factor objective: %.2e\n", memory->obj_factor);

    // scale QP cost
    // print_ocp_qp_in(qp_in);
    if (memory->obj_factor != 1.0)
    {
        ocp_qp_scale_objective(qp_in, memory->obj_factor);
    }

    // printf("AFTER SCALING\n");
    // print_ocp_qp_in(qp_in);
    // printf("ocp_nlp_qpscaling_scale_qp: computed obj_factor = %e\n", memory->obj_factor);
}


// calculate scaling factors for all inequality constraints (except bounds) of the QP.
// The scaling factor is calculated as the maximum of the row norm and the maximum of the bounds.
void ocp_nlp_qpscaling_scale_constraints(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *nb = qp_in->dim->nb;
    int *ng = qp_in->dim->ng;
    int *ns = qp_in->dim->ns;
    int N = dims->N;
    int i, j, s_idx;
    double mask_value_lower, mask_value_upper;
    ocp_nlp_qpscaling_memory *memory = mem_;
    // ocp_nlp_qpscaling_opts *opts = opts_; // option for what norm to be used
    double row_norm, scaling_factor;

    for (i = 0; i <= N; i++)
    {
        for (j = 0; j < ng[i]; j++)
        {
            row_norm = norm_inf_matrix_row(j, nu[i]+nx[i],  &qp_in->DCt[i]);
            mask_value_lower = BLASFEO_DVECEL(qp_in->d_mask+i, nb[i]+j);
            mask_value_upper = BLASFEO_DVECEL(qp_in->d_mask+i, 2*nb[i]+ng[i]+j);

            // calculate scaling factor from row norm
            double bound_max = fmax(fabs(mask_value_lower * BLASFEO_DVECEL(qp_in->d+i, nb[i]+j)),
                                    fabs(mask_value_upper * BLASFEO_DVECEL(qp_in->d+i, 2*nb[i]+ng[i]+j)));
            scaling_factor = 1 / fmax(1.0, fmax(bound_max, row_norm));

            // store scaling factor in memory
            BLASFEO_DVECEL(memory->constraints_scaling_vec+i, j) = scaling_factor;

            // scale the row
            scale_matrix_row(j, nu[i]+nx[i],  &qp_in->DCt[i], scaling_factor);

            s_idx = qp_in->idxs_rev[i][nb[i] + j];  // index of slack corresponding to this constraint
            if (s_idx != -1)
            {
                // printf("Scaling slack %d for constraint %d at stage %d with factor %.2e\n", s_idx, j, i, scaling_factor);
                // scale associated slack cost
                // lower
                BLASFEO_DVECEL(qp_in->rqz+i, nu[i]+nx[i]+s_idx) *= 1/scaling_factor;
                BLASFEO_DVECEL(qp_in->Z+i, s_idx) *= 1/(scaling_factor*scaling_factor);
                // upper
                BLASFEO_DVECEL(qp_in->rqz+i, nu[i]+nx[i]+ns[i]+s_idx) *= 1/scaling_factor;
                BLASFEO_DVECEL(qp_in->Z+i, ns[i]+s_idx) *= 1/(scaling_factor*scaling_factor);
            }

            // scale lower bound
            if (mask_value_lower == 1.0)
            {
                BLASFEO_DVECEL(qp_in->d+i, nb[i]+j) = BLASFEO_DVECEL(qp_in->d+i, nb[i]+j) * scaling_factor;
            }

            // scale upper bound
            if (mask_value_upper == 1.0)
            {
                BLASFEO_DVECEL(qp_in->d+i, 2*nb[i]+ng[i]+j) = BLASFEO_DVECEL(qp_in->d+i, 2*nb[i]+ng[i]+j) * scaling_factor;
            }
        }
    }
}

void ocp_nlp_qpscaling_scale_qp(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    ocp_nlp_qpscaling_opts *opts = opts_;

    // printf("qp_in BEFORE SCALING\n");
    // print_ocp_qp_in(qp_in);
    if (opts->scale_qp_objective)
    {
        ocp_nlp_qpscaling_scale_objective(dims, opts_, mem_, qp_in);
    }
    else
    {
        // set the obj_factor to 1.0, for consinstency
        ocp_nlp_qpscaling_memory *memory = mem_;
        memory->obj_factor = 1.0;
    }
    if (opts->scale_qp_constraints)
    {
        ocp_nlp_qpscaling_scale_constraints(dims, opts_, mem_, qp_in);
    }
//     printf("qp_in AFTER SCALING\n");
//     print_ocp_qp_in(qp_in);
}


void ocp_nlp_qpscaling_rescale_solution(ocp_nlp_qpscaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    ocp_nlp_qpscaling_memory *memory = mem_;
    ocp_nlp_qpscaling_opts *opts = opts_;

    if (memory->obj_factor != 1.0 && opts->scale_qp_objective)
    {
        ocp_qp_out_scale_duals(qp_out, 1/memory->obj_factor);
        ocp_qp_scale_objective(qp_in, 1/memory->obj_factor);
    }
    if (opts->scale_qp_constraints)
    {
        rescale_solution_constraint_scaling(qp_in, qp_out, memory);
    }

    return;
}
