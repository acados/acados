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


#include "acados_c/ocp_nlp_interface.h"

// external
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_external.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_conl.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_disc.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgp.h"
#include "acados/ocp_nlp/ocp_nlp_reg_convexify.h"
#include "acados/ocp_nlp/ocp_nlp_reg_glm.h"
#include "acados/ocp_nlp/ocp_nlp_reg_mirror.h"
#include "acados/ocp_nlp/ocp_nlp_reg_project.h"
#include "acados/ocp_nlp/ocp_nlp_reg_project_reduc_hess.h"
#include "acados/ocp_nlp/ocp_nlp_reg_noreg.h"
#include "acados/ocp_nlp/ocp_nlp_qpscaling.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_fixed_step.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_merit_backtracking.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_funnel.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_nlp/ocp_nlp_sqp_with_feasible_qp.h"
#include "acados/ocp_nlp/ocp_nlp_sqp_rti.h"
#include "acados/ocp_nlp/ocp_nlp_ddp.h"
#include "acados/utils/mem.h"
#include "acados/utils/strsep.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_blas.h"


/************************************************
* plan
************************************************/

static acados_size_t ocp_nlp_plan_calculate_size(int N)
{
    // N - number of shooting nodes
    acados_size_t bytes = sizeof(ocp_nlp_plan_t);
    bytes += N * sizeof(sim_solver_plan_t);
    bytes += (N + 1) * sizeof(ocp_nlp_cost_t);
    bytes += N * sizeof(ocp_nlp_dynamics_t);
    bytes += (N+1) * sizeof(ocp_nlp_constraints_t);
    return bytes;
}



static ocp_nlp_plan_t *ocp_nlp_plan_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_plan_t *plan = (ocp_nlp_plan_t *) c_ptr;
    c_ptr += sizeof(ocp_nlp_plan_t);

    plan->sim_solver_plan = (sim_solver_plan_t *) c_ptr;
    c_ptr += N * sizeof(sim_solver_plan_t);

    plan->nlp_cost = (ocp_nlp_cost_t *) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_cost_t);

    plan->nlp_dynamics = (ocp_nlp_dynamics_t *) c_ptr;
    c_ptr += N * sizeof(ocp_nlp_dynamics_t);

    plan->nlp_constraints = (ocp_nlp_constraints_t *) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_constraints_t);

    plan->N = N;

    return plan;
}



static void ocp_nlp_plan_initialize_default(ocp_nlp_plan_t *plan)
{
    int ii;

    int N = plan->N;

    // initialize to default value !=0 to detect empty plans
    // cost
    for (ii=0; ii <= N; ii++)
    {
        plan->nlp_cost[ii] = INVALID_COST;
    }
    // dynamics
    for (ii=0; ii < N; ii++)
    {
        plan->nlp_dynamics[ii] = INVALID_DYNAMICS;
    }
    // constraints
    for (ii=0; ii <= N; ii++)
    {
        plan->nlp_constraints[ii] = INVALID_CONSTRAINT;
    }
    // nlp solver
    plan->nlp_solver = INVALID_NLP_SOLVER;
    // qp solver
    plan->ocp_qp_solver_plan.qp_solver = INVALID_QP_SOLVER;
    // sim solver
    for (ii = 0; ii < N; ii++)
    {
        plan->sim_solver_plan[ii].sim_solver = INVALID_SIM_SOLVER;
    }

    // regularization: no reg by default
    plan->regularization = NO_REGULARIZE;

    // globalization: fixed step by default
    plan->globalization = FIXED_STEP;

    return;
}



ocp_nlp_plan_t *ocp_nlp_plan_create(int N)
{
    acados_size_t bytes = ocp_nlp_plan_calculate_size(N);
    void *ptr = acados_malloc(bytes, 1);
    assert(ptr != 0);

    ocp_nlp_plan_t *plan = ocp_nlp_plan_assign(N, ptr);

    ocp_nlp_plan_initialize_default(plan);

    return plan;
}



void ocp_nlp_plan_destroy(void* plan_)
{
    free(plan_);
}



/************************************************
* config
************************************************/

ocp_nlp_config *ocp_nlp_config_create(ocp_nlp_plan_t plan)
{
    int N = plan.N;

    /* calculate_size & malloc & assign */

    acados_size_t bytes = ocp_nlp_config_calculate_size(N);
    void *config_mem = acados_calloc(1, bytes);
    assert(config_mem != 0);
    ocp_nlp_config *config = ocp_nlp_config_assign(N, config_mem);

    /* initialize config according plan */

    // NLP solver
    switch (plan.nlp_solver)
    {
        case SQP:
            ocp_nlp_sqp_config_initialize_default(config);
            break;
        case SQP_WITH_FEASIBLE_QP:
            ocp_nlp_sqp_wfqp_config_initialize_default(config);
            break;
        case SQP_RTI:
            ocp_nlp_sqp_rti_config_initialize_default(config);
            break;
        case DDP:
            ocp_nlp_ddp_config_initialize_default(config);
            break;
        case INVALID_NLP_SOLVER:
            printf("\nerror: ocp_nlp_config_create: forgot to initialize plan->nlp_solver\n");
            exit(1);
        default:
            printf("\nerror: ocp_nlp_config_create: unsupported plan->nlp_solver\n");
            exit(1);
    }

    // QP solver
    ocp_qp_xcond_solver_config_initialize_from_plan(plan.ocp_qp_solver_plan.qp_solver,
                                                    config->qp_solver);

    // relaxed QP solver
    ocp_qp_xcond_solver_config_initialize_from_plan(plan.ocp_qp_solver_plan.qp_solver,
                                                    config->relaxed_qp_solver);

    // regularization
    switch (plan.regularization)
    {
        case NO_REGULARIZE:
            ocp_nlp_reg_noreg_config_initialize_default(config->regularize);
            break;
        case MIRROR:
            ocp_nlp_reg_mirror_config_initialize_default(config->regularize);
            break;
        case PROJECT:
            ocp_nlp_reg_project_config_initialize_default(config->regularize);
            break;
        case PROJECT_REDUC_HESS:
            ocp_nlp_reg_project_reduc_hess_config_initialize_default(config->regularize);
            break;
        case CONVEXIFY:
            ocp_nlp_reg_convexify_config_initialize_default(config->regularize);
            break;
        case GERSHGORIN_LEVENBERG_MARQUARDT:
            ocp_nlp_reg_glm_config_initialize_default(config->regularize);
            break;
        default:
            printf("\nerror: ocp_nlp_config_create: unsupported plan->regularization\n");
            exit(1);
    }

    // globalization
    switch (plan.globalization)
    {
        case FIXED_STEP:
            ocp_nlp_globalization_fixed_step_config_initialize_default(config->globalization);
            break;
        case MERIT_BACKTRACKING:
            ocp_nlp_globalization_merit_backtracking_config_initialize_default(config->globalization);
            break;
        case FUNNEL_L1PEN_LINESEARCH:
            ocp_nlp_globalization_funnel_config_initialize_default(config->globalization);
            break;
        default:
            printf("\nerror: ocp_nlp_config_create: unsupported plan->globalization\n");
            exit(1);
    }

    // globalization
    // NLP solver
    if (plan.nlp_solver == DDP && plan.globalization == MERIT_BACKTRACKING)
    {
        config->globalization->find_acceptable_iterate = &ocp_nlp_globalization_merit_backtracking_find_acceptable_iterate_for_ddp;
        config->globalization->needs_qp_objective_value = &ocp_nlp_globalization_merit_backtracking_ddp_needs_qp_objective_value;
    }

    // cost
    for (int i = 0; i <= N; ++i)
    {
        switch (plan.nlp_cost[i])
        {
            case LINEAR_LS:
                ocp_nlp_cost_ls_config_initialize_default(config->cost[i], i);
                break;
            case NONLINEAR_LS:
                ocp_nlp_cost_nls_config_initialize_default(config->cost[i], i);
                break;
            case CONVEX_OVER_NONLINEAR:
                ocp_nlp_cost_conl_config_initialize_default(config->cost[i], i);
                break;
            case EXTERNAL:
                ocp_nlp_cost_external_config_initialize_default(config->cost[i], i);
                break;
            case INVALID_COST:
                printf("\nerror: ocp_nlp_config_create: forgot to initialize plan->nlp_cost\n");
                exit(1);
            default:
                printf("\nerror: ocp_nlp_config_create: unsupported plan->nlp_cost\n");
                exit(1);
        }
    }

    // Dynamics
    for (int i = 0; i < N; ++i)
    {
        switch (plan.nlp_dynamics[i])
        {
            case CONTINUOUS_MODEL:
                ocp_nlp_dynamics_cont_config_initialize_default(config->dynamics[i], i);
                //                config->dynamics[i]->sim_solver =
                //                sim_config_create(plan.sim_solver[i]);
                sim_solver_t solver_name = plan.sim_solver_plan[i].sim_solver;

                switch (solver_name)
                {
                    case ERK:
                        sim_erk_config_initialize_default(config->dynamics[i]->sim_solver);
                        break;
                    case IRK:
                        sim_irk_config_initialize_default(config->dynamics[i]->sim_solver);
                        break;
                    case GNSF:
                        sim_gnsf_config_initialize_default(config->dynamics[i]->sim_solver);
                        break;
                    case LIFTED_IRK:
                        sim_lifted_irk_config_initialize_default(config->dynamics[i]->sim_solver);
                        break;
                    default:
                        printf("\nerror: ocp_nlp_config_create: unsupported plan->sim_solver\n");
                        exit(1);
                }

                break;
            case DISCRETE_MODEL:
                ocp_nlp_dynamics_disc_config_initialize_default(config->dynamics[i], i);
                break;
            case INVALID_DYNAMICS:
                printf("\nerror: ocp_nlp_config_create: forgot to initialize plan->nlp_dynamics\n");
                exit(1);
            default:
                printf("\nerror: ocp_nlp_config_create: unsupported plan->nlp_dynamics\n");
                exit(1);
        }
    }

    // Constraints
    for (int i = 0; i <= N; ++i)
    {
        switch (plan.nlp_constraints[i])
        {
            case BGH:
                ocp_nlp_constraints_bgh_config_initialize_default(config->constraints[i], i);
                break;
            case BGP:
                ocp_nlp_constraints_bgp_config_initialize_default(config->constraints[i], i);
                break;
            case INVALID_CONSTRAINT:
                printf(
                    "\nerror: ocp_nlp_config_create: forgot to initialize plan->nlp_constraints\n");
                exit(1);
            default:
                printf("\nerror: ocp_nlp_config_create: unsupported plan->nlp_constraints\n");
                exit(1);
        }
    }

    return config;
}



void ocp_nlp_config_destroy(void *config_)
{
    free(config_);
}


/************************************************
* dims
************************************************/

ocp_nlp_dims *ocp_nlp_dims_create(void *config_)
{
    ocp_nlp_config *config = config_;

    acados_size_t bytes = ocp_nlp_dims_calculate_size(config);

    void *ptr = acados_calloc(1, bytes);
    assert(ptr != 0);

    ocp_nlp_dims *dims = ocp_nlp_dims_assign(config, ptr);
    dims->raw_memory = ptr;

    return dims;
}



void ocp_nlp_dims_destroy(void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    free(dims->raw_memory);
}



/************************************************
* NLP inputs
************************************************/

ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    acados_size_t bytes = ocp_nlp_in_calculate_size(config, dims);

    void *ptr = acados_calloc(1, bytes);
    assert(ptr != 0);

    ocp_nlp_in *nlp_in = ocp_nlp_in_assign(config, dims, ptr);
    nlp_in->raw_memory = ptr;

    return nlp_in;
}



void ocp_nlp_in_destroy(void *in_)
{
    ocp_nlp_in *in = in_;
    free(in->raw_memory);
}



void ocp_nlp_in_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage,
        const char *field, void *value)
{
    if (!strcmp(field, "Ts"))
    {
        double *Ts_value = value;
        in->Ts[stage] = Ts_value[0];
    }
    else if (!strcmp(field, "parameter_values"))
    {
        double *parameter_values = value;
        for (int ii = 0; ii < dims->np[stage]; ii++)
        {
            in->parameter_values[stage][ii] = parameter_values[ii];
        }
    }
    else
    {
        printf("\nerror: ocp_nlp_in_set: field %s not available\n", field);
        exit(1);
    }
    return;
}


void ocp_nlp_in_set_params_sparse(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage,
        int *idx, double *p, int n_update)
{
    for (int ii = 0; ii < n_update; ii++)
    {
        in->parameter_values[stage][idx[ii]] = p[ii];
    }

    return;
}




void ocp_nlp_in_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage,
        const char *field, void *value)
{
    if (!strcmp(field, "Ts"))
    {
        double *Ts_value = value;
        Ts_value[0] = in->Ts[stage];
    }
    else if (!strcmp(field, "parameter_pointer"))
    {
        double **ptr = value;
        ptr[0] = in->parameter_values[stage];
    }
    else if (!strcmp(field, "p"))
    {
        double *out = (double *) value;
        for (int ii = 0; ii < dims->np[stage]; ii++)
        {
            out[ii] = in->parameter_values[stage][ii];
        }
    }
    else
    {
        printf("\nerror: ocp_nlp_in_get: field %s not available\n", field);
        exit(1);
    }
    return;
}


int ocp_nlp_dynamics_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
        int stage, const char *field, void *value)
{
    ocp_nlp_dynamics_config *dynamics_config = config->dynamics[stage];

    dynamics_config->model_set(dynamics_config, dims->dynamics[stage], in->dynamics[stage], field, value);

    return ACADOS_SUCCESS;
}



int ocp_nlp_cost_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_in *in, int stage, const char *field, void *value)
{
    ocp_nlp_cost_config *cost_config = config->cost[stage];
    return cost_config->model_set(cost_config, dims->cost[stage], in->cost[stage], field, value);
}

int ocp_nlp_cost_model_get(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_in *in, int stage, const char *field, void *value)
{
    ocp_nlp_cost_config *cost_config = config->cost[stage];
    return cost_config->model_get(cost_config, dims->cost[stage], in->cost[stage], field, value);
}

int ocp_nlp_constraints_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_in *in, ocp_nlp_out *out, int stage, const char *field, void *value)
{
    ocp_nlp_constraints_config *constr_config = config->constraints[stage];

    // this updates both the bounds and the mask
    int status = constr_config->model_set(constr_config, dims->constraints[stage],
            in->constraints[stage], field, value);
    // multiply lam with new mask to ensure that multipliers associated with masked constraints are zero.
    blasfeo_dvecmul(2*dims->ni[stage], &in->dmask[stage], 0, &out->lam[stage], 0, &out->lam[stage], 0);

    return status;
}


int ocp_nlp_dynamics_model_set_external_param_fun(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
        int stage, const char *field, void *ext_fun_)
{
    ocp_nlp_dynamics_config *dynamics_config = config->dynamics[stage];
    external_function_external_param_generic * ext_fun = (external_function_external_param_generic *) ext_fun_;

    ext_fun->set_param_pointer(ext_fun, in->parameter_values[stage]);

    if (dims->n_global_data > 0)
        ext_fun->set_global_data_pointer(ext_fun, in->global_data);

    dynamics_config->model_set(dynamics_config, dims->dynamics[stage], in->dynamics[stage], field, ext_fun);

    return ACADOS_SUCCESS;
}



int ocp_nlp_cost_model_set_external_param_fun(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_in *in, int stage, const char *field, void *ext_fun_)
{
    ocp_nlp_cost_config *cost_config = config->cost[stage];
    external_function_external_param_generic * ext_fun = (external_function_external_param_generic *) ext_fun_;

    ext_fun->set_param_pointer(ext_fun, in->parameter_values[stage]);

    if (dims->n_global_data > 0)
        ext_fun->set_global_data_pointer(ext_fun, in->global_data);

    return cost_config->model_set(cost_config, dims->cost[stage], in->cost[stage], field, ext_fun);

}



int ocp_nlp_constraints_model_set_external_param_fun(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_in *in, int stage, const char *field, void *ext_fun_)
{
    ocp_nlp_constraints_config *constr_config = config->constraints[stage];
    external_function_external_param_generic * ext_fun = (external_function_external_param_generic *) ext_fun_;

    ext_fun->set_param_pointer(ext_fun, in->parameter_values[stage]);

    if (dims->n_global_data > 0)
        ext_fun->set_global_data_pointer(ext_fun, in->global_data);

    return constr_config->model_set(constr_config, dims->constraints[stage],
            in->constraints[stage], field, ext_fun);
}



void ocp_nlp_constraints_model_get(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_in *in, int stage, const char *field, void *value)
{
    ocp_nlp_constraints_config *constr_config = config->constraints[stage];

    constr_config->model_get(constr_config, dims->constraints[stage],
            in->constraints[stage], field, value);
    return;
}




/************************************************
* out
************************************************/

ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    acados_size_t bytes = ocp_nlp_out_calculate_size(config, dims);

    void *ptr = acados_calloc(1, bytes);
    assert(ptr != 0);

    ocp_nlp_out *nlp_out = ocp_nlp_out_assign(config, dims, ptr);
    nlp_out->raw_memory = ptr;

    return nlp_out;
}



void ocp_nlp_out_destroy(void *out_)
{
    ocp_nlp_out *out = out_;
    free(out->raw_memory);
}



void ocp_nlp_out_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out, ocp_nlp_in *in,
        int stage, const char *field, void *value)
{
    double *double_values = value;
    if (!strcmp(field, "x"))
    {
        blasfeo_pack_dvec(dims->nx[stage], double_values, 1, &out->ux[stage], dims->nu[stage]);
    }
    else if (!strcmp(field, "u"))
    {
        blasfeo_pack_dvec(dims->nu[stage], double_values, 1, &out->ux[stage], 0);
    }
    else if (!strcmp(field, "sl"))
    {
        blasfeo_pack_dvec(dims->ns[stage], double_values, 1, &out->ux[stage],
                            dims->nu[stage] + dims->nx[stage]);
    }
    else if (!strcmp(field, "su"))
    {
        blasfeo_pack_dvec(dims->ns[stage], double_values, 1, &out->ux[stage],
                            dims->nu[stage] + dims->nx[stage] + dims->ns[stage]);
    }
    else if (!strcmp(field, "pi"))
    {
        blasfeo_pack_dvec(dims->nx[stage+1], double_values, 1, &out->pi[stage], 0);
    }
    else if (!strcmp(field, "lam"))
    {
        blasfeo_pack_dvec(2*dims->ni[stage], double_values, 1, &out->lam[stage], 0);
        // multiply with mask to ensure that multiplier associated with masked constraints are zero
        blasfeo_dvecmul(2*dims->ni[stage], &in->dmask[stage], 0, &out->lam[stage], 0, &out->lam[stage], 0);
    }
    else if (!strcmp(field, "z"))
    {
        blasfeo_pack_dvec(dims->nz[stage], double_values, 1, &out->z[stage], 0);
    }
    else
    {
        printf("\nerror: ocp_nlp_out_set: field %s not available\n", field);
        exit(1);
    }
}


void ocp_nlp_out_set_values_to_zero(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out)
{
    int N = dims->N;
    for (int i = 0; i<=N; i++)
    {
        blasfeo_dvecse(dims->nv[i], 0.0, &out->ux[i], 0);
        blasfeo_dvecse(dims->nz[i], 0.0, &out->z[i], 0);
        blasfeo_dvecse(dims->nx[i+1], 0.0, &out->pi[i], 0);
        blasfeo_dvecse(2*dims->ni[i], 0.0, &out->lam[i], 0);
    }
}


void ocp_nlp_out_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, void *value)
{
    if (!strcmp(field, "x"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->nx[stage], &out->ux[stage], dims->nu[stage], double_values, 1);
    }
    else if (!strcmp(field, "u"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->nu[stage], &out->ux[stage], 0, double_values, 1);
    }
    else if (!strcmp(field, "sl"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->ns[stage], &out->ux[stage],
             dims->nu[stage] + dims->nx[stage], double_values, 1);
    }
    else if (!strcmp(field, "su"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->ns[stage], &out->ux[stage],
             dims->nu[stage] + dims->nx[stage] + dims->ns[stage], double_values, 1);
    }
    else if (!strcmp(field, "z"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->nz[stage], &out->z[stage], 0, double_values, 1);
    }
    else if (!strcmp(field, "pi"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->nx[stage+1], &out->pi[stage], 0, double_values, 1);
    }
    else if (!strcmp(field, "lam"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(2*dims->ni[stage], &out->lam[stage], 0, double_values, 1);
    }
    else if ((!strcmp(field, "kkt_norm_inf")) || (!strcmp(field, "kkt_norm")))
    {
        double *double_values = value;
        double_values[0] = out->inf_norm_res;
    }
    else
    {
        printf("\nerror: ocp_nlp_out_get: field %s not available\n", field);
        exit(1);
    }
}


int ocp_nlp_dims_get_total_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out, const char *field)
{
    if (!strcmp(field, "x"))
    {
        return dims->nx_total;
    }
    else if (!strcmp(field, "u"))
    {
        return dims->nu_total;
    }
    else if (!strcmp(field, "sl") || !strcmp(field, "su") ||
             !strcmp(field, "zl") || !strcmp(field, "zu") ||
             !strcmp(field, "Zl") || !strcmp(field, "Zu") ||
             !strcmp(field, "cost_z") || !strcmp(field, "cost_Z"))
    {
        return dims->ns_total;
    }
    else if (!strcmp(field, "s"))
    {
        return 2*dims->ns_total;
    }
    else if (!strcmp(field, "z"))
    {
        return dims->nz_total;
    }
    else if (!strcmp(field, "pi"))
    {
        return dims->nx_total - dims->nx[0];
    }
    else if (!strcmp(field, "lam"))
    {
        return 2*dims->ni_total;
    }
    else if (!strcmp(field, "p"))
    {
        return dims->np_total;
    }
    else if (!strcmp(field, "lbx") || !strcmp(field, "ubx") || !strcmp(field, "nbx"))
    {
        return dims->nbx_total;
    }
    else if (!strcmp(field, "lbu") || !strcmp(field, "ubu") || !strcmp(field, "nbu"))
    {
        return dims->nbu_total;
    }
    else if (!strcmp(field, "lg") || !strcmp(field, "ug") || !strcmp(field, "ng"))
    {
        return dims->ng_total;
    }
    else if (!strcmp(field, "lh") || !strcmp(field, "uh") || !strcmp(field, "nh"))
    {
        return dims->nh_total;
    }
    else
    {
        printf("\nerror: ocp_nlp_dims_get_total_from_attr: field %s not available\n", field);
        exit(1);
    }
}



int ocp_nlp_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field)
{
    int dims_value = -1;

    // ocp_nlp_dims
    if (!strcmp(field, "x") || !strcmp(field, "nx"))
    {
        return dims->nx[stage];
    }
    else if (!strcmp(field, "pi"))
    {
        return dims->nx[stage+1];
    }
    else if (!strcmp(field, "nv"))
    {
        return dims->nv[stage];
    }
    else if (!strcmp(field, "u") || !strcmp(field, "nu"))
    {
        return dims->nu[stage];
    }
    else if (!strcmp(field, "z") || !strcmp(field, "nz"))
    {
        return dims->nz[stage];
    }
    else if (!strcmp(field, "p") || !strcmp(field, "np"))
    {
        return dims->np[stage];
    }
    else if (!strcmp(field, "p_global") || !strcmp(field, "np_global"))
    {
        return dims->np_global;
    }
    else if (!strcmp(field, "lam"))
    {
        return 2*dims->ni[stage];
    }
    else if (!strcmp(field, "s"))
    {
        return 2*dims->ns[stage];
    }
    else if (!strcmp(field, "sl") || !strcmp(field, "su") ||
             !strcmp(field, "zl") || !strcmp(field, "zu") ||
             !strcmp(field, "Zl") || !strcmp(field, "Zu") ||
             !strcmp(field, "cost_z") || !strcmp(field, "cost_Z"))
    {
        return dims->ns[stage];
    }
    // ocp nlp dynamics
    else if (!strcmp(field, "init_gnsf_phi"))
    {
        config->dynamics[stage]->dims_get(config->dynamics[stage], dims->dynamics[stage],
                                                    "gnsf_nout", &dims_value);
        return dims_value;
    }
    else if (!strcmp(field, "xdot_guess"))
    {
        config->dynamics[stage]->dims_get(config->dynamics[stage], dims->dynamics[stage],
                                                    "nx", &dims_value);
        return dims_value;
    }
    else if (!strcmp(field, "z_guess"))
    {
        config->dynamics[stage]->dims_get(config->dynamics[stage], dims->dynamics[stage],
                                                    "nz", &dims_value);
        return dims_value;
    }
    // ocp_nlp_constraints_dims
    else if (!strcmp(field, "lbx") || !strcmp(field, "ubx") || !strcmp(field, "nbx"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "nbx", &dims_value);
        return dims_value;
    }
    else if (!strcmp(field, "lbu") || !strcmp(field, "ubu") || !strcmp(field, "nbu"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "nbu", &dims_value);
        return dims_value;
    }
    else if (!strcmp(field, "lg") || !strcmp(field, "ug") || !strcmp(field, "ng"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "ng", &dims_value);
        return dims_value;
    }
    else if (!strcmp(field, "lh") || !strcmp(field, "uh") || !strcmp(field, "nh"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "nh", &dims_value);
        return dims_value;
    }
    else if (!strcmp(field, "qpscaling_constr"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "ng", &dims_value);
        dims_value += dims->ni_nl[stage];
        return dims_value;
    }
    // ocp_nlp_cost_dims
    else if (!strcmp(field, "y_ref") || !strcmp(field, "yref") || !strcmp(field, "ny"))
    {
        config->cost[stage]->dims_get(config->cost[stage], dims->cost[stage],
                                            "ny", &dims_value);
        return dims_value;
    }
    else
    {
        printf("\nerror: ocp_nlp_dims_get_from_attr: field %s not available\n", field);
        exit(1);
    }
}


void ocp_nlp_constraint_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, int *dims_out)
{
    // vectors first
    dims_out[1] = 0;
    // ocp_nlp_constraints_dims
    if (!strcmp(field, "lbx") || !strcmp(field, "ubx") || !strcmp(field, "nbx"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "nbx", &dims_out[0]);
        return;
    }
    else if (!strcmp(field, "uphi") || !strcmp(field, "nphi"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "nphi", &dims_out[0]);
        return;
    }
    else if (!strcmp(field, "lbu") || !strcmp(field, "ubu") || !strcmp(field, "nbu"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "nbu", &dims_out[0]);
        return;
    }
    else if (!strcmp(field, "lg") || !strcmp(field, "ug") || !strcmp(field, "ng"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "ng", &dims_out[0]);
        return;
    }
    else if (!strcmp(field, "lh") || !strcmp(field, "uh") || !strcmp(field, "nh"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "nh", &dims_out[0]);
        return;
    }
    else if (!strcmp(field, "ineq_fun"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                                            "ni", &dims_out[0]);
        dims_out[0] *= 2;
    }
    // matrices
    else if (!strcmp(field, "C"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                "ng", &dims_out[0]);

        dims_out[1] = dims->nx[stage];

        return;
    }
    else if (!strcmp(field, "D"))
    {
        config->constraints[stage]->dims_get(config->constraints[stage], dims->constraints[stage],
                "ng", &dims_out[0]);

        dims_out[1] = dims->nu[stage];

        return;
    }
    else
    {
        printf("\nerror: ocp_nlp_constraint_dims_get_from_attr: field %s not available\n", field);
        exit(1);
    }
}


// TODO: ocp_nlp_out not needed?
void ocp_nlp_qp_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, int *dims_out)
{
    // only matrices here matrices
    // dynamics
    if (!strcmp(field, "A") || !strcmp(field, "relaxed_A"))
    {
        dims_out[0] = dims->nx[stage+1];
        dims_out[1] = dims->nx[stage];
    }
    else if (!strcmp(field, "B") || !strcmp(field, "relaxed_B"))
    {
        dims_out[0] = dims->nx[stage+1];
        dims_out[1] = dims->nu[stage];
    }
    else if (!strcmp(field, "b") || !strcmp(field, "relaxed_b"))
    {
        dims_out[0] = 1;
        dims_out[1] = dims->nx[stage+1];
    }
    // cost
    else if (!strcmp(field, "Q") || !strcmp(field, "relaxed_Q") || !strcmp(field, "P") || !strcmp(field, "relaxed_P"))
    {
        dims_out[0] = dims->nx[stage];
        dims_out[1] = dims->nx[stage];
    }
    else if (!strcmp(field, "R") || !strcmp(field, "relaxed_R") || !strcmp(field, "Lr") || !strcmp(field, "relaxed_Lr"))
    {
        dims_out[0] = dims->nu[stage];
        dims_out[1] = dims->nu[stage];
    }
    else if (!strcmp(field, "S") || !strcmp(field, "relaxed_S") || !strcmp(field, "K") || !strcmp(field, "relaxed_K"))
    {
        dims_out[0] = dims->nu[stage];
        dims_out[1] = dims->nx[stage];
    }
    else if (!strcmp(field, "r") || !strcmp(field, "relaxed_r"))
    {
        dims_out[0] = 1;
        dims_out[1] = dims->nu[stage];
    }
    else if (!strcmp(field, "q") || !strcmp(field, "relaxed_q"))
    {
        dims_out[0] = 1;
        dims_out[1] = dims->nx[stage];
    }
    else if (!strcmp(field, "idxs_rev") || !strcmp(field, "relaxed_idxs_rev"))
    {
        dims_out[0] = 1;
        dims_out[1] = dims->nb[stage] + dims->ng[stage] + dims->ni_nl[stage];
    }
    else if (!strcmp(field, "zl") || !strcmp(field, "zu") || !strcmp(field, "Zl") || !strcmp(field, "Zu")  || !strcmp(field, "idxs"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "ns", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "relaxed_zl") || !strcmp(field, "relaxed_zu") || !strcmp(field, "relaxed_Zl") || !strcmp(field, "relaxed_Zu")  || !strcmp(field, "relaxed_idxs"))
    {
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "ns", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "p"))
    {
        dims_out[0] = dims->nx[stage];
        dims_out[1] = 1;
    }
    // constraints
    else if (!strcmp(field, "C"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "ng", &dims_out[0]);
        dims_out[1] = dims->nx[stage];
    }
    else if (!strcmp(field, "D"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "ng", &dims_out[0]);
        dims_out[1] = dims->nu[stage];
    }
    else if (!strcmp(field, "lg") || !strcmp(field, "ug"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "ng", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "lbx") || !strcmp(field, "ubx") || !strcmp(field, "idxbx"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "nbx", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "lbu") || !strcmp(field, "ubu") || !strcmp(field, "idxbu"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "nbu", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "idxb"))
    {
        int tmp_int;
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "nbu", &dims_out[0]);
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "nbx", &tmp_int);
        dims_out[0] += tmp_int;
        dims_out[1] = 1;
    }
    // constraints of relaxed qp
    else if (!strcmp(field, "relaxed_C"))
    {
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "ng", &dims_out[0]);
        dims_out[1] = dims->nx[stage];
    }
    else if (!strcmp(field, "relaxed_D"))
    {
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "ng", &dims_out[0]);
        dims_out[1] = dims->nu[stage];
    }
    else if (!strcmp(field, "relaxed_lg") || !strcmp(field, "relaxed_ug"))
    {
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "ng", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "relaxed_idxb"))
    {
        int tmp_int;
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "nbu", &dims_out[0]);
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "nbx", &tmp_int);
        dims_out[0] += tmp_int;
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "relaxed_lbx") || !strcmp(field, "relaxed_ubx") || !strcmp(field, "relaxed_idxbx"))
    {
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "nbx", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "relaxed_lbu") || !strcmp(field, "relaxed_ubu") || !strcmp(field, "relaxed_idxbu"))
    {
        config->relaxed_qp_solver->dims_get(config->relaxed_qp_solver, dims->relaxed_qp_solver, stage, "nbu", &dims_out[0]);
        dims_out[1] = 1;
    }
    else if (!strcmp(field, "pcond_R"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "pcond_nu", &dims_out[0]);
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "pcond_nu", &dims_out[1]);
    }
    else if (!strcmp(field, "pcond_Q"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "pcond_nx", &dims_out[0]);
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "pcond_nx", &dims_out[1]);
    }
    else if (!strcmp(field, "pcond_S"))
    {
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "pcond_nu", &dims_out[0]);
        config->qp_solver->dims_get(config->qp_solver, dims->qp_solver, stage, "pcond_nx", &dims_out[1]);
    }
    else
    {
        printf("\nerror: ocp_nlp_qp_dims_get_from_attr: field %s not available\n", field);
        exit(1);
    }
}


void ocp_nlp_cost_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
        int stage, const char *field, int *dims_out)
{
    // vectors first
    dims_out[1] = 0;
    if (!strcmp(field, "y_ref") || !strcmp(field, "yref"))
    {
        config->cost[stage]->dims_get(config->cost[stage], dims->cost[stage],
                                            "ny", &dims_out[0]);
    }
    else if (!strcmp(field, "Zl") || !strcmp(field, "Zu") ||
             !strcmp(field, "zl") || !strcmp(field, "zu"))
    {
        dims_out[0] = dims->ns[stage];
    }
    // matrices
    else if (!strcmp(field, "W"))
    {
        config->cost[stage]->dims_get(config->cost[stage], dims->cost[stage],
                                            "ny", &dims_out[0]);

        config->cost[stage]->dims_get(config->cost[stage], dims->cost[stage],
                                            "ny", &dims_out[1]);
    }
    else if (!strcmp(field, "Vx"))
    {
        config->cost[stage]->dims_get(config->cost[stage], dims->cost[stage],
                                            "ny", &dims_out[0]);
        dims_out[1] = dims->nx[stage];
    }
    else if (!strcmp(field, "Vu"))
    {
        config->cost[stage]->dims_get(config->cost[stage], dims->cost[stage],
                                            "ny", &dims_out[0]);
        dims_out[1] = dims->nu[stage];
    }
    else if (!strcmp(field, "Vz"))
    {
        config->cost[stage]->dims_get(config->cost[stage], dims->cost[stage],
                                            "ny", &dims_out[0]);
        dims_out[1] = dims->nz[stage];
    }
    else if (!strcmp(field, "ext_cost_num_hess"))
    {
        dims_out[0] = dims->nx[stage] + dims->nu[stage];
        dims_out[1] = dims->nx[stage] + dims->nu[stage];
    }
    else if (!strcmp(field, "scaling"))
    {
        dims_out[0] = 1;
    }
    else
    {
        printf("\nerror: ocp_nlp_cost_dims_get_from_attr: field %s not available\n", field);
        exit(1);
    }
    return;
}

/************************************************
* opts
************************************************/

void *ocp_nlp_solver_opts_create(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    acados_size_t bytes = config->opts_calculate_size(config, dims);

    void *ptr = acados_calloc(1, bytes);
    assert(ptr != 0);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}



void ocp_nlp_solver_opts_set(ocp_nlp_config *config, void *opts_, const char *field, void *value)
{
    config->opts_set(config, opts_, field, value);
}


// void ocp_nlp_solver_opts_get(ocp_nlp_config *config, void *opts_, const char *field, void *value)
// {
//     config->opts_get(config, opts_, field, value);
// }


void ocp_nlp_solver_opts_set_at_stage(ocp_nlp_config *config, void *opts_, int stage, const char *field, void *value)
{
    config->opts_set_at_stage(config, opts_, stage, field, value);
}


void ocp_nlp_solver_opts_destroy(void *opts)
{
    free(opts);
}



/************************************************
* solver
************************************************/

static acados_size_t ocp_nlp_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_, ocp_nlp_in *nlp_in)
{
    acados_size_t bytes = sizeof(ocp_nlp_solver);

    bytes += config->memory_calculate_size(config, dims, opts_, nlp_in);
    bytes += config->workspace_calculate_size(config, dims, opts_, nlp_in);

    return bytes;
}



static ocp_nlp_solver *ocp_nlp_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                      void *opts_, ocp_nlp_in *nlp_in, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_solver *solver = (ocp_nlp_solver *) c_ptr;
    c_ptr += sizeof(ocp_nlp_solver);

    solver->config = config;
    solver->dims = dims;
    solver->opts = opts_;

    solver->mem = config->memory_assign(config, dims, opts_, nlp_in, c_ptr);
    // printf("\nsolver->mem %p", solver->mem);
    c_ptr += config->memory_calculate_size(config, dims, opts_, nlp_in);

    solver->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(config, dims, opts_, nlp_in);

    assert((char *) raw_memory + ocp_nlp_calculate_size(config, dims, opts_, nlp_in) == c_ptr);

    return solver;
}



ocp_nlp_solver *ocp_nlp_solver_create(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_, ocp_nlp_in *nlp_in)
{
    config->opts_update(config, dims, opts_);

    acados_size_t bytes = ocp_nlp_calculate_size(config, dims, opts_, nlp_in);

    void *ptr = acados_calloc(1, bytes);
    assert(ptr != 0);

    ocp_nlp_solver *solver = ocp_nlp_assign(config, dims, opts_, nlp_in, ptr);

    return solver;
}


void ocp_nlp_solver_destroy(ocp_nlp_solver *solver)
{
    solver->config->terminate(solver->config, solver->mem, solver->work);
    free(solver);
}


void ocp_nlp_solver_reset_qp_memory(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    solver->config->memory_reset_qp_solver(solver->config, solver->dims, nlp_in, nlp_out,
                                    solver->opts, solver->mem, solver->work);
}


int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    return solver->config->evaluate(solver->config, solver->dims, nlp_in, nlp_out,
                                    solver->opts, solver->mem, solver->work);
}


int ocp_nlp_setup_qp_matrices_and_factorize(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    return solver->config->setup_qp_matrices_and_factorize(solver->config, solver->dims, nlp_in, nlp_out,
                                    solver->opts, solver->mem, solver->work);
}


int ocp_nlp_precompute(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    return solver->config->precompute(solver->config, solver->dims, nlp_in, nlp_out,
                                      solver->opts, solver->mem, solver->work);
}



void ocp_nlp_eval_param_sens(ocp_nlp_solver *solver, char *field, int stage, int index,
                             ocp_nlp_out *sens_nlp_out)
{
    solver->config->eval_param_sens(solver->config, solver->dims, solver->opts, solver->mem,
                                    solver->work, field, stage, index, sens_nlp_out);
    return;
}

void ocp_nlp_eval_lagrange_grad_p(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, const char *field, double *out)
{
    solver->config->eval_lagr_grad_p(solver->config, solver->dims, nlp_in, solver->opts, solver->mem, solver->work, field, out);
}


void ocp_nlp_eval_residuals(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    ocp_nlp_config *config = solver->config;
    config->eval_kkt_residual(config, solver->dims, nlp_in, nlp_out, solver->opts, solver->mem, solver->work);
}


void ocp_nlp_eval_cost(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    ocp_nlp_config *config = solver->config;
    ocp_nlp_memory *nlp_mem;
    ocp_nlp_opts *nlp_opts;
    ocp_nlp_workspace *nlp_work;
    ocp_nlp_dims *dims = solver->dims;

    config->get(config, solver->dims, solver->mem, "nlp_mem", &nlp_mem);
    config->opts_get(config, solver->dims, solver->opts, "nlp_opts", &nlp_opts);
    config->work_get(config, solver->dims, solver->work, "nlp_work", &nlp_work);

    ocp_nlp_cost_compute(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
}

void ocp_nlp_eval_constraints(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    ocp_nlp_config *config = solver->config;
    ocp_nlp_memory *nlp_mem;
    ocp_nlp_opts *nlp_opts;
    ocp_nlp_workspace *nlp_work;
    ocp_nlp_dims *dims = solver->dims;

    config->get(config, solver->dims, solver->mem, "nlp_mem", &nlp_mem);
    config->opts_get(config, solver->dims, solver->opts, "nlp_opts", &nlp_opts);
    config->work_get(config, solver->dims, solver->work, "nlp_work", &nlp_work);

    ocp_nlp_eval_constraints_common(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
}


void ocp_nlp_eval_params_jac(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    ocp_nlp_config *config = solver->config;
    ocp_nlp_memory *nlp_mem;
    ocp_nlp_opts *nlp_opts;
    ocp_nlp_workspace *nlp_work;
    ocp_nlp_dims *dims = solver->dims;

    config->get(config, solver->dims, solver->mem, "nlp_mem", &nlp_mem);
    config->opts_get(config, solver->dims, solver->opts, "nlp_opts", &nlp_opts);
    config->work_get(config, solver->dims, solver->work, "nlp_work", &nlp_work);

    ocp_nlp_params_jac_compute(config, dims, nlp_in, nlp_opts, nlp_mem, nlp_work);
}


void ocp_nlp_eval_solution_sens_adj_p(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *sens_nlp_out, const char *field, int stage, double *out)
{
    solver->config->eval_solution_sens_adj_p(solver->config, solver->dims, solver->opts, solver->mem, solver->work, sens_nlp_out, field, stage, out);
}


void ocp_nlp_get(ocp_nlp_solver *solver, const char *field, void *return_value_)
{
    solver->config->get(solver->config, solver->dims, solver->mem, field, return_value_);
}



static void get_from_qp_in(ocp_qp_in *qp_in, int stage, const char *field, void *value)
{
    if (!strcmp(field, "A"))
    {
        double *double_values = value;
        d_ocp_qp_get_A(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "B"))
    {
        double *double_values = value;
        d_ocp_qp_get_B(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "b"))
    {
        double *double_values = value;
        d_ocp_qp_get_b(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "Q"))
    {
        double *double_values = value;
        d_ocp_qp_get_Q(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "R"))
    {
        double *double_values = value;
        d_ocp_qp_get_R(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "S"))
    {
        double *double_values = value;
        d_ocp_qp_get_S(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "r"))
    {
        double *double_values = value;
        d_ocp_qp_get_r(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "q"))
    {
        double *double_values = value;
        d_ocp_qp_get_q(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "lbx"))
    {
        double *double_values = value;
        d_ocp_qp_get_lbx(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "ubx"))
    {
        double *double_values = value;
        d_ocp_qp_get_ubx(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "lbu"))
    {
        double *double_values = value;
        d_ocp_qp_get_lbu(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "ubu"))
    {
        double *double_values = value;
        d_ocp_qp_get_ubu(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "C"))
    {
        double *double_values = value;
        d_ocp_qp_get_C(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "D"))
    {
        double *double_values = value;
        d_ocp_qp_get_D(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "lg"))
    {
        double *double_values = value;
        d_ocp_qp_get_lg(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "ug"))
    {
        double *double_values = value;
        d_ocp_qp_get_ug(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "zl"))
    {
        double *double_values = value;
        d_ocp_qp_get_zl(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "zu"))
    {
        double *double_values = value;
        d_ocp_qp_get_zu(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "Zl"))
    {
        double *double_values = value;
        d_ocp_qp_get_Zl(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "Zu"))
    {
        double *double_values = value;
        d_ocp_qp_get_Zu(stage, qp_in, double_values);
    }
    else if (!strcmp(field, "idxs"))
    {
        int *int_values = value;
        d_ocp_qp_get_idxs(stage, qp_in, int_values);
    }
    else if (!strcmp(field, "idxs_rev"))
    {
        int *int_values = value;
        d_ocp_qp_get_idxs_rev(stage, qp_in, int_values);
    }
    else if (!strcmp(field, "idxb"))
    {
        int *int_values = value;
        d_ocp_qp_get_idxb(stage, qp_in, int_values);
    }
    else
    {
        printf("\nerror: ocp_nlp_get_at_stage: field %s not available\n", field);
        exit(1);
    }
}


void ocp_nlp_get_at_stage(ocp_nlp_solver *solver, int stage, const char *field, void *value)
{
    ocp_nlp_dims *dims = solver->dims;
    ocp_nlp_config *config = solver->config;
    ocp_nlp_memory *nlp_mem;
    config->get(config, dims, solver->mem, "nlp_mem", &nlp_mem);

    if (!strcmp(field, "P") || !strcmp(field, "K") || !strcmp(field, "Lr") || !strcmp(field, "p"))
    {
        ocp_nlp_opts *nlp_opts;
        config->opts_get(config, dims, solver->opts, "nlp_opts", &nlp_opts);

        ocp_qp_xcond_solver_config *xcond_solver_config = config->qp_solver;

        int size1 = 0;
        int size2 = 0;
        if (!strcmp(field, "P"))
        {
            size1 = dims->nx[stage];
            size2 = dims->nx[stage];
        }
        else if (!strcmp(field, "K"))
        {
            size1 = dims->nu[stage];
            size2 = dims->nx[stage];
        }
        else if (!strcmp(field, "p"))
        {
            size1 = dims->nx[stage];
            size2 = 1;
        }
        else if (!strcmp(field, "Lr"))
        {
            size1 = dims->nu[stage];
            size2 = dims->nu[stage];
        }
        xcond_solver_config->solver_get(xcond_solver_config, nlp_mem->qp_in, nlp_mem->qp_out, nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, field, stage, value, size1, size2);
    }
    else if (!strcmp(field, "ineq_fun") || !strcmp(field, "res_stat") || !strcmp(field, "res_eq"))
    {
        ocp_nlp_memory_get_at_stage(config, dims, nlp_mem, stage, field, value);
    }
    else if (!strcmp(field, "pcond_Q"))
    {
        ocp_qp_in *pcond_qp_in;
        ocp_nlp_get(solver, "qp_xcond_in", &pcond_qp_in);
        d_ocp_qp_get_Q(stage, pcond_qp_in, value);
    }
    else if (!strcmp(field, "pcond_R"))
    {
        ocp_qp_in *pcond_qp_in;
        ocp_nlp_get(solver, "qp_xcond_in", &pcond_qp_in);
        d_ocp_qp_get_R(stage, pcond_qp_in, value);
    }
    else if (!strcmp(field, "pcond_S"))
    {
        ocp_qp_in *pcond_qp_in;
        ocp_nlp_get(solver, "qp_xcond_in", &pcond_qp_in);
        d_ocp_qp_get_S(stage, pcond_qp_in, value);
    }
    else
    {
        char *ptr_module = NULL;
        int module_length = 0;
        char module[MAX_STR_LEN];
        extract_module_name(field, module, &module_length, &ptr_module);
        ocp_qp_in *qp_in;
        const char *field_name_getter = field;
        if ( ptr_module!=NULL && (!strcmp(ptr_module, "relaxed")) )
        {
            ocp_nlp_get(solver, "relaxed_qp_in", &qp_in);
            field_name_getter = field+module_length+1;
            get_from_qp_in(qp_in, stage, field_name_getter, value);
        }
        else if ( ptr_module!=NULL && (!strcmp(ptr_module, "qpscaling")) )
        {
            field_name_getter = field+module_length+1;
            ocp_nlp_qpscaling_memory_get(dims->qpscaling, nlp_mem->qpscaling,
                field_name_getter, stage, value);
        }
        else
        {
            qp_in = nlp_mem->qp_in;
            get_from_qp_in(qp_in, stage, field_name_getter, value);
        }
    }
}


void ocp_nlp_get_from_iterate(ocp_nlp_solver *solver, int iter, int stage, const char *field, void *value)
{
    ocp_nlp_config *config = solver->config;
    ocp_nlp_memory *nlp_mem;
    ocp_nlp_dims *dims = solver->dims;

    config->get(config, solver->dims, solver->mem, "nlp_mem", &nlp_mem);

    ocp_nlp_opts *nlp_opts;
    config->opts_get(config, solver->dims, solver->opts, "nlp_opts", &nlp_opts);

    if (!nlp_opts->store_iterates)
    {
        printf("\nerror: ocp_nlp_get_from_iterate: store_iterates needs to be set to true in order to get iterates.\n");
        exit(1);
    }
    ocp_nlp_out_get(config, dims, nlp_mem->iterates[iter], stage, field, value);
}


void ocp_nlp_get_all(ocp_nlp_solver *solver, ocp_nlp_in *in, ocp_nlp_out *out, const char *field, void *value)
{
    ocp_nlp_dims *dims = solver->dims;

    double *double_values = value;
    int tmp_offset = 0;
    int N = dims->N;
    int tmp_int, stage;

    if (!strcmp(field, "x"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = dims->nx[stage];
            blasfeo_unpack_dvec(tmp_int, &out->ux[stage], dims->nu[stage], (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "u"))
    {
        for (stage = 0; stage < N; stage++)
        {
            tmp_int = dims->nu[stage];
            blasfeo_unpack_dvec(tmp_int, &out->ux[stage], 0, (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "sl"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = dims->ns[stage];
            blasfeo_unpack_dvec(tmp_int, &out->ux[stage], dims->nu[stage] + dims->nx[stage], (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "su"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = dims->ns[stage];
            blasfeo_unpack_dvec(tmp_int, &out->ux[stage], dims->nu[stage] + dims->nx[stage] + dims->ns[stage],
                                (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "s"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = 2*dims->ns[stage];
            blasfeo_unpack_dvec(tmp_int, &out->ux[stage], dims->nu[stage] + dims->nx[stage],
                                (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "z"))
    {
        for (stage = 0; stage < N; stage++)
        {
            tmp_int = dims->nz[stage];
            blasfeo_unpack_dvec(tmp_int, &out->z[stage], 0, (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "pi"))
    {
        for (stage = 0; stage < N; stage++)
        {
            tmp_int = dims->nx[stage+1];
            blasfeo_unpack_dvec(tmp_int, &out->pi[stage], 0, (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "lam"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = 2*dims->ni[stage];
            blasfeo_unpack_dvec(tmp_int, &out->lam[stage], 0, (double_values + tmp_offset), 1);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "p"))
    {
        int ii;
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = dims->np[stage];
            for (ii = 0; ii < dims->np[stage]; ii++)
            {
                double_values[tmp_offset + ii] = in->parameter_values[stage][ii];
            }
            tmp_offset += tmp_int;
        }
    }
    else
    {
        printf("\nerror: ocp_nlp_get_all: field %s not available\n", field);
        exit(1);
    }
}


void ocp_nlp_set_all(ocp_nlp_solver *solver, ocp_nlp_in *in, ocp_nlp_out *out, const char *field, void *value)
{
    ocp_nlp_dims *dims = solver->dims;

    double *double_values = value;
    int tmp_offset = 0;
    int N = dims->N;
    int tmp_int, stage;

    if (!strcmp(field, "x"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = dims->nx[stage];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->ux[stage], dims->nu[stage]);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "u"))
    {
        for (stage = 0; stage < N; stage++)
        {
            tmp_int = dims->nu[stage];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->ux[stage], 0);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "sl"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = dims->ns[stage];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->ux[stage], dims->nx[stage] + dims->nu[stage]);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "su"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = dims->ns[stage];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->ux[stage], dims->nx[stage] + dims->nu[stage] + dims->ns[stage]);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "s"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = 2*dims->ns[stage];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->ux[stage], dims->nx[stage] + dims->nu[stage]);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "z"))
    {
        for (stage = 0; stage < N; stage++)
        {
            tmp_int = dims->nz[stage];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->z[stage], 0);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "pi"))
    {
        for (stage = 0; stage < N; stage++)
        {
            tmp_int = dims->nx[stage+1];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->pi[stage], 0);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "lam"))
    {
        for (stage = 0; stage < N+1; stage++)
        {
            tmp_int = 2*dims->ni[stage];
            blasfeo_pack_dvec(tmp_int, double_values + tmp_offset, 1, &out->lam[stage], 0);
            // multiply with mask to ensure that multiplier associated with masked constraints are zero
            blasfeo_dvecmul(2*dims->ni[stage], &in->dmask[stage], 0, &out->lam[stage], 0, &out->lam[stage], 0);
            tmp_offset += tmp_int;
        }
    }
    else if (!strcmp(field, "p"))
    {
        int ii;
        for (stage = 0; stage < N+1; stage++)
        {
            for (ii = 0; ii < dims->np[stage]; ii++)
            {
                in->parameter_values[stage][ii] = double_values[tmp_offset + ii];
            }
            tmp_offset += dims->np[stage];
        }
    }
    else
    {
        printf("\nerror: ocp_nlp_set_all: field %s not available\n", field);
        exit(1);
    }
}


void ocp_nlp_set(ocp_nlp_solver *solver, int stage, const char *field, void *value)
{
    ocp_nlp_memory *mem;
    ocp_nlp_config *config = solver->config;
    config->get(config, solver->dims, solver->mem, "nlp_mem", &mem);
    // printf("called getter: nlp_mem %p\n", mem);

    ocp_nlp_dims *dims = solver->dims;

    if (!strcmp(field, "z_guess"))
    {
        int nz = dims->nz[stage];
        int nx = dims->nx[stage];
        double *double_values = value;
        blasfeo_pack_dvec(nz, double_values, 1, mem->sim_guess + stage, nx);
        if (nz > 0)
            mem->set_sim_guess[stage] = true;
        // printf("set z_guess\n");
        // blasfeo_print_exp_dvec(nz, mem->sim_guess+stage, nx);
    }
    else if (!strcmp(field, "xdot_guess"))
    {
        int nx = dims->nx[stage];
        double *double_values = value;
        blasfeo_pack_dvec(nx, double_values, 1, &mem->sim_guess[stage], 0);
        mem->set_sim_guess[stage] = true;
    }
    else if (!strcmp(field, "gnsf_phi_guess"))
    {
        int nout;
        config->dynamics[stage]->dims_get(config->dynamics[stage], dims->dynamics[stage],
                                            "gnsf_nout", &nout);
        double *double_values = value;
        blasfeo_pack_dvec(nout, double_values, 1, &mem->sim_guess[stage], 0);
        mem->set_sim_guess[stage] = true;
    }
    else
    {
        printf("\nerror: ocp_nlp_set: field %s not available\n", field);
        exit(1);
    }
}
