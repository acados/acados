/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
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
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_disc.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bghp.h"
#include "acados/ocp_nlp/ocp_nlp_reg_conv.h"
#include "acados/ocp_nlp/ocp_nlp_reg_mirror.h"
#include "acados/ocp_nlp/ocp_nlp_sqp.h"
#include "acados/ocp_nlp/ocp_nlp_sqp_rti.h"
#include "acados/utils/mem.h"


/************************************************
* plan
************************************************/

static int ocp_nlp_plan_calculate_size(int N)
{
    // N - number of shooting nodes
    int bytes = sizeof(ocp_nlp_plan);
    bytes += N * sizeof(sim_solver_plan);
    bytes += (N + 1) * sizeof(ocp_nlp_cost_t);
    bytes += N * sizeof(ocp_nlp_dynamics_t);
    bytes += (N+1) * sizeof(ocp_nlp_constraints_t);
    return bytes;
}



static ocp_nlp_plan *ocp_nlp_plan_assign(int N, void *raw_memory)
{
    int ii;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_plan *plan = (ocp_nlp_plan *) c_ptr;
    c_ptr += sizeof(ocp_nlp_plan);

    plan->sim_solver_plan = (sim_solver_plan *) c_ptr;
    c_ptr += N * sizeof(sim_solver_plan);

    plan->nlp_cost = (ocp_nlp_cost_t *) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_cost_t);

    plan->nlp_dynamics = (ocp_nlp_dynamics_t *) c_ptr;
    c_ptr += N * sizeof(ocp_nlp_dynamics_t);

    plan->nlp_constraints = (ocp_nlp_constraints_t *) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_constraints_t);

    // initialize to default value !=0 to detect empty plans
    for (ii=0; ii <= N; ii++)
        plan->nlp_cost[ii] = INVALID_COST;
    for (ii=0; ii < N; ii++)
        plan->nlp_dynamics[ii] = INVALID_MODEL;
    for (ii=0; ii <= N; ii++)
        plan->nlp_constraints[ii] = INVALID_CONSTRAINT;

    return plan;
}



static void ocp_nlp_plan_initialize_default(int N, ocp_nlp_plan *plan)
{
    plan->nlp_solver = SQP;
    plan->regularization = NO_REGULARIZATION;
    plan->N = N;

    for (int ii = 0; ii <= N; ii++)
    {
        plan->nlp_cost[ii] = NONLINEAR_LS;
    }

    for (int ii = 0; ii < N; ii++)
    {
        plan->sim_solver_plan[ii].sim_solver = ERK;
    }
}



ocp_nlp_plan *ocp_nlp_plan_create(int N)
{
    int bytes = ocp_nlp_plan_calculate_size(N);
    void *ptr = acados_malloc(bytes, 1);

    ocp_nlp_plan *plan = ocp_nlp_plan_assign(N, ptr);

    ocp_nlp_plan_initialize_default(N, plan);

    return plan;
}

void ocp_nlp_plan_destroy(void* plan_)
{
    free(plan_);
}


/************************************************
* config
************************************************/

ocp_nlp_config *ocp_nlp_config_create(ocp_nlp_plan plan)
{
    int N = plan.N;

    /* calculate_size & malloc & assign */

    int bytes = ocp_nlp_config_calculate_size(N);
    void *config_mem = acados_calloc(1, bytes);
    ocp_nlp_config *config = ocp_nlp_config_assign(N, config_mem);

    /* initialize config according plan */

    if (plan.nlp_solver == SQP)
    {
        ocp_nlp_sqp_config_initialize_default(config);
    }
    else if (plan.nlp_solver == SQP_RTI)
    {
        ocp_nlp_sqp_rti_config_initialize_default(config);
    }
    else
    {
        printf("Solver not available!\n");
        exit(1);
    }

    // QP solver
    ocp_qp_xcond_solver_config_initialize_default(plan.ocp_qp_solver_plan.qp_solver, config->qp_solver);

    // regularization
    switch (plan.regularization)
    {
        case NO_REGULARIZATION:
            config->regularization = NULL;  // Note(oj): this is maybe a bit dirty
            break;
        case MIRROR:
            ocp_nlp_reg_mirror_config_initialize_default(config->regularization);
            break;
        case CONVEXIFICATION:
            ocp_nlp_reg_conv_config_initialize_default(config->regularization);
            break;
        default:
            printf("Regularization option not available!\n");
            exit(1); 
    }

    // cost
    for (int i = 0; i <= N; ++i)
    {
        switch (plan.nlp_cost[i])
        {
            case LINEAR_LS:
                ocp_nlp_cost_ls_config_initialize_default(config->cost[i]);
                break;
            case NONLINEAR_LS:
                ocp_nlp_cost_nls_config_initialize_default(config->cost[i]);
                break;
            case EXTERNALLY_PROVIDED:
                ocp_nlp_cost_external_config_initialize_default(config->cost[i]);
                break;
            case INVALID_COST:
                printf("\nInvalid cost module type\nForgot to initialize?\n\n");
                exit(1);
            default:
                printf("Cost not available!\n");
                exit(1);
        }
    }

    // Dynamics
    for (int i = 0; i < N; ++i)
    {
        switch (plan.nlp_dynamics[i])
        {
            case CONTINUOUS_MODEL:
                ocp_nlp_dynamics_cont_config_initialize_default(config->dynamics[i]);
//                config->dynamics[i]->sim_solver = sim_config_create(plan.sim_solver_plan[i]);
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
                        printf("\n\nSpecified integrator not available in acados C interface!\n\n");
                        exit(1);
                }

                break;
            case DISCRETE_MODEL:
                ocp_nlp_dynamics_disc_config_initialize_default(config->dynamics[i]);
                break;
            case INVALID_MODEL:
                printf("\nInvalid dynamic module type\nForgot to initialize?\n\n");
                exit(1);
            default:
                printf("Dynamics not available!\n");
                exit(1);
        }
    }

    // Constraints
    for (int i = 0; i <= N; ++i)
    {
        switch (plan.nlp_constraints[i])
        {
            case BGH:
                ocp_nlp_constraints_bgh_config_initialize_default(config->constraints[i]);
                break;
            case BGHP:
                ocp_nlp_constraints_bghp_config_initialize_default(config->constraints[i]);
                break;
            case INVALID_CONSTRAINT:
                printf("\nInvalid constraint module type\nForgot to initialize?\n\n");
                exit(1);
            default:
                printf("\nConstraint not available!\n\n");
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

    int bytes = ocp_nlp_dims_calculate_size(config);

    void *ptr = acados_calloc(1, bytes);

    ocp_nlp_dims *dims = ocp_nlp_dims_assign(config, ptr);

    return dims;
}



void ocp_nlp_dims_destroy(void *dims_)
{
    free(dims_);
}



/************************************************
* input
************************************************/

ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    int bytes = ocp_nlp_in_calculate_size(config, dims);

    void *ptr = acados_calloc(1, bytes);

    ocp_nlp_in *nlp_in = ocp_nlp_in_assign(config, dims, ptr);

    return nlp_in;
}



void ocp_nlp_in_destroy(void *in)
{
    free(in);
}



void ocp_nlp_in_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage,
        const char *field, void *value)
{
    int ii;

    int N = dims->N;

    if (!strcmp(field, "Ts"))
    {
        double *Ts_values = value;
        for (ii=0; ii<N; ii++)
            in->Ts[ii] = *Ts_values;
    }
    else
    {
        printf("\nerror: ocp_nlp_in_set: field %s not available\n", field);
        exit(1);
    }
    return;
}



int ocp_nlp_dynamics_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                           int stage, const char *fun_type, void *fun_ptr)
{
    sim_config *sim_config = config->dynamics[stage]->sim_solver;
    ocp_nlp_dynamics_cont_model *dynamics = in->dynamics[stage];

    int status = sim_config->model_set(dynamics->sim_model, (char *) fun_type, fun_ptr);

    return status;
}



int ocp_nlp_cost_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
                           ocp_nlp_in *in, int stage,
                           const char *field, void *value)
{
    ocp_nlp_cost_config *cost_config = config->cost[stage];
    void *cost_model = in->cost[stage];
    void *cost_dims = dims->cost[stage];

    int status = cost_config->model_set(cost_config, cost_dims, cost_model, field, value);

    return status;
}



// TODO remove and use ocp_nlp_dynamics_model_set instead !!!
int nlp_set_discrete_model_in_stage(ocp_nlp_config *config, ocp_nlp_in *in, int stage,
                                    void *fun_ptr)
{

    ocp_nlp_dynamics_disc_model *dynamics = in->dynamics[stage];
    dynamics->discrete_model = (external_function_generic *) fun_ptr;

    return ACADOS_SUCCESS;
}


int ocp_nlp_constraints_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
             ocp_nlp_in *in, int stage, const char *field, void *value)
{
    ocp_nlp_constraints_config *constr_config = config->constraints[stage];
    void *constr_dims = dims->constraints[stage];

    return constr_config->model_set(constr_config, constr_dims,
                                      in->constraints[stage], field, value);
}

/************************************************
* out
************************************************/

ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    int bytes = ocp_nlp_out_calculate_size(config, dims);

    void *ptr = acados_calloc(1, bytes);

    ocp_nlp_out *nlp_out = ocp_nlp_out_assign(config, dims, ptr);

    // initialize to zeros
    for (int ii = 0; ii <= dims->N; ++ii)
        blasfeo_dvecse(dims->qp_solver->nu[ii] + dims->qp_solver->nx[ii], 0.0, nlp_out->ux + ii, 0);

    return nlp_out;
}



void ocp_nlp_out_destroy(void *out)
{
    free(out);
}



void ocp_nlp_out_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
                     int stage, const char *field, void *value)
{
    if (!strcmp(field, "x"))
    {
        double *double_values = value;
        blasfeo_pack_dvec(dims->nx[stage], double_values, &out->ux[stage], dims->nu[stage]);
    }
    else if (!strcmp(field, "u"))
    {
        double *double_values = value;
        blasfeo_pack_dvec(dims->nu[stage], double_values, &out->ux[stage], 0);
    }
    else
    {
        printf("\nerror: ocp_nlp_out_set: field %s not available\n", field);
        exit(1);
    }
}



void ocp_nlp_out_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
                     int stage, const char *field, void *value)
{
    if (!strcmp(field, "x"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->nx[stage], &out->ux[stage], dims->nu[stage], double_values);
    }
    else if (!strcmp(field, "u"))
    {
        double *double_values = value;
        blasfeo_unpack_dvec(dims->nu[stage], &out->ux[stage], 0, double_values);
    }
    else
    {
        printf("\nerror: ocp_nlp_out_get: field %s not available\n", field);
        exit(1);
    }
}



/************************************************
* opts
************************************************/

void *ocp_nlp_opts_create(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    int bytes = config->opts_calculate_size(config, dims);

    void *ptr = acados_calloc(1, bytes);

    void *opts = config->opts_assign(config, dims, ptr);

    config->opts_initialize_default(config, dims, opts);

    return opts;
}


void ocp_nlp_opts_set(ocp_nlp_config *config, void *opts_,
                      const char *field, const void *value)
{
    config->opts_set(config, opts_, field, value);
}


int ocp_nlp_dynamics_opts_set(ocp_nlp_config *config, void *opts_, int stage,
                                         const char *field, void *value)
{
    return config->dynamics_opts_set(config, opts_, stage, field, value);
}

void ocp_nlp_opts_update(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)
{
    config->opts_update(config, dims, opts_);
}



void ocp_nlp_opts_destroy(void *opts)
{
    free(opts);
}


/************************************************
* solver
************************************************/

static int ocp_nlp_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)
{
    int bytes = sizeof(ocp_nlp_solver);

    bytes += config->memory_calculate_size(config, dims, opts_);
    bytes += config->workspace_calculate_size(config, dims, opts_);

    return bytes;
}



static ocp_nlp_solver *ocp_nlp_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                      void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_solver *solver = (ocp_nlp_solver *) c_ptr;
    c_ptr += sizeof(ocp_nlp_solver);

    solver->config = config;
    solver->dims = dims;
    solver->opts = opts_;

    solver->mem = config->memory_assign(config, dims, opts_, c_ptr);
    c_ptr += config->memory_calculate_size(config, dims, opts_);

    solver->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(config, dims, opts_);

    assert((char *) raw_memory + ocp_nlp_calculate_size(config, dims, opts_) == c_ptr);

    return solver;
}



ocp_nlp_solver *ocp_nlp_solver_create(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_)
{
    config->opts_update(config, dims, opts_);

    int bytes = ocp_nlp_calculate_size(config, dims, opts_);

    void *ptr = acados_calloc(1, bytes);

    ocp_nlp_solver *solver = ocp_nlp_assign(config, dims, opts_, ptr);

    return solver;
}


void ocp_nlp_solver_destroy(void *solver)
{
    free(solver);
}



int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    return solver->config->evaluate(solver->config, solver->dims, nlp_in, nlp_out, solver->opts,
                                    solver->mem, solver->work);
}


int ocp_nlp_precompute(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    return solver->config->precompute(solver->config, solver->dims, nlp_in, nlp_out, solver->opts,
                                    solver->mem, solver->work);
}

void ocp_nlp_get(ocp_nlp_config *config, ocp_nlp_solver *solver,
                 const char *field, void *return_value_)
{
    solver->config->get(solver->config, solver->mem, field, return_value_);
}
