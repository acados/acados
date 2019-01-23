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

#include "acados/ocp_nlp/ocp_nlp_sqp.h"

// external
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if defined(ACADOS_WITH_OPENMP)
#include <omp.h>
#endif

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/ocp_nlp/ocp_nlp_reg_conv.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/math.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



/************************************************
 * options
 ************************************************/

int ocp_nlp_sqp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_opts);

    size += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    if (config->regularization != NULL)
        size += config->regularization->opts_calculate_size();

    // dynamics
    size += (N + 1) * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        size += dynamics[ii]->opts_calculate_size(dynamics[ii], dims->dynamics[ii]);
    }
        
    if (config->dynamics[N] != NULL) 
        size += dynamics[N]->opts_calculate_size(dynamics[N], dims->dynamics[N]);

    // cost
    size += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        size += cost[ii]->opts_calculate_size(cost[ii], dims->cost[ii]);
    }

    // constraints
    size += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        size += constraints[ii]->opts_calculate_size(constraints[ii], dims->constraints[ii]);
    }

    return size;
}



void *ocp_nlp_sqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_opts);

    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    if (config->regularization != NULL)
    {
        opts->reg_opts = config->regularization->opts_assign(c_ptr);
        c_ptr += config->regularization->opts_calculate_size();
    }

    // dynamics
    opts->dynamics = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        opts->dynamics[ii] = dynamics[ii]->opts_assign(dynamics[ii], dims->dynamics[ii], c_ptr);
        c_ptr += dynamics[ii]->opts_calculate_size(dynamics[ii], dims->dynamics[ii]);
    }

    if (config->dynamics[N] != NULL) 
    {
        opts->dynamics[N] = dynamics[N]->opts_assign(dynamics[N], dims->dynamics[N], c_ptr);
        c_ptr += dynamics[N]->opts_calculate_size(dynamics[N], dims->dynamics[N]);
    }

    // cost
    opts->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        opts->cost[ii] = cost[ii]->opts_assign(cost[ii], dims->cost[ii], c_ptr);
        c_ptr += cost[ii]->opts_calculate_size(cost[ii], dims->cost[ii]);
    }

    // constraints
    opts->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        opts->constraints[ii] =
            constraints[ii]->opts_assign(constraints[ii], dims->constraints[ii], c_ptr);
        c_ptr += constraints[ii]->opts_calculate_size(constraints[ii], dims->constraints[ii]);
    }

    assert((char *) raw_memory + ocp_nlp_sqp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_sqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int ii;

    int N = dims->N;

    // SQP opts

    opts->maxIter = 20;
    opts->min_res_g = 1e-8;
    opts->min_res_b = 1e-8;
    opts->min_res_d = 1e-8;
    opts->min_res_m = 1e-8;

    opts->reuse_workspace = 1;
    opts->num_threads = 4;

    // submodules opts

    qp_solver->opts_initialize_default(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    if (config->regularization != NULL)
    {
        opts->reg_opts->delta = 1e-4;
        opts->reg_opts->gamma = 0.0;
    }

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        dynamics[ii]->opts_initialize_default(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
    }

    if (config->dynamics[N] != NULL) 
        dynamics[N]->opts_initialize_default(dynamics[N], dims->dynamics[N], opts->dynamics[N]);

    // cost
    for (ii = 0; ii <= N; ii++)
    {
        cost[ii]->opts_initialize_default(cost[ii], dims->cost[ii], opts->cost[ii]);
    }

    // constraints
    for (ii = 0; ii <= N; ii++)
    {
        constraints[ii]->opts_initialize_default(constraints[ii], dims->constraints[ii],
                                                 opts->constraints[ii]);
    }

    return;
}



void ocp_nlp_sqp_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int ii;

    int N = dims->N;

    qp_solver->opts_update(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        dynamics[ii]->opts_update(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
    }

    if (config->dynamics[N] != NULL) 
        dynamics[N]->opts_update(dynamics[N], dims->dynamics[N], opts->dynamics[N]);

    // cost
    for (ii = 0; ii <= N; ii++)
    {
        cost[ii]->opts_update(cost[ii], dims->cost[ii], opts->cost[ii]);
    }

    // constraints
    for (ii = 0; ii <= N; ii++)
    {
        constraints[ii]->opts_update(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
    }

    return;
}



void ocp_nlp_sqp_opts_set(void *config_, void *opts_, const char *field, const void* value)
{
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;
    ocp_nlp_solver_config *config = config_;

    if (!strcmp(field, "maxIter"))
    {
        int* maxIter = (int *) value;
        opts->maxIter = *maxIter;
    }
    else if (!strcmp(field, "reuse_workspace"))
    {
        int* reuse_workspace = (int *) value;
        opts->reuse_workspace = *reuse_workspace;
    }
    else if (!strcmp(field, "num_threads"))
    {
        int* num_threads = (int *) value;
        opts->num_threads = *num_threads;
    }
    else if (!strcmp(field, "min_res_g"))
    {
        double* min_res_g = (double *) value;
        opts->min_res_g = *min_res_g;
    }
    else if (!strcmp(field, "min_res_b"))
    {
        double* min_res_b = (double *) value;
        opts->min_res_b = *min_res_b;
    }
    else if (!strcmp(field, "min_res_d"))
    {
        double* min_res_d = (double *) value;
        opts->min_res_d = *min_res_d;
    }
    else if (!strcmp(field, "min_res_m"))
    {
        double* min_res_m = (double *) value;
        opts->min_res_m = *min_res_m;
    }
    else
    {
        config->qp_solver->opts_set(config->qp_solver, opts->qp_solver_opts, field, value);
    }
}

int ocp_nlp_sqp_dynamics_opts_set(void *config_, void *opts_, int stage,
                                     const char *field, void *value)
{
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_dynamics_config *dyn_config = config->dynamics[stage];

    return dyn_config->opts_set(dyn_config, opts->dynamics[stage], field, value);

}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_sqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // extract dims
    int N = dims->N;
    // ocp_nlp_cost_dims **cost_dims = dims->cost;
    // int ny;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_memory);

    size += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    if (config->regularization != NULL)
        size += config->regularization->memory_calculate_size(dims->qp_solver);

    // dynamics
    size += (N + 1)* sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        size += dynamics[ii]->memory_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                    opts->dynamics[ii]);
    }

    if (config->dynamics[N] != NULL) 
        size += dynamics[N]->memory_calculate_size(dynamics[N], dims->dynamics[N],
                                                    opts->dynamics[N]);

    // cost
    size += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        size += cost[ii]->memory_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
    }

    // constraints
    size += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        size += constraints[ii]->memory_calculate_size(constraints[ii], dims->constraints[ii],
                                                       opts->constraints[ii]);
    }

    // nlp res
    size += ocp_nlp_res_calculate_size(dims);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims);

    size += 8;  // initial align

    //    make_int_multiple_of(64, &size);

    return size;
}



void *ocp_nlp_sqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int N = dims->N;
    // ocp_nlp_cost_dims **cost_dims = dims->cost;
    // int ny;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_sqp_memory *mem = (ocp_nlp_sqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_memory);

    // QP solver
    mem->qp_solver_mem =
        qp_solver->memory_assign(qp_solver, dims->qp_solver, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    if (config->regularization != NULL)
    {
        mem->reg_mem = config->regularization->memory_assign(dims->qp_solver, c_ptr);
        c_ptr += config->regularization->memory_calculate_size(dims->qp_solver);
    }

    // nlp res
    mem->nlp_res = ocp_nlp_res_assign(dims, c_ptr);
    c_ptr += mem->nlp_res->memsize;

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims);

    // dynamics
    mem->dynamics = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        mem->dynamics[ii] = dynamics[ii]->memory_assign(dynamics[ii], dims->dynamics[ii],
                                                        opts->dynamics[ii], c_ptr);
        c_ptr += dynamics[ii]->memory_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                     opts->dynamics[ii]);
    }

    if (config->dynamics[N] != NULL) 
    {
        mem->dynamics[N] = dynamics[N]->memory_assign(dynamics[N], dims->dynamics[N],
                                                        opts->dynamics[N], c_ptr);
        c_ptr += dynamics[N]->memory_calculate_size(dynamics[N], dims->dynamics[N],
                                                     opts->dynamics[N]);
    }

    // cost
    mem->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        mem->cost[ii] = cost[ii]->memory_assign(cost[ii], dims->cost[ii], opts->cost[ii], c_ptr);
        c_ptr += cost[ii]->memory_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
    }

    // constraints
    mem->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        mem->constraints[ii] = constraints[ii]->memory_assign(
            constraints[ii], dims->constraints[ii], opts->constraints[ii], c_ptr);
        c_ptr += constraints[ii]->memory_calculate_size(constraints[ii], dims->constraints[ii],
                                                        opts->constraints[ii]);
    }

    assert((char *) raw_memory + ocp_nlp_sqp_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

int ocp_nlp_sqp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // loop index
    int ii;

    // extract dims
    int N = dims->N;

    int size = 0;
    int size_tmp = 0;
    int tmp;

    // sqp
    size += sizeof(ocp_nlp_sqp_work);

    // array of pointers
    // cost
    size += (N + 1) * sizeof(void *);
    // dynamics
    size += N * sizeof(void *);
    // constraints
    size += (N + 1) * sizeof(void *);

    // qp in
    size += ocp_qp_in_calculate_size(qp_solver, dims->qp_solver);

    // qp out
    size += ocp_qp_out_calculate_size(qp_solver, dims->qp_solver);

    if (opts->reuse_workspace)
    {

#if defined(ACADOS_WITH_OPENMP)

        // qp solver
        size += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver,
            opts->qp_solver_opts);

        // dynamics
        for (ii = 0; ii < N; ii++)
        {
            size += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                           opts->dynamics[ii]);
        }

        if (config->dynamics[N] != NULL) 
            size += dynamics[N]->workspace_calculate_size(dynamics[N], dims->dynamics[N],
                                                           opts->dynamics[N]);

        // cost
        for (ii = 0; ii <= N; ii++)
        {
            size += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
        }

        // constraints
        for (ii = 0; ii <= N; ii++)
        {
            size += constraints[ii]->workspace_calculate_size(constraints[ii],
                dims->constraints[ii], opts->constraints[ii]);
        }

#else

        // qp solver
        tmp = qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);
        size_tmp = tmp > size_tmp ? tmp : size_tmp;

        // dynamics
        for (ii = 0; ii < N; ii++)
        {
            tmp = dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                           opts->dynamics[ii]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        if (config->dynamics[N] != NULL) 
        {
            tmp = dynamics[N]->workspace_calculate_size(dynamics[N], dims->dynamics[N],
                                                           opts->dynamics[N]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // cost
        for (ii = 0; ii <= N; ii++)
        {
            tmp = cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // constraints
        for (ii = 0; ii <= N; ii++)
        {
            tmp = constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii],
                                                              opts->constraints[ii]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }


        size += size_tmp;

#endif

    }
    else
    {

        // qp solver
        size += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver,
            opts->qp_solver_opts);

        // dynamics
        for (ii = 0; ii < N; ii++)
        {
            size += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                           opts->dynamics[ii]);
        }

        if (config->dynamics[N] != NULL) 
            size += dynamics[N]->workspace_calculate_size(dynamics[N], dims->dynamics[N],
                                                           opts->dynamics[N]);

        // cost
        for (ii = 0; ii <= N; ii++)
        {
            size += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
        }

        // constraints
        for (ii = 0; ii <= N; ii++)
        {
            size += constraints[ii]->workspace_calculate_size(constraints[ii],
                dims->constraints[ii], opts->constraints[ii]);
        }

    }

    return size;
}



// TODO(all): introduce member "memsize" in all structures to make on-line cast cheaper (i.e. avoid
// to calculate size on-line)
static void ocp_nlp_sqp_cast_workspace(void *config_, ocp_nlp_dims *dims, ocp_nlp_sqp_work *work,
                                       ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_opts *opts)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // extract dims
    int N = dims->N;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_work);

    // array of pointers
    //
    work->dynamics = (void **) c_ptr;
    c_ptr += (N + 1)* sizeof(void *);
    //
    work->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    //
    work->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    // qp in
    work->qp_in = ocp_qp_in_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(qp_solver, dims->qp_solver);

    // qp out
    work->qp_out = ocp_qp_out_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(qp_solver, dims->qp_solver);

    if (opts->reuse_workspace)
    {

#if defined(ACADOS_WITH_OPENMP)

        // qp solver
        work->qp_work = (void *) c_ptr;
        c_ptr += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver,
            opts->qp_solver_opts);

        // dynamics
        for (int ii = 0; ii < N; ii++)
        {
            work->dynamics[ii] = c_ptr;
            c_ptr += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                            opts->dynamics[ii]);
        }

        if (config->dynamics[N] != NULL) 
        {
            work->dynamics[N] = c_ptr;
            c_ptr += dynamics[N]->workspace_calculate_size(dynamics[N], dims->dynamics[N],
                                                            opts->dynamics[N]);
        }

        // cost
        for (int ii = 0; ii <= N; ii++)
        {
            work->cost[ii] = c_ptr;
            c_ptr += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
        }

        // constraints
        for (int ii = 0; ii <= N; ii++)
        {
            work->constraints[ii] = c_ptr;
            c_ptr += constraints[ii]->workspace_calculate_size(constraints[ii],
                dims->constraints[ii], opts->constraints[ii]);
        }

#else

        // qp solver
        work->qp_work = (void *) c_ptr;

        // dynamics
        for (int ii = 0; ii < N; ii++)
        {
            work->dynamics[ii] = c_ptr;
        }

        if (config->dynamics[N] != NULL) 
            work->dynamics[N] = c_ptr;

        // cost
        for (int ii = 0; ii <= N; ii++)
        {
            work->cost[ii] = c_ptr;
        }

        // constraints
        for (int ii = 0; ii <= N; ii++)
        {
            work->constraints[ii] = c_ptr;
        }

#endif

    }
    else
    {

        // qp solver
        work->qp_work = (void *) c_ptr;
        c_ptr += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver,
            opts->qp_solver_opts);

        // dynamics
        for (int ii = 0; ii < N; ii++)
        {
            work->dynamics[ii] = c_ptr;
            c_ptr += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                            opts->dynamics[ii]);
        }

        if (config->dynamics[N] != NULL) 
        {
            work->dynamics[N] = c_ptr;
            c_ptr += dynamics[N]->workspace_calculate_size(dynamics[N], dims->dynamics[N],
                                                            opts->dynamics[N]);
        }

        // cost
        for (int ii = 0; ii <= N; ii++)
        {
            work->cost[ii] = c_ptr;
            c_ptr += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
        }

        // constraints
        for (int ii = 0; ii <= N; ii++)
        {
            work->constraints[ii] = c_ptr;
            c_ptr += constraints[ii]->workspace_calculate_size(constraints[ii],
                dims->constraints[ii], opts->constraints[ii]);
        }

    }

    // assert & return
    assert((char *) work + ocp_nlp_sqp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

static void initialize_qp(void *config_, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                          ocp_nlp_out *nlp_out, ocp_nlp_sqp_opts *opts, ocp_nlp_sqp_memory *mem,
                          ocp_nlp_sqp_work *work)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    // loop index
    int ii;

    // extract dims
    int N = dims->N;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (ii = 0; ii <= N; ii++)
    {
        // cost
        config->cost[ii]->initialize(config->cost[ii], dims->cost[ii], nlp_in->cost[ii],
                                     opts->cost[ii], mem->cost[ii], work->cost[ii]);
        // dynamics
        if (ii < N)
            config->dynamics[ii]->initialize(config->dynamics[ii], dims->dynamics[ii],
                                         nlp_in->dynamics[ii], opts->dynamics[ii],
                                         mem->dynamics[ii], work->dynamics[ii]);
        if (config->dynamics[N] != NULL) 
            config->dynamics[N]->initialize(config->dynamics[N], dims->dynamics[N],
                                         nlp_in->dynamics[N], opts->dynamics[N],
                                         mem->dynamics[N], work->dynamics[N]);
        // constraints
        config->constraints[ii]->initialize(config->constraints[ii], dims->constraints[ii],
                                            nlp_in->constraints[ii], opts->constraints[ii],
                                            mem->constraints[ii], work->constraints[ii]);
    }

    return;
}



static void linearize_update_qp_matrices(void *config_, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                                         ocp_nlp_out *nlp_out, ocp_nlp_sqp_opts *opts,
                                         ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_work *work)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    // loop index
    int i;

    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    /* stage-wise multiple shooting lagrangian evaluation */

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i = 0; i <= N; i++)
    {
        // dynamics
        if (i < N)
            config->dynamics[i]->update_qp_matrices(config->dynamics[i], dims->dynamics[i],
                                                nlp_in->dynamics[i], opts->dynamics[i],
                                                mem->dynamics[i], work->dynamics[i]);
        if (config->dynamics[N] != NULL) 
            config->dynamics[N]->update_qp_matrices(config->dynamics[N], dims->dynamics[N],
                                                nlp_in->dynamics[N], opts->dynamics[N],
                                                mem->dynamics[N], work->dynamics[N]);
        // cost
        config->cost[i]->update_qp_matrices(config->cost[i], dims->cost[i], nlp_in->cost[i],
                opts->cost[i], mem->cost[i], work->cost[i]);

        // constraints
        config->constraints[i]->update_qp_matrices(config->constraints[i], dims->constraints[i],
                                                   nlp_in->constraints[i], opts->constraints[i],
                                                   mem->constraints[i], work->constraints[i]);
    }

    /* collect stage-wise evaluations */

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i=0; i <= N; i++)
    {

        // nlp mem: cost_grad
        struct blasfeo_dvec *cost_grad = config->cost[i]->memory_get_grad_ptr(mem->cost[i]);
        blasfeo_dveccp(nv[i], cost_grad, 0, nlp_mem->cost_grad + i, 0);

        // nlp mem: dyn_fun
        if (i < N)
        {
            struct blasfeo_dvec *dyn_fun
                = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
            blasfeo_dveccp(nx[i + 1], dyn_fun, 0, nlp_mem->dyn_fun + i, 0);
        }

        if (config->dynamics[N] != NULL) 
        {
            struct blasfeo_dvec *dyn_fun
                = config->dynamics[N]->memory_get_fun_ptr(mem->dynamics[N]);
            blasfeo_dveccp(nx[N + 1], dyn_fun, 0, nlp_mem->dyn_fun + N, 0);
        }

        // nlp mem: dyn_adj
        if (i < N)
        {
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[i]->memory_get_adj_ptr(mem->dynamics[i]);
            blasfeo_dveccp(nu[i] + nx[i], dyn_adj, 0, nlp_mem->dyn_adj + i, 0);
        }
        if (config->dynamics[N] != NULL) {
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[N]->memory_get_adj_ptr(mem->dynamics[N]);
            blasfeo_dveccp(nu[N] + nx[N], dyn_adj, 0, nlp_mem->dyn_adj + N, 0);
        } 
        else
        {
            blasfeo_dvecse(nu[N] + nx[N], 0.0, nlp_mem->dyn_adj + N, 0);
        }

        if (i > 0)
        {
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[i-1]->memory_get_adj_ptr(mem->dynamics[i-1]);
            blasfeo_daxpy(nx[i], 1.0, dyn_adj, nu[i-1]+nx[i-1], nlp_mem->dyn_adj+i, nu[i],
                nlp_mem->dyn_adj+i, nu[i]);
        }

        if (config->dynamics[N] != NULL) {
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[N-1]->memory_get_adj_ptr(mem->dynamics[N-1]);
            blasfeo_daxpy(nx[N], 1.0, dyn_adj, nu[N-1]+nx[N-1], nlp_mem->dyn_adj+N, nu[N],
                nlp_mem->dyn_adj+N, nu[N]);
        }

        // nlp mem: ineq_fun
        struct blasfeo_dvec *ineq_fun =
            config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, nlp_mem->ineq_fun + i, 0);
        // nlp mem: ineq_adj

        struct blasfeo_dvec *ineq_adj =
            config->constraints[i]->memory_get_adj_ptr(mem->constraints[i]);
        blasfeo_dveccp(nv[i], ineq_adj, 0, nlp_mem->ineq_adj + i, 0);

    }

    // TODO(all): still to clean !!!!!!!!!!!!!

    for (i = 0; i <= N; i++)
    {
        // TODO(rien) where should the update happen??? move to qp update ???
        // TODO(all): fix and move where appropriate
        //  if(i<N)
        //  {
        //   ocp_nlp_dynamics_opts *dynamics_opts = opts->dynamics[i];
        //   sim_rk_opts *opts = dynamics_opts->sim_solver;
        //   if (opts->scheme != NULL && opts->scheme->type != exact)
        //   {
        //    for (int_t j = 0; j < nx; j++)
        //     BLASFEO_DVECEL(nlp_mem->cost_grad+i, nu+j) += work->sim_out[i]->grad[j];
        //    for (int_t j = 0; j < nu; j++)
        //     BLASFEO_DVECEL(nlp_mem->cost_grad+i, j) += work->sim_out[i]->grad[nx+j];
        //   }
        //  }
    }

    return;
}



static void regularize_hessian(void *config_, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                               ocp_nlp_out *nlp_out, ocp_nlp_sqp_opts *opts,
                               ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_work *work)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    if (config->regularization == NULL)
        return;

    config->regularization->evaluate(config->regularization, dims->qp_solver, work->qp_in,
                                     work->qp_out, opts->reg_opts, mem->reg_mem);
}



// update QP rhs for SQP (step prim var, abs dual var)
// TODO(all): move in dynamics, cost, constraints modules ???
static void sqp_update_qp_vectors(void *config_, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                                  ocp_nlp_out *nlp_out, ocp_nlp_sqp_opts *opts,
                                  ocp_nlp_sqp_memory *mem, ocp_nlp_sqp_work *work)
{
    // loop index
    int i;

    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i = 0; i <= N; i++)
    {
        // g
        blasfeo_dveccp(nv[i], nlp_mem->cost_grad + i, 0, work->qp_in->rqz + i, 0);

        // b
        if (i < N)
            blasfeo_dveccp(nx[i + 1], nlp_mem->dyn_fun + i, 0, work->qp_in->b + i, 0);

        // d
        blasfeo_dveccp(2 * ni[i], nlp_mem->ineq_fun + i, 0, work->qp_in->d + i, 0);
    }

    return;
}



static void sqp_update_variables(void *config_, ocp_nlp_dims *dims, ocp_nlp_out *nlp_out,
                                 ocp_nlp_sqp_opts *opts, ocp_nlp_sqp_memory *mem,
                                 ocp_nlp_sqp_work *work)
{
    // loop index
    int i;

    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

    // ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    // TODO(all): fix and move where appropriate
    //    for (i = 0; i < N; i++)
    //    {
    //  nx1 = dims->constraints[i+1]->nx;
    //        for (j = 0; j < nx1; j++)
    //        {
    //            work->sim_in[i]->S_adj[j] = -BLASFEO_DVECEL(&work->qp_out->pi[i], j);
    //        }
    //    }

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i = 0; i <= N; i++)
    {
        // (full) step in primal variables

        blasfeo_daxpy(nv[i], 1.0, work->qp_out->ux + i, 0, nlp_out->ux + i, 0, nlp_out->ux + i, 0);

        // absolute in dual variables

        if (i < N)
            blasfeo_dveccp(nx[i + 1], work->qp_out->pi + i, 0, nlp_out->pi + i, 0);

        blasfeo_dveccp(2 * ni[i], work->qp_out->lam + i, 0, nlp_out->lam + i, 0);

        blasfeo_dveccp(2 * ni[i], work->qp_out->t + i, 0, nlp_out->t + i, 0);

    }

    return;
}



// Simple fixed-step Gauss-Newton based SQP routine
int ocp_nlp_sqp(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{

    // acados timer
    acados_timer timer0, timer1;

    // start timer
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_work *work = work_;

    ocp_nlp_sqp_cast_workspace(config, dims, work, mem, opts);

    // zero timers
    double total_time = 0.0;
    mem->time_qp_sol = 0.0;
    mem->time_lin = 0.0;
    mem->time_tot = 0.0;

    // extract dims
    int N = dims->N;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->num_threads);
    #pragma omp parallel
    { // beginning of parallel region
#endif

    // alias to dynamics_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (int ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->memory_set_ux_ptr(nlp_out->ux + ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_ux1_ptr(nlp_out->ux + ii + 1, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_pi_ptr(nlp_out->pi + ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_BAbt_ptr(work->qp_in->BAbt + ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_RSQrq_ptr(work->qp_in->RSQrq + ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_z_ptr(nlp_out->z, mem->dynamics[ii]);
    }

    if (config->dynamics[N] != NULL)
    {
        config->dynamics[N]->memory_set_ux_ptr(nlp_out->ux + N, mem->dynamics[N]);
        config->dynamics[N]->memory_set_ux1_ptr(nlp_out->ux + N + 1, mem->dynamics[N]);
        config->dynamics[N]->memory_set_pi_ptr(nlp_out->pi + N, mem->dynamics[N]);
        config->dynamics[N]->memory_set_BAbt_ptr(work->qp_in->BAbt + N, mem->dynamics[N]);
        config->dynamics[N]->memory_set_RSQrq_ptr(work->qp_in->RSQrq + N, mem->dynamics[N]);
        config->dynamics[N]->memory_set_z_ptr(nlp_out->z, mem->dynamics[N]);
    }

    // alias to cost_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (int ii = 0; ii < N; ii++)
    {
        config->cost[ii]->memory_set_ux_ptr(nlp_out->ux + ii, mem->cost[ii]);
        if (dims->nz[ii] > 0) {
            config->cost[ii]->memory_set_z_ptr(
                    &(((ocp_nlp_dynamics_cont_memory *) mem->dynamics[ii])->z_out),
                    mem->cost[ii]);
            config->cost[ii]->memory_set_dzdxu_tran_ptr(
                    &(((ocp_nlp_dynamics_cont_memory *) mem->dynamics[ii])->dzdxu_tran),
                    mem->cost[ii]);
        }
        config->cost[ii]->memory_set_RSQrq_ptr(work->qp_in->RSQrq + ii, mem->cost[ii]);
        config->cost[ii]->memory_set_Z_ptr(work->qp_in->Z + ii, mem->cost[ii]);
    }

    if (config->dynamics[N] != NULL)
    {
        config->cost[N]->memory_set_ux_ptr(nlp_out->ux + N, mem->cost[N]);
        if (dims->nz[N] > 0) {
            config->cost[N]->memory_set_z_ptr(
                    &(((ocp_nlp_dynamics_cont_memory *) mem->dynamics[N])->z_out),
                    mem->cost[N]);
            config->cost[N]->memory_set_dzdxu_tran_ptr(
                    &(((ocp_nlp_dynamics_cont_memory *) mem->dynamics[N])->dzdxu_tran),
                    mem->cost[N]);
        }
    }

    config->cost[N]->memory_set_ux_ptr(nlp_out->ux + N, mem->cost[N]);
    config->cost[N]->memory_set_RSQrq_ptr(work->qp_in->RSQrq + N, mem->cost[N]);
    config->cost[N]->memory_set_Z_ptr(work->qp_in->Z + N, mem->cost[N]);

    // alias to constraints_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (int ii = 0; ii <= N; ii++)
    {
        config->constraints[ii]->memory_set_ux_ptr(nlp_out->ux + ii, mem->constraints[ii]);
        config->constraints[ii]->memory_set_lam_ptr(nlp_out->lam + ii, mem->constraints[ii]);
        config->constraints[ii]->memory_set_DCt_ptr(work->qp_in->DCt + ii, mem->constraints[ii]);
        config->constraints[ii]->memory_set_RSQrq_ptr(work->qp_in->RSQrq + ii,
                                                      mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxb_ptr(work->qp_in->idxb[ii], mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxs_ptr(work->qp_in->idxs[ii], mem->constraints[ii]);
    }

    // copy sampling times into dynamics model
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (int ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->model_set_T(nlp_in->Ts[ii], nlp_in->dynamics[ii]);
    }

#if defined(ACADOS_WITH_OPENMP)
    } // end of parallel region
#endif

    // initialize QP
    initialize_qp(config, dims, nlp_in, nlp_out, opts, mem, work);

    // main sqp loop
    int max_sqp_iterations = opts->maxIter;
    int sqp_iter = 0;
    for (; sqp_iter < max_sqp_iterations; sqp_iter++)
    {

        // start timer
        acados_tic(&timer1);

        // linearizate NLP and update QP matrices
        linearize_update_qp_matrices(config, dims, nlp_in, nlp_out, opts, mem, work);

        // stop timer
        mem->time_lin += acados_toc(&timer1);

        // regularize Hessian
        regularize_hessian(config, dims, nlp_in, nlp_out, opts, mem, work);

        // update QP rhs for SQP (step prim var, abs dual var)
        sqp_update_qp_vectors(config, dims, nlp_in, nlp_out, opts, mem, work);

        // compute nlp residuals
        ocp_nlp_res_compute(dims, nlp_in, nlp_out, mem->nlp_res, mem->nlp_mem);

        nlp_out->inf_norm_res = mem->nlp_res->inf_norm_res_g;
        nlp_out->inf_norm_res = (mem->nlp_res->inf_norm_res_b > nlp_out->inf_norm_res) ?
                                    mem->nlp_res->inf_norm_res_b :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (mem->nlp_res->inf_norm_res_d > nlp_out->inf_norm_res) ?
                                    mem->nlp_res->inf_norm_res_d :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (mem->nlp_res->inf_norm_res_m > nlp_out->inf_norm_res) ?
                                    mem->nlp_res->inf_norm_res_m :
                                    nlp_out->inf_norm_res;

        // exit conditions on residuals
        if ((mem->nlp_res->inf_norm_res_g < opts->min_res_g) &
            (mem->nlp_res->inf_norm_res_b < opts->min_res_b) &
            (mem->nlp_res->inf_norm_res_d < opts->min_res_d) &
            (mem->nlp_res->inf_norm_res_m < opts->min_res_m))
        {
            // printf("%d sqp iterations\n", sqp_iter);
            // print_ocp_qp_in(work->qp_in);

            // save sqp iterations number
            mem->sqp_iter = sqp_iter;
            nlp_out->sqp_iter = sqp_iter;

            // stop timer
            total_time += acados_toc(&timer0);
            nlp_out->total_time = total_time;
            mem->time_tot = total_time;

#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
            return 0;
        }

        // printf("\n------- qp_in (sqp iter %d) --------\n", sqp_iter);
        // print_ocp_qp_in(work->qp_in);
        // exit(1);

        // start timer
        acados_tic(&timer1);

        int qp_status =
            qp_solver->evaluate(qp_solver, work->qp_in, work->qp_out, opts->qp_solver_opts,
                                mem->qp_solver_mem, work->qp_work);

        // stop timer
        mem->time_qp_sol += acados_toc(&timer1);

        nlp_out->qp_iter = ((ocp_qp_info *) work->qp_out->misc)->num_iter;

        // printf("\n------- qp_out (sqp iter %d) ---------\n", sqp_iter);
        // print_ocp_qp_out(work->qp_out);
        //  if(sqp_iter==1)
        //  exit(1);

        if (qp_status != 0)
        {
            //   print_ocp_qp_in(work->qp_in);

            // stop timer
            total_time += acados_toc(&timer0);

            mem->time_tot = total_time;
            nlp_out->total_time = total_time;

            printf("QP solver returned error status %d in iteration %d\n", qp_status, sqp_iter);
#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
            return -1;
        }

        sqp_update_variables(config, dims, nlp_out, opts, mem, work);

        // ocp_nlp_dims_print(nlp_out->dims);
        // ocp_nlp_out_print(nlp_out);
        // exit(1);

        // ??? @rien
        //        for (int_t i = 0; i < N; i++)
        //        {
        //   ocp_nlp_dynamics_opts *dynamics_opts = opts->dynamics[i];
        //            sim_rk_opts *rk_opts = dynamics_opts->sim_solver;
        //            if (rk_opts->scheme == NULL)
        //                continue;
        //            rk_opts->sens_adj = (rk_opts->scheme->type != exact);
        //            if (nlp_in->freezeSens) {
        //                // freeze inexact sensitivities after first SQP iteration !!
        //                rk_opts->scheme->freeze = true;
        //            }
        //        }
    }

    // stop timer
    total_time += acados_toc(&timer0);

    // ocp_nlp_out_print(nlp_out);

    // save sqp iterations number
    mem->sqp_iter = sqp_iter;
    mem->time_tot = total_time;

    nlp_out->sqp_iter = sqp_iter;
    nlp_out->total_time = total_time;

    // printf("%d sqp iterations\n", sqp_iter);
    // print_ocp_qp_in(work->qp_in);

    // maximum number of iterations reached
#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
    return 1;
}



int ocp_nlp_sqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_opts *opts = opts_;
    ocp_nlp_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    // ocp_nlp_out *nlp_out = nlp_out_;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_work *work = work_;

    ocp_nlp_sqp_cast_workspace(config, dims, work, mem, opts);

    // extract dims
    int N = dims->N;
    int status = ACADOS_SUCCESS;

    int ii;

    // TODO(fuck_lint) checks
    // TODO(fuck_lint) flag to enable/disable checks
    for (ii = 0; ii <= N; ii++)
    {
        // TODO(fuck_lint) check that ns in opt_var == ns in constraints
    }

    // precompute
    for (ii = 0; ii < N; ii++)
    {
        // set T
        config->dynamics[ii]->model_set_T(nlp_in->Ts[ii], nlp_in->dynamics[ii]);
        // dynamics precompute
        status = config->dynamics[ii]->precompute(config->dynamics[ii], dims->dynamics[ii],
                                            nlp_in->dynamics[ii], opts->dynamics[ii],
                                            mem->dynamics[ii], work->dynamics[ii]);
        if (status != ACADOS_SUCCESS) return status;
    }

    if (config->dynamics[N] != NULL)
    {
        // set T
        config->dynamics[N]->model_set_T(nlp_in->Ts[N], nlp_in->dynamics[N]);
        // dynamics precompute
        status = config->dynamics[N]->precompute(config->dynamics[N], dims->dynamics[N],
                                            nlp_in->dynamics[N], opts->dynamics[N],
                                            mem->dynamics[N], work->dynamics[N]);
        if (status != ACADOS_SUCCESS) return status;
    }
    return status;
}



void ocp_nlp_sqp_get(void *config_, void *mem_, const char *field, void *return_value_)
{
    // ocp_nlp_solver_config *config = config_;
    ocp_nlp_sqp_memory *mem = mem_;

    if (!strcmp("sqp_iter", field))
    {
        int *value = return_value_;
        *value = mem->sqp_iter;
    }
    else if (!strcmp("time_tot", field) || !strcmp("tot_time", field))
    {
        double *value = return_value_;
        *value = mem->time_tot;
    }
    else if (!strcmp("time_qp_sol", field) || !strcmp("time_qp", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_sol;
    }
    else if (!strcmp("time_lin", field))
    {
        double *value = return_value_;
        *value = mem->time_lin;
    }
    else if (!strcmp("nlp_res", field))
    {
        ocp_nlp_res **value = return_value_;
        *value = mem->nlp_res;
    }
    else
    {
        printf("\nerror: output type %s not available in ocp_nlp_sqp module\n", field);
        exit(1);
    }

}



void ocp_nlp_sqp_config_initialize_default(void *config_)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    config->opts_calculate_size = &ocp_nlp_sqp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_sqp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_sqp_opts_initialize_default;
    config->opts_update = &ocp_nlp_sqp_opts_update;
    config->opts_set = &ocp_nlp_sqp_opts_set;
    config->dynamics_opts_set = &ocp_nlp_sqp_dynamics_opts_set;
    config->memory_calculate_size = &ocp_nlp_sqp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_sqp_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_sqp_workspace_calculate_size;
    config->evaluate = &ocp_nlp_sqp;
    config->config_initialize_default = &ocp_nlp_sqp_config_initialize_default;
    config->regularization = NULL;
    config->precompute = &ocp_nlp_sqp_precompute;
    config->get = &ocp_nlp_sqp_get;

    return;
}
