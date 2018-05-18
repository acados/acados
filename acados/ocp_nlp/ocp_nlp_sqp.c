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
#include <stdlib.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/sim/sim_collocation_utils.h"  // TODO(all): remove ???
#include "acados/sim/sim_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

// static int get_max_sim_workspace_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims,
// ocp_nlp_sqp_opts *opts)
// {
//  /* ocp_qp_xcond_solver_config *qp_solver = config->qp_solver; */
//  ocp_nlp_dynamics_config **dynamics = config->dynamics;

//     int sim_work_size;

//     int max_sim_work_size = 0;

//     for (int ii = 0; ii < dims->N; ii++)
//     {
//         // sim_in_size = sim_in_calculate_size(dims->sim[ii]);
//         // if (sim_in_size > *max_sim_in_size) *max_sim_in_size = sim_in_size;
//         // sim_out_size = sim_out_calculate_size(dims->sim[ii]);
//         // if (sim_out_size > *max_sim_out_size) *max_sim_out_size = sim_out_size;
//   ocp_nlp_dynamics_opts *dynamics_opts = opts->dynamics[ii];
//         sim_work_size =
//         dynamics[ii]->sim_solver->workspace_calculate_size(dynamics[ii]->sim_solver,
//         dims->dynamics[ii]->sim, dynamics_opts->sim_solver); if (sim_work_size >
//         max_sim_work_size) max_sim_work_size = sim_work_size;
//     }
//     return max_sim_work_size;
// }



/************************************************
 * options
 ************************************************/

int ocp_nlp_sqp_opts_calculate_size(void *config_, ocp_nlp_dims *dims)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_opts);

    size += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    // dynamics
    size += N * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        size += dynamics[ii]->opts_calculate_size(dynamics[ii], dims->dynamics[ii]);
    }

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



void *ocp_nlp_sqp_opts_assign(void *config_, ocp_nlp_dims *dims, void *raw_memory)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

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

    // dynamics
    opts->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        opts->dynamics[ii] = dynamics[ii]->opts_assign(dynamics[ii], dims->dynamics[ii], c_ptr);
        c_ptr += dynamics[ii]->opts_calculate_size(dynamics[ii], dims->dynamics[ii]);
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



void ocp_nlp_sqp_opts_initialize_default(void *config_, ocp_nlp_dims *dims, void *opts_)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int ii;

    int N = dims->N;

    opts->maxIter = 20;
    opts->min_res_g = 1e-12;
    opts->min_res_b = 1e-12;
    opts->min_res_d = 1e-12;
    opts->min_res_m = 1e-12;

    qp_solver->opts_initialize_default(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        dynamics[ii]->opts_initialize_default(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
    }

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



void ocp_nlp_sqp_opts_update(void *config_, ocp_nlp_dims *dims, void *opts_)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;

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



/************************************************
 * memory
 ************************************************/

int ocp_nlp_sqp_memory_calculate_size(void *config_, ocp_nlp_dims *dims, void *opts_)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;

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

    // dynamics
    size += N * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        size += dynamics[ii]->memory_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                    opts->dynamics[ii]);
    }

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



void *ocp_nlp_sqp_memory_assign(void *config_, ocp_nlp_dims *dims, void *opts_, void *raw_memory)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;

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

    // nlp res
    mem->nlp_res = ocp_nlp_res_assign(dims, c_ptr);
    c_ptr += mem->nlp_res->memsize;

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims);

    // dynamics
    mem->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        mem->dynamics[ii] = dynamics[ii]->memory_assign(dynamics[ii], dims->dynamics[ii],
                                                        opts->dynamics[ii], c_ptr);
        c_ptr += dynamics[ii]->memory_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                     opts->dynamics[ii]);
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

int ocp_nlp_sqp_workspace_calculate_size(void *config_, ocp_nlp_dims *dims, void *opts_)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // loop index
    int ii;

    // extract dims
    int N = dims->N;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_work);

    // qp in
    size += ocp_qp_in_calculate_size(qp_solver, dims->qp_solver);

    // qp out
    size += ocp_qp_out_calculate_size(qp_solver, dims->qp_solver);

    // qp solver
    size += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // dynamics
    size += N * sizeof(void *);
    for (ii = 0; ii < N; ii++)
    {
        size += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                       opts->dynamics[ii]);
    }

    // cost
    size += (N + 1) * sizeof(void *);
    for (ii = 0; ii <= N; ii++)
    {
        size += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
    }

    // constraints
    size += (N + 1) * sizeof(void *);
    for (ii = 0; ii <= N; ii++)
    {
        size += constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii],
                                                          opts->constraints[ii]);
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

    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_work);

    // extract dims
    int N = dims->N;

    // qp in
    work->qp_in = ocp_qp_in_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(qp_solver, dims->qp_solver);

    // qp out
    work->qp_out = ocp_qp_out_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(qp_solver, dims->qp_solver);

    // qp solver
    work->qp_work = (void *) c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // dynamics
    work->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        work->dynamics[ii] = c_ptr;
        c_ptr += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii],
                                                        opts->dynamics[ii]);
    }

    // cost
    work->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        work->cost[ii] = c_ptr;
        c_ptr += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
    }

    // constraints
    work->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);
    for (int ii = 0; ii <= N; ii++)
    {
        work->constraints[ii] = c_ptr;
        c_ptr += constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii],
                                                           opts->constraints[ii]);
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

    for (ii = 0; ii < N; ii++)
    {
        // cost
        config->cost[ii]->initialize(config->cost[ii], dims->cost[ii], nlp_in->cost[ii],
                                     opts->cost[ii], mem->cost[ii], work->cost[ii]);
        // dynamics
        config->dynamics[ii]->initialize(config->dynamics[ii], dims->dynamics[ii],
                                         nlp_in->dynamics[ii], opts->dynamics[ii],
                                         mem->dynamics[ii], work->dynamics[ii]);
        // constraints
        config->constraints[ii]->initialize(config->constraints[ii], dims->constraints[ii],
                                            nlp_in->constraints[ii], opts->constraints[ii],
                                            mem->constraints[ii], work->constraints[ii]);
    }
    ii = N;
    // cost
    config->cost[ii]->initialize(config->cost[ii], dims->cost[ii], nlp_in->cost[ii], opts->cost[ii],
                                 mem->cost[ii], work->cost[ii]);
    // constraints
    config->constraints[ii]->initialize(config->constraints[ii], dims->constraints[ii],
                                        nlp_in->constraints[ii], opts->constraints[ii],
                                        mem->constraints[ii], work->constraints[ii]);

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

    for (i = 0; i < N; i++)
    {
        // cost
        config->cost[i]->update_qp_matrices(config->cost[i], dims->cost[i], nlp_in->cost[i],
                                            opts->cost[i], mem->cost[i], work->cost[i]);
        // dynamics
        config->dynamics[i]->update_qp_matrices(config->dynamics[i], dims->dynamics[i],
                                                nlp_in->dynamics[i], opts->dynamics[i],
                                                mem->dynamics[i], work->dynamics[i]);
        // constraints
        config->constraints[i]->update_qp_matrices(config->constraints[i], dims->constraints[i],
                                                   nlp_in->constraints[i], opts->constraints[i],
                                                   mem->constraints[i], work->constraints[i]);
    }
    i = N;
    // cost
    config->cost[i]->update_qp_matrices(config->cost[i], dims->cost[i], nlp_in->cost[i],
                                        opts->cost[i], mem->cost[i], work->cost[i]);
    // constraints
    config->constraints[i]->update_qp_matrices(config->constraints[i], dims->constraints[i],
                                               nlp_in->constraints[i], opts->constraints[i],
                                               mem->constraints[i], work->constraints[i]);

    /* collect stage-wise evaluations */

    // nlp mem: cost_grad
    for (i = 0; i <= N; i++)
    {
        struct blasfeo_dvec *cost_grad = config->cost[i]->memory_get_grad_ptr(mem->cost[i]);
        blasfeo_dveccp(nv[i], cost_grad, 0, nlp_mem->cost_grad + i, 0);
    }

    // nlp mem: dyn_fun
    for (i = 0; i < N; i++)
    {
        struct blasfeo_dvec *dyn_fun = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        blasfeo_dveccp(nx[i + 1], dyn_fun, 0, nlp_mem->dyn_fun + i, 0);
    }

    // nlp mem: dyn_adj
    for (i = 0; i < N; i++)
    {
        struct blasfeo_dvec *dyn_adj = config->dynamics[i]->memory_get_adj_ptr(mem->dynamics[i]);
        blasfeo_dveccp(nu[i] + nx[i], dyn_adj, 0, nlp_mem->dyn_adj + i, 0);
    }

    blasfeo_dvecse(nu[N] + nx[N], 0.0, nlp_mem->dyn_adj + N, 0);

    for (i = 0; i < N; i++)
    {
        struct blasfeo_dvec *dyn_adj = config->dynamics[i]->memory_get_adj_ptr(mem->dynamics[i]);
        blasfeo_daxpy(nx[i + 1], 1.0, dyn_adj, nu[i] + nx[i], nlp_mem->dyn_adj + i + 1, nu[i + 1],
                      nlp_mem->dyn_adj + i + 1, nu[i + 1]);
    }

    // nlp mem: ineq_fun
    for (i = 0; i <= N; i++)
    {
        struct blasfeo_dvec *ineq_fun =
            config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, nlp_mem->ineq_fun + i, 0);
    }

    // nlp mem: ineq_adj
    for (i = 0; i <= N; i++)
    {
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

    // g
    for (i = 0; i <= N; i++)
    {
        blasfeo_dveccp(nv[i], nlp_mem->cost_grad + i, 0, work->qp_in->rqz + i, 0);
    }

    // b
    for (i = 0; i < N; i++)
    {
        blasfeo_dveccp(nx[i + 1], nlp_mem->dyn_fun + i, 0, work->qp_in->b + i, 0);
    }

    // d
    for (i = 0; i <= N; i++)
    {
        blasfeo_dveccp(2 * ni[i], nlp_mem->ineq_fun + i, 0, work->qp_in->d + i, 0);
    }

    return;
}



static void sqp_update_variables(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out, ocp_nlp_sqp_opts *opts,
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

    // TODO(all): fix and move where appropriate
    //    for (i = 0; i < N; i++)
    //    {
    //  nx1 = dims->constraints[i+1]->nx;
    //        for (j = 0; j < nx1; j++)
    //        {
    //            work->sim_in[i]->S_adj[j] = -BLASFEO_DVECEL(&work->qp_out->pi[i], j);
    //        }
    //    }

    // (full) step in primal variables
    for (i = 0; i <= N; i++)
    {
        blasfeo_daxpy(nv[i], 1.0, work->qp_out->ux + i, 0, nlp_out->ux + i, 0, nlp_out->ux + i, 0);
    }

    // absolute in dual variables
    for (i = 0; i < N; i++)
    {
        blasfeo_dveccp(nx[i + 1], work->qp_out->pi + i, 0, nlp_out->pi + i, 0);
    }

    for (i = 0; i <= N; i++)
    {
        blasfeo_dveccp(2 * ni[i], work->qp_out->lam + i, 0, nlp_out->lam + i, 0);
    }

    for (i = 0; i <= N; i++)
    {
        blasfeo_dveccp(2 * ni[i], work->qp_out->t + i, 0, nlp_out->t + i, 0);
    }

    return;
}



// Simple fixed-step Gauss-Newton based SQP routine
int ocp_nlp_sqp(void *config_, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;
    ocp_nlp_sqp_opts *opts = (ocp_nlp_sqp_opts *) opts_;
    ocp_nlp_sqp_memory *mem = (ocp_nlp_sqp_memory *) mem_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_work *work = (ocp_nlp_sqp_work *) work_;

    ocp_nlp_sqp_cast_workspace(config, dims, work, mem, opts);

    int N = dims->N;

    // alias to dynamics_memory
    for (int ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->memory_set_ux_ptr(nlp_out->ux + ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_ux1_ptr(nlp_out->ux + ii + 1, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_pi_ptr(nlp_out->pi + ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_BAbt_ptr(work->qp_in->BAbt + ii, mem->dynamics[ii]);
    }

    // alias to cost_memory
    for (int ii = 0; ii <= N; ii++)
    {
        config->cost[ii]->memory_set_ux_ptr(nlp_out->ux + ii, mem->cost[ii]);
        config->cost[ii]->memory_set_RSQrq_ptr(work->qp_in->RSQrq + ii, mem->cost[ii]);
        config->cost[ii]->memory_set_Z_ptr(work->qp_in->Z + ii, mem->cost[ii]);
    }

    // alias to constraints_memory
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
    for (int ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->model_set_T(nlp_in->Ts[ii], nlp_in->dynamics[ii]);
    }

    // initialize QP
    initialize_qp(config, dims, nlp_in, nlp_out, opts, mem, work);

    // start timer
    acados_timer timer;
    double total_time = 0;
    acados_tic(&timer);

    // main sqp loop
    int max_sqp_iterations = opts->maxIter;
    int sqp_iter = 0;
    for (; sqp_iter < max_sqp_iterations; sqp_iter++)
    {
        // linearizate NLP and update QP matrices
        linearize_update_qp_matrices(config, dims, nlp_in, nlp_out, opts, mem, work);

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
            total_time += acados_toc(&timer);

            return 0;
        }

        // printf("\n------- qp_in (sqp iter %d) --------\n", sqp_iter);
        //  print_ocp_qp_in(work->qp_in);

        int qp_status =
            qp_solver->evaluate(qp_solver, work->qp_in, work->qp_out, opts->qp_solver_opts,
                                mem->qp_solver_mem, work->qp_work);

        // printf("\n------- qp_out (sqp iter %d) ---------\n", sqp_iter);
        //  print_ocp_qp_out(work->qp_out);
        //  if(sqp_iter==1)
        //  exit(1);

        if (qp_status != 0)
        {
            //   print_ocp_qp_in(work->qp_in);

            printf("QP solver returned error status %d in iteration %d\n", qp_status, sqp_iter);
            return -1;
        }

        sqp_update_variables(dims, nlp_out, opts, mem, work);

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
    total_time += acados_toc(&timer);

    // ocp_nlp_out_print(nlp_out);

    // save sqp iterations number
    mem->sqp_iter = sqp_iter;
    nlp_out->sqp_iter = sqp_iter;

    // printf("%d sqp iterations\n", sqp_iter);
    // print_ocp_qp_in(work->qp_in);

    // maximum number of iterations reached
    return 1;
}



void ocp_nlp_sqp_config_initialize_default(void *config_)
{
    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) config_;

    config->opts_calculate_size = &ocp_nlp_sqp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_sqp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_sqp_opts_initialize_default;
    config->opts_update = &ocp_nlp_sqp_opts_update;
    config->memory_calculate_size = &ocp_nlp_sqp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_sqp_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_sqp_workspace_calculate_size;
    config->evaluate = &ocp_nlp_sqp;
    config->config_initialize_default = &ocp_nlp_sqp_config_initialize_default;

    return;
}
