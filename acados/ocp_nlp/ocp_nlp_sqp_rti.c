/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


#include "acados/ocp_nlp/ocp_nlp_sqp_rti.h"

// external
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
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
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



/************************************************
 * options
 ************************************************/

int ocp_nlp_sqp_rti_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_rti_opts);

    size += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    size += config->regularize->opts_calculate_size();

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



void *ocp_nlp_sqp_rti_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_sqp_rti_opts *opts = (ocp_nlp_sqp_rti_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_rti_opts);

    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    opts->regularize = config->regularize->opts_assign(c_ptr);
    c_ptr += config->regularize->opts_calculate_size();

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

    assert((char *) raw_memory + ocp_nlp_sqp_rti_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_sqp_rti_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;
    ocp_nlp_reg_config *regularize = config->regularize;

    int ii;

    int N = dims->N;

    // SQP RTI opts

//    opts->compute_dual_sol = 1;

    opts->reuse_workspace = 1;
#if defined(ACADOS_WITH_OPENMP)
    opts->num_threads = ACADOS_NUM_THREADS;
#endif

	opts->ext_qp_res = 0;

    // submodules opts

    // do not compute adjoint in dynamics and constraints
    int compute_adj = 0;

    // qp solver
    qp_solver->opts_initialize_default(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // regularization
    regularize->opts_initialize_default(regularize, dims->regularize, opts->regularize);

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        dynamics[ii]->opts_initialize_default(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
        dynamics[ii]->opts_set(dynamics[ii], opts->dynamics[ii], "compute_adj", &compute_adj);
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
        constraints[ii]->opts_set(constraints[ii], opts->constraints[ii], "compute_adj", &compute_adj);
    }

    return;
}



void ocp_nlp_sqp_rti_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;

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



void ocp_nlp_sqp_rti_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_sqp_rti_opts *opts = (ocp_nlp_sqp_rti_opts *) opts_;
    ocp_nlp_config *config = config_;

	int ii;

	char module[MAX_STR_LEN];
	char *ptr_module = NULL;
	int module_length = 0;

	// extract module name
	char *char_ = strchr(field, '_');
	if(char_!=NULL)
	{
		module_length = char_-field;
		for(ii=0; ii<module_length; ii++)
			module[ii] = field[ii];
		module[module_length] = '\0'; // add end of string
		ptr_module = module;
	}

	// pass options to QP module
	if(!strcmp(ptr_module, "qp"))
	{
		config->qp_solver->opts_set(config->qp_solver, opts->qp_solver_opts, field+module_length+1, value);

		if(!strcmp(field, "qp_warm_start"))
		{
			int* i_ptr = (int *) value;
			opts->qp_warm_start = *i_ptr;
		}
	}
	else // nlp opts
	{
		if (!strcmp(field, "num_threads"))
		{
			int* num_threads = (int *) value;
			opts->num_threads = *num_threads;
		}
		else if (!strcmp(field, "exact_hess"))
		{
			int N = config->N;
			// cost
			for (ii=0; ii<=N; ii++)
				config->cost[ii]->opts_set(config->cost[ii], opts->cost[ii], "exact_hess", value);
			// dynamics
			for (ii=0; ii<N; ii++)
				config->dynamics[ii]->opts_set(config->dynamics[ii], opts->dynamics[ii], "compute_hess", value);
//			// constraints TODO disabled for now as prevents convergence !!!
//			for (ii=0; ii<=N; ii++)
//				config->constraints[ii]->opts_set(config->constraints[ii], opts->constraints[ii], "compute_hess", value);
		}
		else if (!strcmp(field, "ext_qp_res"))
		{
			int* ext_qp_res = (int *) value;
			opts->ext_qp_res = *ext_qp_res;
		}
		else
		{
			printf("\nerror: ocp_nlp_sqp_rti_opts_set: wrong field: %s\n", field);
			exit(1);
		}
	}

	return;

}



void ocp_nlp_sqp_rti_dynamics_opts_set(void *config_, void *opts_, int stage,
        const char *field, void *value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_dynamics_config *dyn_config = config->dynamics[stage];

    dyn_config->opts_set(dyn_config, opts->dynamics[stage], field, value);

    return;

}



void ocp_nlp_sqp_rti_cost_opts_set(void *config_, void *opts_, int stage,
        const char *field, void *value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_cost_config *cost_config = config->cost[stage];

    cost_config->opts_set(cost_config, opts->cost[stage], field, value);

    return;

}



void ocp_nlp_sqp_rti_constraints_opts_set(void *config_, void *opts_, int stage,
        const char *field, void *value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_constraints_config *constraints_config = config->constraints[stage];

    constraints_config->opts_set(constraints_config, opts->constraints[stage], (char *) field, value);

    return;

}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_sqp_rti_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // extract dims
    int N = dims->N;
    // ocp_nlp_cost_dims **cost_dims = dims->cost;
    // int ny;

    int size = 0;

    size += sizeof(ocp_nlp_sqp_rti_memory);

    size += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    size += config->regularize->memory_calculate_size(config->regularize, dims->regularize, opts->regularize);

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

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims);

	// stat
	int stat_m = 1+1;
	int stat_n = 2;
	if(opts->ext_qp_res)
		stat_n += 4;
	size += stat_n*stat_m*sizeof(double);

    size += 8;  // initial align

    //    make_int_multiple_of(64, &size);

    return size;
}



void *ocp_nlp_sqp_rti_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;

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

    ocp_nlp_sqp_rti_memory *mem = (ocp_nlp_sqp_rti_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_rti_memory);

    // QP solver
    mem->qp_solver_mem =
        qp_solver->memory_assign(qp_solver, dims->qp_solver, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // regularization
    mem->regularize_mem = config->regularize->memory_assign(config->regularize, dims->regularize, opts->regularize, c_ptr);
    c_ptr += config->regularize->memory_calculate_size(config->regularize, dims->regularize, opts->regularize);

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

	// stat
	mem->stat = (double *) c_ptr;
	mem->stat_m = 1+1;
	mem->stat_n = 2;
	if(opts->ext_qp_res)
		mem->stat_n += 4;
	c_ptr += mem->stat_m*mem->stat_n*sizeof(double);

    mem->status = ACADOS_READY;

    assert((char *) raw_memory+ocp_nlp_sqp_rti_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

int ocp_nlp_sqp_rti_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // loop index
    int ii;

    // extract dims
    int N = dims->N;
	int *nx = dims->nx;
	int *nu = dims->nu;
	int *nz = dims->nz;

    int size = 0;
    int size_tmp = 0;
    int tmp;

    // sqp
    size += sizeof(ocp_nlp_sqp_rti_work);

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

	if(opts->ext_qp_res)
	{
		// qp res
		size += ocp_qp_res_calculate_size(dims->qp_solver);

		// qp res ws
		size += ocp_qp_res_workspace_calculate_size(dims->qp_solver);
	}

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

	// dzduxt
	size += (N+1)*sizeof(struct blasfeo_dmat);
	for(ii=0; ii<=N; ii++)
		size += blasfeo_memsize_dmat(nu[ii]+nx[ii], nz[ii]);
	// z_alg
	size += (N+1)*sizeof(struct blasfeo_dvec);
	for(ii=0; ii<=N; ii++)
		size += blasfeo_memsize_dvec(nz[ii]);

    size += 1*8;  // blasfeo_str align
    size += 1*64;  // blasfeo_mem align

    return size;
}



// TODO(all): introduce member "memsize" in all structures to make on-line cast cheaper (i.e. avoid
// to calculate size on-line)
static void ocp_nlp_sqp_rti_cast_workspace(void *config_, ocp_nlp_dims *dims,
                                           ocp_nlp_sqp_rti_work *work,
                                           ocp_nlp_sqp_rti_memory *mem, ocp_nlp_sqp_rti_opts *opts)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // loop index
	int ii;

    // extract dims
    int N = dims->N;
	int *nx = dims->nx;
	int *nu = dims->nu;
	int *nz = dims->nz;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_sqp_rti_work);

    // array of pointers
    //
    work->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);
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

	if(opts->ext_qp_res)
	{
		// qp res
		work->qp_res = ocp_qp_res_assign(dims->qp_solver, c_ptr);
		c_ptr += ocp_qp_res_calculate_size(dims->qp_solver);

		// qp res ws
		work->qp_res_ws = ocp_qp_res_workspace_assign(dims->qp_solver, c_ptr);
		c_ptr += ocp_qp_res_workspace_calculate_size(dims->qp_solver);
	}

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

		int size_tmp = 0;
		int tmp;

        // qp solver
        work->qp_work = (void *) c_ptr;
        tmp = qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);
        size_tmp = tmp > size_tmp ? tmp : size_tmp;

        // dynamics
        for (int ii = 0; ii < N; ii++)
        {
            work->dynamics[ii] = c_ptr;
            tmp = dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // cost
        for (int ii = 0; ii <= N; ii++)
        {
            work->cost[ii] = c_ptr;
            tmp = cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // constraints
        for (int ii = 0; ii <= N; ii++)
        {
            work->constraints[ii] = c_ptr;
            tmp = constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        c_ptr += size_tmp;

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

    // blasfeo_str align
    align_char_to(8, &c_ptr);

	// dzduxt
	work->dzduxt = (struct blasfeo_dmat *) c_ptr;
	c_ptr += (N+1)*sizeof(struct blasfeo_dmat);
	// z_alg
	work->z_alg = (struct blasfeo_dvec *) c_ptr;
	c_ptr += (N+1)*sizeof(struct blasfeo_dvec);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

	// dzduxt
	for(ii=0; ii<=N; ii++)
		{
		blasfeo_create_dmat(nu[ii]+nx[ii], nz[ii], work->dzduxt+ii, c_ptr);
		c_ptr += blasfeo_memsize_dmat(nu[ii]+nx[ii], nz[ii]);
		}
	// z_alg
	for(ii=0; ii<=N; ii++)
		{
		blasfeo_create_dvec(nz[ii], work->z_alg+ii, c_ptr);
		c_ptr += blasfeo_memsize_dvec(nz[ii]);
		}

    // assert & return
    assert((char *) work + ocp_nlp_sqp_rti_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

static void initialize_qp(void *config_, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                          ocp_nlp_out *nlp_out, ocp_nlp_sqp_rti_opts *opts,
                          ocp_nlp_sqp_rti_memory *mem, ocp_nlp_sqp_rti_work *work)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

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
        // constraints
        config->constraints[ii]->initialize(config->constraints[ii], dims->constraints[ii],
                                            nlp_in->constraints[ii], opts->constraints[ii],
                                            mem->constraints[ii], work->constraints[ii]);
    }

    return;
}



static void linearize_update_qp_matrices(void *config_, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                                         ocp_nlp_out *nlp_out, ocp_nlp_sqp_rti_opts *opts,
                                         ocp_nlp_sqp_rti_memory *mem, ocp_nlp_sqp_rti_work *work)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

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
        // init Hessian to 0 
        blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, work->qp_in->RSQrq+i, 0, 0);
        // dynamics
        if (i < N)
            config->dynamics[i]->update_qp_matrices(config->dynamics[i], dims->dynamics[i],
                    nlp_in->dynamics[i], opts->dynamics[i],
                    mem->dynamics[i], work->dynamics[i]);
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

        // nlp mem: dyn_adj
        if (i < N)
        {
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[i]->memory_get_adj_ptr(mem->dynamics[i]);
            blasfeo_dveccp(nu[i] + nx[i], dyn_adj, 0, nlp_mem->dyn_adj + i, 0);
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
        //   sim_opts *opts = dynamics_opts->sim_solver;
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
                                  ocp_nlp_out *nlp_out, ocp_nlp_sqp_rti_opts *opts,
                                  ocp_nlp_sqp_rti_memory *mem, ocp_nlp_sqp_rti_work *work)
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



static void sqp_update_variables(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out,
                                 ocp_nlp_sqp_rti_opts *opts, ocp_nlp_sqp_rti_memory *mem,
                                 ocp_nlp_sqp_rti_work *work)
{
    // loop index
    int i;

    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;

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

        if (i < N)
			blasfeo_dveccp(nz[i], work->z_alg+i, 0, nlp_out->z+i, 0);

    }

    return;
}



// Simple fixed-step Gauss-Newton based SQP routine
int ocp_nlp_sqp_rti(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
 

    // acados timer
    acados_timer timer0, timer1;

    // start timer
    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_sqp_rti_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_rti_work *work = work_;

    ocp_nlp_sqp_rti_cast_workspace(config, dims, work, mem, opts);

    // zero timers
    double total_time = 0.0;
    mem->time_qp_sol = 0.0;
    mem->time_lin = 0.0;
    mem->time_reg = 0.0;
    mem->time_tot = 0.0;

    // extract dims
    int N = dims->N;

    int ii;

	int qp_iter = 0;
	int qp_status = 0;

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
    #pragma omp for nowait
#endif
    for (ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->memory_set_ux_ptr(nlp_out->ux+ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_ux1_ptr(nlp_out->ux+ii+1, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_pi_ptr(nlp_out->pi+ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_BAbt_ptr(work->qp_in->BAbt+ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_RSQrq_ptr(work->qp_in->RSQrq+ii, mem->dynamics[ii]);
//        config->dynamics[ii]->memory_set_z_alg_ptr(nlp_out->z+ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_dzduxt_ptr(work->dzduxt+ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_z_guess_ptr(nlp_out->z+ii, mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_z_alg_ptr(work->z_alg+ii, mem->dynamics[ii]);
    }

    // alias to cost_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (ii = 0; ii <= N; ii++)
    {
        config->cost[ii]->memory_set_ux_ptr(nlp_out->ux + ii, mem->cost[ii]);
		config->cost[ii]->memory_set_z_alg_ptr(work->z_alg+ii, mem->cost[ii]);
		config->cost[ii]->memory_set_dzdux_tran_ptr(work->dzduxt+ii, mem->cost[ii]);
        config->cost[ii]->memory_set_RSQrq_ptr(work->qp_in->RSQrq + ii, mem->cost[ii]);
        config->cost[ii]->memory_set_Z_ptr(work->qp_in->Z + ii, mem->cost[ii]);
    }

    // alias to constraints_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (ii = 0; ii <= N; ii++)
    {
        config->constraints[ii]->memory_set_ux_ptr(nlp_out->ux+ii, mem->constraints[ii]);
        config->constraints[ii]->memory_set_lam_ptr(nlp_out->lam+ii, mem->constraints[ii]);
        config->constraints[ii]->memory_set_DCt_ptr(work->qp_in->DCt+ii, mem->constraints[ii]);
        config->constraints[ii]->memory_set_RSQrq_ptr(work->qp_in->RSQrq+ii, mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxb_ptr(work->qp_in->idxb[ii], mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxs_ptr(work->qp_in->idxs[ii], mem->constraints[ii]);
    }

    // alias to regularize memory
    config->regularize->memory_set_RSQrq_ptr(dims->regularize, work->qp_in->RSQrq, mem->regularize_mem);
    config->regularize->memory_set_rq_ptr(dims->regularize, work->qp_in->rqz, mem->regularize_mem);
    config->regularize->memory_set_BAbt_ptr(dims->regularize, work->qp_in->BAbt, mem->regularize_mem);
    config->regularize->memory_set_b_ptr(dims->regularize, work->qp_in->b, mem->regularize_mem);
	config->regularize->memory_set_idxb_ptr(dims->regularize, work->qp_in->idxb, mem->regularize_mem);
    config->regularize->memory_set_DCt_ptr(dims->regularize, work->qp_in->DCt, mem->regularize_mem);
    config->regularize->memory_set_ux_ptr(dims->regularize, work->qp_out->ux, mem->regularize_mem);
    config->regularize->memory_set_pi_ptr(dims->regularize, work->qp_out->pi, mem->regularize_mem);
    config->regularize->memory_set_lam_ptr(dims->regularize, work->qp_out->lam, mem->regularize_mem);

    // copy sampling times into dynamics model
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (int ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->model_set(config->dynamics[ii], dims->dynamics[ii], nlp_in->dynamics[ii], "T", nlp_in->Ts+ii);
    }

#if defined(ACADOS_WITH_OPENMP)
    } // end of parallel region
#endif

    // initialize QP
    initialize_qp(config, dims, nlp_in, nlp_out, opts, mem, work);

        

    // SQP body

    // start timer
    acados_tic(&timer1);

    // linearizate NLP and update QP matrices
    linearize_update_qp_matrices(config, dims, nlp_in, nlp_out, opts, mem, work);

    // stop timer
    mem->time_lin += acados_toc(&timer1);

    // update QP rhs for SQP (step prim var, abs dual var)
    sqp_update_qp_vectors(config, dims, nlp_in, nlp_out, opts, mem, work);

	// save statistics
//	mem->stat[mem->stat_n*1+0] = qp_status;
//	mem->stat[mem->stat_n*1+1] = qp_iter;

	// start timer
	acados_tic(&timer1);
    // regularize Hessian
    config->regularize->regularize_hessian(config->regularize, dims->regularize, opts->regularize, mem->regularize_mem);
	// stop timer
	mem->time_reg += acados_toc(&timer1);

    // printf("\n------- qp_in (sqp iter %d) --------\n", sqp_iter);
    // print_ocp_qp_in(work->qp_in);
    // exit(1);

	// TODO no warm start across NLP solutions (yet)
	int tmp_int = 0;
	config->qp_solver->opts_set(config->qp_solver, opts->qp_solver_opts, "warm_start", &tmp_int);

    // start timer
    acados_tic(&timer1);
	// TODO move qp_out in memory !!!!! (it has to be preserved to do warm start)
    qp_status = qp_solver->evaluate(qp_solver, work->qp_in, work->qp_out, opts->qp_solver_opts, mem->qp_solver_mem, work->qp_work);
    // stop timer
    mem->time_qp_sol += acados_toc(&timer1);

	// start timer
	acados_tic(&timer1);
    // compute correct dual solution in case of Hessian regularization
    config->regularize->correct_dual_sol(config->regularize, dims->regularize, opts->regularize, mem->regularize_mem);
	// stop timer
	mem->time_reg += acados_toc(&timer1);

	// TODO move into QP solver memory ???
	nlp_out->qp_iter = ((ocp_qp_info *) work->qp_out->misc)->num_iter;
	qp_iter = ((ocp_qp_info *) work->qp_out->misc)->num_iter;

	// compute external QP residuals (for debugging)
	if(opts->ext_qp_res)
	{
		ocp_qp_res_compute(work->qp_in, work->qp_out, work->qp_res, work->qp_res_ws);
		ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*1+2));
//		printf("\nsqp_iter %d, res %e %e %e %e\n", sqp_iter, inf_norm_qp_res[0], inf_norm_qp_res[1], inf_norm_qp_res[2], inf_norm_qp_res[3]);
	}

    // printf("\n------- qp_out (sqp iter %d) ---------\n", sqp_iter);
    //  print_ocp_qp_out(work->qp_out);
    //  if(sqp_iter==1)
    //  exit(1);

	// save statistics
	mem->stat[mem->stat_n*1+0] = qp_status;
	mem->stat[mem->stat_n*1+1] = qp_iter;

	if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
    {
        //   print_ocp_qp_in(work->qp_in);

        // stop timer
        total_time += acados_toc(&timer0);

        mem->time_tot = total_time;
        nlp_out->total_time = total_time;

        printf("QP solver returned error status %d\n", qp_status);
#if defined(ACADOS_WITH_OPENMP)
        // restore number of threads
        omp_set_num_threads(num_threads_bkp);
#endif
        mem->status = ACADOS_QP_FAILURE;
        return mem->status;
    }

    sqp_update_variables(dims, nlp_out, opts, mem, work);

    // ocp_nlp_dims_print(nlp_out->dims);
    // ocp_nlp_out_print(nlp_out);
    // exit(1);

    // stop timer
    total_time += acados_toc(&timer0);

    mem->time_tot = total_time;
    nlp_out->total_time = total_time;

    // ocp_nlp_out_print(nlp_out);

    // print_ocp_qp_in(work->qp_in);

#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
    mem->status = ACADOS_SUCCESS;
    return mem->status;
}



int ocp_nlp_sqp_rti_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_opts *opts = opts_;
    ocp_nlp_sqp_rti_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    // ocp_nlp_out *nlp_out = nlp_out_;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_sqp_rti_work *work = work_;

    ocp_nlp_sqp_rti_cast_workspace(config, dims, work, mem, opts);

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
        config->dynamics[ii]->model_set(config->dynamics[ii], dims->dynamics[ii], nlp_in->dynamics[ii], "T", nlp_in->Ts+ii);
        // dynamics precompute
        status = config->dynamics[ii]->precompute(config->dynamics[ii], dims->dynamics[ii],
                                            nlp_in->dynamics[ii], opts->dynamics[ii],
                                            mem->dynamics[ii], work->dynamics[ii]);
        if (status != ACADOS_SUCCESS) return status;
    }

    return status;
}



void ocp_nlp_sqp_rti_get(void *config_, void *mem_, const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_sqp_rti_memory *mem = mem_;

    if (!strcmp("sqp_iter", field))
    {
        int *value = return_value_;
        *value = 1;
    }
    else if (!strcmp("status", field))
    {
        int *value = return_value_;
        *value = mem->status;
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
    else if (!strcmp("time_reg", field))
    {
        double *value = return_value_;
        *value = mem->time_reg;
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("stat_m", field))
    {
        int *value = return_value_;
        *value = mem->stat_m;
    }
    else if (!strcmp("stat_n", field))
    {
        int *value = return_value_;
        *value = mem->stat_n;
    }
    else
    {
        printf("\nerror: output type %s not available in ocp_nlp_sqp_rti module\n", field);
        exit(1);
    }
}


void ocp_nlp_sqp_rti_config_initialize_default(void *config_)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    config->opts_calculate_size = &ocp_nlp_sqp_rti_opts_calculate_size;
    config->opts_assign = &ocp_nlp_sqp_rti_opts_assign;
    config->opts_initialize_default = &ocp_nlp_sqp_rti_opts_initialize_default;
    config->opts_update = &ocp_nlp_sqp_rti_opts_update;
    config->opts_set = &ocp_nlp_sqp_rti_opts_set;
    config->dynamics_opts_set = &ocp_nlp_sqp_rti_dynamics_opts_set;
    config->cost_opts_set = &ocp_nlp_sqp_rti_cost_opts_set;
    config->constraints_opts_set = &ocp_nlp_sqp_rti_constraints_opts_set;
    config->memory_calculate_size = &ocp_nlp_sqp_rti_memory_calculate_size;
    config->memory_assign = &ocp_nlp_sqp_rti_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_sqp_rti_workspace_calculate_size;
    config->evaluate = &ocp_nlp_sqp_rti;
    config->config_initialize_default = &ocp_nlp_sqp_rti_config_initialize_default;
    config->precompute = &ocp_nlp_sqp_rti_precompute;
    config->get = &ocp_nlp_sqp_rti_get;

    return;
}

