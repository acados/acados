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


#include "acados/ocp_nlp/ocp_nlp_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// blasfeo
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
// acados
#include "acados/utils/mem.h"



/************************************************
 * config
 ************************************************/

int ocp_nlp_config_calculate_size(int N)
{
    int ii;

    int size = 0;

    // self
    size += sizeof(ocp_nlp_config);

    // qp solver
    size += 1 * ocp_qp_xcond_solver_config_calculate_size();

    // regularization
    size += ocp_nlp_reg_config_calculate_size();


    // dynamics
    size += N * sizeof(ocp_nlp_dynamics_config *);

    for (ii = 0; ii < N; ii++) size += ocp_nlp_dynamics_config_calculate_size();

    // cost
    size += (N + 1) * sizeof(ocp_nlp_cost_config *);

    for (ii = 0; ii <= N; ii++) size += ocp_nlp_cost_config_calculate_size();

    // constraints
    size += (N + 1) * sizeof(ocp_nlp_constraints_config *);

    for (ii = 0; ii <= N; ii++) size += ocp_nlp_constraints_config_calculate_size();

    return size;
}



ocp_nlp_config *ocp_nlp_config_assign(int N, void *raw_memory)
{
    int ii;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_config *config = (ocp_nlp_config *) c_ptr;
    c_ptr += sizeof(ocp_nlp_config);

    config->N = N;

    // qp solver
    config->qp_solver = ocp_qp_xcond_solver_config_assign(c_ptr);
    c_ptr += ocp_qp_xcond_solver_config_calculate_size();

    // regularization
    config->regularize = ocp_nlp_reg_config_assign(c_ptr);
    c_ptr += ocp_nlp_reg_config_calculate_size();

    // dynamics
    config->dynamics = (ocp_nlp_dynamics_config **) c_ptr;
    c_ptr += N * sizeof(ocp_nlp_dynamics_config *);

    for (ii = 0; ii < N; ii++)
    {
        config->dynamics[ii] = ocp_nlp_dynamics_config_assign(c_ptr);
        c_ptr += ocp_nlp_dynamics_config_calculate_size();
    }

    // cost
    config->cost = (ocp_nlp_cost_config **) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_cost_config *);

    for (ii = 0; ii <= N; ii++)
    {
        config->cost[ii] = ocp_nlp_cost_config_assign(c_ptr);
        c_ptr += ocp_nlp_cost_config_calculate_size();
    }

    // constraints
    config->constraints = (ocp_nlp_constraints_config **) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_constraints_config *);

    for (ii = 0; ii <= N; ii++)
    {
        config->constraints[ii] = ocp_nlp_constraints_config_assign(c_ptr);
        c_ptr += ocp_nlp_constraints_config_calculate_size();
    }

    return config;
}



/************************************************
 * dims
 ************************************************/

static int ocp_nlp_dims_calculate_size_self(int N)
{
    int size = 0;

    size += sizeof(ocp_nlp_dims);

    // nlp sizes
    size += 6 * (N + 1) * sizeof(int);  // nv, nx, nu, ni, nz, ns

    // dynamics
    size += N * sizeof(void *);

    // cost
    size += (N + 1) * sizeof(void *);

    // constraints
    size += (N + 1) * sizeof(void *);

    // regularization
    size += ocp_nlp_reg_dims_calculate_size(N);

    size += sizeof(ocp_nlp_reg_dims);

    size += 8;  // initial align

    return size;
}



int ocp_nlp_dims_calculate_size(void *config_)
{
    ocp_nlp_config *config = config_;

    int N = config->N;

    int ii;

    int size = 0;

    // self
    size += ocp_nlp_dims_calculate_size_self(N);

    // dynamics
    for (ii = 0; ii < N; ii++)
        size += config->dynamics[ii]->dims_calculate_size(config->dynamics[ii]);

    // cost
    for (ii = 0; ii <= N; ii++) size += config->cost[ii]->dims_calculate_size(config->cost[ii]);

    // constraints
    for (ii = 0; ii <= N; ii++)
        size += config->constraints[ii]->dims_calculate_size(config->constraints[ii]);

    // qp solver
    size += config->qp_solver->dims_calculate_size(config->qp_solver, N);

    return size;
}



static ocp_nlp_dims *ocp_nlp_dims_assign_self(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    int ii;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_nlp_dims *dims = (ocp_nlp_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dims);

    // nv
    assign_and_advance_int(N + 1, &dims->nv, &c_ptr);
    // nx
    assign_and_advance_int(N + 1, &dims->nx, &c_ptr);
    // nu
    assign_and_advance_int(N + 1, &dims->nu, &c_ptr);
    // ni
    assign_and_advance_int(N + 1, &dims->ni, &c_ptr);
    // nz
    assign_and_advance_int(N + 1, &dims->nz, &c_ptr);
    // ns
    assign_and_advance_int(N + 1, &dims->ns, &c_ptr);

    // dynamics
    dims->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);

    // cost
    dims->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    // constraints
    dims->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    // regularization
    dims->regularize = ocp_nlp_reg_dims_assign(N, c_ptr);
    c_ptr += ocp_nlp_reg_dims_calculate_size(N);

    /* initialize qp_solver dimensions */
//    dims->qp_solver->N = N;
//    for (ii = 0; ii <= N; ii++)
//    {
        // TODO(dimitris): values below are needed for reformulation of QP when soft constraints
        //   are not supported. Make this a bit more transparent as it clushes with nbx/nbu above.
//        dims->qp_solver->nsbx[ii] = 0;
//        dims->qp_solver->nsbu[ii] = 0;
//        dims->qp_solver->nsg[ii] = 0;
//    }

    // N
    dims->N = N;

	// initialize dimensions to zero by default
	// nv
	for(ii=0; ii<=N; ii++)
		dims->nv[ii] = 0;
	// nx
	for(ii=0; ii<=N; ii++)
		dims->nx[ii] = 0;
	// nu
	for(ii=0; ii<=N; ii++)
		dims->nu[ii] = 0;
	// ni
	for(ii=0; ii<=N; ii++)
		dims->ni[ii] = 0;
	// nz
	for(ii=0; ii<=N; ii++)
		dims->nz[ii] = 0;
	// ns
	for(ii=0; ii<=N; ii++)
		dims->ns[ii] = 0;
	// TODO initialize dims to zero by default also in modules !!!!!!!

    // assert
    assert((char *) raw_memory + ocp_nlp_dims_calculate_size_self(N) >= c_ptr);

    return dims;
}



ocp_nlp_dims *ocp_nlp_dims_assign(void *config_, void *raw_memory)
{
    ocp_nlp_config *config = config_;

    int N = config->N;

    int ii;

    char *c_ptr = (char *) raw_memory;

    // self
    ocp_nlp_dims *dims = ocp_nlp_dims_assign_self(N, c_ptr);
    c_ptr += ocp_nlp_dims_calculate_size_self(N);

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        dims->dynamics[ii] = config->dynamics[ii]->dims_assign(config->dynamics[ii], c_ptr);
        c_ptr += config->dynamics[ii]->dims_calculate_size(config->dynamics[ii]);
    }

    // cost
    for (ii = 0; ii <= N; ii++)
    {
        dims->cost[ii] = config->cost[ii]->dims_assign(config->cost[ii], c_ptr);
        c_ptr += config->cost[ii]->dims_calculate_size(config->cost[ii]);
    }

    // constraints
    for (ii = 0; ii <= N; ii++)
    {
        dims->constraints[ii] =
            config->constraints[ii]->dims_assign(config->constraints[ii], c_ptr);
        c_ptr += config->constraints[ii]->dims_calculate_size(config->constraints[ii]);
    }

    // qp solver
    dims->qp_solver = config->qp_solver->dims_assign(config->qp_solver, N, c_ptr);
    c_ptr += config->qp_solver->dims_calculate_size(config->qp_solver, N);

    // assert
    assert((char *) raw_memory + ocp_nlp_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void ocp_nlp_dims_set_opt_vars(void *config_, void *dims_, const char *field,
                                    const void* value_array)
{
    // to set dimension nx, nu, nz, ns (number of slacks = number of soft constraints)
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;

    int ii;

    int N = config->N;
    int *int_array = (int *) value_array;

    /* set ocp_nlp dimension */
    if (!strcmp(field, "nx"))
    {
        // opt var
        for (ii = 0; ii <= N; ii++)
        {
            // set nx
            dims->nx[ii] = int_array[ii];
            // update nv
            dims->nv[ii] = dims->nu[ii] + dims->nx[ii] + 2 * dims->ns[ii];
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            config->cost[i]->dims_set(config->cost[i],
                                      dims->cost[i], "nx", &int_array[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                          dims->dynamics[i], "nx", &int_array[i]);
        }
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                           dims->dynamics[i], "nx1", &int_array[i+1]);
        }
        // constraints
        for (int i = 0; i <= N; i++)
        {
            config->constraints[i]->dims_set(config->constraints[i], dims->constraints[i],
                                                 "nx", &int_array[i]);
        }
        // qp solver
        for (int i = 0; i <= N; i++)
        {
            config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "nx", &int_array[i]);
        }
        // regularization
        for (ii = 0; ii <= N; ii++)
        {
            config->regularize->dims_set(config->regularize, dims->regularize, ii, "nx", &int_array[ii]);
        }
    }
    else if (!strcmp(field, "nu"))
    {
        // nlp opt var
        for (int ii = 0; ii <= N; ii++)
        {
            // set nu
            dims->nu[ii] = int_array[ii];
            // update nv
            dims->nv[ii] = dims->nu[ii] + dims->nx[ii] + 2 * dims->ns[ii];
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            config->cost[i]->dims_set(config->cost[i],
                                      dims->cost[i], "nu", &int_array[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                          dims->dynamics[i], "nu", &int_array[i]);
        }
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                           dims->dynamics[i], "nu1", &int_array[i+1]);
        }
        // constraints
        for (int i = 0; i <= N; i++)
        {
            config->constraints[i]->dims_set(config->constraints[i], dims->constraints[i],
                                                 "nu", &int_array[i]);
        }
        // qp solver
        for (int i = 0; i <= N; i++)
        {
            config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "nu", &int_array[i]);
        }
        // regularization
        for (ii = 0; ii <= N; ii++)
        {
            config->regularize->dims_set(config->regularize, dims->regularize, ii, "nu", &int_array[ii]);
        }
    }
    else if (!strcmp(field, "nz"))
    {
        // nlp opt var
        for (int ii = 0; ii <= N; ii++)
        {
            // set nz
            dims->nz[ii] = int_array[ii];
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            config->cost[i]->dims_set(config->cost[i],
                                      dims->cost[i], "nz", &int_array[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                          dims->dynamics[i], "nz", &int_array[i]);
        }
        // constraints
        for (int i = 0; i <= N; i++)
        {
            config->constraints[i]->dims_set(config->constraints[i], dims->constraints[i],
                                                 "nz", &int_array[i]);
        }
    }
    else if (!strcmp(field, "ns"))
    {
        // nlp opt var
        for (int ii = 0; ii <= N; ii++)
        {
            // set ns
            dims->ns[ii] = int_array[ii];
            // update nv
            dims->nv[ii] = dims->nu[ii] + dims->nx[ii] + 2 * dims->ns[ii];
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            config->cost[i]->dims_set(config->cost[i],
                                      dims->cost[i], "ns", &int_array[i]);
        }
        // qp solver
        for (int i = 0; i <= N; i++)
        {
            config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "ns",
                                        &int_array[i]);
        }
    }
    else
    {
        printf("error: dims type not available in module ocp_nlp: %s", field);
        exit(1);
    }


#if 0
    /* set ocp_nlp submodule dimensions */
    if (strcmp(field, "ns"))  //  dynamics do not contain slack/soft constraints
    {
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                          dims->dynamics[i], field, &int_array[i]);
        }
    }

    if (!strcmp(field, "nu"))
    {
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                           dims->dynamics[i], "nu1", &int_array[i+1]);
        }
    }
    if (!strcmp(field, "nx"))
    {
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i],
                                           dims->dynamics[i], "nx1", &int_array[i+1]);
        }
    }

    for (int i = 0; i <= N; i++)  // cost
    {
        config->cost[i]->dims_set(config->cost[i],
                                  dims->cost[i], field, &int_array[i]);
    }

    for (int i = 0; i <= N; i++)  // constraints
    {
        config->constraints[i]->dims_set(config->constraints[i], dims->constraints[i],
                                             field, &int_array[i]);
    }

    if (strcmp(field, "nz"))  //  qp_solver does not contain nz
    {
        for (int i = 0; i <= N; i++)  // qp_solver
        {
            config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, field,
                                        &int_array[i]);
        }
    }
#endif

    return;

}



void ocp_nlp_dims_set_constraints(void *config_, void *dims_, int stage, const char *field,
                                  const void* value_)
{
    // to set dimension nbx, nbu, ng, nh, nq (quadratic over nonlinear)
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;

    int *int_value = (int *) value_;
    int i = stage;

    // set in constraint module
    config->constraints[i]->dims_set(config->constraints[i], dims->constraints[i],
                                        field, int_value);
    // update ni in ocp_nlp dimensions
    config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                                        "ni", &dims->ni[i]);

    // update qp_solver dims
    if ( (!strcmp(field, "nbx")) || (!strcmp(field, "nbu")) )
    {
		// qp solver
        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, field, int_value);

		// regularization
        config->regularize->dims_set(config->regularize, dims->regularize, i, (char *) field, int_value);
    }
    else if ( (!strcmp(field, "nsbx")) || (!strcmp(field, "nsbu")) )
    {
		// qp solver
        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, field, int_value);
    }
    else if ( (!strcmp(field, "ng")) || (!strcmp(field, "nh")) || (!strcmp(field, "nphi")))
    {
        // update ng_qp_solver in qp_solver
        int ng_qp_solver;
        config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                                         "ng_qp_solver", &ng_qp_solver);

		// qp solver
        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "ng", &ng_qp_solver);

		// regularization
        config->regularize->dims_set(config->regularize, dims->regularize, i, "ng", &ng_qp_solver);
    }
    else if ( (!strcmp(field, "nsg")) || (!strcmp(field, "nsh")) || (!strcmp(field, "nsphi")))
    {
        // update ng_qp_solver in qp_solver
        int nsg_qp_solver;
        config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i], "nsg_qp_solver", &nsg_qp_solver);

		// qp solver
        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "nsg", &nsg_qp_solver);
    }

	return;
}



void ocp_nlp_dims_set_cost(void *config_, void *dims_, int stage,
                           const char *field, const void* value_)
{
    // to set dimension ny (output)
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;

    int *int_value = (int *) value_;

    config->cost[stage]->dims_set(config->cost[stage], dims->cost[stage], field, int_value);
}



void ocp_nlp_dims_set_dynamics(void *config_, void *dims_, int stage,
                               const char *field, const void* value)
{
    // mainly for gnsf dimensions
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;

    int *int_value = (int *) value;

    config->dynamics[stage]->dims_set(config->dynamics[stage], dims->dynamics[stage], field, int_value);
}




/************************************************
 * in
 ************************************************/

int ocp_nlp_in_calculate_size_self(int N)
{
    int size = sizeof(ocp_nlp_in);

    size += N * sizeof(double);  // Ts

    size += N * sizeof(void *);  // dynamics

    size += (N + 1) * sizeof(void *);  // cost

    size += (N + 1) * sizeof(void *);  // constraints

    return size;
}



int ocp_nlp_in_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    int ii;

    int N = dims->N;

    int size = ocp_nlp_in_calculate_size_self(N);

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        size +=
            config->dynamics[ii]->model_calculate_size(config->dynamics[ii], dims->dynamics[ii]);
    }

    // cost
    for (ii = 0; ii <= N; ii++)
    {
        size += config->cost[ii]->model_calculate_size(config->cost[ii], dims->cost[ii]);
    }

    // constraints
    for (ii = 0; ii <= N; ii++)
    {
        size += config->constraints[ii]->model_calculate_size(config->constraints[ii],
                                                              dims->constraints[ii]);
    }

    size += 8;  // initial align

    //  make_int_multiple_of(64, &size);

    return size;
}



ocp_nlp_in *ocp_nlp_in_assign_self(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_nlp_in *in = (ocp_nlp_in *) c_ptr;
    c_ptr += sizeof(ocp_nlp_in);

    // Ts
    in->Ts = (double *) c_ptr;
    c_ptr += N * sizeof(double);

    // dynamics
    in->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);

    // cost
    in->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    // constraints
    in->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    return in;
}



ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
    int ii;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    // struct
    ocp_nlp_in *in = ocp_nlp_in_assign_self(N, c_ptr);
    c_ptr += ocp_nlp_in_calculate_size_self(N);

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        in->dynamics[ii] =
            config->dynamics[ii]->model_assign(config->dynamics[ii], dims->dynamics[ii], c_ptr);
        c_ptr +=
            config->dynamics[ii]->model_calculate_size(config->dynamics[ii], dims->dynamics[ii]);
    }

    // cost
    for (ii = 0; ii <= N; ii++)
    {
        in->cost[ii] = config->cost[ii]->model_assign(config->cost[ii], dims->cost[ii], c_ptr);
        c_ptr += config->cost[ii]->model_calculate_size(config->cost[ii], dims->cost[ii]);
    }

    // constraints
    for (ii = 0; ii <= N; ii++)
    {
        in->constraints[ii] = config->constraints[ii]->model_assign(config->constraints[ii],
                                                                    dims->constraints[ii], c_ptr);
        c_ptr += config->constraints[ii]->model_calculate_size(config->constraints[ii],
                                                               dims->constraints[ii]);
    }

    assert((char *) raw_memory + ocp_nlp_in_calculate_size(config, dims) >= c_ptr);

    return in;
}



/************************************************
 * out
 ************************************************/

int ocp_nlp_out_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;

    int size = sizeof(ocp_nlp_out);

    size += 4 * (N + 1) * sizeof(struct blasfeo_dvec);  // ux, lam, t, z
    size += 1 * N * sizeof(struct blasfeo_dvec);        // pi

    for (int ii = 0; ii < N; ii++)
    {
        size += 1 * blasfeo_memsize_dvec(nv[ii]);      // ux
        size += 1 * blasfeo_memsize_dvec(nz[ii]);      // z
        size += 2 * blasfeo_memsize_dvec(2 * ni[ii]);  // lam, t
        size += 1 * blasfeo_memsize_dvec(nx[ii + 1]);  // pi
    }
    size += 1 * blasfeo_memsize_dvec(nv[N]);      // ux
    size += 1 * blasfeo_memsize_dvec(nz[N]);     // z
    size += 2 * blasfeo_memsize_dvec(2 * ni[N]);  // lam, t

    size += 8;   // initial align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    //  make_int_multiple_of(64, &size);

    return size;
}



ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
	// loop index
	int ii;

    // extract sizes
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;

    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_out *out = (ocp_nlp_out *) c_ptr;
    c_ptr += sizeof(ocp_nlp_out);

    // blasfeo_struct align
    align_char_to(8, &c_ptr);

    // blasfeo_dvec_struct
    // ux
    assign_and_advance_blasfeo_dvec_structs(N + 1, &out->ux, &c_ptr);
    // z
    assign_and_advance_blasfeo_dvec_structs(N + 1, &out->z, &c_ptr);
    // pi
    assign_and_advance_blasfeo_dvec_structs(N, &out->pi, &c_ptr);
    // lam
    assign_and_advance_blasfeo_dvec_structs(N + 1, &out->lam, &c_ptr);
    // t
    assign_and_advance_blasfeo_dvec_structs(N + 1, &out->t, &c_ptr);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dvec
    // ux
    for (int ii = 0; ii <= N; ++ii)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[ii], out->ux + ii, &c_ptr);
    }
    // z
    for (int ii = 0; ii <= N; ++ii)
    {
        assign_and_advance_blasfeo_dvec_mem(nz[ii], out->z + ii, &c_ptr);
    }
    // pi
    for (int ii = 0; ii < N; ++ii)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[ii + 1], out->pi + ii, &c_ptr);
    }
    // lam
    for (int ii = 0; ii <= N; ++ii)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[ii], out->lam + ii, &c_ptr);
    }
    // t
    for (int ii = 0; ii <= N; ++ii)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[ii], out->t + ii, &c_ptr);
    }

	// zero solution
	for(ii=0; ii<N; ii++)
	{
		blasfeo_dvecse(nv[ii], 0.0, out->ux+ii, 0);
		blasfeo_dvecse(nz[ii], 0.0, out->z+ii, 0);
		blasfeo_dvecse(nx[ii+1], 0.0, out->pi+ii, 0);
		blasfeo_dvecse(2*ni[ii], 0.0, out->lam+ii, 0);
		blasfeo_dvecse(2*ni[ii], 0.0, out->t+ii, 0);
	}
	ii = N;
	blasfeo_dvecse(nv[ii], 0.0, out->ux+ii, 0);
	blasfeo_dvecse(nz[ii], 0.0, out->z+ii, 0);
	blasfeo_dvecse(2*ni[ii], 0.0, out->lam+ii, 0);
	blasfeo_dvecse(2*ni[ii], 0.0, out->t+ii, 0);

    assert((char *) raw_memory + ocp_nlp_out_calculate_size(config, dims) >= c_ptr);

    return out;
}



/************************************************
 * options
 ************************************************/

int ocp_nlp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    int size = 0;

    size += sizeof(ocp_nlp_opts);

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



void *ocp_nlp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_opts *opts = (ocp_nlp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_opts);

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

    assert((char *) raw_memory + ocp_nlp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_opts *opts = opts_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;
    ocp_nlp_reg_config *regularize = config->regularize;

    int ii;

    int N = dims->N;

    opts->reuse_workspace = 1;
#if defined(ACADOS_WITH_OPENMP)
    opts->num_threads = ACADOS_NUM_THREADS;
#endif

    opts->step_length = 1.0;

    // submodules opts

    // qp solver
    qp_solver->opts_initialize_default(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // regularization
    regularize->opts_initialize_default(regularize, dims->regularize, opts->regularize);

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
        constraints[ii]->opts_initialize_default(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
    }

    return;
}



void ocp_nlp_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_opts *opts = opts_;

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



void ocp_nlp_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_opts *opts = (ocp_nlp_opts *) opts_;
    ocp_nlp_config *config = config_;

    int ii;

    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name, i.e. substring in field before '_'
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (ii=0; ii<module_length; ii++)
            module[ii] = field[ii];
        module[module_length] = '\0'; // add end of string
        ptr_module = module;
    }

    // pass options to QP module
    if ( ptr_module!=NULL && (!strcmp(ptr_module, "qp")) )
    {
        config->qp_solver->opts_set(config->qp_solver, opts->qp_solver_opts,
                                    field+module_length+1, value);
    }
    // pass options to dynamics module
    else // nlp opts
    {
        if (!strcmp(field, "reuse_workspace"))
        {
            int* reuse_workspace = (int *) value;
            opts->reuse_workspace = *reuse_workspace;
        }
        else if (!strcmp(field, "num_threads"))
        {
            int* num_threads = (int *) value;
            opts->num_threads = *num_threads;
        }
        else if (!strcmp(field, "step_length"))
        {
            double* step_length = (double *) value;
            opts->step_length = *step_length;
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
            // constraints TODO disabled for now as prevents convergence !!!
//            for (ii=0; ii<=N; ii++)
//                config->constraints[ii]->opts_set(config->constraints[ii], opts->constraints[ii], "compute_hess", value);
        }
        else
        {
            printf("\nerror: ocp_nlp_opts_set: wrong field: %s\n", field);
            exit(1);
        }
    }

    return;

}



void ocp_nlp_opts_set_at_stage(void *config_, void *opts_, int stage, const char *field, void* value)
{
    ocp_nlp_opts *opts = (ocp_nlp_opts *) opts_;
    ocp_nlp_config *config = config_;

    int ii;

    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (ii=0; ii<module_length; ii++)
            module[ii] = field[ii];
        module[module_length] = '\0'; // add end of string
        ptr_module = module;
    }

    // pass options to dynamics module
    if ( ptr_module!=NULL && (!strcmp(ptr_module, "dynamics")) )
    {
        config->dynamics[stage]->opts_set( config->dynamics[stage], opts->dynamics[stage],
                                           field+module_length+1, value );
    }
    // pass options to cost module
    else if ( ptr_module!=NULL && (!strcmp(ptr_module, "cost")) )
    {
        config->cost[stage]->opts_set( config->cost[stage], opts->cost[stage],
                                                 field+module_length+1, value);
    }
    // pass options to constraint module
    else if ( ptr_module!=NULL && (!strcmp(ptr_module, "constraints")) )
    {
        config->constraints[stage]->opts_set( config->constraints[stage], opts->constraints[stage],
                                              (char *) field+module_length+1, value);
    }
	else
	{
		printf("\nerror: ocp_nlp_opts_set_at_stage: wrong field: %s\n", field);
		exit(1);
    }

    return;

}



/************************************************
 * memory
 ************************************************/

int ocp_nlp_memory_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nz = dims->nz;
    int *nu = dims->nu;
    int *ni = dims->ni;

    int size = sizeof(ocp_nlp_memory);

    // qp in
    size += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // qp out
    size += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    // qp solver
    size += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // regularization
    size += config->regularize->memory_calculate_size(config->regularize, dims->regularize, opts->regularize);

    // dynamics
    size += N * sizeof(void *);
    for (int ii = 0; ii < N; ii++)
    {
        size += dynamics[ii]->memory_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
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
        size += constraints[ii]->memory_calculate_size(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
    }

    size += (N+1)*sizeof(bool); // set_sim_guess

    size += (N+1)*sizeof(struct blasfeo_dmat); // dzduxt
    size += 6*(N+1)*sizeof(struct blasfeo_dvec);  // cost_grad ineq_fun ineq_adj dyn_adj sim_guess z_alg
    size += 1*N*sizeof(struct blasfeo_dvec);        // dyn_fun

    for (int ii = 0; ii < N; ii++)
    {
		size += 1*blasfeo_memsize_dmat(nu[ii]+nx[ii], nz[ii]); // dzduxt
		size += 1*blasfeo_memsize_dvec(nz[ii]); // z_alg
        size += 2*blasfeo_memsize_dvec(nv[ii]);           // cost_grad ineq_adj
        size += 1*blasfeo_memsize_dvec(nu[ii] + nx[ii]);  // dyn_adj
        size += 1*blasfeo_memsize_dvec(nx[ii + 1]);       // dyn_fun
        size += 1*blasfeo_memsize_dvec(2 * ni[ii]);       // ineq_fun
        size += 1*blasfeo_memsize_dvec(nx[ii] + nz[ii]); // sim_guess
    }
	size += 1*blasfeo_memsize_dmat(nu[N]+nx[N], nz[N]); // dzduxt
	size += 1*blasfeo_memsize_dvec(nz[N]); // z_alg
    size += 2*blasfeo_memsize_dvec(nv[N]);          // cost_grad ineq_adj
    size += 1*blasfeo_memsize_dvec(nu[N] + nx[N]);  // dyn_adj
    size += 1*blasfeo_memsize_dvec(2 * ni[N]);      // ineq_fun
    size += 1*blasfeo_memsize_dvec(nx[N] + nz[N]);  // sim_guess

    size += 8;   // initial align
    size += 8;   // middle align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    make_int_multiple_of(8, &size);

    return size;
}



ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts, void *raw_memory)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // extract sizes
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nz = dims->nz;
    int *nu = dims->nu;
    int *ni = dims->ni;

    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_nlp_memory *mem = (ocp_nlp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_memory);

    // dynamics
    mem->dynamics = (void **) c_ptr;
    c_ptr += N*sizeof(void *);

    // cost
    mem->cost = (void **) c_ptr;
    c_ptr += (N+1)*sizeof(void *);

    // constraints
    mem->constraints = (void **) c_ptr;
    c_ptr += (N+1)*sizeof(void *);

    // middle align
    align_char_to(8, &c_ptr);

    // qp in
    mem->qp_in = ocp_qp_in_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // qp out
    mem->qp_out = ocp_qp_out_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    // QP solver
    mem->qp_solver_mem = qp_solver->memory_assign(qp_solver, dims->qp_solver, opts->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // regularization
    mem->regularize_mem = config->regularize->memory_assign(config->regularize, dims->regularize,
                                                            opts->regularize, c_ptr);
    c_ptr += config->regularize->memory_calculate_size(config->regularize, dims->regularize,
                                                       opts->regularize);

    // dynamics
    for (int ii = 0; ii < N; ii++)
    {
        mem->dynamics[ii] = dynamics[ii]->memory_assign(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii], c_ptr);
        c_ptr += dynamics[ii]->memory_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
    }

    // cost
    for (int ii = 0; ii <= N; ii++)
    {
        mem->cost[ii] = cost[ii]->memory_assign(cost[ii], dims->cost[ii], opts->cost[ii], c_ptr);
        c_ptr += cost[ii]->memory_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
    }

    // constraints
    for (int ii = 0; ii <= N; ii++)
    {
        mem->constraints[ii] = constraints[ii]->memory_assign(constraints[ii],
                                            dims->constraints[ii], opts->constraints[ii], c_ptr);
        c_ptr += constraints[ii]->memory_calculate_size( constraints[ii], dims->constraints[ii],
                                                                 opts->constraints[ii]);
    }

    // set_sim_guess
    assign_and_advance_bool(N+1, &mem->set_sim_guess, &c_ptr);
    for (int ii = 0; ii <= N; ++ii)
    {
        mem->set_sim_guess[ii] = false;
    }

    // blasfeo_struct align
    align_char_to(8, &c_ptr);

    // dzduxt
    mem->dzduxt = (struct blasfeo_dmat *) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dmat);
    // z_alg
    mem->z_alg = (struct blasfeo_dvec *) c_ptr;
    c_ptr += (N+1)*sizeof(struct blasfeo_dvec);

    // cost_grad
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->cost_grad, &c_ptr);
    // ineq_fun
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->ineq_fun, &c_ptr);
    // ineq_adj
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->ineq_adj, &c_ptr);
    // dyn_fun
    assign_and_advance_blasfeo_dvec_structs(N, &mem->dyn_fun, &c_ptr);
    // dyn_adj
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->dyn_adj, &c_ptr);
    // sim_guess
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->sim_guess, &c_ptr);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // dzduxt
    for (int ii=0; ii<=N; ii++)
    {
        blasfeo_create_dmat(nu[ii]+nx[ii], nz[ii], mem->dzduxt+ii, c_ptr);
        c_ptr += blasfeo_memsize_dmat(nu[ii]+nx[ii], nz[ii]);
    }
    // z_alg
    for (int ii=0; ii<=N; ii++)
    {
        blasfeo_create_dvec(nz[ii], mem->z_alg+ii, c_ptr);
        c_ptr += blasfeo_memsize_dvec(nz[ii]);
    }

    // cost_grad
    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[ii], mem->cost_grad + ii, &c_ptr);
    }
    // ineq_fun
    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[ii], mem->ineq_fun + ii, &c_ptr);
    }
    // ineq_adj
    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[ii], mem->ineq_adj + ii, &c_ptr);
    }
    // dyn_fun
    for (int ii = 0; ii < N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[ii + 1], mem->dyn_fun + ii, &c_ptr);
    }
    // dyn_adj
    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nu[ii] + nx[ii], mem->dyn_adj + ii, &c_ptr);
    }
    // sim_guess
    for (int ii = 0; ii <= N; ++ii)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[ii] + nz[ii], mem->sim_guess + ii, &c_ptr);
        // set to 0;
		blasfeo_dvecse(nx[ii] + nz[ii], 0.0, mem->sim_guess+ii, 0);
        // printf("sim_guess ii %d: %p\n", ii, mem->sim_guess+ii);
    }
    // printf("created memory %p\n", mem);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

int ocp_nlp_workspace_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int ii;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    int size = 0;
    int size_tmp = 0;
    int tmp;

    // nlp
    size += sizeof(ocp_nlp_workspace);

	// tmp_nlp_out
	size += ocp_nlp_out_calculate_size(config, dims);

	// weights_nlp_out
	size += ocp_nlp_out_calculate_size(config, dims);

    // array of pointers
    // cost
    size += (N+1)*sizeof(void *);
    // dynamics
    size += N*sizeof(void *);
    // constraints
    size += (N+1)*sizeof(void *);

	// module workspace
    if (opts->reuse_workspace)
    {

#if defined(ACADOS_WITH_OPENMP)

        // qp solver
        size += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver,
            opts->qp_solver_opts);

        // dynamics
        for (ii = 0; ii < N; ii++)
        {
            size += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
        }

        // cost
        for (ii = 0; ii <= N; ii++)
        {
            size += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
        }

        // constraints
        for (ii = 0; ii <= N; ii++)
        {
            size += constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
        }

#else

        // qp solver
        tmp = qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);
        size_tmp = tmp > size_tmp ? tmp : size_tmp;

        // dynamics
        for (ii = 0; ii < N; ii++)
        {
            tmp = dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
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
            tmp = constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
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
            size += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
        }

        // cost
        for (ii = 0; ii <= N; ii++)
        {
            size += cost[ii]->workspace_calculate_size(cost[ii], dims->cost[ii], opts->cost[ii]);
        }

        // constraints
        for (ii = 0; ii <= N; ii++)
        {
            size += constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
        }

    }

    return size;
}



ocp_nlp_workspace *ocp_nlp_workspace_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts, ocp_nlp_memory *mem, void *raw_memory)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    char *c_ptr = (char *) raw_memory;

	ocp_nlp_workspace *work = (ocp_nlp_workspace *) c_ptr;
    c_ptr += sizeof(ocp_nlp_workspace);

	// tmp_nlp_out
	work->tmp_nlp_out = ocp_nlp_out_assign(config, dims, c_ptr);
	c_ptr += ocp_nlp_out_calculate_size(config, dims);

	// weights_nlp_out
	work->weights_nlp_out = ocp_nlp_out_assign(config, dims, c_ptr);
	c_ptr += ocp_nlp_out_calculate_size(config, dims);

    // array of pointers
    //
    work->dynamics = (void **) c_ptr;
    c_ptr += N*sizeof(void *);
    //
    work->cost = (void **) c_ptr;
    c_ptr += (N+1)*sizeof(void *);
    //
    work->constraints = (void **) c_ptr;
    c_ptr += (N+1)*sizeof(void *);

    if (opts->reuse_workspace)
    {

#if defined(ACADOS_WITH_OPENMP)

        // qp solver
        work->qp_work = (void *) c_ptr;
        c_ptr += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

        // dynamics
        for (int ii = 0; ii < N; ii++)
        {
            work->dynamics[ii] = c_ptr;
            c_ptr += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
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
            c_ptr += constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
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
            c_ptr += dynamics[ii]->workspace_calculate_size(dynamics[ii], dims->dynamics[ii], opts->dynamics[ii]);
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
            c_ptr += constraints[ii]->workspace_calculate_size(constraints[ii], dims->constraints[ii], opts->constraints[ii]);
        }

    }

    assert((char *) work + ocp_nlp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return work;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_initialize_qp(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
         ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{

    int ii;

    int N = dims->N;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (ii = 0; ii <= N; ii++)
    {
        // cost
        config->cost[ii]->initialize(config->cost[ii], dims->cost[ii], in->cost[ii],
        		opts->cost[ii], mem->cost[ii], work->cost[ii]);
        // dynamics
        if (ii < N)
            config->dynamics[ii]->initialize(config->dynamics[ii], dims->dynamics[ii],
            		in->dynamics[ii], opts->dynamics[ii], mem->dynamics[ii], work->dynamics[ii]);
        // constraints
        config->constraints[ii]->initialize(config->constraints[ii], dims->constraints[ii],
        		in->constraints[ii], opts->constraints[ii], mem->constraints[ii], work->constraints[ii]);
    }

    return;
}



void ocp_nlp_approximate_qp_matrices(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem,
    ocp_nlp_workspace *work)
{

    int i;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    /* stage-wise multiple shooting lagrangian evaluation */

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i = 0; i <= N; i++)
    {
        // init Hessian to 0 
        blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, mem->qp_in->RSQrq+i, 0, 0);

        // dynamics
        if (i < N)
            config->dynamics[i]->update_qp_matrices(config->dynamics[i], dims->dynamics[i],
                    in->dynamics[i], opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);

        // cost
        config->cost[i]->update_qp_matrices(config->cost[i], dims->cost[i], in->cost[i],
                opts->cost[i], mem->cost[i], work->cost[i]);

        // constraints
        config->constraints[i]->update_qp_matrices(config->constraints[i], dims->constraints[i],
                in->constraints[i], opts->constraints[i], mem->constraints[i], work->constraints[i]);
    }

    /* collect stage-wise evaluations */

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i=0; i <= N; i++)
    {

        // nlp mem: cost_grad
        struct blasfeo_dvec *cost_grad = config->cost[i]->memory_get_grad_ptr(mem->cost[i]);
        blasfeo_dveccp(nv[i], cost_grad, 0, mem->cost_grad + i, 0);

        // nlp mem: dyn_fun
        if (i < N)
        {
            struct blasfeo_dvec *dyn_fun
                = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
            blasfeo_dveccp(nx[i + 1], dyn_fun, 0, mem->dyn_fun + i, 0);
        }

        // nlp mem: dyn_adj
        if (i < N)
        {
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[i]->memory_get_adj_ptr(mem->dynamics[i]);
            blasfeo_dveccp(nu[i] + nx[i], dyn_adj, 0, mem->dyn_adj + i, 0);
        }
        else
        {
            blasfeo_dvecse(nu[N] + nx[N], 0.0, mem->dyn_adj + N, 0);
        }
        if (i > 0)
        {
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[i-1]->memory_get_adj_ptr(mem->dynamics[i-1]);
            blasfeo_daxpy(nx[i], 1.0, dyn_adj, nu[i-1]+nx[i-1], mem->dyn_adj+i, nu[i],
                mem->dyn_adj+i, nu[i]);
        }

        // nlp mem: ineq_fun
        struct blasfeo_dvec *ineq_fun =
            config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, mem->ineq_fun + i, 0);

        // nlp mem: ineq_adj
        struct blasfeo_dvec *ineq_adj =
            config->constraints[i]->memory_get_adj_ptr(mem->constraints[i]);
        blasfeo_dveccp(nv[i], ineq_adj, 0, mem->ineq_adj + i, 0);

    }

    for (i = 0; i <= N; i++)
    {
        // TODO(rien) where should the update happen??? move to qp update ???
        // TODO(all): fix and move where appropriate
        //  if (i<N)
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
void ocp_nlp_approximate_qp_vectors_sqp(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int i;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i = 0; i <= N; i++)
    {
        // g
        blasfeo_dveccp(nv[i], mem->cost_grad + i, 0, mem->qp_in->rqz + i, 0);

        // b
        if (i < N)
            blasfeo_dveccp(nx[i + 1], mem->dyn_fun + i, 0, mem->qp_in->b + i, 0);

        // d
        blasfeo_dveccp(2 * ni[i], mem->ineq_fun + i, 0, mem->qp_in->d + i, 0);
    }

    return;
}




void ocp_nlp_embed_initial_value(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int *ni = dims->ni;

    // constraints
    config->constraints[0]->bounds_update(config->constraints[0], dims->constraints[0],
            in->constraints[0], opts->constraints[0], mem->constraints[0], work->constraints[0]);

    // nlp mem: ineq_fun
    struct blasfeo_dvec *ineq_fun =
        config->constraints[0]->memory_get_fun_ptr(mem->constraints[0]);
    blasfeo_dveccp(2 * ni[0], ineq_fun, 0, mem->ineq_fun, 0);

    // d
    blasfeo_dveccp(2 * ni[0], mem->ineq_fun, 0, mem->qp_in->d, 0);

    return;
}



double ocp_nlp_evaluate_merit_fun(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{

	int i, j;

	int N = dims->N;
	int *nx = dims->nx;
	int *ni = dims->ni;

	double merit_fun = 0.0;

	// compute fun value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i=0; i<=N; i++)
    {
        // cost
        config->cost[i]->compute_fun(config->cost[i], dims->cost[i], in->cost[i], opts->cost[i], mem->cost[i], work->cost[i]);
    }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i=0; i<N; i++)
    {
        // cost
        config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i], in->dynamics[i], opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
    }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i=0; i<=N; i++)
    {
        // constr
        config->constraints[i]->compute_fun(config->constraints[i], dims->constraints[i],
                                            in->constraints[i], opts->constraints[i],
                                            mem->constraints[i], work->constraints[i]);
    }

	double *tmp_fun;
	double tmp;
	struct blasfeo_dvec *tmp_fun_vec;

	double cost_fun = 0.0;
	for(i=0; i<=N; i++)
	{
		tmp_fun = config->cost[i]->memory_get_fun_ptr(mem->cost[i]);
		cost_fun += *tmp_fun;
	}

	double dyn_fun = 0.0;
	for(i=0; i<N; i++)
	{
//		printf("\ni %d\n", i);
		tmp_fun_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
//		blasfeo_print_exp_tran_dvec(nx[i+1], tmp_fun_vec, 0);
//		blasfeo_print_exp_tran_dvec(nx[i+1], work->weights_nlp_out->pi+i, 0);
		for(j=0; j<nx[i+1]; j++)
		{
//			printf("\n%e %e\n", fabs(BLASFEO_DVECEL(work->weights_nlp_out->pi+i, j)), fabs(BLASFEO_DVECEL(tmp_fun_vec, j)));
			dyn_fun += fabs(BLASFEO_DVECEL(work->weights_nlp_out->pi+i, j)) * fabs(BLASFEO_DVECEL(tmp_fun_vec, j));
		}
	}

	double constr_fun = 0.0;
	for(i=0; i<=N; i++)
	{
//		printf("\ni %d\n", i);
		tmp_fun_vec = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
//		blasfeo_print_exp_tran_dvec(2*ni[i], tmp_fun_vec, 0);
//		blasfeo_print_exp_tran_dvec(2*ni[i], work->weights_nlp_out->lam+i, 0);
		for(j=0; j<2*ni[i]; j++)
		{
			tmp = BLASFEO_DVECEL(tmp_fun_vec, j);
			tmp = tmp>0.0 ? tmp : 0.0;
//			printf("\n%e %e\n", fabs(BLASFEO_DVECEL(work->weights_nlp_out->pi+i, j)), fabs(BLASFEO_DVECEL(tmp_fun_vec, j)));
			constr_fun += fabs(BLASFEO_DVECEL(work->weights_nlp_out->lam+i, j)) * tmp;
		}
	}

	merit_fun = cost_fun + dyn_fun + constr_fun;

	printf("\n%e %e %e %e\n", merit_fun, cost_fun, dyn_fun, constr_fun);

	return merit_fun;
}



void ocp_nlp_update_variables_sqp(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int i;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;

    // ocp_nlp_config *config = (ocp_nlp_config *) config_;

	// (fixed) step length
    double alpha = opts->step_length;

#if 0 // XXX test piece of code
	double tmp0, tmp1;

	// current point
    for (i = 0; i <= N; i++)
        blasfeo_dveccp(nv[i], out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);
    
	for (i = 0; i < N; i++)
		blasfeo_dveccp(nx[i+1], out->pi+i, 0, work->tmp_nlp_out->pi+i, 0);

	for (i = 0; i <= N; i++)
		blasfeo_dveccp(2*ni[i], out->lam+i, 0, work->tmp_nlp_out->lam+i, 0);

        // linear update of algebraic variables using state and input sensitivity
//        if (i < N)
//        {
//            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, mem->dzduxt+i, 0, 0, mem->qp_out->ux+i, 0, 1.0, mem->z_alg+i, 0, out->z+i, 0); 
//        }

	// initialize weights
	if(mem->sqp_iter[0]==0)
	{
		for (i = 0; i < N; i++)
			blasfeo_dveccp(nx[i+1], out->pi+i, 0, work->weights_nlp_out->pi+i, 0);

		for (i = 0; i <= N; i++)
			blasfeo_dveccp(2*ni[i], out->lam+i, 0, work->weights_nlp_out->lam+i, 0);
	}

	// update weigths
	for (i = 0; i < N; i++)
	{
		for(j=0; j<nx[i+1]; j++)
		{
			tmp0 = fabs(BLASFEO_DVECEL(work->weights_nlp_out->pi+i, j));
			tmp1 = 0.5 * (tmp0 + fabs(BLASFEO_DVECEL(mem->qp_out->pi+i, j)));
			BLASFEO_DVECEL(work->weights_nlp_out->pi+i, j) = tmp0>tmp1 ? tmp0 : tmp1;
		}
	}
	for (i = 0; i <= N; i++)
	{
		for(j=0; j<2*ni[i]; j++)
		{
			tmp0 = fabs(BLASFEO_DVECEL(work->weights_nlp_out->lam+i, j));
			tmp1 = 0.5 * (tmp0 + fabs(BLASFEO_DVECEL(mem->qp_out->lam+i, j)));
			BLASFEO_DVECEL(work->weights_nlp_out->lam+i, j) = tmp0>tmp1 ? tmp0 : tmp1;
		}
	}

	printf("\n\nmerit fun value\n");
	double merit_fun0 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);

	double alpha_min = 0.2;

	for (j=0; j<10 & alpha>alpha_min; j++)
	{

		for (i = 0; i <= N; i++)
			blasfeo_daxpy(nv[i], alpha, mem->qp_out->ux+i, 0, out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);

		printf("\n%d tmp merit fun value\n", j);
		double merit_fun1 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);

		if(merit_fun1 < merit_fun0)
		{
			break;
		}
		else
		{
			alpha *= 0.7;
		}
	
	}

	printf("\nalpha %f\n", alpha);

#endif






#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i = 0; i <= N; i++)
    {
        // (full) step in primal variables
        blasfeo_daxpy(nv[i], alpha, mem->qp_out->ux + i, 0, out->ux + i, 0, out->ux + i, 0);
    
        // update dual variables
        if (i < N)
        {
            blasfeo_dvecsc(nx[i+1], 1.0-alpha, out->pi+i, 0);
            blasfeo_daxpy(nx[i+1], alpha, mem->qp_out->pi+i, 0, out->pi+i, 0, out->pi+i, 0);
        }

        blasfeo_dvecsc(2*ni[i], 1.0-alpha, out->lam+i, 0);
        blasfeo_daxpy(2*ni[i], alpha, mem->qp_out->lam+i, 0, out->lam+i, 0, out->lam+i, 0);

        // update slack values
        blasfeo_dvecsc(2*ni[i], 1.0-alpha, out->t+i, 0);
        blasfeo_daxpy(2*ni[i], alpha, mem->qp_out->t+i, 0, out->t+i, 0, out->t+i, 0);

        // linear update of algebraic variables using state and input sensitivity
        if (i < N)
        {
            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, mem->dzduxt+i, 0, 0, mem->qp_out->ux+i, 0, 1.0, mem->z_alg+i, 0, out->z+i, 0); 
        }
    }

    return;
}



/************************************************
 * residuals
 ************************************************/

int ocp_nlp_res_calculate_size(ocp_nlp_dims *dims)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

    int size = sizeof(ocp_nlp_res);

    size += 3 * (N + 1) * sizeof(struct blasfeo_dvec);  // res_g res_d res_m
    size += 1 * N * sizeof(struct blasfeo_dvec);        // res_b

    for (int ii = 0; ii < N; ii++)
    {
        size += 1 * blasfeo_memsize_dvec(nv[ii]);      // res_g
        size += 1 * blasfeo_memsize_dvec(nx[ii + 1]);  // res_b
        size += 2 * blasfeo_memsize_dvec(2 * ni[ii]);  // res_d res_m
    }
    size += 1 * blasfeo_memsize_dvec(nv[N]);      // res_g
    size += 2 * blasfeo_memsize_dvec(2 * ni[N]);  // res_d res_m

    size += 8;   // initial align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    //  make_int_multiple_of(64, &size);

    return size;
}

ocp_nlp_res *ocp_nlp_res_assign(ocp_nlp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // extract sizes
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_nlp_res *res = (ocp_nlp_res *) c_ptr;
    c_ptr += sizeof(ocp_nlp_res);

    // blasfeo_struct align
    align_char_to(8, &c_ptr);

    // res_g
    assign_and_advance_blasfeo_dvec_structs(N + 1, &res->res_g, &c_ptr);
    // res_b
    assign_and_advance_blasfeo_dvec_structs(N, &res->res_b, &c_ptr);
    // res_d
    assign_and_advance_blasfeo_dvec_structs(N + 1, &res->res_d, &c_ptr);
    // res_m
    assign_and_advance_blasfeo_dvec_structs(N + 1, &res->res_m, &c_ptr);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // res_g
    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[ii], res->res_g + ii, &c_ptr);
    }
    // res_b
    for (int ii = 0; ii < N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[ii + 1], res->res_b + ii, &c_ptr);
    }
    // res_d
    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[ii], res->res_d + ii, &c_ptr);
    }
    // res_m
    for (int ii = 0; ii <= N; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[ii], res->res_m + ii, &c_ptr);
    }

    res->memsize = ocp_nlp_res_calculate_size(dims);

    return res;
}



void ocp_nlp_res_compute(ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_res *res,
                         ocp_nlp_memory *mem)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    double tmp_res;

    // res_g
    res->inf_norm_res_g = 0.0;
    for (int ii = 0; ii <= N; ii++)
    {
        blasfeo_daxpy(nv[ii], -1.0, mem->ineq_adj + ii, 0, mem->cost_grad + ii, 0, res->res_g + ii,
                      0);
        blasfeo_daxpy(nu[ii] + nx[ii], -1.0, mem->dyn_adj + ii, 0, res->res_g + ii, 0,
                      res->res_g + ii, 0);
        blasfeo_dvecnrm_inf(nv[ii], res->res_g + ii, 0, &tmp_res);
        res->inf_norm_res_g = tmp_res > res->inf_norm_res_g ? tmp_res : res->inf_norm_res_g;
    }

    // res_b
    res->inf_norm_res_b = 0.0;
    for (int ii = 0; ii < N; ii++)
    {
        blasfeo_dveccp(nx[ii + 1], mem->dyn_fun + ii, 0, res->res_b + ii, 0);
        blasfeo_dvecnrm_inf(nx[ii + 1], res->res_b + ii, 0, &tmp_res);
        res->inf_norm_res_b = tmp_res > res->inf_norm_res_b ? tmp_res : res->inf_norm_res_b;
    }

    // res_d
    res->inf_norm_res_d = 0.0;
    for (int ii = 0; ii <= N; ii++)
    {
        blasfeo_daxpy(2 * ni[ii], 1.0, out->t + ii, 0, mem->ineq_fun + ii, 0, res->res_d + ii, 0);
        blasfeo_dvecnrm_inf(2 * ni[ii], res->res_d + ii, 0, &tmp_res);
        res->inf_norm_res_d = tmp_res > res->inf_norm_res_d ? tmp_res : res->inf_norm_res_d;
    }

    // res_m
    res->inf_norm_res_m = 0.0;
    for (int ii = 0; ii <= N; ii++)
    {
        blasfeo_dvecmul(2 * ni[ii], out->lam + ii, 0, out->t + ii, 0, res->res_m + ii, 0);
        blasfeo_dvecnrm_inf(2 * ni[ii], res->res_m + ii, 0, &tmp_res);
        res->inf_norm_res_m = tmp_res > res->inf_norm_res_m ? tmp_res : res->inf_norm_res_m;
    }

    return;
}


