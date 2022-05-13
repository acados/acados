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
#include "acados/utils/print.h"
// openmp
#if defined(ACADOS_WITH_OPENMP)
#include <omp.h>
#endif


/************************************************
 * config
 ************************************************/

acados_size_t ocp_nlp_config_calculate_size(int N)
{
    acados_size_t size = 0;

    // self
    size += sizeof(ocp_nlp_config);

    // qp solver
    size += 1 * ocp_qp_xcond_solver_config_calculate_size();

    // regularization
    size += ocp_nlp_reg_config_calculate_size();


    // dynamics
    size += N * sizeof(ocp_nlp_dynamics_config *);

    for (int i = 0; i < N; i++) size += ocp_nlp_dynamics_config_calculate_size();

    // cost
    size += (N + 1) * sizeof(ocp_nlp_cost_config *);

    for (int i = 0; i <= N; i++) size += ocp_nlp_cost_config_calculate_size();

    // constraints
    size += (N + 1) * sizeof(ocp_nlp_constraints_config *);

    for (int i = 0; i <= N; i++) size += ocp_nlp_constraints_config_calculate_size();

    return size;
}



ocp_nlp_config *ocp_nlp_config_assign(int N, void *raw_memory)
{
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

    for (int i = 0; i < N; i++)
    {
        config->dynamics[i] = ocp_nlp_dynamics_config_assign(c_ptr);
        c_ptr += ocp_nlp_dynamics_config_calculate_size();
    }

    // cost
    config->cost = (ocp_nlp_cost_config **) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_cost_config *);

    for (int i = 0; i <= N; i++)
    {
        config->cost[i] = ocp_nlp_cost_config_assign(c_ptr);
        c_ptr += ocp_nlp_cost_config_calculate_size();
    }

    // constraints
    config->constraints = (ocp_nlp_constraints_config **) c_ptr;
    c_ptr += (N + 1) * sizeof(ocp_nlp_constraints_config *);

    for (int i = 0; i <= N; i++)
    {
        config->constraints[i] = ocp_nlp_constraints_config_assign(c_ptr);
        c_ptr += ocp_nlp_constraints_config_calculate_size();
    }

    return config;
}



/************************************************
 * dims
 ************************************************/

static acados_size_t ocp_nlp_dims_calculate_size_self(int N)
{
    acados_size_t size = 0;

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
    size += 8;  // intermediate align
    make_int_multiple_of(8, &size);

    return size;
}



acados_size_t ocp_nlp_dims_calculate_size(void *config_)
{
    ocp_nlp_config *config = config_;

    int N = config->N;

    acados_size_t size = 0;

    // self
    size += ocp_nlp_dims_calculate_size_self(N);

    // dynamics
    for (int i = 0; i < N; i++)
        size += config->dynamics[i]->dims_calculate_size(config->dynamics[i]);

    // cost
    for (int i = 0; i <= N; i++) size += config->cost[i]->dims_calculate_size(config->cost[i]);

    // constraints
    for (int i = 0; i <= N; i++)
        size += config->constraints[i]->dims_calculate_size(config->constraints[i]);

    // qp solver
    size += config->qp_solver->dims_calculate_size(config->qp_solver, N);

    return size;
}



static ocp_nlp_dims *ocp_nlp_dims_assign_self(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_nlp_dims *dims = (ocp_nlp_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dims);

    // dynamics
    dims->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);

    // cost
    dims->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    // constraints
    dims->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

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

    // intermediate align
    align_char_to(8, &c_ptr);

    // regularization
    dims->regularize = ocp_nlp_reg_dims_assign(N, c_ptr);
    c_ptr += ocp_nlp_reg_dims_calculate_size(N);

    /* initialize qp_solver dimensions */
//    dims->qp_solver->N = N;
//    for (int i = 0; i <= N; i++)
//    {
        // TODO(dimitris): values below are needed for reformulation of QP when soft constraints
        //   are not supported. Make this a bit more transparent as it clushes with nbx/nbu above.
//        dims->qp_solver->nsbx[i] = 0;
//        dims->qp_solver->nsbu[i] = 0;
//        dims->qp_solver->nsg[i] = 0;
//    }

    // N
    dims->N = N;

    // initialize dimensions to zero by default
    // nv
    for(int i=0; i<=N; i++)
        dims->nv[i] = 0;
    // nx
    for(int i=0; i<=N; i++)
        dims->nx[i] = 0;
    // nu
    for(int i=0; i<=N; i++)
        dims->nu[i] = 0;
    // ni
    for(int i=0; i<=N; i++)
        dims->ni[i] = 0;
    // nz
    for(int i=0; i<=N; i++)
        dims->nz[i] = 0;
    // ns
    for(int i=0; i<=N; i++)
        dims->ns[i] = 0;
    // TODO initialize dims to zero by default also in modules !!!!!!!

    // assert
    assert((char *) raw_memory + ocp_nlp_dims_calculate_size_self(N) >= c_ptr);

    return dims;
}



ocp_nlp_dims *ocp_nlp_dims_assign(void *config_, void *raw_memory)
{
    ocp_nlp_config *config = config_;

    int N = config->N;

    char *c_ptr = (char *) raw_memory;

    // self
    ocp_nlp_dims *dims = ocp_nlp_dims_assign_self(N, c_ptr);
    c_ptr += ocp_nlp_dims_calculate_size_self(N);

    // dynamics
    for (int i = 0; i < N; i++)
    {
        dims->dynamics[i] = config->dynamics[i]->dims_assign(config->dynamics[i], c_ptr);
        c_ptr += config->dynamics[i]->dims_calculate_size(config->dynamics[i]);
    }

    // cost
    for (int i = 0; i <= N; i++)
    {
        dims->cost[i] = config->cost[i]->dims_assign(config->cost[i], c_ptr);
        c_ptr += config->cost[i]->dims_calculate_size(config->cost[i]);
    }

    // constraints
    for (int i = 0; i <= N; i++)
    {
        dims->constraints[i] =
            config->constraints[i]->dims_assign(config->constraints[i], c_ptr);
        c_ptr += config->constraints[i]->dims_calculate_size(config->constraints[i]);
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

    int N = config->N;
    int *int_array = (int *) value_array;

    /* set ocp_nlp dimension */
    if (!strcmp(field, "nx"))
    {
        // opt var
        for (int i = 0; i <= N; i++)
        {
            // set nx
            dims->nx[i] = int_array[i];
            // update nv
            dims->nv[i] = dims->nu[i] + dims->nx[i] + 2 * dims->ns[i];
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
        for (int i = 0; i <= N; i++)
        {
            config->regularize->dims_set(config->regularize, dims->regularize, i, "nx", &int_array[i]);
        }
    }
    else if (!strcmp(field, "nu"))
    {
        // nlp opt var
        for (int i = 0; i <= N; i++)
        {
            // set nu
            dims->nu[i] = int_array[i];
            // update nv
            dims->nv[i] = dims->nu[i] + dims->nx[i] + 2 * dims->ns[i];
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
        for (int i = 0; i <= N; i++)
        {
            config->regularize->dims_set(config->regularize, dims->regularize, i, "nu", &int_array[i]);
        }
    }
    else if (!strcmp(field, "nz"))
    {
        // nlp opt var
        for (int i = 0; i <= N; i++)
        {
            // set nz
            dims->nz[i] = int_array[i];
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
        for (int i = 0; i <= N; i++)
        {
            // set ns
            dims->ns[i] = int_array[i];
            // update nv
            dims->nv[i] = dims->nu[i] + dims->nx[i] + 2 * dims->ns[i];
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
    else if ( (!strcmp(field, "nbxe")) || (!strcmp(field, "nbue")) )
    {
        // qp solver
        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, field, int_value);
    }
    else if ( (!strcmp(field, "nge")) || (!strcmp(field, "nhe")) || (!strcmp(field, "nphie")))
    {
        // update ng_qp_solver in qp_solver
        int ng_qp_solver;
        config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                                         "nge_qp_solver", &ng_qp_solver);

        // qp solver
        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "nge", &ng_qp_solver);
    }
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

acados_size_t ocp_nlp_in_calculate_size_self(int N)
{
    acados_size_t size = sizeof(ocp_nlp_in);

    size += N * sizeof(double);  // Ts

    size += N * sizeof(void *);  // dynamics

    size += (N + 1) * sizeof(void *);  // cost

    size += (N + 1) * sizeof(void *);  // constraints

    size += 3*8;  // aligns
    return size;
}



acados_size_t ocp_nlp_in_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    int N = dims->N;

    acados_size_t size = ocp_nlp_in_calculate_size_self(N);

    // dynamics
    for (int i = 0; i < N; i++)
    {
        size +=
            config->dynamics[i]->model_calculate_size(config->dynamics[i], dims->dynamics[i]);
    }

    // cost
    for (int i = 0; i <= N; i++)
    {
        size += config->cost[i]->model_calculate_size(config->cost[i], dims->cost[i]);
    }

    // constraints
    for (int i = 0; i <= N; i++)
    {
        size += config->constraints[i]->model_calculate_size(config->constraints[i],
                                                              dims->constraints[i]);
    }
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

    // align
    align_char_to(8, &c_ptr);

    // Ts
    assign_and_advance_double(N, &in->Ts, &c_ptr);

    // dynamics
    in->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);

    // cost
    in->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    // constraints
    in->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_in_calculate_size_self(N) >= c_ptr);

    return in;
}



ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    // struct
    ocp_nlp_in *in = ocp_nlp_in_assign_self(N, c_ptr);
    c_ptr += ocp_nlp_in_calculate_size_self(N);

    // dynamics
    for (int i = 0; i < N; i++)
    {
        in->dynamics[i] =
            config->dynamics[i]->model_assign(config->dynamics[i], dims->dynamics[i], c_ptr);
        c_ptr +=
            config->dynamics[i]->model_calculate_size(config->dynamics[i], dims->dynamics[i]);
    }

    // cost
    for (int i = 0; i <= N; i++)
    {
        in->cost[i] = config->cost[i]->model_assign(config->cost[i], dims->cost[i], c_ptr);
        c_ptr += config->cost[i]->model_calculate_size(config->cost[i], dims->cost[i]);
    }

    // constraints
    for (int i = 0; i <= N; i++)
    {
        in->constraints[i] = config->constraints[i]->model_assign(config->constraints[i],
                                                                    dims->constraints[i], c_ptr);
        c_ptr += config->constraints[i]->model_calculate_size(config->constraints[i],
                                                               dims->constraints[i]);
    }

    assert((char *) raw_memory + ocp_nlp_in_calculate_size(config, dims) >= c_ptr);

    return in;
}



/************************************************
 * out
 ************************************************/

acados_size_t ocp_nlp_out_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;

    acados_size_t size = sizeof(ocp_nlp_out);

    size += 4 * (N + 1) * sizeof(struct blasfeo_dvec);  // ux, lam, t, z
    size += 1 * N * sizeof(struct blasfeo_dvec);        // pi

    for (int i = 0; i < N; i++)
    {
        size += 1 * blasfeo_memsize_dvec(nv[i]);      // ux
        size += 1 * blasfeo_memsize_dvec(nz[i]);      // z
        size += 2 * blasfeo_memsize_dvec(2 * ni[i]);  // lam, t
        size += 1 * blasfeo_memsize_dvec(nx[i + 1]);  // pi
    }
    size += 1 * blasfeo_memsize_dvec(nv[N]);      // ux
    size += 1 * blasfeo_memsize_dvec(nz[N]);     // z
    size += 2 * blasfeo_memsize_dvec(2 * ni[N]);  // lam, t

    size += 8;   // initial align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    make_int_multiple_of(8, &size);

    return size;
}



ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
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
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[i], out->ux + i, &c_ptr);
    }
    // z
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(nz[i], out->z + i, &c_ptr);
    }
    // pi
    for (int i = 0; i < N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[i + 1], out->pi + i, &c_ptr);
    }
    // lam
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[i], out->lam + i, &c_ptr);
    }
    // t
    for (int i = 0; i <= N; ++i)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[i], out->t + i, &c_ptr);
    }

    // zero solution
    for(int i=0; i<N; i++)
    {
        blasfeo_dvecse(nv[i], 0.0, out->ux+i, 0);
        blasfeo_dvecse(nz[i], 0.0, out->z+i, 0);
        blasfeo_dvecse(nx[i+1], 0.0, out->pi+i, 0);
        blasfeo_dvecse(2*ni[i], 0.0, out->lam+i, 0);
        blasfeo_dvecse(2*ni[i], 0.0, out->t+i, 0);
    }
    blasfeo_dvecse(nv[N], 0.0, out->ux+N, 0);
    blasfeo_dvecse(nz[N], 0.0, out->z+N, 0);
    blasfeo_dvecse(2*ni[N], 0.0, out->lam+N, 0);
    blasfeo_dvecse(2*ni[N], 0.0, out->t+N, 0);

    assert((char *) raw_memory + ocp_nlp_out_calculate_size(config, dims) >= c_ptr);

    return out;
}



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_opts);

    size += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    size += config->regularize->opts_calculate_size();

    // dynamics
    size += N * sizeof(void *);
    for (int i = 0; i < N; i++)
    {
        size += dynamics[i]->opts_calculate_size(dynamics[i], dims->dynamics[i]);
    }

    // cost
    size += (N + 1) * sizeof(void *);
    for (int i = 0; i <= N; i++)
    {
        size += cost[i]->opts_calculate_size(cost[i], dims->cost[i]);
    }

    // constraints
    size += (N + 1) * sizeof(void *);
    for (int i = 0; i <= N; i++)
    {
        size += constraints[i]->opts_calculate_size(constraints[i], dims->constraints[i]);
    }

    size += 2*8;  // 2 aligns

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
    align_char_to(8, &c_ptr);

    ocp_nlp_opts *opts = (ocp_nlp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_opts);

    /* pointers to substructures */
    opts->dynamics = (void **) c_ptr;
    c_ptr += N * sizeof(void *);

    opts->cost = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    opts->constraints = (void **) c_ptr;
    c_ptr += (N + 1) * sizeof(void *);

    align_char_to(8, &c_ptr);

    /* substructures */
    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    opts->regularize = config->regularize->opts_assign(c_ptr);
    c_ptr += config->regularize->opts_calculate_size();

    // dynamics
    for (int i = 0; i < N; i++)
    {
        opts->dynamics[i] = dynamics[i]->opts_assign(dynamics[i], dims->dynamics[i], c_ptr);
        c_ptr += dynamics[i]->opts_calculate_size(dynamics[i], dims->dynamics[i]);
    }

    // cost
    for (int i = 0; i <= N; i++)
    {
        opts->cost[i] = cost[i]->opts_assign(cost[i], dims->cost[i], c_ptr);
        c_ptr += cost[i]->opts_calculate_size(cost[i], dims->cost[i]);
    }

    // constraints
    for (int i = 0; i <= N; i++)
    {
        opts->constraints[i] =
            constraints[i]->opts_assign(constraints[i], dims->constraints[i], c_ptr);
        c_ptr += constraints[i]->opts_calculate_size(constraints[i], dims->constraints[i]);
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

    int N = dims->N;

    opts->reuse_workspace = 1;
#if defined(ACADOS_WITH_OPENMP)
    #if defined(ACADOS_NUM_THREADS)
    opts->num_threads = ACADOS_NUM_THREADS;
    // printf("\nocp_nlp: openmp threads from macro = %d\n", opts->num_threads);
    #else
    opts->num_threads = omp_get_max_threads();
    // printf("\nocp_nlp: omp_get_max_threads %d", omp_get_max_threads());
    #endif
#endif
    // printf("\nocp_nlp: openmp threads = %d\n", opts->num_threads);

    opts->globalization = FIXED_STEP;
    opts->print_level = 0;
    opts->step_length = 1.0;
    opts->levenberg_marquardt = 0.0;


    /* submodules opts */
    // qp solver
    qp_solver->opts_initialize_default(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // regularization
    regularize->opts_initialize_default(regularize, dims->regularize, opts->regularize);

    // dynamics
    for (int i = 0; i < N; i++)
    {
        dynamics[i]->opts_initialize_default(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
    }

    // cost
    for (int i = 0; i <= N; i++)
    {
        cost[i]->opts_initialize_default(cost[i], dims->cost[i], opts->cost[i]);
    }

    // constraints
    for (int i = 0; i <= N; i++)
    {
        constraints[i]->opts_initialize_default(constraints[i], dims->constraints[i], opts->constraints[i]);
    }

    // globalization
    opts->alpha_min = 0.05;
    opts->alpha_reduction = 0.7;
    opts->full_step_dual = 0;
    opts->line_search_use_sufficient_descent = 0;
    opts->globalization_use_SOC = 0;
    opts->eps_sufficient_descent = 1e-4; // Leineweber1999: MUSCOD-I eps_T = 1e-4 (p.89); Note: eps_T = 0.1 originally proposed by Powell 1978 (Leineweber 1999, p. 53)

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

    int N = dims->N;

    qp_solver->opts_update(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // dynamics
    for (int i = 0; i < N; i++)
    {
        dynamics[i]->opts_update(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
    }

    // cost
    for (int i = 0; i <= N; i++)
    {
        cost[i]->opts_update(cost[i], dims->cost[i], opts->cost[i]);
    }

    // constraints
    for (int i = 0; i <= N; i++)
    {
        constraints[i]->opts_update(constraints[i], dims->constraints[i], opts->constraints[i]);
    }

    return;
}



void ocp_nlp_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_opts *opts = (ocp_nlp_opts *) opts_;
    ocp_nlp_config *config = config_;

    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name, i.e. substring in field before '_'
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (int i=0; i<module_length; i++)
            module[i] = field[i];
        module[module_length] = '\0'; // add end of string
        ptr_module = module;
    }

    // pass options to QP module
    if ( ptr_module!=NULL && (!strcmp(ptr_module, "qp")) )
    {
        config->qp_solver->opts_set(config->qp_solver, opts->qp_solver_opts,
                                    field+module_length+1, value);
    }
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
        else if (!strcmp(field, "alpha_reduction"))
        {
            double* alpha_reduction = (double *) value;
            opts->alpha_reduction = *alpha_reduction;
        }
        else if (!strcmp(field, "alpha_min"))
        {
            double* alpha_min = (double *) value;
            opts->alpha_min = *alpha_min;
        }
        else if (!strcmp(field, "eps_sufficient_descent"))
        {
            double* eps_sufficient_descent = (double *) value;
            opts->eps_sufficient_descent = *eps_sufficient_descent;
        }
        else if (!strcmp(field, "full_step_dual"))
        {
            int* full_step_dual = (int *) value;
            opts->full_step_dual = *full_step_dual;
        }
        else if (!strcmp(field, "line_search_use_sufficient_descent"))
        {
            int* line_search_use_sufficient_descent = (int *) value;
            opts->line_search_use_sufficient_descent = *line_search_use_sufficient_descent;
        }
        else if (!strcmp(field, "globalization_use_SOC"))
        {
            int* globalization_use_SOC = (int *) value;
            opts->globalization_use_SOC = *globalization_use_SOC;
        }
        else if (!strcmp(field, "globalization"))
        {
            char* globalization = (char *) value;
            if (!strcmp(globalization, "fixed_step"))
            {
                opts->globalization = FIXED_STEP;
            }
            else if (!strcmp(globalization, "merit_backtracking"))
            {
                opts->globalization = MERIT_BACKTRACKING;
            }
            else
            {
                printf("\nerror: ocp_nlp_opts_set: not supported value for globalization, got: %s\n",
                       globalization);
                exit(1);
            }
        }
        else if (!strcmp(field, "levenberg_marquardt"))
        {
            double* levenberg_marquardt = (double *) value;
            opts->levenberg_marquardt = *levenberg_marquardt;
        }
        else if (!strcmp(field, "exact_hess"))
        {
            int N = config->N;
            // cost
            for (int i=0; i<=N; i++)
                config->cost[i]->opts_set(config->cost[i], opts->cost[i], "exact_hess", value);
            // dynamics
            for (int i=0; i<N; i++)
                config->dynamics[i]->opts_set(config->dynamics[i], opts->dynamics[i],
                                               "compute_hess", value);
            // constraints
            for (int i=0; i<=N; i++)
                config->constraints[i]->opts_set(config->constraints[i], opts->constraints[i],
                                                  "compute_hess", value);
        }
        // selectively turn on exact hessian contributions
        else if (!strcmp(field, "exact_hess_cost"))
        {
            int N = config->N;
            for (int i=0; i<=N; i++)
                config->cost[i]->opts_set(config->cost[i], opts->cost[i], "exact_hess", value);
        }
        else if (!strcmp(field, "exact_hess_dyn"))
        {
            int N = config->N;
            for (int i=0; i<N; i++)
                config->dynamics[i]->opts_set(config->dynamics[i], opts->dynamics[i],
                                               "compute_hess", value);
        }
        else if (!strcmp(field, "exact_hess_constr"))
        {
            int N = config->N;
            for (int i=0; i<=N; i++)
                config->constraints[i]->opts_set(config->constraints[i], opts->constraints[i],
                                                  "compute_hess", value);
        }
        else if (!strcmp(field, "print_level"))
        {
            int* print_level = (int *) value;
            if (*print_level < 0)
            {
                printf("\nerror: ocp_nlp_opts_set: invalid value for print_level field, need int >=0, got %d.\n", *print_level);
                exit(1);
            }
            opts->print_level = *print_level;
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

    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (int i=0; i<module_length; i++)
            module[i] = field[i];
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

acados_size_t ocp_nlp_memory_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts)
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

    acados_size_t size = sizeof(ocp_nlp_memory);

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
    for (int i = 0; i < N; i++)
    {
        size += dynamics[i]->memory_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
    }

    // cost
    size += (N + 1) * sizeof(void *);
    for (int i = 0; i <= N; i++)
    {
        size += cost[i]->memory_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
    }

    // constraints
    size += (N + 1) * sizeof(void *);
    for (int i = 0; i <= N; i++)
    {
        size += constraints[i]->memory_calculate_size(constraints[i], dims->constraints[i], opts->constraints[i]);
    }

    // nlp res
    size += ocp_nlp_res_calculate_size(dims);

    size += (N+1)*sizeof(bool); // set_sim_guess

    size += (N+1)*sizeof(struct blasfeo_dmat); // dzduxt
    size += 6*(N+1)*sizeof(struct blasfeo_dvec);  // cost_grad ineq_fun ineq_adj dyn_adj sim_guess z_alg
    size += 1*N*sizeof(struct blasfeo_dvec);        // dyn_fun

    for (int i = 0; i < N; i++)
    {
        size += 1*blasfeo_memsize_dmat(nu[i]+nx[i], nz[i]); // dzduxt
        size += 1*blasfeo_memsize_dvec(nz[i]); // z_alg
        size += 2*blasfeo_memsize_dvec(nv[i]);           // cost_grad ineq_adj
        size += 1*blasfeo_memsize_dvec(nu[i] + nx[i]);  // dyn_adj
        size += 1*blasfeo_memsize_dvec(nx[i + 1]);       // dyn_fun
        size += 1*blasfeo_memsize_dvec(2 * ni[i]);       // ineq_fun
        size += 1*blasfeo_memsize_dvec(nx[i] + nz[i]); // sim_guess
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

    /* pointers to substructures */
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

    /* substructures */
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
    for (int i = 0; i < N; i++)
    {
        mem->dynamics[i] = dynamics[i]->memory_assign(dynamics[i], dims->dynamics[i], opts->dynamics[i], c_ptr);
        c_ptr += dynamics[i]->memory_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
    }

    // cost
    for (int i = 0; i <= N; i++)
    {
        mem->cost[i] = cost[i]->memory_assign(cost[i], dims->cost[i], opts->cost[i], c_ptr);
        c_ptr += cost[i]->memory_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
    }

    // constraints
    for (int i = 0; i <= N; i++)
    {
        mem->constraints[i] = constraints[i]->memory_assign(constraints[i],
                                            dims->constraints[i], opts->constraints[i], c_ptr);
        c_ptr += constraints[i]->memory_calculate_size( constraints[i], dims->constraints[i],
                                                                 opts->constraints[i]);
    }

    // nlp res
    mem->nlp_res = ocp_nlp_res_assign(dims, c_ptr);
    c_ptr += mem->nlp_res->memsize;

    // blasfeo_struct align
    align_char_to(8, &c_ptr);

    // dzduxt
    assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->dzduxt, &c_ptr);

    // z_alg
    assign_and_advance_blasfeo_dvec_structs(N + 1, &mem->z_alg, &c_ptr);
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

    // set_sim_guess
    assign_and_advance_bool(N+1, &mem->set_sim_guess, &c_ptr);
    for (int i = 0; i <= N; ++i)
    {
        mem->set_sim_guess[i] = false;
    }

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // dzduxt
    for (int i=0; i<=N; i++)
    {
        assign_and_advance_blasfeo_dmat_mem(nu[i]+nx[i], nz[i], mem->dzduxt+i, &c_ptr);
    }
    // z_alg
    for (int i=0; i<=N; i++)
    {
        blasfeo_create_dvec(nz[i], mem->z_alg+i, c_ptr);
        c_ptr += blasfeo_memsize_dvec(nz[i]);
    }

    // cost_grad
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[i], mem->cost_grad + i, &c_ptr);
    }
    // ineq_fun
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[i], mem->ineq_fun + i, &c_ptr);
    }
    // ineq_adj
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[i], mem->ineq_adj + i, &c_ptr);
    }
    // dyn_fun
    for (int i = 0; i < N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[i + 1], mem->dyn_fun + i, &c_ptr);
    }
    // dyn_adj
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nu[i] + nx[i], mem->dyn_adj + i, &c_ptr);
    }
    // sim_guess
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[i] + nz[i], mem->sim_guess + i, &c_ptr);
        // set to 0;
        blasfeo_dvecse(nx[i] + nz[i], 0.0, mem->sim_guess+i, 0);
        // printf("sim_guess i %d: %p\n", i, mem->sim_guess+i);
    }
    // printf("created memory %p\n", mem);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_workspace_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;
    // int *nz = dims->nz;

    acados_size_t size = 0;

    // nlp
    size += sizeof(ocp_nlp_workspace);

    // tmp_nlp_out
    size += ocp_nlp_out_calculate_size(config, dims);

    // weight_merit_fun
    size += ocp_nlp_out_calculate_size(config, dims);

    // blasfeo_dvec
    int nxu_max = 0;
    int nx_max = 0;
    int ni_max = 0;
    for (int i = 0; i <= N; i++)
    {
        nx_max = nx_max > nx[i] ? nx_max : nx[i];
        nxu_max = nxu_max > (nx[i]+nu[i]) ? nxu_max : (nx[i]+nu[i]);
        ni_max = ni_max > ni[i] ? ni_max : ni[i];
    }
    size += 1 * blasfeo_memsize_dvec(nx_max);
    size += 1 * blasfeo_memsize_dvec(nxu_max);
    size += 1 * blasfeo_memsize_dvec(ni_max);

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
        for (int i = 0; i < N; i++)
        {
            size += dynamics[i]->workspace_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
        }

        // cost
        for (int i = 0; i <= N; i++)
        {
            size += cost[i]->workspace_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
        }

        // constraints
        for (int i = 0; i <= N; i++)
        {
            size += constraints[i]->workspace_calculate_size(constraints[i], dims->constraints[i], opts->constraints[i]);
        }

#else
        acados_size_t size_tmp = 0;
        int tmp;

        // qp solver
        tmp = qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);
        size_tmp = tmp > size_tmp ? tmp : size_tmp;

        // dynamics
        for (int i = 0; i < N; i++)
        {
            tmp = dynamics[i]->workspace_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // cost
        for (int i = 0; i <= N; i++)
        {
            tmp = cost[i]->workspace_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // constraints
        for (int i = 0; i <= N; i++)
        {
            tmp = constraints[i]->workspace_calculate_size(constraints[i], dims->constraints[i], opts->constraints[i]);
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
        for (int i = 0; i < N; i++)
        {
            size += dynamics[i]->workspace_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
        }

        // cost
        for (int i = 0; i <= N; i++)
        {
            size += cost[i]->workspace_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
        }

        // constraints
        for (int i = 0; i <= N; i++)
        {
            size += constraints[i]->workspace_calculate_size(constraints[i], dims->constraints[i], opts->constraints[i]);
        }

    }

    size += 8; // struct align
    return size;
}



ocp_nlp_workspace *ocp_nlp_workspace_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                             ocp_nlp_opts *opts, ocp_nlp_memory *mem, void *raw_memory)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;
    int *nx = dims->nx;
    // int *nv = dims->nv;
    int *nu = dims->nu;
    int *ni = dims->ni;
    // int *nz = dims->nz;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_workspace *work = (ocp_nlp_workspace *) c_ptr;
    c_ptr += sizeof(ocp_nlp_workspace);

    /* pointers to substructures */
    //
    work->dynamics = (void **) c_ptr;
    c_ptr += N*sizeof(void *);
    //
    work->cost = (void **) c_ptr;
    c_ptr += (N+1)*sizeof(void *);
    //
    work->constraints = (void **) c_ptr;
    c_ptr += (N+1)*sizeof(void *);

    align_char_to(8, &c_ptr);

    /* substructures */
    // tmp_nlp_out
    work->tmp_nlp_out = ocp_nlp_out_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_out_calculate_size(config, dims);

    // weight_merit_fun
    work->weight_merit_fun = ocp_nlp_out_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_out_calculate_size(config, dims);

    // blasfeo_dvec
    int nxu_max = 0;
    int nx_max = 0;
    int ni_max = 0;
    for (int i = 0; i <= N; i++)
    {
        nx_max = nx_max > nx[i] ? nx_max : nx[i];
        nxu_max = nxu_max > (nx[i]+nu[i]) ? nxu_max : (nx[i]+nu[i]);
        ni_max = ni_max > ni[i] ? ni_max : ni[i];
    }
    assign_and_advance_blasfeo_dvec_mem(nxu_max, &work->tmp_nxu, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(ni_max, &work->tmp_ni, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx_max, &work->dxnext_dy, &c_ptr);

    if (opts->reuse_workspace)
    {

#if defined(ACADOS_WITH_OPENMP)

        // qp solver
        work->qp_work = (void *) c_ptr;
        c_ptr += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);

        // dynamics
        for (int i = 0; i < N; i++)
        {
            work->dynamics[i] = c_ptr;
            c_ptr += dynamics[i]->workspace_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
        }

        // cost
        for (int i = 0; i <= N; i++)
        {
            work->cost[i] = c_ptr;
            c_ptr += cost[i]->workspace_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
        }

        // constraints
        for (int i = 0; i <= N; i++)
        {
            work->constraints[i] = c_ptr;
            c_ptr += constraints[i]->workspace_calculate_size(constraints[i], dims->constraints[i], opts->constraints[i]);
        }

#else

        acados_size_t size_tmp = 0;
        int tmp;

        // qp solver
        work->qp_work = (void *) c_ptr;
        tmp = qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, opts->qp_solver_opts);
        size_tmp = tmp > size_tmp ? tmp : size_tmp;

        // dynamics
        for (int i = 0; i < N; i++)
        {
            work->dynamics[i] = c_ptr;
            tmp = dynamics[i]->workspace_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // cost
        for (int i = 0; i <= N; i++)
        {
            work->cost[i] = c_ptr;
            tmp = cost[i]->workspace_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
            size_tmp = tmp > size_tmp ? tmp : size_tmp;
        }

        // constraints
        for (int i = 0; i <= N; i++)
        {
            work->constraints[i] = c_ptr;
            tmp = constraints[i]->workspace_calculate_size(constraints[i], dims->constraints[i], opts->constraints[i]);
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
        for (int i = 0; i < N; i++)
        {
            work->dynamics[i] = c_ptr;
            c_ptr += dynamics[i]->workspace_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
        }

        // cost
        for (int i = 0; i <= N; i++)
        {
            work->cost[i] = c_ptr;
            c_ptr += cost[i]->workspace_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
        }

        // constraints
        for (int i = 0; i <= N; i++)
        {
            work->constraints[i] = c_ptr;
            c_ptr += constraints[i]->workspace_calculate_size(constraints[i], dims->constraints[i], opts->constraints[i]);
        }

    }

    assert((char *) work + ocp_nlp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return work;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_alias_memory_to_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
         ocp_nlp_out *nlp_out, ocp_nlp_opts *opts, ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work)
{
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel
    { // beginning of parallel region
#endif

    int N = dims->N;

    // alias to dynamics_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (int i = 0; i < N; i++)
    {
        config->dynamics[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_ux1_ptr(nlp_out->ux+i+1, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_tmp_ux1_ptr(nlp_work->tmp_nlp_out->ux+i+1, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_pi_ptr(nlp_out->pi+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_tmp_pi_ptr(nlp_work->tmp_nlp_out->pi+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_BAbt_ptr(nlp_mem->qp_in->BAbt+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_dzduxt_ptr(nlp_mem->dzduxt+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_sim_guess_ptr(nlp_mem->sim_guess+i, nlp_mem->set_sim_guess+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_z_alg_ptr(nlp_mem->z_alg+i, nlp_mem->dynamics[i]);
    }

    // alias to cost_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (int i = 0; i <= N; i++)
    {
        config->cost[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->cost[i]);
        config->cost[i]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+i, nlp_mem->cost[i]);
        config->cost[i]->memory_set_z_alg_ptr(nlp_mem->z_alg+i, nlp_mem->cost[i]);
        config->cost[i]->memory_set_dzdux_tran_ptr(nlp_mem->dzduxt+i, nlp_mem->cost[i]);
        config->cost[i]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+i, nlp_mem->cost[i]);
        config->cost[i]->memory_set_Z_ptr(nlp_mem->qp_in->Z+i, nlp_mem->cost[i]);
    }

    // alias to constraints_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (int i = 0; i <= N; i++)
    {
        config->constraints[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_lam_ptr(nlp_out->lam+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_tmp_lam_ptr(nlp_work->tmp_nlp_out->lam+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_z_alg_ptr(nlp_mem->z_alg+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_dzdux_tran_ptr(nlp_mem->dzduxt+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_DCt_ptr(nlp_mem->qp_in->DCt+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_idxb_ptr(nlp_mem->qp_in->idxb[i], nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_idxs_rev_ptr(nlp_mem->qp_in->idxs_rev[i], nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_idxe_ptr(nlp_mem->qp_in->idxe[i], nlp_mem->constraints[i]);
    }

    // alias to regularize memory
    config->regularize->memory_set_RSQrq_ptr(dims->regularize, nlp_mem->qp_in->RSQrq, nlp_mem->regularize_mem);
    config->regularize->memory_set_rq_ptr(dims->regularize, nlp_mem->qp_in->rqz, nlp_mem->regularize_mem);
    config->regularize->memory_set_BAbt_ptr(dims->regularize, nlp_mem->qp_in->BAbt, nlp_mem->regularize_mem);
    config->regularize->memory_set_b_ptr(dims->regularize, nlp_mem->qp_in->b, nlp_mem->regularize_mem);
    config->regularize->memory_set_idxb_ptr(dims->regularize, nlp_mem->qp_in->idxb, nlp_mem->regularize_mem);
    config->regularize->memory_set_DCt_ptr(dims->regularize, nlp_mem->qp_in->DCt, nlp_mem->regularize_mem);
    config->regularize->memory_set_ux_ptr(dims->regularize, nlp_mem->qp_out->ux, nlp_mem->regularize_mem);
    config->regularize->memory_set_pi_ptr(dims->regularize, nlp_mem->qp_out->pi, nlp_mem->regularize_mem);
    config->regularize->memory_set_lam_ptr(dims->regularize, nlp_mem->qp_out->lam, nlp_mem->regularize_mem);

    // copy sampling times into dynamics model
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    // NOTE(oj): this will lead in an error for irk_gnsf, T must be set in precompute;
    //    -> remove here and make sure precompute is called everywhere (e.g. Python interface).
    for (int i = 0; i < N; i++)
    {
        config->dynamics[i]->model_set(config->dynamics[i], dims->dynamics[i],
                                         nlp_in->dynamics[i], "T", nlp_in->Ts+i);
    }

#if defined(ACADOS_WITH_OPENMP)
    } // end of parallel region
#endif

    return;
}


void ocp_nlp_initialize_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
         ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int N = dims->N;

    // NOTE: initialize is called at the start of every NLP solver call.
    // It computes things in submodules based on stuff that can be changed by the user between
    // subsequent solver calls, e.g. factorization of weight matrix.
    // IN CONTRAST: precompute is only called once after solver creation
    //  -> computes things that are not expected to change between subsequent solver calls
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // cost
        config->cost[i]->initialize(config->cost[i], dims->cost[i], in->cost[i],
                opts->cost[i], mem->cost[i], work->cost[i]);
        // dynamics
        if (i < N)
            config->dynamics[i]->initialize(config->dynamics[i], dims->dynamics[i],
                    in->dynamics[i], opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
        // constraints
        config->constraints[i]->initialize(config->constraints[i], dims->constraints[i],
                in->constraints[i], opts->constraints[i], mem->constraints[i], work->constraints[i]);
    }

    return;
}


void ocp_nlp_initialize_t_slacks(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
         ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    struct blasfeo_dvec *ineq_fun;
    int N = dims->N;
    int *ni = dims->ni;
    int *ns = dims->ns;
    int *nx = dims->nx;
    int *nu = dims->nu;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // copy out->ux to tmp_nlp_out->ux, since this is used in compute_fun
        blasfeo_dveccp(nx[i]+nu[i]+2*ns[i], out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);

        // evaluate inequalities
        config->constraints[i]->compute_fun(config->constraints[i], dims->constraints[i],
                                             in->constraints[i], opts->constraints[i],
                                             mem->constraints[i], work->constraints[i]);
        ineq_fun = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        // t = -ineq_fun
        blasfeo_dveccpsc(2 * ni[i], -1.0, ineq_fun, 0, out->t + i, 0);
    }

    return;
}



void ocp_nlp_approximate_qp_matrices(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem,
    ocp_nlp_workspace *work)
{
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    /* stage-wise multiple shooting lagrangian evaluation */

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // init Hessian to 0 
        blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, mem->qp_in->RSQrq+i, 0, 0);


        if (i < N)
        {
            // Levenberg Marquardt term: Ts[i] * levenberg_marquardt * eye()
            if (opts->levenberg_marquardt > 0.0)
                blasfeo_ddiare(nu[i] + nx[i], in->Ts[i] * opts->levenberg_marquardt,
                               mem->qp_in->RSQrq+i, 0, 0);

            // dynamics
            config->dynamics[i]->update_qp_matrices(config->dynamics[i], dims->dynamics[i],
                    in->dynamics[i], opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
        }
        else
        {
            // Levenberg Marquardt term: 1.0 * levenberg_marquardt * eye()
            if (opts->levenberg_marquardt > 0.0)
                blasfeo_ddiare(nu[i] + nx[i], opts->levenberg_marquardt,
                               mem->qp_in->RSQrq+i, 0, 0);
        }

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
    for (int i=0; i <= N; i++)
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

    for (int i = 0; i <= N; i++)
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
}



// update QP rhs for SQP (step prim var, abs dual var)
// TODO(all): move in dynamics, cost, constraints modules ???
void ocp_nlp_approximate_qp_vectors_sqp(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // g
        blasfeo_dveccp(nv[i], mem->cost_grad + i, 0, mem->qp_in->rqz + i, 0);

        // b
        if (i < N)
            blasfeo_dveccp(nx[i + 1], mem->dyn_fun + i, 0, mem->qp_in->b + i, 0);

        // d
        blasfeo_dveccp(2 * ni[i], mem->ineq_fun + i, 0, mem->qp_in->d + i, 0);
    }
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
}



double ocp_nlp_compute_merit_gradient(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    /* computes merit function gradient at iterate: out -- using already evaluated gradients of submodules
       with weights: work->weight_merit_fun */
    int i, j;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    double merit_grad = 0.0;
    double weight;

    // NOTE: step is in: mem->qp_out->ux
    struct blasfeo_dvec *tmp_vec; // size nv
    struct blasfeo_dvec tmp_vec_nxu = work->tmp_nxu;  // size nxu
    struct blasfeo_dvec dxnext_dy = work->dxnext_dy;  // size nx

    // cost
    for (i=0; i<=N; i++)
    {
        tmp_vec = config->cost[i]->memory_get_grad_ptr(mem->cost[i]);
        merit_grad += blasfeo_ddot(nv[i], tmp_vec, 0, mem->qp_out->ux + i, 0);
    }
    double merit_grad_cost = merit_grad;

    /* dynamics */
    double merit_grad_dyn = 0.0;
    for (i=0; i<N; i++)
    {
        // get shooting node gap x_next(x_n, u_n) - x_{n+1};
        tmp_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);

        /* compute directional derivative of xnext with direction y -> dxnext_dy */
        blasfeo_dgemv_t(nx[i]+nu[i], nx[i+1], 1.0, mem->qp_in->BAbt+i, 0, 0, mem->qp_out->ux+i, 0,
                        0.0, &dxnext_dy, 0, &dxnext_dy, 0);

        /* add merit gradient contributions depending on sign of shooting gap */
        for (j = 0; j < nx[i+1]; j++)
        {
            weight = BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j);
            double deqj_dy = BLASFEO_DVECEL(&dxnext_dy, j) - BLASFEO_DVECEL(mem->qp_out->ux+(i+1), nu[i+1]+j);
            {
                if (BLASFEO_DVECEL(tmp_vec, j) > 0)
                {
                    merit_grad_dyn += weight * deqj_dy;
                    // printf("\ndyn_contribution +%e, weight %e, deqj_dy %e, i %d, j %d", weight * deqj_dy, weight, deqj_dy, i, j);
                }
                else
                {
                    merit_grad_dyn -= weight * deqj_dy;
                    // printf("\ndyn_contribution %e, weight %e, deqj_dy %e, i %d, j %d", -weight * deqj_dy, weight, deqj_dy, i, j);
                }
            }
        }
    }

    /* inequality contributions */
    // NOTE: slack bound inequalities are not considered here.
    // They should never be infeasible. Only if explicitly initialized infeasible from outside.
    int constr_index, slack_index_in_ux, slack_index;
    ocp_qp_dims* qp_dims = mem->qp_in->dim;
    int *nb = qp_dims->nb;
    int *ng = qp_dims->ng;
    int *ns = qp_dims->ns;
    double merit_grad_ineq = 0.0;
    double slack_step;

    for (i=0; i<=N; i++)
    {
        tmp_vec = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        int *idxb = mem->qp_in->idxb[i];
        if (ni[i] > 0)
        {
            // NOTE: loop could be simplified handling lower and upper constraints together.
            for (j = 0; j < 2 * (nb[i] + ng[i]); j++) // 2 * ni
            {
                double constraint_val = BLASFEO_DVECEL(tmp_vec, j);
                if (constraint_val > 0)
                {
                    weight = BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j);

                    // find corresponding slack value
                    constr_index = j < nb[i]+ng[i] ? j : j-(nb[i]+ng[i]);
                    slack_index = mem->qp_in->idxs_rev[i][constr_index];
                    // if softened: add slack contribution
                    if (slack_index >= 0)
                    {
                        slack_index_in_ux = j < (nb[i]+ng[i]) ? nx[i] + nu[i] + slack_index
                                                              : nx[i] + nu[i] + slack_index + ns[i];
                        slack_step = BLASFEO_DVECEL(mem->qp_out->ux+i, slack_index_in_ux);
                        merit_grad_ineq -= weight * slack_step;
                        // printf("at node %d, ineq %d, idxs_rev[%d] = %d\n", i, j, constr_index, slack_index);
                        // printf("slack contribution: uxs[%d] = %e\n", slack_index_in_ux, slack_step);
                    }


                    // NOTE: the inequalities are internally organized in the following order:
                    //     [ lbu lbx lg lh lphi ubu ubx ug uh uphi;
                    //     lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]
                    // printf("constraint %d %d is active with value %e", i, j, constraint_val);
                    if (j < nb[i])
                    {
                        // printf("lower idxb[%d] = %d dir %f, constraint_val %f, nb = %d\n", j, idxb[j], BLASFEO_DVECEL(mem->qp_out->ux, idxb[j]), constraint_val, nb[i]);
                        merit_grad_ineq += weight * BLASFEO_DVECEL(mem->qp_out->ux+i, idxb[j]);
                    }
                    else if (j < nb[i] + ng[i])
                    {
                        // merit_grad_ineq += weight * mem->qp_in->DCt_j * dux
                        blasfeo_dcolex(nx[i] + nu[i], mem->qp_in->DCt+i, j - nb[i], 0, &tmp_vec_nxu, 0);
                        merit_grad_ineq += weight * blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0);
                        // printf("general linear constraint lower contribution = %e, val = %e\n", blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0), constraint_val);
                    }
                    else if (j < 2*nb[i] + ng[i])
                    {
                        // printf("upper idxb[%d] = %d dir %f, constraint_val %f, nb = %d\n", j-nb[i]-ng[i], idxb[j-nb[i]-ng[i]], BLASFEO_DVECEL(mem->qp_out->ux, idxb[j-nb[i]-ng[i]]), constraint_val, nb[i]);
                        merit_grad_ineq += weight * BLASFEO_DVECEL(mem->qp_out->ux+i, idxb[j-nb[i]-ng[i]]);
                    }
                    else if (j < 2*nb[i] + 2*ng[i])
                    {
                        blasfeo_dcolex(nx[i] + nu[i], mem->qp_in->DCt+i, j - 2*nb[i] - ng[i], 0, &tmp_vec_nxu, 0);
                        merit_grad_ineq += weight * blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0);
                        // printf("general linear constraint upper contribution = %e, val = %e\n", blasfeo_ddot(nx[i] + nu[i], &tmp_vec_nxu, 0, mem->qp_out->ux+i, 0), constraint_val);
                    }
                }
            }
        }
    }
    // print_ocp_qp_dims(qp_dims);
    // print_ocp_qp_in(mem->qp_in);

    merit_grad = merit_grad_cost + merit_grad_dyn + merit_grad_ineq;
    if (opts->print_level > 1)
        printf("computed merit_grad = %e, merit_grad_cost = %e, merit_grad_dyn = %e, merit_grad_ineq = %e\n", merit_grad, merit_grad_cost, merit_grad_dyn, merit_grad_ineq);

    return merit_grad;
}



static double ocp_nlp_get_violation(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    // computes constraint violation infinity norm
    // assumes constraint functions are evaluated before, e.g. done in ocp_nlp_evaluate_merit_fun
    int i, j;
    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    struct blasfeo_dvec *tmp_fun_vec;
    double violation = 0.0;
    double tmp;
    for (i=0; i<N; i++)
    {
        tmp_fun_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        for (j=0; j<nx[i+1]; j++)
        {
            tmp = fabs(BLASFEO_DVECEL(tmp_fun_vec, j));
            violation = tmp > violation ? tmp : violation;
        }
    }

    for (i=0; i<=N; i++)
    {
        tmp_fun_vec = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        for (j=0; j<2*ni[i]; j++)
        {
            // Note constraint violation corresponds to > 0
            tmp = BLASFEO_DVECEL(tmp_fun_vec, j);
            violation = tmp > violation ? tmp : violation;
        }
    }

    return violation;
}



double ocp_nlp_evaluate_merit_fun(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                  ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
                                  ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    /* computes merit function value at iterate: tmp_nlp_out, with weights: work->weight_merit_fun */
    //int j;

    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;

    double merit_fun = 0.0;

    // compute fun value
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<=N; i++)
    {
        // cost
        config->cost[i]->compute_fun(config->cost[i], dims->cost[i], in->cost[i], opts->cost[i],
                                    mem->cost[i], work->cost[i]);
    }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<N; i++)
    {
        // dynamics
        config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i], in->dynamics[i],
                                         opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
    }
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<=N; i++)
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
    for(int i=0; i<=N; i++)
    {
        tmp_fun = config->cost[i]->memory_get_fun_ptr(mem->cost[i]);
        cost_fun += *tmp_fun;
    }

    double dyn_fun = 0.0;
    for(int i=0; i<N; i++)
    {
        tmp_fun_vec = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        // printf("\nMerit: dyn will multiply tmp_fun, weights %d\n", i);
        // blasfeo_print_exp_tran_dvec(nx[i+1], tmp_fun_vec, 0);
        // blasfeo_print_exp_tran_dvec(nx[i+1], work->weight_merit_fun->pi+i, 0);
        for(int j=0; j<nx[i+1]; j++)
        {
//            printf("\n%e %e\n", fabs(BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j)), fabs(BLASFEO_DVECEL(tmp_fun_vec, j)));
            dyn_fun += fabs(BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j)) * fabs(BLASFEO_DVECEL(tmp_fun_vec, j));
        }
    }

    double constr_fun = 0.0;
    for(int i=0; i<=N; i++)
    {
//        printf("\ni %d\n", i);
        tmp_fun_vec = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
//        blasfeo_print_exp_tran_dvec(2*ni[i], tmp_fun_vec, 0);
//        blasfeo_print_exp_tran_dvec(2*ni[i], work->weight_merit_fun->lam+i, 0);
        for (int j=0; j<2*ni[i]; j++)
        {
            tmp = BLASFEO_DVECEL(tmp_fun_vec, j);
            if (tmp > 0.0)
            {
                // tmp = constraint violation
                // printf("IN merit fun: ineq i %d, j %d tmp_fun %e, multiplier %e\n", i, j, tmp, BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j));
                constr_fun += fabs(BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j)) * tmp;
            }
        }
    }

    merit_fun = cost_fun + dyn_fun + constr_fun;

	// printf("Merit fun: %e cost: %e dyn: %e constr: %e\n", merit_fun, cost_fun, dyn_fun, constr_fun);

    return merit_fun;
}



double ocp_nlp_line_search(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
            int check_early_termination)
{
    int i, j;

    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *ni = dims->ni;

    double alpha = opts->step_length;
    double tmp0, tmp1, merit_fun1;
    ocp_qp_out *qp_out = mem->qp_out;

    // Line search version Jonathan
    // Following Leineweber1999, Section "3.5.1 Line Search Globalization"
    // TODO: check out more advanced step search Leineweber1995

    if (opts->globalization == MERIT_BACKTRACKING)
    {
        // copy out (current iterate) to work->tmp_nlp_out
        for (i = 0; i <= N; i++)
            blasfeo_dveccp(nv[i], out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);

        // for (i = 0; i < N; i++)
        //     blasfeo_dveccp(nx[i+1], out->pi+i, 0, work->tmp_nlp_out->pi+i, 0);

        // for (i = 0; i <= N; i++)
        //     blasfeo_dveccp(2*ni[i], out->lam+i, 0, work->tmp_nlp_out->lam+i, 0);

            // linear update of algebraic variables using state and input sensitivity
    //        if (i < N)
    //        {
    //            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, mem->dzduxt+i, 0, 0, mem->qp_out->ux+i, 0, 1.0, mem->z_alg+i, 0, out->z+i, 0);
    //        }

        /* modify/initialize merit function weights (Leineweber1999 M5.1, p.89) */
        if (mem->sqp_iter[0]==0)
        {
            // initialize weights
            // equality merit weights = abs( eq multipliers of qp_sol )
            for (i = 0; i < N; i++)
            {
                for (j=0; j<nx[i+1]; j++)
                {
                    // tmp0 = fabs(BLASFEO_DVECEL(out->pi+i, j));
                    tmp0 = fabs(BLASFEO_DVECEL(qp_out->pi+i, j));
                }
            }

            for (i = 0; i <= N; i++)
            {
                blasfeo_dveccp(2*ni[i], qp_out->lam+i, 0, work->weight_merit_fun->lam+i, 0);
            }
        }
        else
        {
            // update weights
            // printf("merit fun: update weights, sqp_iter = %d\n", mem->sqp_iter[0]);
            for (i = 0; i < N; i++)
            {
                for(j=0; j<nx[i+1]; j++)
                {
                    // abs(lambda) (LW)
                    tmp0 = fabs(BLASFEO_DVECEL(qp_out->pi+i, j));
                    // .5 * (abs(lambda) + sigma)
                    tmp1 = 0.5 * (tmp0 + BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j));
                    BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j) = tmp0 > tmp1 ? tmp0 : tmp1;
                }
            }
            for (i = 0; i <= N; i++)
            {
                for(j=0; j<2*ni[i]; j++)
                {
                    // mu (LW)
                    tmp0 = BLASFEO_DVECEL(qp_out->lam+i, j);
                    // .5 * (mu + tau)
                    tmp1 = 0.5 * (tmp0 + BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j));
                    BLASFEO_DVECEL(work->weight_merit_fun->lam+i, j) = tmp0>tmp1 ? tmp0 : tmp1;
                }
            }
        }

        if (1) // (mem->sqp_iter[0]!=0) // TODO: why does Leineweber do full step in first SQP iter?
        {
            double merit_fun0 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);

            double reduction_factor = opts->alpha_reduction;
            double max_next_merit_fun_val = merit_fun0;
            double eps_sufficient_descent = opts->eps_sufficient_descent;
            double dmerit_dy = 0.0;
            alpha = 1.0;

            // to avoid armijo evaluation and loop when checking if SOC should be done
            if (check_early_termination)
            {
                // TMP:
                // printf("tmp: merit_grad eval in early termination\n");
                // dmerit_dy = ocp_nlp_compute_merit_gradient(config, dims, in, out, opts, mem, work);

                // TODO(oj): should the merit weight update be undone in case of early termination?
                double violation_current = ocp_nlp_get_violation(config, dims, in, out, opts, mem, work);

                // tmp_nlp_out = out + alpha * qp_out
                for (i = 0; i <= N; i++)
                    blasfeo_daxpy(nv[i], alpha, qp_out->ux+i, 0, out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);
                merit_fun1 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);

                double violation_step = ocp_nlp_get_violation(config, dims, in, out, opts, mem, work);
                if (opts->print_level > 0)
                {
                    printf("\npreliminary line_search: merit0 %e, merit1 %e; viol_current %e, viol_step %e\n", merit_fun0, merit_fun1, violation_current, violation_step);
                }

                if (merit_fun1 < merit_fun0 && violation_step < violation_current)
                {
                    // full step if merit and constraint violation improves
                    // TODO: check armijo in this case?
                    return alpha;
                }
                else // this implies SOC will be done
                {
                    return reduction_factor * reduction_factor;
                }
            }

            /* actual Line Search*/
            if (opts->line_search_use_sufficient_descent)
            {
                // check Armijo-type sufficient descent condition Leinweber1999 (2.35);
                dmerit_dy = ocp_nlp_compute_merit_gradient(config, dims, in, out, opts, mem, work);
                if (dmerit_dy > 0.0)
                {
                    if (dmerit_dy > 1e-6 && opts->print_level > 0)
                    {
                        printf("\nacados line search: found dmerit_dy = %e > 0. Setting it to 0.0 instead\n", dmerit_dy);
                    }
                    dmerit_dy = 0.0;
                }
            }

            // From Leineweber1999: eq (3.64) -> only relevant for adaptive integrators looking at Remark 3.2.
            // "It is noteworthy that our practical implementation takes into account the potential nonsmoothness introduced by the fact that certain components of the penalty function - namely the continuity condition residuals - are evaluated only within integration tolerance."
            // double sum_pi = 0.0;
            // for (i = 0; i < N; i++)
            // {
            //     for (j = 0; j < dims->nx[i+1]; j++)
            //         sum_pi += BLASFEO_DVECEL(work->weight_merit_fun->pi+i, j);
            // }
            // double relaxed_val = 2.0 * 1e-6 * sum_pi;
            // if (abs(merit_fun0 - merit_fun1) < relaxed_val)
            // {
            //     printf("\nexiting because of relaxed_val.");
            //     break;
            // }

            for (j=0; alpha*reduction_factor > opts->alpha_min; j++)
            {
                // tmp_nlp_out = out + alpha * qp_out
                for (i = 0; i <= N; i++)
                    blasfeo_daxpy(nv[i], alpha, qp_out->ux+i, 0, out->ux+i, 0, work->tmp_nlp_out->ux+i, 0);

                merit_fun1 = ocp_nlp_evaluate_merit_fun(config, dims, in, out, opts, mem, work);
                if (opts->print_level > 1)
                {
                    printf("backtracking %d, merit_fun1 = %e, merit_fun0 %e\n", j, merit_fun1, merit_fun0);
                }

                // if (merit_fun1 < merit_fun0 && merit_fun1 > max_next_merit_fun_val)
                // {
                //     printf("\nalpha %f would be accepted without sufficient descent condition", alpha);
                // }

                max_next_merit_fun_val = merit_fun0 + eps_sufficient_descent * dmerit_dy * alpha;
                if (merit_fun1 < max_next_merit_fun_val)
                {
                    break;
                }
                else
                {
                    alpha *= reduction_factor;
                }
            }
        }
    }

    return alpha;
}


void ocp_nlp_update_variables_sqp(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work, double alpha)
{
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;


#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // step in primal variables
        blasfeo_daxpy(nv[i], alpha, mem->qp_out->ux + i, 0, out->ux + i, 0, out->ux + i, 0);

        // update dual variables
        if (opts->full_step_dual)
        {
            blasfeo_dveccp(2*ni[i], mem->qp_out->lam+i, 0, out->lam+i, 0);
            if (i < N)
            {
                blasfeo_dveccp(nx[i+1], mem->qp_out->pi+i, 0, out->pi+i, 0);
            }
        }
        else
        {
            // update duals with alpha step
            blasfeo_dvecsc(2*ni[i], 1.0-alpha, out->lam+i, 0);
            blasfeo_daxpy(2*ni[i], alpha, mem->qp_out->lam+i, 0, out->lam+i, 0, out->lam+i, 0);
            if (i < N)
            {
                blasfeo_dvecsc(nx[i+1], 1.0-alpha, out->pi+i, 0);
                blasfeo_daxpy(nx[i+1], alpha, mem->qp_out->pi+i, 0, out->pi+i, 0, out->pi+i, 0);
            }
        }

        // update slack values
        blasfeo_dvecsc(2*ni[i], 1.0-alpha, out->t+i, 0);
        blasfeo_daxpy(2*ni[i], alpha, mem->qp_out->t+i, 0, out->t+i, 0, out->t+i, 0);

        // linear update of algebraic variables using state and input sensitivity
        if (i < N)
        {
            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, mem->dzduxt+i, 0, 0,
                    mem->qp_out->ux+i, 0, 1.0, mem->z_alg+i, 0, out->z+i, 0);
        }
    }
}



/************************************************
 * residuals
 ************************************************/

acados_size_t ocp_nlp_res_calculate_size(ocp_nlp_dims *dims)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

    acados_size_t size = sizeof(ocp_nlp_res);

    size += 3 * (N + 1) * sizeof(struct blasfeo_dvec);  // res_stat res_ineq res_comp
    size += 1 * N * sizeof(struct blasfeo_dvec);        // res_eq

    for (int i = 0; i < N; i++)
    {
        size += 1 * blasfeo_memsize_dvec(nv[i]);      // res_stat
        size += 1 * blasfeo_memsize_dvec(nx[i + 1]);  // res_eq
        size += 2 * blasfeo_memsize_dvec(2 * ni[i]);  // res_ineq res_comp
    }
    size += 1 * blasfeo_memsize_dvec(nv[N]);      // res_stat
    size += 2 * blasfeo_memsize_dvec(2 * ni[N]);  // res_ineq res_comp

    size += 1 * blasfeo_memsize_dvec(N);      // tmp

    size += 8;   // initial align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    make_int_multiple_of(8, &size);

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

    // res_stat
    assign_and_advance_blasfeo_dvec_structs(N + 1, &res->res_stat, &c_ptr);
    // res_eq
    assign_and_advance_blasfeo_dvec_structs(N, &res->res_eq, &c_ptr);
    // res_ineq
    assign_and_advance_blasfeo_dvec_structs(N + 1, &res->res_ineq, &c_ptr);
    // res_comp
    assign_and_advance_blasfeo_dvec_structs(N + 1, &res->res_comp, &c_ptr);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // res_stat
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[i], res->res_stat + i, &c_ptr);
    }
    // res_eq
    for (int i = 0; i < N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[i + 1], res->res_eq + i, &c_ptr);
    }
    // res_ineq
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[i], res->res_ineq + i, &c_ptr);
    }
    // res_comp
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[i], res->res_comp + i, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem(N, &res->tmp, &c_ptr);

    res->memsize = ocp_nlp_res_calculate_size(dims);

    assert((char *) raw_memory + res->memsize >= c_ptr);

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

    // res_stat
    for (int i = 0; i <= N; i++)
    {
        blasfeo_daxpy(nv[i], -1.0, mem->ineq_adj + i, 0, mem->cost_grad + i, 0, res->res_stat + i,
                      0);
        blasfeo_daxpy(nu[i] + nx[i], -1.0, mem->dyn_adj + i, 0, res->res_stat + i, 0,
                      res->res_stat + i, 0);
        blasfeo_dvecnrm_inf(nv[i], res->res_stat + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N+1, &res->tmp, 0, &res->inf_norm_res_stat);

    // res_eq
    for (int i = 0; i < N; i++)
    {
        blasfeo_dveccp(nx[i + 1], mem->dyn_fun + i, 0, res->res_eq + i, 0);
        blasfeo_dvecnrm_inf(nx[i + 1], res->res_eq + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N, &res->tmp, 0, &res->inf_norm_res_eq);

    // res_ineq
    res->inf_norm_res_ineq = 0.0;
    for (int i = 0; i <= N; i++)
    {
        blasfeo_daxpy(2 * ni[i], 1.0, out->t + i, 0, mem->ineq_fun + i, 0, res->res_ineq + i, 0);
        blasfeo_dvecnrm_inf(2 * ni[i], res->res_ineq + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N+1, &res->tmp, 0, &res->inf_norm_res_ineq);

    // res_comp
    res->inf_norm_res_comp = 0.0;
    for (int i = 0; i <= N; i++)
    {
        blasfeo_dvecmul(2 * ni[i], out->lam + i, 0, out->t + i, 0, res->res_comp + i, 0);
        blasfeo_dvecnrm_inf(2 * ni[i], res->res_comp + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N+1, &res->tmp, 0, &res->inf_norm_res_comp);


    // printf("computed residuals stat: %e, eq: %e, ineq: %e, comp: %e\n", res->inf_norm_res_stat, res->inf_norm_res_eq,
    //        res->inf_norm_res_ineq, res->inf_norm_res_comp);
}


void ocp_nlp_res_get_inf_norm(ocp_nlp_res *res, double *out)
{
    double norm = res->inf_norm_res_stat;
    norm = (res->inf_norm_res_eq > norm) ? res->inf_norm_res_eq : norm;
    norm = (res->inf_norm_res_ineq > norm) ? res->inf_norm_res_ineq : norm;
    norm = (res->inf_norm_res_comp > norm) ? res->inf_norm_res_comp : norm;
    *out = norm;
    return;
}


void ocp_nlp_cost_compute(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    // extract dims
    int N = dims->N;

    double* tmp_cost = NULL;
    double total_cost = 0.0;

    for (int i = 0; i <= N; i++)
    {
        // set pointers
        // NOTE(oj): the cost compute function takes the tmp_ux_ptr as input,
        //  since it is also used for globalization,
        //  especially with primal variables that are NOT current SQP iterates.
        config->cost[i]->memory_set_tmp_ux_ptr(out->ux+i, mem->cost[i]);

        config->cost[i]->compute_fun(config->cost[i], dims->cost[i], in->cost[i],
                    opts->cost[i], mem->cost[i], work->cost[i]);
        tmp_cost = config->cost[i]->memory_get_fun_ptr(mem->cost[i]);
        // printf("cost at stage %d = %e, total = %e\n", i, *tmp_cost, total_cost);
        total_cost += *tmp_cost;
    }
    mem->cost_value = total_cost;

    // printf("\ncomputed total cost: %e\n", total_cost);
}

