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


#include "acados/ocp_nlp/ocp_nlp_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// blasfeo
#include "blasfeo_common.h"
#include "blasfeo_d_blas.h"
// hpipm
#include "hpipm/include/hpipm_d_ocp_qp_dim.h"
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/strsep.h"
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

    // globalization
    size += ocp_nlp_globalization_config_calculate_size();

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

    // globalization
    config->globalization = ocp_nlp_globalization_config_assign(c_ptr);
    c_ptr += ocp_nlp_globalization_config_calculate_size();

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
    size += 10 * (N + 1) * sizeof(int);  // nv, nx, nu, ni, nz, ns, np, ng, nb, ni_nl

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
    // nz
    assign_and_advance_int(N + 1, &dims->nz, &c_ptr);
    // ns
    assign_and_advance_int(N + 1, &dims->ns, &c_ptr);
    // np
    assign_and_advance_int(N + 1, &dims->np, &c_ptr);
    // ni
    assign_and_advance_int(N + 1, &dims->ni, &c_ptr);
    // nb
    assign_and_advance_int(N + 1, &dims->nb, &c_ptr);
    // ng
    assign_and_advance_int(N + 1, &dims->ng, &c_ptr);
    // ni_nl
    assign_and_advance_int(N + 1, &dims->ni_nl, &c_ptr);

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
    // nz
    for(int i=0; i<=N; i++)
        dims->nz[i] = 0;
    // ns
    for(int i=0; i<=N; i++)
        dims->ns[i] = 0;
    // np
    for(int i=0; i<=N; i++)
        dims->np[i] = 0;
    // ni
    for(int i=0; i<=N; i++)
        dims->ni[i] = 0;
    // nb
    for(int i=0; i<=N; i++)
        dims->nb[i] = 0;
    // ng
    for(int i=0; i<=N; i++)
        dims->ng[i] = 0;
    // ni_nl
    for(int i=0; i<=N; i++)
        dims->ni_nl[i] = 0;

    dims->n_global_data = 0;
    dims->np_global = 0;

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


void ocp_nlp_dims_set_global(void *config_, void *dims_, const char *field, int value_field)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;
    int N = dims->N;

    if (!strcmp(field, "np_global"))
    {
        dims->np_global = value_field;
        // cost
        for (int i = 0; i <= N; i++)
        {
            config->cost[i]->dims_set(config->cost[i], dims->cost[i], "np_global", &value_field);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i], dims->dynamics[i], "np_global", &value_field);
        }
        // constraints
        for (int i = 0; i <= N; i++)
        {
            config->constraints[i]->dims_set(config->constraints[i], dims->constraints[i], "np_global", &value_field);
        }
    }
    else if (!strcmp(field, "n_global_data"))
    {
        dims->n_global_data = value_field;
    }
    else
    {
        printf("ocp_nlp_dims_set_global: field %s not supported.\n", field);
        exit(1);
    }
}


void ocp_nlp_dims_set_opt_vars(void *config_, void *dims_, const char *field,
                                    const void* value_array)
{
    // to set dimension nx, nu, nz, ns (number of slacks = number of soft constraints), np
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
    else if (!strcmp(field, "np"))
    {
        // nlp opt var -- TODO: not really an opt var, but we need to set it here
        // -> rename function to ocp_nlp_dims_set_common?
        for (int i = 0; i <= N; i++)
        {
            // set np
            dims->np[i] = int_array[i];
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            config->cost[i]->dims_set(config->cost[i], dims->cost[i], "np", &int_array[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            config->dynamics[i]->dims_set(config->dynamics[i], dims->dynamics[i], "np", &int_array[i]);
        }
        // TODO: implement np for constraints
        // // constraints
        // for (int i = 0; i <= N; i++)
        // {
        //     config->constraints[i]->dims_set(config->constraints[i], dims->constraints[i], "np", &int_array[i]);
        // }

    }
    else
    {
        printf("error: dims type not available in module ocp_nlp: %s", field);
        exit(1);
    }
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
    // update ocp_nlp dimensions
    config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                                        "ni", &dims->ni[i]);
    config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                                        "nb", &dims->nb[i]);
    config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                                        "ng", &dims->ng[i]);
    config->constraints[i]->dims_get(config->constraints[i], dims->constraints[i],
                                        "ni_nl", &dims->ni_nl[i]);

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

static acados_size_t ocp_nlp_in_calculate_size_self(ocp_nlp_dims *dims)
{
    int N = dims->N;
    acados_size_t size = sizeof(ocp_nlp_in);

    size += N * sizeof(double);  // Ts
    // parameter values
    for (int i = 0; i <= N; i++)
    {
        size += dims->np[i] * sizeof(double);
    }
    // global_data
    size += dims->n_global_data * sizeof(double);

    size += (N + 1) * sizeof(double *);

    size += N * sizeof(void *);  // dynamics

    size += (N + 1) * sizeof(void *);  // cost

    size += (N + 1) * sizeof(void *);  // constraints

    size += 4*8;  // aligns
    return size;
}



acados_size_t ocp_nlp_in_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    int N = dims->N;

    acados_size_t size = ocp_nlp_in_calculate_size_self(dims);

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

    make_int_multiple_of(8, &size);

    return size;
}



static ocp_nlp_in *ocp_nlp_in_assign_self(ocp_nlp_dims *dims, void *raw_memory)
{
    int N = dims->N;
    char *c_ptr = (char *) raw_memory;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_nlp_in *in = (ocp_nlp_in *) c_ptr;
    c_ptr += sizeof(ocp_nlp_in);

    // align
    align_char_to(8, &c_ptr);

    // double pointers
    assign_and_advance_double_ptrs(N+1, &in->parameter_values, &c_ptr);

    align_char_to(8, &c_ptr);

    // doubles
    // Ts
    assign_and_advance_double(N, &in->Ts, &c_ptr);

    // parameter values
    for (int i = 0; i <= N; i++)
    {
        assign_and_advance_double(dims->np[i], &in->parameter_values[i], &c_ptr);
        for (int ip = 0; ip < dims->np[i]; ip++)
        {
            in->parameter_values[i][ip] = 0.0;
        }
    }
    assign_and_advance_double(dims->n_global_data, &in->global_data, &c_ptr);


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

    assert((char *) raw_memory + ocp_nlp_in_calculate_size_self(dims) >= c_ptr);

    return in;
}



ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
    int N = dims->N;

    char *c_ptr = (char *) raw_memory;

    // struct
    ocp_nlp_in *in = ocp_nlp_in_assign_self(dims, c_ptr);
    c_ptr += ocp_nlp_in_calculate_size_self(dims);

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

    size += 3 * (N + 1) * sizeof(struct blasfeo_dvec);  // ux, lam, z
    size += 1 * N * sizeof(struct blasfeo_dvec);        // pi

    for (int i = 0; i < N; i++)
    {
        size += 1 * blasfeo_memsize_dvec(nv[i]);      // ux
        size += 1 * blasfeo_memsize_dvec(nz[i]);      // z
        size += 1 * blasfeo_memsize_dvec(2 * ni[i]);  // lam
        size += 1 * blasfeo_memsize_dvec(nx[i + 1]);  // pi
    }
    size += 1 * blasfeo_memsize_dvec(nv[N]);      // ux
    size += 1 * blasfeo_memsize_dvec(nz[N]);     // z
    size += 1 * blasfeo_memsize_dvec(2 * ni[N]);  // lam

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

    // zero solution
    for(int i=0; i<N; i++)
    {
        blasfeo_dvecse(nv[i], 0.0, out->ux+i, 0);
        blasfeo_dvecse(nz[i], 0.0, out->z+i, 0);
        blasfeo_dvecse(nx[i+1], 0.0, out->pi+i, 0);
        blasfeo_dvecse(2*ni[i], 0.0, out->lam+i, 0);
    }
    blasfeo_dvecse(nv[N], 0.0, out->ux+N, 0);
    blasfeo_dvecse(nz[N], 0.0, out->z+N, 0);
    blasfeo_dvecse(2*ni[N], 0.0, out->lam+N, 0);

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

    size += config->globalization->opts_calculate_size(config, dims);

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

    opts->globalization = config->globalization->opts_assign(config, dims, c_ptr);
    c_ptr += config->globalization->opts_calculate_size(config, dims);

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
    ocp_nlp_globalization_config *globalization = config->globalization;

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

    opts->print_level = 0;
    opts->levenberg_marquardt = 0.0;
    opts->log_primal_step_norm = 0;
    opts->max_iter = 1;

    /* submodules opts */
    // qp solver
    qp_solver->opts_initialize_default(qp_solver, dims->qp_solver, opts->qp_solver_opts);

    // regularization
    regularize->opts_initialize_default(regularize, dims->regularize, opts->regularize);

    // globalization
    globalization->opts_initialize_default(globalization, dims, opts->globalization);

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

    opts->with_solution_sens_wrt_params = 0;
    opts->with_value_sens_wrt_params = 0;

    // adaptive Levenberg-Marquardt options
    opts->adaptive_levenberg_marquardt_mu_min = 1e-16;
    opts->adaptive_levenberg_marquardt_lam = 5.0;
    opts->with_adaptive_levenberg_marquardt = false;

    opts->ext_qp_res = 0;
    opts->store_iterates = false;

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

    char *ptr_module = NULL;
    int module_length = 0;
    char module[MAX_STR_LEN];
    extract_module_name(field, module, &module_length, &ptr_module);

    // pass options to QP module
    if ( ptr_module!=NULL && (!strcmp(ptr_module, "qp")) )
    {
        config->qp_solver->opts_set(config->qp_solver, opts->qp_solver_opts,
                                    field+module_length+1, value);
    }
    else if ( ptr_module!=NULL && (!strcmp(ptr_module, "reg")) )
    {
        config->regularize->opts_set(config->regularize, opts->regularize,
                                    field+module_length+1, value);
    }
    else if ( ptr_module!=NULL && (!strcmp(ptr_module, "globalization")) )
    {
        config->globalization->opts_set(config->globalization, opts->globalization,
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
        else if (!strcmp(field, "ext_qp_res"))
        {
            int* ext_qp_res = (int *) value;
            opts->ext_qp_res = *ext_qp_res;
        }
        else if (!strcmp(field, "store_iterates"))
        {
            bool* store_iterates = (bool *) value;

            if (*store_iterates && config->is_real_time_algorithm())
            {
                printf("Warning: Can not store intermediate iterates for real-time solvers.\n");
            }
            else
            {
                opts->store_iterates = *store_iterates;
            }
        }
        else if (!strcmp(field, "levenberg_marquardt"))
        {
            double* levenberg_marquardt = (double *) value;
            opts->levenberg_marquardt = *levenberg_marquardt;
        }
        // newly added options for DDP and SQP
        else if (!strcmp(field, "with_adaptive_levenberg_marquardt"))
        {
            bool* with_adaptive_levenberg_marquardt = (bool *) value;
            opts->with_adaptive_levenberg_marquardt = *with_adaptive_levenberg_marquardt;
        }
        else if (!strcmp(field, "adaptive_levenberg_marquardt_lam"))
        {
            double* adaptive_levenberg_marquardt_lam = (double *) value;
            opts->adaptive_levenberg_marquardt_lam = *adaptive_levenberg_marquardt_lam;
        }
        else if (!strcmp(field, "adaptive_levenberg_marquardt_mu_min"))
        {
            double* adaptive_levenberg_marquardt_mu_min = (double *) value;
            opts->adaptive_levenberg_marquardt_mu_min = *adaptive_levenberg_marquardt_mu_min;
        }
        else if (!strcmp(field, "adaptive_levenberg_marquardt_mu0"))
        {
            double* adaptive_levenberg_marquardt_mu0 = (double *) value;
            opts->adaptive_levenberg_marquardt_mu0 = *adaptive_levenberg_marquardt_mu0;
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
        else if (!strcmp(field, "log_primal_step_norm"))
        {
            int* log_primal_step_norm = (int *) value;
            opts->log_primal_step_norm = *log_primal_step_norm;
        }
        else if (!strcmp(field, "max_iter"))
        {
            int* max_iter = (int *) value;

            if (*max_iter > 0 && config->is_real_time_algorithm())
            {
                printf("Warning: can not set max_iter > 1 for real-time solvers.");
            }
            else
            {
                opts->max_iter = *max_iter;
            }
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
        else if (!strcmp(field, "fixed_hess"))
        {
            int* fixed_hess = (int *) value;
            opts->fixed_hess = *fixed_hess;
        }
        else if (!strcmp(field, "with_solution_sens_wrt_params"))
        {
            int N = config->N;

            int* with_solution_sens_wrt_params = (int *) value;
            opts->with_solution_sens_wrt_params = *with_solution_sens_wrt_params;
            // cost
            for (int i=0; i<=N; i++)
                config->cost[i]->opts_set(config->cost[i], opts->cost[i], "with_solution_sens_wrt_params", value);
            // dynamics
            for (int i=0; i<N; i++)
                config->dynamics[i]->opts_set(config->dynamics[i], opts->dynamics[i],
                                               "with_solution_sens_wrt_params", value);
            // constraints
            for (int i=0; i<=N; i++)
                config->constraints[i]->opts_set(config->constraints[i], opts->constraints[i],
                                                  "with_solution_sens_wrt_params", value);
        }
        else if (!strcmp(field, "with_value_sens_wrt_params"))
        {
            int* with_value_sens_wrt_params = (int *) value;
            opts->with_value_sens_wrt_params = *with_value_sens_wrt_params;
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

    char *ptr_module = NULL;
    int module_length = 0;
    char module[MAX_STR_LEN];
    extract_module_name(field, module, &module_length, &ptr_module);

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

acados_size_t ocp_nlp_memory_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts, ocp_nlp_in *nlp_in)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // extract dims
    int N = dims->N;
    int np_global = dims->np_global;

    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nz = dims->nz;
    int *nu = dims->nu;
    int *ni = dims->ni;
    int *ni_nl = dims->ni_nl;

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

    // intermediate iterates
    if (opts->store_iterates)
    {
        size += (opts->max_iter + 1) * sizeof(struct ocp_nlp_out *);

        for (int i = 0; i <= opts->max_iter; i++)
        {
            size += ocp_nlp_out_calculate_size(config, dims);
        }
    }

    if (opts->with_solution_sens_wrt_params)
    {
        size += 2*(N+1)*sizeof(struct blasfeo_dmat); // jac_lag_stat_p_global, jac_ineq_p_global
        size += N * sizeof(struct blasfeo_dmat);  // jac_dyn_p_global
        for (int i = 0; i <= N; i++)
        {
            size += blasfeo_memsize_dmat(nv[i], np_global);  // jac_lag_stat_p_global
            size += blasfeo_memsize_dmat(ni_nl[i], np_global);  // jac_ineq_p_global
        }
        for (int i = 0; i < N; i++)
        {
            size += blasfeo_memsize_dmat(nx[i+1], np_global);  // jac_dyn_p_global
        }
    }

    // nlp res
    size += ocp_nlp_res_calculate_size(dims);

    // timings
    size += sizeof(struct ocp_nlp_timings);

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
    size += 1 * blasfeo_memsize_dvec(np_global); //  out_np_global;

    size += 8;   // initial align
    size += 8;   // middle align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    make_int_multiple_of(8, &size);

    return size;
}



ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts, ocp_nlp_in *in, void *raw_memory)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    // extract sizes
    int N = dims->N;
    int np_global = dims->np_global;

    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nz = dims->nz;
    int *nu = dims->nu;
    int *ni = dims->ni;
    int *ni_nl = dims->ni_nl;

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

    // intermediate iterates
    if (opts->store_iterates)
    {
        mem->iterates = (struct ocp_nlp_out **) c_ptr;
        c_ptr += (opts->max_iter + 1)*sizeof(struct ocp_nlp_out *);
    }

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

    // globalization
    mem->globalization = config->globalization->memory_assign(config->globalization, dims, c_ptr);
    c_ptr += config->globalization->memory_calculate_size(config->globalization, dims);
    // ->memory_calculate_size(config->globalization, dims);

    int i;
    // dynamics
    for (i = 0; i < N; i++)
    {
        mem->dynamics[i] = dynamics[i]->memory_assign(dynamics[i], dims->dynamics[i], opts->dynamics[i], c_ptr);
        c_ptr += dynamics[i]->memory_calculate_size(dynamics[i], dims->dynamics[i], opts->dynamics[i]);
    }

    // cost
    for (i = 0; i <= N; i++)
    {
        mem->cost[i] = cost[i]->memory_assign(cost[i], dims->cost[i], opts->cost[i], c_ptr);
        c_ptr += cost[i]->memory_calculate_size(cost[i], dims->cost[i], opts->cost[i]);
    }

    // constraints
    for (i = 0; i <= N; i++)
    {
        mem->constraints[i] = constraints[i]->memory_assign(constraints[i],
                                            dims->constraints[i], opts->constraints[i], c_ptr);
        c_ptr += constraints[i]->memory_calculate_size( constraints[i], dims->constraints[i],
                                                                 opts->constraints[i]);
    }

    // intermediate iterates
    if (opts->store_iterates)
    {
        for (i = 0; i <= opts->max_iter; i++)
        {
            mem->iterates[i] = ocp_nlp_out_assign(config, dims, c_ptr);
            c_ptr += ocp_nlp_out_calculate_size(config, dims);
        }
    }

    // nlp res
    mem->nlp_res = ocp_nlp_res_assign(dims, c_ptr);
    c_ptr += mem->nlp_res->memsize;

    // timings
    mem->nlp_timings = (ocp_nlp_timings*) c_ptr;
    c_ptr += sizeof(ocp_nlp_timings);


    // zero timings
    ocp_nlp_timings_reset(mem->nlp_timings);
    mem->nlp_timings->time_feedback = 0;
    mem->nlp_timings->time_preparation = 0;
    mem->nlp_timings->time_solution_sensitivities = 0;

    // blasfeo_struct align
    align_char_to(8, &c_ptr);

    if (opts->with_solution_sens_wrt_params)
    {
        assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->jac_lag_stat_p_global, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(N + 1, &mem->jac_ineq_p_global, &c_ptr);
        assign_and_advance_blasfeo_dmat_structs(N, &mem->jac_dyn_p_global, &c_ptr);
    }

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
    for (i = 0; i <= N; ++i)
    {
        mem->set_sim_guess[i] = false;
    }

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dmat
    if (opts->with_solution_sens_wrt_params)
    {
        for (i = 0; i <= N; i++)
        {
            assign_and_advance_blasfeo_dmat_mem(nv[i], np_global, mem->jac_lag_stat_p_global+i, &c_ptr);
            assign_and_advance_blasfeo_dmat_mem(ni_nl[i], np_global, mem->jac_ineq_p_global+i, &c_ptr);
        }
        for (i = 0; i < N; i++)
        {
            assign_and_advance_blasfeo_dmat_mem(nx[i+1], np_global, mem->jac_dyn_p_global+i, &c_ptr);
        }
    }

    // dzduxt
    for (i=0; i<=N; i++)
    {
        assign_and_advance_blasfeo_dmat_mem(nu[i]+nx[i], nz[i], mem->dzduxt+i, &c_ptr);
    }
    // z_alg
    for (i=0; i<=N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nz[i], mem->z_alg + i, &c_ptr);
    }
    // cost_grad
    for (i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[i], mem->cost_grad + i, &c_ptr);
    }
    // ineq_fun
    for (i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(2 * ni[i], mem->ineq_fun + i, &c_ptr);
    }
    // ineq_adj
    for (i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nv[i], mem->ineq_adj + i, &c_ptr);
    }
    // dyn_fun
    for (i = 0; i < N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[i + 1], mem->dyn_fun + i, &c_ptr);
    }
    // dyn_adj
    for (i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nu[i] + nx[i], mem->dyn_adj + i, &c_ptr);
    }
    // sim_guess
    for (i = 0; i <= N; i++)
    {
        assign_and_advance_blasfeo_dvec_mem(nx[i] + nz[i], mem->sim_guess + i, &c_ptr);
        // set to 0;
        blasfeo_dvecse(nx[i] + nz[i], 0.0, mem->sim_guess+i, 0);
        // printf("sim_guess i %d: %p\n", i, mem->sim_guess+i);
    }
    assign_and_advance_blasfeo_dvec_mem(np_global, &mem->out_np_global, &c_ptr);

    mem->compute_hess = 1;

    return mem;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_workspace_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts, ocp_nlp_in *in)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;
    int np_global = dims->np_global;

    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;
    // int *np = dims->np;
    // int *ns = dims->ns;
    int *nv = dims->nv;

    // int *nz = dims->nz;

    acados_size_t size = 0;

    // nlp
    size += sizeof(ocp_nlp_workspace);

    // tmp_nlp_out
    size += ocp_nlp_out_calculate_size(config, dims);

    // weight_merit_fun
    size += ocp_nlp_out_calculate_size(config, dims);

    // tmp_qp_in
    size += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // tmp_qp_out
    size += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    // blasfeo_dvec
    int nv_max = 0;
    int nx_max = 0;
    int ni_max = 0;
    // int np_max = 0;

    for (int i = 0; i <= N; i++)
    {
        nx_max = nx_max > nx[i] ? nx_max : nx[i];
        nv_max = nv_max > nv[i] ? nv_max : nv[i];
        ni_max = ni_max > ni[i] ? ni_max : ni[i];
        // np_max = np_max > np[i] ? np_max : np[i];
    }
    size += 1 * blasfeo_memsize_dvec(nx_max);
    size += 1 * blasfeo_memsize_dvec(nv_max);
    size += 1 * blasfeo_memsize_dvec(ni_max);

    size += 1 * blasfeo_memsize_dvec(np_global); //  tmp_np_global;

    // array of pointers
    // cost
    size += (N+1)*sizeof(void *);
    // dynamics
    size += N*sizeof(void *);
    // constraints
    size += (N+1)*sizeof(void *);

    // doubles
    size += nv_max * sizeof(double); // tmp_nv_double

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

    size_t ext_fun_workspace_size = 0;
    if (opts->reuse_workspace)
    {
#if defined(ACADOS_WITH_OPENMP)
        // constraints
        for (int i = 0; i <= N; i++)
        {
            ext_fun_workspace_size += constraints[i]->get_external_fun_workspace_requirement(constraints[i], dims->constraints[i], opts->constraints[i], in->constraints[i]);
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            ext_fun_workspace_size += cost[i]->get_external_fun_workspace_requirement(cost[i], dims->cost[i], opts->cost[i], in->cost[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            ext_fun_workspace_size += dynamics[i]->get_external_fun_workspace_requirement(dynamics[i], dims->dynamics[i], opts->dynamics[i], in->dynamics[i]);
        }
#else
        size_t tmp_size;
        // constraints
        for (int i = 0; i <= N; i++)
        {
            tmp_size = constraints[i]->get_external_fun_workspace_requirement(constraints[i], dims->constraints[i], opts->constraints[i], in->constraints[i]);
            ext_fun_workspace_size = tmp_size > ext_fun_workspace_size ? tmp_size : ext_fun_workspace_size;
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            tmp_size = cost[i]->get_external_fun_workspace_requirement(cost[i], dims->cost[i], opts->cost[i], in->cost[i]);
            ext_fun_workspace_size = tmp_size > ext_fun_workspace_size ? tmp_size : ext_fun_workspace_size;
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            tmp_size = dynamics[i]->get_external_fun_workspace_requirement(dynamics[i], dims->dynamics[i], opts->dynamics[i], in->dynamics[i]);
            ext_fun_workspace_size = tmp_size > ext_fun_workspace_size ? tmp_size : ext_fun_workspace_size;
        }
#endif
    }
    else
    {
        // constraints
        for (int i = 0; i <= N; i++)
        {
            ext_fun_workspace_size += constraints[i]->get_external_fun_workspace_requirement(constraints[i], dims->constraints[i], opts->constraints[i], in->constraints[i]);
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            ext_fun_workspace_size += cost[i]->get_external_fun_workspace_requirement(cost[i], dims->cost[i], opts->cost[i], in->cost[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            ext_fun_workspace_size += dynamics[i]->get_external_fun_workspace_requirement(dynamics[i], dims->dynamics[i], opts->dynamics[i], in->dynamics[i]);
        }
    }
    size += 64; // ext_fun_workspace_size align
    size += ext_fun_workspace_size;

    size += 8; // struct align
    size += 64; // blasfeo align
    return size;
}



ocp_nlp_workspace *ocp_nlp_workspace_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                             ocp_nlp_opts *opts, ocp_nlp_in *nlp_in, ocp_nlp_memory *mem, void *raw_memory)
{
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    ocp_nlp_dynamics_config **dynamics = config->dynamics;
    ocp_nlp_cost_config **cost = config->cost;
    ocp_nlp_constraints_config **constraints = config->constraints;

    int N = dims->N;
    int np_global = dims->np_global;
    int *nx = dims->nx;
    int *nv = dims->nv;
    // int *ns = dims->ns;
    // int *nu = dims->nu;
    int *ni = dims->ni;
    // int *nz = dims->nz;
    // int *np = dims->np;
    int nv_max = 0;
    int nx_max = 0;
    int ni_max = 0;
    // int np_max = 0;

    for (int i = 0; i <= N; i++)
    {
        nx_max = nx_max > nx[i] ? nx_max : nx[i];
        nv_max = nv_max > nv[i] ? nv_max : nv[i];
        ni_max = ni_max > ni[i] ? ni_max : ni[i];
        // np_max = np_max > np[i] ? np_max : np[i];
    }

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

    // tmp qp in
    work->tmp_qp_in = ocp_qp_in_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // tmp qp out
    work->tmp_qp_out = ocp_qp_out_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    // weight_merit_fun
    work->weight_merit_fun = ocp_nlp_out_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_out_calculate_size(config, dims);

    assign_and_advance_double(nv_max, &work->tmp_nv_double, &c_ptr);

    // align for blasfeo mem
    align_char_to(64, &c_ptr);

    // blasfeo_dvec
    assign_and_advance_blasfeo_dvec_mem(nv_max, &work->tmp_nv, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(ni_max, &work->tmp_ni, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx_max, &work->dxnext_dy, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(np_global, &work->tmp_np_global, &c_ptr);

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
    }

    // align for external_function workspace
    align_char_to(64, &c_ptr);

    if (opts->reuse_workspace)
    {
#if defined(ACADOS_WITH_OPENMP)
        /* dont reuse workspace */
        // constraints
        for (int i = 0; i <= N; i++)
        {
            constraints[i]->set_external_fun_workspaces(constraints[i], dims->constraints[i], opts->constraints[i], nlp_in->constraints[i], c_ptr);
            c_ptr += constraints[i]->get_external_fun_workspace_requirement(constraints[i], dims->constraints[i], opts->constraints[i], nlp_in->constraints[i]);
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            cost[i]->set_external_fun_workspaces(cost[i], dims->cost[i], opts->cost[i], nlp_in->cost[i], c_ptr);
            c_ptr += cost[i]->get_external_fun_workspace_requirement(cost[i], dims->cost[i], opts->cost[i], nlp_in->cost[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            dynamics[i]->set_external_fun_workspaces(dynamics[i], dims->dynamics[i], opts->dynamics[i], nlp_in->dynamics[i], c_ptr);
            c_ptr += dynamics[i]->get_external_fun_workspace_requirement(dynamics[i], dims->dynamics[i], opts->dynamics[i], nlp_in->dynamics[i]);
        }
#else
        /* Reuse workspace */
        // constraints
        for (int i = 0; i <= N; i++)
        {
            constraints[i]->set_external_fun_workspaces(constraints[i], dims->constraints[i], opts->constraints[i], nlp_in->constraints[i], c_ptr);
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            cost[i]->set_external_fun_workspaces(cost[i], dims->cost[i], opts->cost[i], nlp_in->cost[i], c_ptr);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            dynamics[i]->set_external_fun_workspaces(dynamics[i], dims->dynamics[i], opts->dynamics[i], nlp_in->dynamics[i], c_ptr);
        }
#endif
    }
    else
    {
        /* dont reuse workspace */
        // constraints
        for (int i = 0; i <= N; i++)
        {
            constraints[i]->set_external_fun_workspaces(constraints[i], dims->constraints[i], opts->constraints[i], nlp_in->constraints[i], c_ptr);
            c_ptr += constraints[i]->get_external_fun_workspace_requirement(constraints[i], dims->constraints[i], opts->constraints[i], nlp_in->constraints[i]);
        }
        // cost
        for (int i = 0; i <= N; i++)
        {
            cost[i]->set_external_fun_workspaces(cost[i], dims->cost[i], opts->cost[i], nlp_in->cost[i], c_ptr);
            c_ptr += cost[i]->get_external_fun_workspace_requirement(cost[i], dims->cost[i], opts->cost[i], nlp_in->cost[i]);
        }
        // dynamics
        for (int i = 0; i < N; i++)
        {
            dynamics[i]->set_external_fun_workspaces(dynamics[i], dims->dynamics[i], opts->dynamics[i], nlp_in->dynamics[i], c_ptr);
            c_ptr += dynamics[i]->get_external_fun_workspace_requirement(dynamics[i], dims->dynamics[i], opts->dynamics[i], nlp_in->dynamics[i]);
        }
    }

    assert((char *) work + mem->workspace_size >= c_ptr);

    return work;
}



/************************************************
 * functions
 ************************************************/
double ocp_nlp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_nlp_workspace *nlp_work)
{
    // Compute the QP objective function value
    double qp_cost = 0.0;
    int i, nux, ns;
    int N = dims->N;
    // Sum over stages 0 to N
    for (i = 0; i <= N; i++)
    {
        nux = dims->nx[i] + dims->nu[i];
        ns = dims->ns[i];
        // Calculate 0.5 * d.T H d
        blasfeo_dsymv_l(nux, 0.5, &qp_in->RSQrq[i], 0, 0, &qp_out->ux[i], 0,
                        0.0, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &nlp_work->tmp_nv, 0);

        // slack QP objective value, compare to computation in cost modules;
        // tmp_nv = 2 * z + Z .* slack;
        blasfeo_dveccpsc(2*ns, 2.0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
        blasfeo_dvecmulacc(2*ns, &qp_in->Z[i], 0, &qp_out->ux[i], nux, &nlp_work->tmp_nv, 0);
        // qp_cost += .5 * (tmp_nv .* slack)
        qp_cost += 0.5 * blasfeo_ddot(2*ns, &nlp_work->tmp_nv, 0, &qp_out->ux[i], nux);
        // Calculate g.T d
        qp_cost += blasfeo_ddot(nux, &qp_out->ux[i], 0, &qp_in->rqz[i], 0);

        // Calculate gradient of slacks
        qp_cost += blasfeo_ddot(2 * ns, &qp_out->ux[i], nux, &qp_in->rqz[i], nux);
    }
    return qp_cost;
}


double ocp_nlp_get_l1_infeasibility(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_memory *nlp_mem)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *ni = dims->ni;
    int i;
    int j;

    // compute current l1 infeasibility
    double tmp;
    struct blasfeo_dvec *tmp_fun_vec;
    double dyn_l1_infeasibility = 0.0;
    for (i=0; i<N; i++)
    {
        tmp_fun_vec = config->dynamics[i]->memory_get_fun_ptr(nlp_mem->dynamics[i]);
        for (j=0; j<nx[i+1]; j++)
        {
            dyn_l1_infeasibility += fabs(BLASFEO_DVECEL(tmp_fun_vec, j));
        }
    }

    double constraint_l1_infeasibility = 0.0;
    for(i=0; i<=N; i++)
    {
        tmp_fun_vec = config->constraints[i]->memory_get_fun_ptr(nlp_mem->constraints[i]);
        // tmp_fun_vec = out->t+i;
        for (j=0; j<2*ni[i]; j++)
        {
            tmp = BLASFEO_DVECEL(tmp_fun_vec, j);
            if (tmp > 0.0)
            {
                constraint_l1_infeasibility += tmp;
            }
        }
    }
    return dyn_l1_infeasibility + constraint_l1_infeasibility;
}

void ocp_nlp_set_primal_variable_pointers_in_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                                                       ocp_nlp_out *nlp_out, ocp_nlp_memory *nlp_mem)
{
    int N = dims->N;
    for (int i = 0; i < N; i++)
    {
        config->dynamics[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_ux1_ptr(nlp_out->ux+i+1, nlp_mem->dynamics[i]);
    }
    for (int i = 0; i <= N; i++)
    {
        config->cost[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->cost[i]);
        config->constraints[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->constraints[i]);
    }
    return;
}


static void ocp_nlp_regularize_set_qp_in_ptrs(ocp_nlp_reg_config *reg_config, ocp_nlp_reg_dims *reg_dims, void *reg_mem, ocp_qp_in *qp_in)
{
    reg_config->memory_set_RSQrq_ptr(reg_dims, qp_in->RSQrq, reg_mem);
    reg_config->memory_set_rq_ptr(reg_dims, qp_in->rqz, reg_mem);
    reg_config->memory_set_BAbt_ptr(reg_dims, qp_in->BAbt, reg_mem);
    reg_config->memory_set_b_ptr(reg_dims, qp_in->b, reg_mem);
    reg_config->memory_set_idxb_ptr(reg_dims, qp_in->idxb, reg_mem);
    reg_config->memory_set_DCt_ptr(reg_dims, qp_in->DCt, reg_mem);
}

static void ocp_nlp_regularize_set_qp_out_ptrs(ocp_nlp_reg_config *reg_config, ocp_nlp_reg_dims *reg_dims, void *reg_mem, ocp_qp_out *qp_out)
{
    reg_config->memory_set_ux_ptr(reg_dims, qp_out->ux, reg_mem);
    reg_config->memory_set_pi_ptr(reg_dims, qp_out->pi, reg_mem);
    reg_config->memory_set_lam_ptr(reg_dims, qp_out->lam, reg_mem);
}


void ocp_nlp_alias_memory_to_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
         ocp_nlp_out *nlp_out, ocp_nlp_opts *opts, ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work)
{
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel
    { // beginning of parallel region
#endif

    int N = dims->N;
    // TODO: For z, why dont we use nlp_out->z+i instead of nlp_mem->z_alg+i? as is done for ux.
    //  - z_alg contains values from integrator, used in cost and constraint linearization.
    //  - nlp_out->z is updated as nlp_out->z = mem->z_alg + alpha * dzdux * qp_out->ux
    // Probably, this can also be achieved without mem->z_alg.
    // Would it work to initialize integrator always with z_out? Probably no, e.g. for lifted IRK.


    // alias to dynamics_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (int i = 0; i < N; i++)
    {
        config->dynamics[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_ux1_ptr(nlp_out->ux+i+1, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_pi_ptr(nlp_out->pi+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_BAbt_ptr(nlp_mem->qp_in->BAbt+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_dzduxt_ptr(nlp_mem->dzduxt+i, nlp_mem->dynamics[i]);
        config->dynamics[i]->memory_set_sim_guess_ptr(nlp_mem->sim_guess+i, nlp_mem->set_sim_guess+i, nlp_mem->dynamics[i]);
        // NOTE: no z at terminal stage, since dynamics modules dont compute it.
        config->dynamics[i]->memory_set_z_alg_ptr(nlp_mem->z_alg+i, nlp_mem->dynamics[i]);

        if (opts->with_solution_sens_wrt_params)
        {
            config->dynamics[i]->memory_set_dyn_jac_p_global_ptr(nlp_mem->jac_dyn_p_global+i, nlp_mem->dynamics[i]);
            config->dynamics[i]->memory_set_jac_lag_stat_p_global_ptr(nlp_mem->jac_lag_stat_p_global+i, nlp_mem->dynamics[i]);
        }

        int cost_integration;
        config->dynamics[i]->opts_get(config->dynamics[i], opts->dynamics[i],
                                    "cost_computation", &cost_integration);
        if (cost_integration)
        {
            // set pointers to cost function & gradient in integrator
            double *cost_fun = config->cost[i]->memory_get_fun_ptr(nlp_mem->cost[i]);
            struct blasfeo_dvec *cost_grad = config->cost[i]->memory_get_grad_ptr(nlp_mem->cost[i]);
            struct blasfeo_dvec *y_ref = config->cost[i]->model_get_y_ref_ptr(nlp_in->cost[i]);
            struct blasfeo_dmat *W_chol = config->cost[i]->memory_get_W_chol_ptr(nlp_mem->cost[i]);
            struct blasfeo_dvec *W_chol_diag = config->cost[i]->memory_get_W_chol_diag_ptr(nlp_mem->cost[i]);
            double *outer_hess_is_diag = config->cost[i]->get_outer_hess_is_diag_ptr(nlp_mem->cost[i], nlp_in->cost[i]);

            config->dynamics[i]->memory_set(config->dynamics[i], dims->dynamics[i], nlp_mem->dynamics[i], "cost_grad", cost_grad);
            config->dynamics[i]->memory_set(config->dynamics[i], dims->dynamics[i], nlp_mem->dynamics[i], "cost_fun", cost_fun);
            config->dynamics[i]->memory_set(config->dynamics[i], dims->dynamics[i], nlp_mem->dynamics[i], "y_ref", y_ref);
            config->dynamics[i]->memory_set(config->dynamics[i], dims->dynamics[i], nlp_mem->dynamics[i], "W_chol", W_chol);
            config->dynamics[i]->memory_set(config->dynamics[i], dims->dynamics[i], nlp_mem->dynamics[i], "W_chol_diag", W_chol_diag);
            config->dynamics[i]->memory_set(config->dynamics[i], dims->dynamics[i], nlp_mem->dynamics[i], "outer_hess_is_diag", outer_hess_is_diag);
        }
    }

    // alias to cost_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for nowait
#endif
    for (int i = 0; i <= N; i++)
    {
        if (opts->with_solution_sens_wrt_params)
        {
            config->cost[i]->memory_set_jac_lag_stat_p_global_ptr(nlp_mem->jac_lag_stat_p_global+i, nlp_mem->cost[i]);
        }
        config->cost[i]->memory_set_ux_ptr(nlp_out->ux+i, nlp_mem->cost[i]);
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
        config->constraints[i]->memory_set_lam_ptr(nlp_out->lam+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_z_alg_ptr(nlp_mem->z_alg+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_dzdux_tran_ptr(nlp_mem->dzduxt+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_DCt_ptr(nlp_mem->qp_in->DCt+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+i, nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_idxb_ptr(nlp_mem->qp_in->idxb[i], nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_idxs_rev_ptr(nlp_mem->qp_in->idxs_rev[i], nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_idxe_ptr(nlp_mem->qp_in->idxe[i], nlp_mem->constraints[i]);
        config->constraints[i]->memory_set_dmask_ptr(nlp_mem->qp_in->d_mask+i, nlp_mem->constraints[i]);
        if (opts->with_solution_sens_wrt_params)
        {
            config->constraints[i]->memory_set_jac_lag_stat_p_global_ptr(nlp_mem->jac_lag_stat_p_global+i, nlp_mem->constraints[i]);
            config->constraints[i]->memory_set_jac_ineq_p_global_ptr(nlp_mem->jac_ineq_p_global+i, nlp_mem->constraints[i]);
        }
    }

    // alias to regularize memory
    ocp_nlp_regularize_set_qp_in_ptrs(config->regularize, dims->regularize, nlp_mem->regularize_mem, nlp_mem->qp_in);
    ocp_nlp_regularize_set_qp_out_ptrs(config->regularize, dims->regularize, nlp_mem->regularize_mem, nlp_mem->qp_out);

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



static void adaptive_levenberg_marquardt_update_mu(double iter, double step_size, ocp_nlp_opts *opts, ocp_nlp_memory *mem)
{
    if (iter == 0)
    {
        mem->adaptive_levenberg_marquardt_mu = opts->adaptive_levenberg_marquardt_mu0;
        mem->adaptive_levenberg_marquardt_mu_bar = opts->adaptive_levenberg_marquardt_mu0;
    }
    else
    {
        if (step_size == 1.0)
        {
            double mu_tmp = mem->adaptive_levenberg_marquardt_mu;
            mem->adaptive_levenberg_marquardt_mu = fmax(opts->adaptive_levenberg_marquardt_mu_min,
                            mem->adaptive_levenberg_marquardt_mu_bar/(opts->adaptive_levenberg_marquardt_lam));
            mem->adaptive_levenberg_marquardt_mu_bar = mu_tmp;
        }
        else
        {
            mem->adaptive_levenberg_marquardt_mu = fmin(opts->adaptive_levenberg_marquardt_lam * mem->adaptive_levenberg_marquardt_mu, 1.0);
        }
    }
}

void ocp_nlp_add_levenberg_marquardt_term(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem,
    ocp_nlp_workspace *work, double alpha, int iter)
{
    if (opts->with_adaptive_levenberg_marquardt)
    {
        adaptive_levenberg_marquardt_update_mu(iter, alpha, opts, mem);
        double reg_param = 2*mem->cost_value*mem->adaptive_levenberg_marquardt_mu;
        opts->levenberg_marquardt = reg_param;
    }
    // Only add the Levenberg-Marquardt term when it is bigger than zero
    if (mem->compute_hess && opts->levenberg_marquardt > 0.0)
    {
        int N = dims->N;
        int *nx = dims->nx;
        int *nu = dims->nu;
        for (int i = 0; i <= N; i++)
        {
            if (i < N)
            {
                // Levenberg Marquardt term: Ts[i] * levenberg_marquardt * eye()
                blasfeo_ddiare(nu[i] + nx[i], in->Ts[i] * opts->levenberg_marquardt,
                                mem->qp_in->RSQrq+i, 0, 0);
            }
            else
            {
                // Levenberg Marquardt term: 1.0 * levenberg_marquardt * eye()
                blasfeo_ddiare(nu[i] + nx[i], opts->levenberg_marquardt,
                                mem->qp_in->RSQrq+i, 0, 0);
            }
        }
    } // else: do nothing
}



static void collect_integrator_timings(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_memory *mem)
{
    /* collect stage-wise timings */
    ocp_nlp_timings *nlp_timings = mem->nlp_timings;
    for (int ii=0; ii < dims->N; ii++)
    {
        double tmp_time;
        config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->dynamics[ii], "time_sim", &tmp_time);
        nlp_timings->time_sim += tmp_time;
        config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->dynamics[ii], "time_sim_la", &tmp_time);
        nlp_timings->time_sim_la += tmp_time;
        config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->dynamics[ii], "time_sim_ad", &tmp_time);
        nlp_timings->time_sim_ad += tmp_time;
    }
}

void ocp_nlp_approximate_qp_matrices(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem,
    ocp_nlp_workspace *work)
{
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;

    /* stage-wise multiple shooting lagrangian evaluation */
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // init Hessian to 0
        if (mem->compute_hess)
        {
            blasfeo_dgese(nu[i] + nx[i], nu[i] + nx[i], 0.0, mem->qp_in->RSQrq+i, 0, 0);
        }

        if (i < N)
        {
            // dynamics
            config->dynamics[i]->update_qp_matrices(config->dynamics[i], dims->dynamics[i],
                    in->dynamics[i], opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
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
            // TODO: this could be simplified by not copying pi in the dynamics module.
            struct blasfeo_dvec *dyn_adj
                = config->dynamics[i-1]->memory_get_adj_ptr(mem->dynamics[i-1]);
            blasfeo_daxpy(nx[i], 1.0, dyn_adj, nu[i-1]+nx[i-1], mem->dyn_adj+i, nu[i],
                mem->dyn_adj+i, nu[i]);
        }

        // nlp mem: ineq_adj
        struct blasfeo_dvec *ineq_adj =
            config->constraints[i]->memory_get_adj_ptr(mem->constraints[i]);
        blasfeo_dveccp(nv[i], ineq_adj, 0, mem->ineq_adj + i, 0);
    }

    collect_integrator_timings(config, dims, mem);
}



// update QP rhs for SQP (step prim var, abs dual var)
// - use cost gradient and dynamics residual from memory
// - evaluate constraints wrt bounds -> allows to update all bounds between preparation and feedback phase.
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

        // evaluate constraint residuals
        config->constraints[i]->update_qp_vectors(config->constraints[i], dims->constraints[i],
            in->constraints[i], opts->constraints[i], mem->constraints[i], work->constraints[i]);

        // copy ineq function value into nlp mem, then into QP
        struct blasfeo_dvec *ineq_fun = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, mem->ineq_fun + i, 0);

        // d
        blasfeo_dveccp(2 * ni[i], mem->ineq_fun + i, 0, mem->qp_in->d + i, 0);
    }
}

// zero order update QP: Update all constraint evaluations in QP
void ocp_nlp_zero_order_qp_update(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int N = dims->N;
    // int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // evaluate constraint residuals
        config->constraints[i]->compute_fun(config->constraints[i], dims->constraints[i],
            in->constraints[i], opts->constraints[i], mem->constraints[i], work->constraints[i]);
        // copy ineq function value into QP
        struct blasfeo_dvec *ineq_fun = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, mem->qp_in->d + i, 0);
        // copy into nlp_mem
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, mem->ineq_fun + i, 0);
    }

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<N; i++)
    {
        // dynamics
        config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i], in->dynamics[i],
                                         opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);

        struct blasfeo_dvec *dyn_fun = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        blasfeo_dveccp(nx[i + 1], dyn_fun, 0, mem->qp_in->b + i, 0);
        blasfeo_dveccp(nx[i + 1], dyn_fun, 0, mem->dyn_fun + i, 0);
    }

    // add gradient correction
    // rqz += Hess * last_step = RQ * qp_out
    for (int i = 0; i <= N; i++)
    {
        // NOTE: only lower triagonal of RSQ is stored
        blasfeo_dsymv_l_mn(nx[i]+nu[i], nx[i]+nu[i], 1.0, mem->qp_in->RSQrq+i, 0, 0,
                        mem->qp_out->ux+i, 0, 1.0, mem->qp_in->rqz+i, 0, mem->qp_in->rqz+i, 0);
        // TODO: fix for ns > 0.
    }
}


// Level C iterations Update all constraint evaluations in QP and Lagrange gradient
void ocp_nlp_level_c_update(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i = 0; i <= N; i++)
    {
        // evaluate constraint residuals
        config->constraints[i]->compute_fun(config->constraints[i], dims->constraints[i],
            in->constraints[i], opts->constraints[i], mem->constraints[i], work->constraints[i]);
        // copy ineq function value into QP
        struct blasfeo_dvec *ineq_fun = config->constraints[i]->memory_get_fun_ptr(mem->constraints[i]);
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, mem->qp_in->d + i, 0);
        // copy into nlp_mem
        blasfeo_dveccp(2 * ni[i], ineq_fun, 0, mem->ineq_fun + i, 0);
    }

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<=N; i++)
    {
        // nlp mem: cost_grad
        config->cost[i]->compute_gradient(config->cost[i], dims->cost[i], in->cost[i], opts->cost[i], mem->cost[i], work->cost[i]);
        struct blasfeo_dvec *cost_grad = config->cost[i]->memory_get_grad_ptr(mem->cost[i]);
        blasfeo_dveccp(nv[i], cost_grad, 0, mem->cost_grad + i, 0);
        blasfeo_dveccp(nv[i], mem->cost_grad + i, 0, mem->qp_in->rqz + i, 0);
    }


#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (int i=0; i<N; i++)
    {
        // dynamics
        // config->dynamics[i]->update_qp_matrices(config->dynamics[i], dims->dynamics[i], in->dynamics[i],
        config->dynamics[i]->compute_fun_and_adj(config->dynamics[i], dims->dynamics[i], in->dynamics[i],
                                         opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);

        struct blasfeo_dvec *dyn_fun = config->dynamics[i]->memory_get_fun_ptr(mem->dynamics[i]);
        blasfeo_dveccp(nx[i + 1], dyn_fun, 0, mem->qp_in->b + i, 0);
        blasfeo_dveccp(nx[i + 1], dyn_fun, 0, mem->dyn_fun + i, 0);

        // add adjoint contribution to gradient
        struct blasfeo_dvec *dyn_adj = config->dynamics[i]->memory_get_adj_ptr(mem->dynamics[i]);
        blasfeo_dvecad(nu[i] + nx[i], -1.0, dyn_adj, 0, mem->qp_in->rqz+i, 0);
        // add adjoint contribution C * lambda_k
        blasfeo_dgemv_n(nu[i] + nx[i], nx[i+1], -1.0, mem->qp_in->BAbt+i, 0, 0, out->pi+i, 0, 1.0, mem->qp_in->rqz+i, 0, mem->qp_in->rqz+i, 0);

        // - I part is linear, so dont need to add that!
        // blasfeo_dvecad(nx[i+1], 1.0, out->pi+i, 0, mem->qp_in->rqz+i, 0)

        // DEBUG:
        // printf("\ndyn_adj i %d\n", i);
        // blasfeo_print_exp_tran_dvec(nu[i] + nx[i], dyn_adj, 0);
        // blasfeo_dgemv_n(nu[i] + nx[i], nx[i+1], 1.0, mem->qp_in->BAbt+i, 0, 0, out->pi+i, 0, 0.0, &work->tmp_nv, 0, &work->tmp_nv, 0);
        // printf("C * lam\n");
        // blasfeo_print_exp_tran_dvec(nu[i] + nx[i], &work->tmp_nv, 0);
    }

    // TODO:
    // - adjoint call for inequalities as for dynamics
}


/*
calculates new iterate or trial iterate in 'out_destination' with step 'mem->qp_out',
step size 'alpha', and current iterate 'out_start'.
 */
void ocp_nlp_update_variables_sqp(void *config_, void *dims_,
            void *in_, void *out_, void *opts_, void *mem_,
            void *work_, void *out_destination_,
            void *solver_mem, double alpha, bool full_step_dual)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_out *out_start = out_;
    ocp_nlp_memory *mem = mem_;
    ocp_nlp_out *out_destination = out_destination_;
    // solver_mem is not used in this function, but needed for DDP
    // the function is used in the config->globalization->step_update
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
        blasfeo_daxpy(nv[i], alpha, mem->qp_out->ux + i, 0, out_start->ux + i, 0, out_destination->ux + i, 0);

        // update dual variables
        if (full_step_dual)
        {
            blasfeo_dveccp(2*ni[i], mem->qp_out->lam+i, 0, out_destination->lam+i, 0);
            if (i < N)
            {
                blasfeo_dveccp(nx[i+1], mem->qp_out->pi+i, 0, out_destination->pi+i, 0);
            }
        }
        else
        {
            // update duals with alpha step
            blasfeo_daxpby(2*ni[i], 1.0-alpha, out_start->lam+i, 0, alpha, mem->qp_out->lam+i, 0, out_destination->lam+i, 0);
            // blasfeo_dvecsc(2*ni[i], 1.0-alpha, out->lam+i, 0);
            // blasfeo_daxpy(2*ni[i], alpha, mem->qp_out->lam+i, 0, out->lam+i, 0, out->lam+i, 0);
            if (i < N)
            {
                // blasfeo_dvecsc(nx[i+1], 1.0-alpha, out->pi+i, 0);
                // blasfeo_daxpy(nx[i+1], alpha, mem->qp_out->pi+i, 0, out->pi+i, 0, out->pi+i, 0);
                blasfeo_daxpby(nx[i+1], 1.0-alpha, out_start->pi+i, 0, alpha, mem->qp_out->pi+i, 0, out_destination->pi+i, 0);
            }
        }

        // linear update of algebraic variables using state and input sensitivity
        if (i < N)
        {
            // out->z = mem->z_alg + alpha * dzdux * qp_out->ux
            blasfeo_dgemv_t(nu[i]+nx[i], nz[i], alpha, mem->dzduxt+i, 0, 0,
                    mem->qp_out->ux+i, 0, 1.0, mem->z_alg+i, 0, out_destination->z+i, 0);
        }
    }
}

void ocp_nlp_initialize_qp_from_nlp(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_qp_in *qp_in,
            ocp_nlp_out *out, ocp_qp_out *qp_out)
{
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *ni = dims->ni;

    for (int i = 0; i <= N; i++)
    {
        // set primal variables to zero
        blasfeo_dvecse(nv[i], 0.0, qp_out->ux+i, 0);

        // copy multipliers from ocp_nlp_out to ocp_qp_out
        blasfeo_dveccp(2*ni[i], out->lam+i, 0, qp_out->lam+i, 0);
        if (i < N)
            blasfeo_dveccp(nx[i+1], out->pi+i, 0, qp_out->pi+i, 0);
    }
    // compute t
    ocp_qp_compute_t(qp_in, qp_out);
}


int ocp_nlp_precompute_common(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int N = dims->N;
    int status = ACADOS_SUCCESS;
    int ii;

    for (ii = 0; ii <= N; ii++)
    {
        int module_val;
        config->constraints[ii]->dims_get(config->constraints[ii], dims->constraints[ii], "ns", &module_val);
        if (dims->ns[ii] != module_val)
        {
            printf("ocp_nlp_sqp_precompute: inconsistent dimension ns for stage %d with constraint module, got %d, module: %d.",
                   ii, dims->ns[ii], module_val);
            exit(1);
        }
    }

    // precompute
    for (ii = 0; ii < N; ii++)
    {
        // set T
        config->dynamics[ii]->model_set(config->dynamics[ii], dims->dynamics[ii],
                                        in->dynamics[ii], "T", in->Ts+ii);
        // dynamics precompute
        status = config->dynamics[ii]->precompute(config->dynamics[ii], dims->dynamics[ii],
                                                in->dynamics[ii], opts->dynamics[ii],
                                                mem->dynamics[ii], work->dynamics[ii]);
        if (status != ACADOS_SUCCESS)
            return status;
    }
    for (ii = 0; ii <= N; ii++)
    {
        // cost precompute
        config->cost[ii]->precompute(config->cost[ii], dims->cost[ii], in->cost[ii],
                                     opts->cost[ii], mem->cost[ii], work->cost[ii]);
    }

    ocp_nlp_alias_memory_to_submodules(config, dims, in, out, opts, mem, work);
    if (opts->fixed_hess)
    {
        mem->compute_hess = 1;
        ocp_nlp_approximate_qp_matrices(config, dims, in, out, opts, mem, work);
        mem->compute_hess = 0;
        for (ii=0; ii<=N; ii++)
            config->cost[ii]->opts_set(config->cost[ii], opts->cost[ii], "compute_hess", &mem->compute_hess);
    }

    return status;
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
    double tmp;

    // res_stat
    for (int i = 0; i <= N; i++)
    {
        blasfeo_daxpy(nv[i], -1.0, mem->ineq_adj + i, 0, mem->cost_grad + i, 0,
                      res->res_stat + i, 0);
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
        for (int j=0; j<2*ni[i]; j++)
        {
            tmp = BLASFEO_DVECEL(mem->ineq_fun+i, j);
            if (tmp > res->inf_norm_res_ineq)
            {
                res->inf_norm_res_ineq = tmp;
            }
        }
    }

    // res_comp
    res->inf_norm_res_comp = 0.0;
    for (int i = 0; i <= N; i++)
    {
        blasfeo_dvecmul(2 * ni[i], out->lam + i, 0, mem->ineq_fun+i, 0, res->res_comp + i, 0);
        blasfeo_dvecnrm_inf(2 * ni[i], res->res_comp + i, 0, &tmp_res);
        blasfeo_dvecse(1, tmp_res, &res->tmp, i);
    }
    blasfeo_dvecnrm_inf(N+1, &res->tmp, 0, &res->inf_norm_res_comp);
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


void copy_ocp_nlp_out(ocp_nlp_dims *dims, ocp_nlp_out *from, ocp_nlp_out *to)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;
    for (int i = 0; i <= N; i++)
    {
        blasfeo_dveccp(nv[i], from->ux+i, 0, to->ux+i, 0);
        blasfeo_dveccp(nz[i], from->z+i, 0, to->z+i, 0);
        blasfeo_dveccp(2*ni[i], from->lam+i, 0, to->lam+i, 0);
    }

    for (int i = 0; i < N; i++)
        blasfeo_dveccp(nx[i+1], from->pi+i, 0, to->pi+i, 0);

    return;
}

void ocp_nlp_get_cost_value_from_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int N = dims->N;

    double* tmp_cost = NULL;
    double total_cost = 0.0;

    for (int i = 0; i <= N; i++)
    {
        tmp_cost = config->cost[i]->memory_get_fun_ptr(mem->cost[i]);
        total_cost += *tmp_cost;
    }
    mem->cost_value = total_cost;
}


void ocp_nlp_cost_compute(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    int N = dims->N;

    double* tmp_cost = NULL;
    double total_cost = 0.0;

    int cost_integration;

    for (int i = 0; i <= N; i++)
    {
        if (i < N)
        {
            config->dynamics[i]->opts_get(config->dynamics[i], opts->dynamics[i], "cost_computation", &cost_integration);

            if (cost_integration)
            {
                config->dynamics[i]->compute_fun(config->dynamics[i], dims->dynamics[i],
                        in->dynamics[i], opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
            }
        }

        config->cost[i]->compute_fun(config->cost[i], dims->cost[i], in->cost[i],
                    opts->cost[i], mem->cost[i], work->cost[i]);
        tmp_cost = config->cost[i]->memory_get_fun_ptr(mem->cost[i]);
        // printf("cost at stage %d = %e, total = %e\n", i, *tmp_cost, total_cost);
        total_cost += *tmp_cost;
    }
    mem->cost_value = total_cost;

    // printf("\ncomputed total cost: %e\n", total_cost);
}



void ocp_nlp_params_jac_compute(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work)
{
    // This function sets up: jac_lag_stat_p_global, jac_ineq_p_global, jac_dyn_p_global
    // - jac_lag_stat_p_global: first dynamics writes its contribution, then cost and constraints modules add their contribution.
    // - jac_dyn_p_global is computed in dynamics module
    // - jac_ineq_p_global is computed in constraints module

    if (!opts->with_solution_sens_wrt_params)
    {
        printf("ocp_nlp_params_jac_compute: option with_solution_sens_wrt_params has to be true to evaluate solution sensitivities wrt. global parameters.\n");
        exit(1);
    }

    int N = dims->N;
    int np_global = dims->np_global;
    int i;

    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    struct blasfeo_dmat *jac_lag_stat_p_global = mem->jac_lag_stat_p_global;
    // struct blasfeo_dmat *jac_ineq_p_global = mem->jac_ineq_p_global;
    // struct blasfeo_dmat *jac_dyn_p_global = mem->jac_dyn_p_global;

#if defined(ACADOS_WITH_OPENMP)
    #pragma omp parallel for
#endif
    for (i = 0; i <= N; i++)
    {
        if (i < N)
        {
            // first nx+nu rows are overwritten by dynamics -> initialize ns part
            blasfeo_dgese(2*ns[i], np_global, 0., &jac_lag_stat_p_global[i], nx[i]+nu[i], 0);
            config->dynamics[i]->compute_jac_hess_p(config->dynamics[i], dims->dynamics[i], in->dynamics[i],
                        opts->dynamics[i], mem->dynamics[i], work->dynamics[i]);
        }
        else
        {
            // initialize jac_lag_stat_p_global = 0 as dynamics dont contribute
            blasfeo_dgese(nv[i], np_global, 0., &jac_lag_stat_p_global[i], 0, 0);
        }
        config->cost[i]->compute_jac_p(config->cost[i], dims->cost[i], in->cost[i],
                            opts->cost[i], mem->cost[i], work->cost[i]);
        config->constraints[i]->compute_jac_hess_p(config->constraints[i], dims->constraints[i],
                    in->constraints[i], opts->constraints[i], mem->constraints[i], work->constraints[i]);
    }
}


void ocp_nlp_common_eval_param_sens(ocp_nlp_config *config, ocp_nlp_dims *dims,
                        ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
                        char *field, int stage, int index, ocp_nlp_out *sens_nlp_out)
{
    int i;

    int N = dims->N;
    int *nv = dims->nv;
    int *ni = dims->ni;
    int *nx = dims->nx;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ni_nl = dims->ni_nl;

    struct blasfeo_dmat *jac_lag_stat_p_global = mem->jac_lag_stat_p_global;
    struct blasfeo_dmat *jac_ineq_p_global = mem->jac_ineq_p_global;
    struct blasfeo_dmat *jac_dyn_p_global = mem->jac_dyn_p_global;

    ocp_qp_in *tmp_qp_in = work->tmp_qp_in;
    ocp_qp_out *tmp_qp_out = work->tmp_qp_out;
    d_ocp_qp_copy_all(mem->qp_in, tmp_qp_in);
    d_ocp_qp_set_rhs_zero(tmp_qp_in);

    if ((!strcmp("ex", field)) && (stage==0))
    {
        double one = 1.0;
        d_ocp_qp_set_el("lbx", stage, index, &one, tmp_qp_in);
        d_ocp_qp_set_el("ubx", stage, index, &one, tmp_qp_in);
    }
    else if (!strcmp("p_global", field))
    {
        for (i = 0; i <= N; i++)
        {
            // stationarity
            blasfeo_dcolex(nv[i], &jac_lag_stat_p_global[i], 0, index, &tmp_qp_in->rqz[i], 0);
            // dynamics
            if (i < N)
                blasfeo_dcolex(nx[i+1], &jac_dyn_p_global[i], 0, index, &tmp_qp_in->b[i], 0);
            // inequalities
            blasfeo_dcolex(ni_nl[i], &jac_ineq_p_global[i], 0, index, &tmp_qp_in->d[i], nb[i]+ng[i]);
            blasfeo_dvecsc(ni_nl[i], -1.0, &tmp_qp_in->d[i], nb[i]+ng[i]);
            blasfeo_daxpy(ni_nl[i], -1.0, &tmp_qp_in->d[i], nb[i]+ng[i], &tmp_qp_in->d[i], 2*(nb[i]+ng[i])+ni_nl[i],
                                                                         &tmp_qp_in->d[i], 2*(nb[i]+ng[i])+ni_nl[i]);
        }
    }
    else
    {
        printf("\nerror: field %s at stage %d not available in ocp_nlp_common_eval_param_sens\n", field, stage);
        exit(1);
    }

    // d_ocp_qp_print(tmp_qp_in->dim, tmp_qp_in);
    config->qp_solver->eval_sens(config->qp_solver, dims->qp_solver, tmp_qp_in, tmp_qp_out,
                            opts->qp_solver_opts, mem->qp_solver_mem, work->qp_work);
    // d_ocp_qp_sol_print(tmp_qp_out->dim, tmp_qp_out);

    /* copy tmp_qp_out into sens_nlp_out */
    for (i = 0; i <= N; i++)
    {
        blasfeo_dveccp(nv[i], tmp_qp_out->ux + i, 0, sens_nlp_out->ux + i, 0);

        if (i < N)
            blasfeo_dveccp(nx[i + 1], tmp_qp_out->pi + i, 0, sens_nlp_out->pi + i, 0);

        blasfeo_dveccp(2 * ni[i], tmp_qp_out->lam + i, 0, sens_nlp_out->lam + i, 0);
    }
}


void ocp_nlp_common_eval_solution_sens_adj_p(ocp_nlp_config *config, ocp_nlp_dims *dims,
                        ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
                        ocp_nlp_out *sens_nlp_out, const char *field, int stage, void *grad_p)
{
    if (!opts->with_solution_sens_wrt_params)
    {
        printf("ocp_nlp_common_eval_solution_sens_adj_p: option with_solution_sens_wrt_params has to be true to evaluate solution sensitivities wrt. global parameters.\n");
        exit(1);
    }
    int i;
    int N = dims->N;
    int np_global = dims->np_global;

    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ni_nl = dims->ni_nl;

    struct blasfeo_dmat *jac_lag_stat_p_global = mem->jac_lag_stat_p_global;
    struct blasfeo_dmat *jac_ineq_p_global = mem->jac_ineq_p_global;
    struct blasfeo_dmat *jac_dyn_p_global = mem->jac_dyn_p_global;

    ocp_qp_in *tmp_qp_in = work->tmp_qp_in;
    ocp_qp_out *tmp_qp_out = work->tmp_qp_out;

    d_ocp_qp_copy_all(mem->qp_in, tmp_qp_in);
    d_ocp_qp_set_rhs_zero(tmp_qp_in);

    /* copy sens_nlp_out to tmp_qp_in */
    for (i = 0; i <= N; i++)
    {
        blasfeo_dveccp(nv[i], sens_nlp_out->ux + i, 0, tmp_qp_in->rqz + i, 0);
        // NOTE: noone needs sensitivities in adj dir pi, lam, t wrt. p
        // if (i < N)
        //     blasfeo_dveccp(nx[i + 1], sens_nlp_out->pi + i, 0, tmp_qp_in->b + i, 0);
        // blasfeo_dveccp(2 * ni[i], sens_nlp_out->lam + i, ?);
        // blasfeo_dveccp(2 * ni[i], sens_nlp_out->t + i, ?);
    }

    config->qp_solver->eval_sens(config->qp_solver, dims->qp_solver, tmp_qp_in, tmp_qp_out,
                            opts->qp_solver_opts, mem->qp_solver_mem, work->qp_work);

    if (!strcmp("p_global", field))
    {
        blasfeo_dvecse(np_global, 0., &mem->out_np_global, 0);
        for (i = 0; i <= N; i++)
        {
            /* multiply J.T with result of backsolve and add to in mem->out_np_global */
            // stationarity
            blasfeo_dgemv_t(nv[i], np_global, 1.0, &jac_lag_stat_p_global[i], 0, 0, tmp_qp_out->ux+i, 0, 1.0, &mem->out_np_global, 0, &mem->out_np_global, 0);
            // inequalities: upper
            blasfeo_dgemv_t(ni_nl[i], np_global, -1.0, &jac_ineq_p_global[i], 0, 0, tmp_qp_out->lam+i, nb[i]+ng[i], 1.0, &mem->out_np_global, 0, &mem->out_np_global, 0);
            // inequalities: lower
            blasfeo_dgemv_t(ni_nl[i], np_global, 1.0, &jac_ineq_p_global[i], 0, 0, tmp_qp_out->lam+i, 2*(nb[i]+ng[i])+ni_nl[i], 1.0, &mem->out_np_global, 0, &mem->out_np_global, 0);
            // dynamics
            if (i < N)
            {
                blasfeo_dgemv_t(nx[i+1], np_global, 1.0, &jac_dyn_p_global[i], 0, 0, tmp_qp_out->pi+i, 0, 1.0, &mem->out_np_global, 0, &mem->out_np_global, 0);
            }
        }

        // unpack
        blasfeo_unpack_dvec(np_global, &mem->out_np_global, 0, grad_p, 1);
    }
    else
    {
        printf("\nerror: field %s at stage %d not available in ocp_nlp_common_eval_solution_sens_adj_p\n", field, stage);
        exit(1);
    }
}


void ocp_nlp_common_eval_lagr_grad_p(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                        ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
                        const char *field, void *grad_p)
{
    int i;
    int N = dims->N;
    int np_global = dims->np_global;

    if (!strcmp("p_global", field))
    {
        // initialize to zero
        blasfeo_dvecse(np_global, 0., &mem->out_np_global, 0);

        for (i = 0; i < N; i++)
        {
            // dynamics contribution
            config->dynamics[i]->compute_adj_p(config->dynamics[i], dims->dynamics[i], in->dynamics[i], opts->dynamics[i],
                                    mem->dynamics[i], &work->tmp_np_global);
            blasfeo_dvecad(np_global, 1., &work->tmp_np_global, 0, &mem->out_np_global, 0);

            // cost contribution
            config->cost[i]->eval_grad_p(config->cost[i], dims->cost[i], in->cost[i], opts->cost[i],
                                    mem->cost[i], work->cost[i], &work->tmp_np_global);
            blasfeo_dvecad(np_global, 1., &work->tmp_np_global, 0, &mem->out_np_global, 0);

            // constraints contribution
            config->constraints[i]->compute_adj_p(config->constraints[i], dims->constraints[i], in->constraints[i], opts,
                                    mem->constraints[i], work->constraints[i], &work->tmp_np_global);

            blasfeo_dvecad(np_global, 1., &work->tmp_np_global, 0, &mem->out_np_global, 0);
        }

        // terminal cost contribution
        config->cost[N]->eval_grad_p(config->cost[N], dims->cost[N], in->cost[N], opts->cost[N],
                                    mem->cost[N], work->cost[N], &work->tmp_np_global);
        blasfeo_dvecad(np_global, 1., &work->tmp_np_global, 0, &mem->out_np_global, 0);

        blasfeo_unpack_dvec(np_global, &mem->out_np_global, 0, grad_p, 1);
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_common_eval_lagr_grad_p\n", field);
        exit(1);
    }
}


int ocp_nlp_solve_qp_and_correct_dual(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *nlp_opts,
                     ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work,
                     bool precondensed_lhs, ocp_qp_in *qp_in_, ocp_qp_out *qp_out_,
                     ocp_qp_xcond_solver *xcond_solver)
{
    acados_timer timer;

    // xcond_solver is "optional", if NULL is given use stuff from nlp_dims, mem etc.
    ocp_qp_xcond_solver_config *qp_solver;
    ocp_qp_xcond_solver_dims *qp_dims;
    ocp_qp_xcond_solver_opts *qp_opts;
    ocp_qp_xcond_solver_memory *qp_mem;
    ocp_qp_xcond_solver_workspace *qp_work;
    if (xcond_solver == NULL)
    {
        qp_solver = config->qp_solver;
        qp_dims = dims->qp_solver;
        qp_opts = nlp_opts->qp_solver_opts;
        qp_mem = nlp_mem->qp_solver_mem;
        qp_work = nlp_work->qp_work;
    }
    else
    {
        qp_solver = xcond_solver->config;
        qp_dims = xcond_solver->dims;
        qp_opts = xcond_solver->opts;
        qp_mem = xcond_solver->mem;
        qp_work = xcond_solver->work;
    }

    // qp_in_, qp_out_ are "optional", if NULL is given use nlp_mem->qp_in, nlp_mem->qp_out
    ocp_qp_in *qp_in;
    if (qp_in_ == NULL)
    {
        qp_in = nlp_mem->qp_in;
    }
    else
    {
        qp_in = qp_in_;
        ocp_nlp_regularize_set_qp_in_ptrs(config->regularize, dims->regularize, nlp_mem->regularize_mem, qp_in);
    }

    ocp_qp_out *qp_out = nlp_mem->qp_out;
    if (qp_out_ == NULL)
    {
        qp_out = nlp_mem->qp_out;
    }
    else
    {
        qp_out = qp_out_;
        ocp_nlp_regularize_set_qp_out_ptrs(config->regularize, dims->regularize, nlp_mem->regularize_mem, qp_out);
    }

    ocp_nlp_timings *nlp_timings = nlp_mem->nlp_timings;

    double tmp_time;
    int qp_status;

    // solve qp
    acados_tic(&timer);
    if (precondensed_lhs)
    {
        qp_status = qp_solver->condense_rhs_and_solve(qp_solver, qp_dims,
                qp_in, qp_out, qp_opts, qp_mem, qp_work);
    }
    else
    {
        qp_status = qp_solver->evaluate(qp_solver, qp_dims,
                qp_in, qp_out, qp_opts, qp_mem, qp_work);
    }
    // add qp timings
    nlp_timings->time_qp_sol += acados_toc(&timer);
    // NOTE: timings within qp solver are added internally (lhs+rhs)
    qp_solver->memory_get(qp_solver, qp_mem, "time_qp_solver_call", &tmp_time);
    nlp_timings->time_qp_solver_call += tmp_time;
    qp_solver->memory_get(qp_solver, qp_mem, "time_qp_xcond", &tmp_time);
    nlp_timings->time_qp_xcond += tmp_time;

    // compute correct dual solution in case of Hessian regularization
    acados_tic(&timer);
    config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                            nlp_opts->regularize, nlp_mem->regularize_mem);
    nlp_timings->time_reg += acados_toc(&timer);

    // reset regularize pointers if necessary
    if (qp_in_ != NULL)
    {
        ocp_nlp_regularize_set_qp_in_ptrs(config->regularize, dims->regularize, nlp_mem->regularize_mem, nlp_mem->qp_in);
    }
    if (qp_out_ != NULL)
    {
        ocp_nlp_regularize_set_qp_out_ptrs(config->regularize, dims->regularize, nlp_mem->regularize_mem, nlp_mem->qp_out);
    }

    return qp_status;
}


void ocp_nlp_dump_qp_in_to_file(ocp_qp_in *qp_in, int sqp_iter, int soc)
{
    char filename[100];
    if (soc > 0)
        sprintf(filename, "soc_qp_in_%d.txt", sqp_iter);
    else
        sprintf(filename, "qp_in_%d.txt", sqp_iter);
    FILE *out_file = fopen(filename, "w");
    print_ocp_qp_in_to_file(out_file, qp_in);
    fclose(out_file);
}


void ocp_nlp_dump_qp_out_to_file(ocp_qp_out *qp_out, int sqp_iter, int soc)
{
    char filename[100];
    if (soc > 0)
        sprintf(filename, "soc_qp_out_%d.txt", sqp_iter);
    else
        sprintf(filename, "qp_out_%d.txt", sqp_iter);
    FILE *out_file = fopen(filename, "w");
    print_ocp_qp_out_to_file(out_file, qp_out);
    fclose(out_file);
}


void ocp_nlp_timings_get(ocp_nlp_config *config, ocp_nlp_timings *timings, const char *field, void *return_value_)
{
    double *value = return_value_;
    if (!strcmp("time_tot", field))
    {
        *value = timings->time_tot;
    }
    else if (!strcmp("time_qp_sol", field) || !strcmp("time_qp", field))
    {
        *value = timings->time_qp_sol;
    }
    else if (!strcmp("time_qp_solver", field) || !strcmp("time_qp_solver_call", field))
    {
        *value = timings->time_qp_solver_call;
    }
    else if (!strcmp("time_qp_xcond", field))
    {
        *value = timings->time_qp_xcond;
    }
    else if (!strcmp("time_lin", field))
    {
        *value = timings->time_lin;
    }
    else if (!strcmp("time_reg", field))
    {
        *value = timings->time_reg;
    }
    else if (!strcmp("time_glob", field))
    {
        *value = timings->time_glob;
    }
    else if (!strcmp("time_solution_sensitivities", field))
    {
        *value = timings->time_solution_sensitivities;
    }
    else if (!strcmp("time_sim", field))
    {
        *value = timings->time_sim;
    }
    else if (!strcmp("time_sim_la", field))
    {
        *value = timings->time_sim_la;
    }
    else if (!strcmp("time_sim_ad", field))
    {
        *value = timings->time_sim_ad;
    }
    else if (!strcmp("time_preparation", field))
    {
        *value = timings->time_preparation;
    }
    else if (!strcmp("time_feedback", field))
    {
        if (config->is_real_time_algorithm())
        {
            *value = timings->time_feedback;
        }
        else
        {
            *value = timings->time_tot;
        }
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_timings_get\n", field);
        exit(1);
    }
}

void ocp_nlp_memory_get(ocp_nlp_config *config, ocp_nlp_memory *nlp_mem, const char *field, void *return_value_)
{
    if (!strcmp("sqp_iter", field) || !strcmp("nlp_iter", field) || !strcmp("ddp_iter", field))
    {
        int *value = return_value_;
        *value = nlp_mem->iter;
    }
    else if (!strcmp("status", field))
    {
        int *value = return_value_;
        *value = nlp_mem->status;
    }
    else if (!strcmp("nlp_mem", field))
    {
        void **value = return_value_;
        *value = nlp_mem;
    }
    else if (!strcmp("nlp_res", field))
    {
        ocp_nlp_res **value = return_value_;
        *value = nlp_mem->nlp_res;
    }
    else if (!strcmp("qp_xcond_in", field))
    {
        void **value = return_value_;
        *value = nlp_mem->qp_solver_mem->xcond_qp_in;
    }
    else if (!strcmp("qp_xcond_out", field))
    {
        void **value = return_value_;
        *value = nlp_mem->qp_solver_mem->xcond_qp_out;
    }
    else if (!strcmp("qp_in", field))
    {
        void **value = return_value_;
        *value = nlp_mem->qp_in;
    }
    else if (!strcmp("qp_out", field))
    {
        void **value = return_value_;
        *value = nlp_mem->qp_out;
    }
    else if (!strcmp("qp_iter", field))
    {
        config->qp_solver->memory_get(config->qp_solver,
            nlp_mem->qp_solver_mem, "iter", return_value_);
    }
    else if (!strcmp("qp_status", field))
    {
        config->qp_solver->memory_get(config->qp_solver,
            nlp_mem->qp_solver_mem, "status", return_value_);
    }
    else if (!strcmp("res_stat", field))
    {
        double *value = return_value_;
        *value = nlp_mem->nlp_res->inf_norm_res_stat;
    }
    else if (!strcmp("res_eq", field))
    {
        double *value = return_value_;
        *value = nlp_mem->nlp_res->inf_norm_res_eq;
    }
    else if (!strcmp("res_ineq", field))
    {
        double *value = return_value_;
        *value = nlp_mem->nlp_res->inf_norm_res_ineq;
    }
    else if (!strcmp("res_comp", field))
    {
        double *value = return_value_;
        *value = nlp_mem->nlp_res->inf_norm_res_comp;
    }
    else if (!strcmp("cost_value", field))
    {
        double *value = return_value_;
        *value = nlp_mem->cost_value;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_memory_get\n", field);
        exit(1);
    }
}


void ocp_nlp_timings_reset(ocp_nlp_timings *timings)
{
    timings->time_qp_sol = 0.0;
    timings->time_qp_solver_call = 0.0;
    timings->time_qp_xcond = 0.0;
    timings->time_lin = 0.0;
    timings->time_reg = 0.0;
    timings->time_glob = 0.0;
    timings->time_sim = 0.0;
    timings->time_sim_la = 0.0;
    timings->time_sim_ad = 0.0;
}
