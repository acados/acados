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

#include "acados/ocp_nlp/ocp_nlp_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

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
    config->regularization = ocp_nlp_reg_config_assign(c_ptr);
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

// TODO(oj): should be static, but used by current ocp_nlp c++ interface
int ocp_nlp_dims_calculate_size_self(int N)
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

    // qp solver
    size += ocp_qp_dims_calculate_size(N);

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

    return size;
}



// TODO(oj): should be static, but used by current ocp_nlp c++ interface
ocp_nlp_dims *ocp_nlp_dims_assign_self(int N, void *raw_memory)
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

    // qp solver
    dims->qp_solver = ocp_qp_dims_assign(N, c_ptr);
    c_ptr += ocp_qp_dims_calculate_size(N);

    /* initialize qp_solver dimensions */
    dims->qp_solver->N = N;
    for (ii = 0; ii <= N; ii++)
    {
        // TODO(dimitris): values below are needed for reformulation of QP when soft constraints
        //   are not supported. Make this a bit more transparent as it clushes with nbx/nbu above.
        dims->qp_solver->nsbx[ii] = 0;
        dims->qp_solver->nsbu[ii] = 0;
        dims->qp_solver->nsg[ii] = 0;
    }

    // N
    dims->N = N;

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
        for (int ii = 0; ii <= N; ii++)
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
            config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "nx",
                                        &int_array[i]);
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
            config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "nu",
                                        &int_array[i]);
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
    config->constraints[i]->get_dims(config->constraints[i], dims->constraints[i],
                                        "ni", &dims->ni[i]);

    // update qp_solver dims
    if ( (!strcmp(field, "nbx")) || (!strcmp(field, "nbu")) )
    {
        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i,
                                    field, int_value);
        // Note(oj): how to decide if ocp_qp or dense? plan not available here..
        // Note(giaf): the NLP solver always calls a xcond_solver,
        //             so only its dimensions should be set now
        // if( is_ocp_qp_solver )
        // {
            // set nbx, respectively nbu in qp_solver
            // nb in qp_solver gets updated in the setter
        // config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i,field, &value[i]);
        // }
        // else  // dense_qp
        // {
        //     // update nb in qp_solver
        //     int nb;
        //     config->constraints[i]->get_dims(config->constraints[i],
        //                            dims->constraints[i], "nb", &nb);
        //     config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "nb", &nb);
        // }
    }
    if ( (!strcmp(field, "ng")) || (!strcmp(field, "nh")) )
    {
        // update ng_qp_solver in qp_solver
        int ng, nh, ng_qp_solver;
        config->constraints[i]->get_dims(config->constraints[i],
                                            dims->constraints[i], "ng", &ng);
        config->constraints[i]->get_dims(config->constraints[i],
                                            dims->constraints[i], "nh", &nh);

        ng_qp_solver = ng + nh;

        config->qp_solver->dims_set(config->qp_solver, dims->qp_solver, i, "ng", &ng_qp_solver);
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

    config->dynamics[stage]->dims_set(config->dynamics[stage],
                                  dims->dynamics[stage], field, int_value);
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

    size += 4 * (N + 1) * sizeof(struct blasfeo_dvec);  // ux lam z
    size += 1 * N * sizeof(struct blasfeo_dvec);        // pi

    for (int ii = 0; ii < N; ii++)
    {
        size += 1 * blasfeo_memsize_dvec(nv[ii]);      // ux
        size += 1 * blasfeo_memsize_dvec(nz[ii]);      // z
        size += 1 * blasfeo_memsize_dvec(nx[ii + 1]);  // pi
        size += 2 * blasfeo_memsize_dvec(2 * ni[ii]);  // lam t
    }
    size += 1 * blasfeo_memsize_dvec(nv[N]);      // ux
    size += 1 * blasfeo_memsize_dvec(nz[N]);     // z
    size += 2 * blasfeo_memsize_dvec(2 * ni[N]);  // lam t

    size += 8;   // initial align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    //  make_int_multiple_of(64, &size);

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

    assert((char *) raw_memory + ocp_nlp_out_calculate_size(config, dims) >= c_ptr);

    return out;
}

/************************************************
 * memory
 ************************************************/

int ocp_nlp_memory_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    int size = sizeof(ocp_nlp_memory);

    size += 4 * (N + 1) * sizeof(struct blasfeo_dvec);  // cost_grad ineq_fun ineq_adj dyn_adj
    size += 1 * N * sizeof(struct blasfeo_dvec);        // dyn_fun

    for (int ii = 0; ii < N; ii++)
    {
        size += 2 * blasfeo_memsize_dvec(nv[ii]);           // cost_grad ineq_adj
        size += 1 * blasfeo_memsize_dvec(nu[ii] + nx[ii]);  // dyn_adj
        size += 1 * blasfeo_memsize_dvec(nx[ii + 1]);       // dyn_fun
        size += 1 * blasfeo_memsize_dvec(2 * ni[ii]);       // ineq_fun
    }
    size += 2 * blasfeo_memsize_dvec(nv[N]);          // cost_grad ineq_adj
    size += 1 * blasfeo_memsize_dvec(nu[N] + nx[N]);  // dyn_adj
    size += 1 * blasfeo_memsize_dvec(2 * ni[N]);      // ineq_fun

    size += 8;   // initial align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    //  make_int_multiple_of(64, &size);

    return size;
}

ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                      void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // extract sizes
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    ocp_nlp_memory *mem = (ocp_nlp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_memory);

    // blasfeo_struct align
    align_char_to(8, &c_ptr);

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

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

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

    return mem;
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
