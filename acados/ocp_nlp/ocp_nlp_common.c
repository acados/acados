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

int ocp_nlp_solver_config_calculate_size(int N)
{
    int ii;

    int size = 0;

    // qp solver
    size += sizeof(ocp_nlp_solver_config);

    size += 1 * ocp_qp_xcond_solver_config_calculate_size();

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



ocp_nlp_solver_config *ocp_nlp_solver_config_assign(int N, void *raw_memory)
{
    int ii;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_solver_config *config = (ocp_nlp_solver_config *) c_ptr;
    c_ptr += sizeof(ocp_nlp_solver_config);

    config->N = N;

    // qp solver
    config->qp_solver = ocp_qp_xcond_solver_config_assign(c_ptr);
    c_ptr += ocp_qp_xcond_solver_config_calculate_size();

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

int ocp_nlp_dims_calculate_size_self(int N)
{
    int size = 0;

    size += sizeof(ocp_nlp_dims);

    // nlp sizes
    size += 4 * (N + 1) * sizeof(int);  // nv, nx, nu, ni

    // dynamics
    size += N * sizeof(void *);

    // cost
    size += (N + 1) * sizeof(void *);

    // constraints
    size += (N + 1) * sizeof(void *);

    // qp solver
    size += ocp_qp_dims_calculate_size(N);

    size += 8;  // initial align

    return size;
}



int ocp_nlp_dims_calculate_size(void *config_)
{
    ocp_nlp_solver_config *config = config_;

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



ocp_nlp_dims *ocp_nlp_dims_assign_self(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

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

    // N
    dims->N = N;

    // assert
    assert((char *) raw_memory + ocp_nlp_dims_calculate_size_self(N) >= c_ptr);

    return dims;
}



ocp_nlp_dims *ocp_nlp_dims_assign(void *config_, void *raw_memory)
{
    ocp_nlp_solver_config *config = config_;

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



void ocp_nlp_dims_initialize(void *config_, int *nx, int *nu, int *ny, int *nbx, int *nbu, int *ng,
                             int *nh, int *nq, int *ns, int *nz, ocp_nlp_dims *dims)
{
    ocp_nlp_solver_config *config = config_;

    int ii;

    int N = dims->N;

    // nlp dims
    for (ii = 0; ii <= N; ii++)
    {
        dims->nv[ii] = nu[ii] + nx[ii] + 2 * ns[ii];
        dims->nx[ii] = nx[ii];
        dims->nu[ii] = nu[ii];
        dims->ni[ii] = nbx[ii] + nbu[ii] + ng[ii] + nh[ii] + ns[ii];
    }

    // dynamics
    for (ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->dims_initialize(config->dynamics[ii], dims->dynamics[ii], nx[ii],
                                              nu[ii], nx[ii + 1], nu[ii + 1], nz[ii]);
    }

    for (ii = 0; ii <= N; ii++)
    {
        config->cost[ii]->dims_initialize(config->cost[ii], dims->cost[ii], nx[ii], nu[ii], ny[ii],
                                          ns[ii]);
    }

    for (ii = 0; ii <= N; ii++)
    {
        config->constraints[ii]->dims_initialize(config->constraints[ii], dims->constraints[ii],
                                                 nx[ii], nu[ii], nbx[ii], nbu[ii], ng[ii], nh[ii],
                                                 nq[ii], ns[ii]);
    }

    dims->qp_solver->N = N;
    for (ii = 0; ii <= N; ii++)
    {
        dims->qp_solver->nx[ii] = nx[ii];
        dims->qp_solver->nu[ii] = nu[ii];
        dims->qp_solver->nbx[ii] = nbx[ii];
        dims->qp_solver->nbu[ii] = nbu[ii];
        dims->qp_solver->nb[ii] = nbx[ii] + nbu[ii];
        dims->qp_solver->ng[ii] = ng[ii] + nh[ii];
        dims->qp_solver->ns[ii] = ns[ii];
    }

    return;
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

int ocp_nlp_in_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
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

ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory)
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

int ocp_nlp_out_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{
    // extract dims
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

    int size = sizeof(ocp_nlp_out);

    size += 3 * (N + 1) * sizeof(struct blasfeo_dvec);  // ux lam
    size += 1 * N * sizeof(struct blasfeo_dvec);        // pi

    for (int ii = 0; ii < N; ii++)
    {
        size += 1 * blasfeo_memsize_dvec(nv[ii]);      // ux
        size += 1 * blasfeo_memsize_dvec(nx[ii + 1]);  // pi
        size += 2 * blasfeo_memsize_dvec(2 * ni[ii]);  // lam t
    }
    size += 1 * blasfeo_memsize_dvec(nv[N]);      // ux
    size += 2 * blasfeo_memsize_dvec(2 * ni[N]);  // lam t

    size += 8;   // initial align
    size += 8;   // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    //  make_int_multiple_of(64, &size);

    return size;
}

ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
    // extract sizes
    int N = dims->N;
    int *nv = dims->nv;
    int *nx = dims->nx;
    // int *nu = dims->nu;
    int *ni = dims->ni;

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

int ocp_nlp_memory_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
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

ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims,
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
