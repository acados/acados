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
#include <assert.h>

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
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

	size += 1*ocp_qp_xcond_solver_config_calculate_size();

	// dynamics
	size += N*sizeof(ocp_nlp_dynamics_config *);

	for (ii=0; ii<N; ii++)
		size += ocp_nlp_dynamics_config_calculate_size();

	// cost
	size += (N+1)*sizeof(ocp_nlp_cost_config *);

	for (ii=0; ii<=N; ii++)
		size += ocp_nlp_cost_config_calculate_size();

	// constraints
	size += (N+1)*sizeof(ocp_nlp_constraints_config *);

	for (ii=0; ii<=N; ii++)
		size += ocp_nlp_constraints_config_calculate_size();

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
	c_ptr += N*sizeof(ocp_nlp_dynamics_config *);

	for (ii=0; ii<N; ii++)
	{
		config->dynamics[ii] = ocp_nlp_dynamics_config_assign(c_ptr);
		c_ptr += ocp_nlp_dynamics_config_calculate_size();
	}

	// cost
	config->cost = (ocp_nlp_cost_config **) c_ptr;
	c_ptr += (N+1)*sizeof(ocp_nlp_cost_config *);

	for (ii=0; ii<=N; ii++)
	{
		config->cost[ii] = ocp_nlp_cost_config_assign(c_ptr);
		c_ptr += ocp_nlp_cost_config_calculate_size();
	}

	// constraints
	config->constraints = (ocp_nlp_constraints_config **) c_ptr;
	c_ptr += (N+1)*sizeof(ocp_nlp_constraints_config *);

	for (ii=0; ii<=N; ii++)
	{
		config->constraints[ii] = ocp_nlp_constraints_config_assign(c_ptr);
		c_ptr += ocp_nlp_constraints_config_calculate_size();
	}

	return config;
}



/************************************************
* dims
************************************************/

int ocp_nlp_dims_calculate_size(int N)
{
	int size = 0;

    size += sizeof(ocp_nlp_dims);

	// nlp sizes
	size += 3*(N+1)*sizeof(int); // nx, nu, ni

	// dynamics_dims
	size += N*sizeof(ocp_nlp_dynamics_dims *);
	size += N*ocp_nlp_dynamics_dims_calculate_size();

	// cost_dims
	size += (N+1)*sizeof(ocp_nlp_cost_dims *);
	size += (N+1)*ocp_nlp_cost_dims_calculate_size();

	// constraints_dims
	size += (N+1)*sizeof(ocp_nlp_constraints_dims *);
	size += (N+1)*ocp_nlp_constraints_dims_calculate_size();

	// qp solver
	size += ocp_qp_dims_calculate_size(N);

	size += 8; // initial align

	return size;
}



ocp_nlp_dims *ocp_nlp_dims_assign(int N, void *raw_memory)
{
	int ii;

    char *c_ptr = (char *) raw_memory;

	// initial align
	align_char_to(8, &c_ptr);

	// struct
	ocp_nlp_dims *dims = (ocp_nlp_dims *) c_ptr;
	c_ptr += sizeof(ocp_nlp_dims);

	// nx
	assign_and_advance_int(N+1, &dims->nx, &c_ptr);
	// nu
	assign_and_advance_int(N+1, &dims->nu, &c_ptr);
	// ni
	assign_and_advance_int(N+1, &dims->ni, &c_ptr);

	// dynamics
	dims->dynamics = (ocp_nlp_dynamics_dims **) c_ptr;
	c_ptr += N*sizeof(ocp_nlp_dynamics_dims *);

	for (ii=0; ii<N; ii++)
	{
		dims->dynamics[ii] = ocp_nlp_dynamics_dims_assign(c_ptr);
		c_ptr += ocp_nlp_dynamics_dims_calculate_size();
	}

	// cost dims
	dims->cost = (ocp_nlp_cost_dims **) c_ptr;
	c_ptr += (N+1)*sizeof(ocp_nlp_cost_dims *);

	for (ii=0; ii<=N; ii++)
	{
		dims->cost[ii] = ocp_nlp_cost_dims_assign(c_ptr);
		c_ptr += ocp_nlp_cost_dims_calculate_size();
	}

	// constraints
	dims->constraints = (ocp_nlp_constraints_dims **) c_ptr;
	c_ptr += (N+1)*sizeof(ocp_nlp_constraints_dims *);

	for (ii=0; ii<=N; ii++)
	{
		dims->constraints[ii] = ocp_nlp_constraints_dims_assign(c_ptr);
		c_ptr += ocp_nlp_constraints_dims_calculate_size();
	}

	// qp solver
	dims->qp_solver = ocp_qp_dims_assign(N, c_ptr);
	c_ptr += ocp_qp_dims_calculate_size(N);

	// N
	dims->N = N;

	// assert
    assert((char *) raw_memory + ocp_nlp_dims_calculate_size(N) >= c_ptr);

	return dims;
}



void ocp_nlp_dims_initialize(int *nx, int *nu, int *ny, int *nbx, int *nbu, int *ng, int *nh, int *nq, int *ns, ocp_nlp_dims *dims)
{
	int ii;

	int N = dims->N;

	// nlp dims
	for (ii=0; ii<=N; ii++)
	{
		dims->nx[ii] = nx[ii];
		dims->nu[ii] = nu[ii];
		dims->ni[ii] = nbx[ii]+nbu[ii]+ng[ii]+nh[ii];
	}

	// TODO dynamics and sim dims initialize ???
	for (ii=0; ii<N; ii++)
	{
		dims->dynamics[ii]->nx = nx[ii];
		dims->dynamics[ii]->nu = nu[ii];
		dims->dynamics[ii]->nx1 = nx[ii+1];
		dims->dynamics[ii]->nu1 = nu[ii+1];
		dims->dynamics[ii]->sim->nx = nx[ii];
		dims->dynamics[ii]->sim->nu = nu[ii];
	}

	// TODO cost dims initialize ???
	for (ii=0; ii<=N; ii++)
	{
		dims->cost[ii]->nx = nx[ii];
		dims->cost[ii]->nu = nu[ii];
		dims->cost[ii]->ny = ny[ii];
	}

	// TODO constraints dims initialize ???
	for (ii=0; ii<=N; ii++)
	{
		dims->constraints[ii]->nx = nx[ii];
		dims->constraints[ii]->nu = nu[ii];
		dims->constraints[ii]->nbx = nbx[ii];
		dims->constraints[ii]->nbu = nbu[ii];
		dims->constraints[ii]->nb = nbx[ii]+nbu[ii];
		dims->constraints[ii]->ng = ng[ii];
		dims->constraints[ii]->nh = nh[ii];
		dims->constraints[ii]->nq = nq[ii];
		dims->constraints[ii]->ns = ns[ii];
	}

	dims->qp_solver->N = N;
	for (ii=0; ii<=N; ii++)
	{
		dims->qp_solver->nx[ii] = nx[ii];
		dims->qp_solver->nu[ii] = nu[ii];
		dims->qp_solver->nbx[ii] = nbx[ii];
		dims->qp_solver->nbu[ii] = nbu[ii];
		dims->qp_solver->nb[ii] = nbx[ii]+nbu[ii];
		dims->qp_solver->ng[ii] = ng[ii] + nh[ii];
		dims->qp_solver->ns[ii] = ns[ii];
	}

	return;
}



/************************************************
* in
************************************************/

int ocp_nlp_in_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{

	int ii;

	int N = dims->N;

    int size = sizeof(ocp_nlp_in);

	size += N*sizeof(double); // Ts

	// dynamics
	size += N*sizeof(void *);
	for (ii=0; ii<N; ii++)
	{
		size += config->dynamics[ii]->model_calculate_size(config->dynamics[ii], dims->dynamics[ii]);
	}

	// cost
	size += (N+1)*sizeof(void *);
	for (ii=0; ii<=N; ii++)
	{
		size += config->cost[ii]->model_calculate_size(config->cost[ii], dims->cost[ii]);
	}

	// constraints
	size += (N+1)*sizeof(void *);
	for (ii=0; ii<=N; ii++)
	{
		size += config->constraints[ii]->model_calculate_size(config->constraints[ii], dims->constraints[ii]);
	}

	size += 8; // initial align

//	make_int_multiple_of(64, &size);

    return size;
}



ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory)
{

	int ii;

	int N = dims->N;

    char *c_ptr = (char *) raw_memory;

	// initial align
	align_char_to(8, &c_ptr);

	// struct
    ocp_nlp_in *in = (ocp_nlp_in *) c_ptr;
    c_ptr += sizeof(ocp_nlp_in);

	in->dims = dims;

	// Ts
	in->Ts = (double *) c_ptr;
	c_ptr += N*sizeof(double);

	// dynamics
	in->dynamics = (void **) c_ptr;
	c_ptr += N*sizeof(void *);
	for (ii=0; ii<N; ii++)
	{
		in->dynamics[ii] = config->dynamics[ii]->model_assign(config->dynamics[ii], dims->dynamics[ii], c_ptr);
		c_ptr += config->dynamics[ii]->model_calculate_size(config->dynamics[ii], dims->dynamics[ii]);
	}

	// cost
	in->cost = (void **) c_ptr;
	c_ptr += (N+1)*sizeof(void *);
	for (ii=0; ii<=N; ii++)
	{
		in->cost[ii] = config->cost[ii]->model_assign(config->cost[ii], dims->cost[ii], c_ptr);
		c_ptr += config->cost[ii]->model_calculate_size(config->cost[ii], dims->cost[ii]);
	}

	// constraints
	in->constraints = (void **) c_ptr;
	c_ptr += (N+1)*sizeof(void *);
	for (ii=0; ii<=N; ii++)
	{
		in->constraints[ii] = config->constraints[ii]->model_assign(config->constraints[ii], dims->constraints[ii], c_ptr);
		c_ptr += config->constraints[ii]->model_calculate_size(config->constraints[ii], dims->constraints[ii]);
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
	int nx, nu, nb, ng, nh, nx1;

    int size = sizeof(ocp_nlp_out);

	size += 3*(N+1)*sizeof(struct blasfeo_dvec); // ux lam
	size += 1*N*sizeof(struct blasfeo_dvec); // pi

    for (int ii = 0; ii < N; ii++)
    {
		nx = dims->constraints[ii]->nx;
		nx1 = dims->constraints[ii+1]->nx;
		nu = dims->constraints[ii]->nu;
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		size += 1*blasfeo_memsize_dvec(nu+nx); // ux
		size += 1*blasfeo_memsize_dvec(nx1); // pi
		size += 2*blasfeo_memsize_dvec(2*nb+2*ng+2*nh); // lam t
    }
	nx = dims->constraints[N]->nx;
	nu = dims->constraints[N]->nu;
	nb = dims->constraints[N]->nb;
	ng = dims->constraints[N]->ng;
	nh = dims->constraints[N]->nh;
	size += 1*blasfeo_memsize_dvec(nu+nx); // ux
	size += 2*blasfeo_memsize_dvec(2*nb+2*ng+2*nh); // lam t

	size += 8; // initial align
	size += 8; // blasfeo_struct align
	size += 64; // blasfeo_mem align

//	make_int_multiple_of(64, &size);

    return size;
}



ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory)
{

	// extract sizes
    int N = dims->N;
	int nx, nu, nb, ng, nh, nx1;

    char *c_ptr = (char *) raw_memory;

	// initial align
	align_char_to(8, &c_ptr);

    ocp_nlp_out *out = (ocp_nlp_out *)c_ptr;
    c_ptr += sizeof(ocp_nlp_out);

	// blasfeo_struct align
	align_char_to(8, &c_ptr);

	// blasfeo_dvec_struct
	// ux
	assign_and_advance_blasfeo_dvec_structs(N+1, &out->ux, &c_ptr);
	// pi
	assign_and_advance_blasfeo_dvec_structs(N, &out->pi, &c_ptr);
	// lam
	assign_and_advance_blasfeo_dvec_structs(N+1, &out->lam, &c_ptr);
	// t
	assign_and_advance_blasfeo_dvec_structs(N+1, &out->t, &c_ptr);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// blasfeo_dvec
	// ux
    for (int ii = 0; ii <= N; ++ii)
	{
		nx = dims->constraints[ii]->nx;
		nu = dims->constraints[ii]->nu;
		assign_and_advance_blasfeo_dvec_mem(nu+nx, out->ux+ii, &c_ptr);
	}
	// pi
    for (int ii = 0; ii < N; ++ii)
	{
		nx1 = dims->constraints[ii+1]->nx;
		assign_and_advance_blasfeo_dvec_mem(nx1, out->pi+ii, &c_ptr);
	}
	// lam
    for (int ii = 0; ii <= N; ++ii)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, out->lam+ii, &c_ptr);
	}
	// t
    for (int ii = 0; ii <= N; ++ii)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, out->t+ii, &c_ptr);
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
	int nx, nu, nb, ng, nh, nx1;

    int size = sizeof(ocp_nlp_memory);

	size += 4*(N+1)*sizeof(struct blasfeo_dvec); // cost_grad ineq_fun ineq_adj dyn_adj
	size += 1*N*sizeof(struct blasfeo_dvec); // dyn_fun

	for (int ii = 0; ii < N; ii++)
	{
		nx = dims->constraints[ii]->nx;
		nx1 = dims->constraints[ii+1]->nx;
		nu = dims->constraints[ii]->nu;
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		size += 3*blasfeo_memsize_dvec(nu+nx); // cost_grad ineq_adj dyn_adj
		size += 1*blasfeo_memsize_dvec(nx1); // dyn_fun
		size += 1*blasfeo_memsize_dvec(2*nb+2*ng+2*nh); // ineq_fun
	}
	nx = dims->constraints[N]->nx;
	nu = dims->constraints[N]->nu;
	nb = dims->constraints[N]->nb;
	ng = dims->constraints[N]->ng;
	nh = dims->constraints[N]->nh;
	size += 3*blasfeo_memsize_dvec(nu+nx); // cost_grad ineq_adj dyn_adj
	size += 1*blasfeo_memsize_dvec(2*nb+2*ng+2*nh); // ineq_fun


	size += 8; // initial align
	size += 8; // blasfeo_struct align
	size += 64; // blasfeo_mem align

//	make_int_multiple_of(64, &size);

	return size;

}



ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

	// extract sizes
    int N = dims->N;
	int nx, nu, nb, ng, nh, nx1;

	// initial align
	align_char_to(8, &c_ptr);

	// struct
    ocp_nlp_memory *mem = (ocp_nlp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_memory);

	// blasfeo_struct align
	align_char_to(8, &c_ptr);

	// cost_grad
	assign_and_advance_blasfeo_dvec_structs(N+1, &mem->cost_grad, &c_ptr);
	// ineq_fun
	assign_and_advance_blasfeo_dvec_structs(N+1, &mem->ineq_fun, &c_ptr);
	// ineq_adj
	assign_and_advance_blasfeo_dvec_structs(N+1, &mem->ineq_adj, &c_ptr);
	// dyn_fun
	assign_and_advance_blasfeo_dvec_structs(N, &mem->dyn_fun, &c_ptr);
	// dyn_adj
	assign_and_advance_blasfeo_dvec_structs(N+1, &mem->dyn_adj, &c_ptr);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// cost_grad
	for (int ii = 0; ii <= N; ii++)
	{
		nx = dims->constraints[ii]->nx;
		nu = dims->constraints[ii]->nu;
		assign_and_advance_blasfeo_dvec_mem(nu+nx, mem->cost_grad+ii, &c_ptr);
	}
	// ineq_fun
	for (int ii = 0; ii <= N; ii++)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, mem->ineq_fun+ii, &c_ptr);
	}
	// ineq_adj
	for (int ii = 0; ii <= N; ii++)
	{
		nx = dims->constraints[ii]->nx;
		nu = dims->constraints[ii]->nu;
		assign_and_advance_blasfeo_dvec_mem(nu+nx, mem->ineq_adj+ii, &c_ptr);
	}
	// dyn_fun
	for (int ii = 0; ii < N; ii++)
	{
		nx1 = dims->dynamics[ii]->nx1;
		assign_and_advance_blasfeo_dvec_mem(nx1, mem->dyn_fun+ii, &c_ptr);
	}
	// dyn_adj
	for (int ii = 0; ii <= N; ii++)
	{
		nx = dims->constraints[ii]->nx;
		nu = dims->constraints[ii]->nu;
//		nx = dims->dynamics[ii]->nx; // XXX not defined for N
//		nu = dims->dynamics[ii]->nu; // XXX not defined for N
		assign_and_advance_blasfeo_dvec_mem(nu+nx, mem->dyn_adj+ii, &c_ptr);
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
	int nx, nu, nb, ng, nh, nx1;

    int size = sizeof(ocp_nlp_res);

	size += 3*(N+1)*sizeof(struct blasfeo_dvec); // res_g res_d res_m
	size += 1*N*sizeof(struct blasfeo_dvec); // res_b

	for(int ii = 0; ii < N; ii++)
	{
		nx = dims->constraints[ii]->nx;
		nx1 = dims->constraints[ii+1]->nx;
		nu = dims->constraints[ii]->nu;
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		size += 1*blasfeo_memsize_dvec(nu+nx); // res_g
		size += 1*blasfeo_memsize_dvec(nx1); // res_b
		size += 2*blasfeo_memsize_dvec(2*nb+2*ng+2*nh); // res_d res_m
	}
	nx = dims->constraints[N]->nx;
	nu = dims->constraints[N]->nu;
	nb = dims->constraints[N]->nb;
	ng = dims->constraints[N]->ng;
	nh = dims->constraints[N]->nh;
	size += 1*blasfeo_memsize_dvec(nu+nx); // res_g
	size += 2*blasfeo_memsize_dvec(2*nb+2*ng+2*nh); // res_d res_m

	size += 8; // initial align
	size += 8; // blasfeo_struct align
	size += 64; // blasfeo_mem align

//	make_int_multiple_of(64, &size);

	return size;

}



ocp_nlp_res *ocp_nlp_res_assign(ocp_nlp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

	// extract sizes
    int N = dims->N;
	int nx, nu, nb, ng, nh, nx1;

	// initial align
	align_char_to(8, &c_ptr);

	// struct
    ocp_nlp_res *res = (ocp_nlp_res *) c_ptr;
    c_ptr += sizeof(ocp_nlp_res);

	// blasfeo_struct align
	align_char_to(8, &c_ptr);

	// res_g
	assign_and_advance_blasfeo_dvec_structs(N+1, &res->res_g, &c_ptr);
	// res_b
	assign_and_advance_blasfeo_dvec_structs(N, &res->res_b, &c_ptr);
	// res_d
	assign_and_advance_blasfeo_dvec_structs(N+1, &res->res_d, &c_ptr);
	// res_m
	assign_and_advance_blasfeo_dvec_structs(N+1, &res->res_m, &c_ptr);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// res_g
	for (int ii = 0; ii <= N; ii++)
	{
		nx = dims->cost[ii]->nx;
		nu = dims->cost[ii]->nu;
		assign_and_advance_blasfeo_dvec_mem(nu+nx, res->res_g+ii, &c_ptr);
	}
	// res_b
	for (int ii = 0; ii < N; ii++)
	{
		nx1 = dims->dynamics[ii]->nx1;
		assign_and_advance_blasfeo_dvec_mem(nx1, res->res_b+ii, &c_ptr);
	}
	// res_d
	for (int ii = 0; ii <= N; ii++)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, res->res_d+ii, &c_ptr);
	}
	// res_m
	for (int ii = 0; ii <= N; ii++)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, res->res_m+ii, &c_ptr);
	}

	res->memsize = ocp_nlp_res_calculate_size(dims);

	return res;

}



void ocp_nlp_res_compute(ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_res *res, ocp_nlp_memory *mem) //, ocp_nlp_res_workspace *work)
{

	// extract dims
    int N = dims->N;
	int nx, nu, nb, ng, nh, nx1;

	double tmp_res;

	// res_g
	res->inf_norm_res_g = 0.0;
	for (int ii = 0; ii <= N; ii++)
	{
		nx = dims->cost[ii]->nx;
		nu = dims->cost[ii]->nu;
		blasfeo_daxpy(nu+nx, -1.0, mem->dyn_adj+ii, 0, mem->cost_grad+ii, 0, res->res_g+ii, 0);
		blasfeo_daxpy(nu+nx, -1.0, mem->ineq_adj+ii, 0, res->res_g+ii, 0, res->res_g+ii, 0);
		blasfeo_dvecnrm_inf(nu+nx, res->res_g+ii, 0, &tmp_res);
		res->inf_norm_res_g = tmp_res>res->inf_norm_res_g ? tmp_res : res->inf_norm_res_g;
	}

	// res_b
	res->inf_norm_res_b = 0.0;
	for (int ii = 0; ii < N; ii++)
	{
		nx1 = dims->dynamics[ii]->nx1;
		blasfeo_dveccp(nx1, mem->dyn_fun+ii, 0, res->res_b+ii, 0);
		blasfeo_dvecnrm_inf(nx1, res->res_b+ii, 0, &tmp_res);
		res->inf_norm_res_b = tmp_res>res->inf_norm_res_b ? tmp_res : res->inf_norm_res_b;
	}

	// res_d
	res->inf_norm_res_d = 0.0;
	for (int ii = 0; ii <= N; ii++)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		blasfeo_daxpy(2*nb+2*ng+2*nh, 1.0, out->t+ii, 0, mem->ineq_fun+ii, 0, res->res_d+ii, 0);
		blasfeo_dvecnrm_inf(2*nb+2*ng+2*nh, res->res_d+ii, 0, &tmp_res);
		res->inf_norm_res_d = tmp_res>res->inf_norm_res_d ? tmp_res : res->inf_norm_res_d;
	}

	// res_m
	res->inf_norm_res_m = 0.0;
	for (int ii = 0; ii <= N; ii++)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nh = dims->constraints[ii]->nh;
		blasfeo_dvecmul(2*nb+2*ng+2*nh, out->lam+ii, 0, out->t+ii, 0, res->res_m+ii, 0);
		blasfeo_dvecnrm_inf(2*nb+2*ng+2*nh, res->res_m+ii, 0, &tmp_res);
		res->inf_norm_res_m = tmp_res>res->inf_norm_res_m ? tmp_res : res->inf_norm_res_m;
	}

	return;

}



