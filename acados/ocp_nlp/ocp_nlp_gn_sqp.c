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

#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"

// external
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_casadi_wrapper.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_collocation.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"
#include "acados/utils/casadi_wrapper.h"



static int get_max_sim_workspace_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
    sim_dims sim_dims;
    int sim_work_size;

    int max_sim_work_size = 0;

    for (int ii = 0; ii < dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        // sim_in_size = sim_in_calculate_size(&sim_dims);
        // if (sim_in_size > *max_sim_in_size) *max_sim_in_size = sim_in_size;
        // sim_out_size = sim_out_calculate_size(&sim_dims);
        // if (sim_out_size > *max_sim_out_size) *max_sim_out_size = sim_out_size;
        sim_work_size = args->sim_solvers[ii]->calculate_workspace_size(&sim_dims, args->sim_solvers_args[ii]);
        if (sim_work_size > max_sim_work_size) max_sim_work_size = sim_work_size;
    }
    return max_sim_work_size;
}



/************************************************
* arguments
************************************************/



int ocp_nlp_gn_sqp_calculate_args_size(ocp_nlp_dims *dims, ocp_qp_xcond_solver_fcn_ptrs *qp_solver, sim_solver_fcn_ptrs *sim_solvers)
{
    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);
    size += sizeof(ocp_nlp_gn_sqp_args);
    size += sizeof(ocp_qp_xcond_solver_fcn_ptrs);

    size += qp_solver->calculate_args_size(&qp_dims, qp_solver->qp_solver);

    sim_dims sim_dims;

    size += dims->N*sizeof(sim_solver_fcn_ptrs *);
    size += dims->N*sizeof(void *);  //sim_solvers_args

    for (int ii = 0; ii < dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        size += sizeof(sim_solver_fcn_ptrs);
        size += sim_solvers[ii].calculate_args_size(&sim_dims);
    }

    return size;
}



ocp_nlp_gn_sqp_args *ocp_nlp_gn_sqp_assign_args(ocp_nlp_dims *dims, ocp_qp_xcond_solver_fcn_ptrs *qp_solver, sim_solver_fcn_ptrs *sim_solvers, void *raw_memory)
{
    ocp_nlp_gn_sqp_args *args;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    char *c_ptr = (char *) raw_memory;

    args = (ocp_nlp_gn_sqp_args *) c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_args);

    args->qp_solver = (ocp_qp_xcond_solver_fcn_ptrs*) c_ptr;
    c_ptr += sizeof(ocp_qp_xcond_solver_fcn_ptrs);

    // copy function pointers
    *args->qp_solver = *qp_solver;

    args->qp_solver_args = args->qp_solver->assign_args(&qp_dims, qp_solver->qp_solver, c_ptr);
    c_ptr += args->qp_solver->calculate_args_size(&qp_dims, qp_solver->qp_solver);

    sim_dims sim_dims;

    args->sim_solvers = (sim_solver_fcn_ptrs **) c_ptr;
    c_ptr += dims->N*sizeof(sim_solver_fcn_ptrs *);

    args->sim_solvers_args = (void **) c_ptr;
    c_ptr += dims->N*sizeof(void *);

    for (int ii = 0; ii < dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);

        args->sim_solvers[ii] = (sim_solver_fcn_ptrs *) c_ptr;
        c_ptr += sizeof(sim_solver_fcn_ptrs);

        // copy function pointers
        *args->sim_solvers[ii] = sim_solvers[ii];

        args->sim_solvers_args[ii] = args->sim_solvers[ii]->assign_args(&sim_dims, c_ptr);
        c_ptr += args->sim_solvers[ii]->calculate_args_size(&sim_dims);
    }

	// default arguments
	args->maxIter = 20;
	args->min_res_g = 1e-12;
	args->min_res_b = 1e-12;
	args->min_res_d = 1e-12;
	args->min_res_m = 1e-12;

    assert((char*)raw_memory + ocp_nlp_gn_sqp_calculate_args_size(dims, qp_solver, sim_solvers) == c_ptr);

    return args;
}



/************************************************
* memory
************************************************/



int ocp_nlp_gn_sqp_calculate_memory_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
	// loop index
	int ii;

	// extract dims
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;

    int size = 0;

	// ???
    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    size += sizeof(ocp_nlp_gn_sqp_memory);

    size += args->qp_solver->calculate_memory_size(&qp_dims, args->qp_solver_args);

    sim_dims sim_dims;

    size += N*sizeof(void *);  // sim_solvers_mem

    for (int ii = 0; ii < N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        size += args->sim_solvers[ii]->calculate_memory_size(&sim_dims, args->sim_solvers_args[ii]);
    }

	// nlp res
	size += ocp_nlp_res_calculate_size(dims);

	// nlp mem
	size += ocp_nlp_mem_calculate_size(dims);

    size += 8; // initial align

    return size;
}



ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_assign_memory(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

	// loop index
	int ii;

	// extract dims
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

	// initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_memory);

    // QP solver
    mem->qp_solver_mem = args->qp_solver->assign_memory(&qp_dims, args->qp_solver_args, c_ptr);
    c_ptr += args->qp_solver->calculate_memory_size(&qp_dims, args->qp_solver_args);

    // integrators
    sim_dims sim_dims;

    mem->sim_solvers_mem = (void **) c_ptr;
    c_ptr += N*sizeof(void *);

    for (ii = 0; ii < N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        mem->sim_solvers_mem[ii] = args->sim_solvers[ii]->assign_memory(&sim_dims, args->sim_solvers_args[ii], c_ptr);
        c_ptr += args->sim_solvers[ii]->calculate_memory_size(&sim_dims, args->sim_solvers_args[ii]);
    }

	// nlp res
	mem->nlp_res = ocp_nlp_res_assign(dims, c_ptr);
	c_ptr += mem->nlp_res->memsize;

	// nlp mem
	mem->nlp_mem = ocp_nlp_mem_assign(dims, c_ptr);
	c_ptr += mem->nlp_mem->memsize;

	// dims
    mem->dims = dims;

    assert((char *)raw_memory + ocp_nlp_gn_sqp_calculate_memory_size(dims, args) >= c_ptr);

    return mem;
}



/************************************************
* workspace
************************************************/



int ocp_nlp_gn_sqp_calculate_workspace_size(ocp_nlp_dims *dims, ocp_nlp_gn_sqp_args *args)
{
	// loop index
	int ii;

	// extract dims
	int N = dims->N;
	int *nx = dims->nx;
	int *nu = dims->nu;
	int *nb = dims->nb;
	int *ng = dims->ng;
	int *ns = dims->ns;
	ocp_nlp_cost_ls_dims *cost_dims = dims->cost_dims;
	int *nv = cost_dims->nv;
	int *ny = cost_dims->ny;

    int size = 0;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, dims);

    size += sizeof(ocp_nlp_gn_sqp_work);

    size += ocp_qp_in_calculate_size(&qp_dims);
    size += ocp_qp_out_calculate_size(&qp_dims);
    size += args->qp_solver->calculate_workspace_size(&qp_dims, args->qp_solver_args);

    sim_dims sim_dims;

    size += N*sizeof(sim_in *);
    size += N*sizeof(sim_out *);
    size += N*sizeof(void *);  // sim_work

    size += get_max_sim_workspace_size(dims, args);

    for (ii = 0; ii < N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, dims, ii);
        size += sim_in_calculate_size(&sim_dims);
        size += sim_out_calculate_size(&sim_dims);
    }


	// temporary stuff
    size += 2*(N+1)*sizeof(struct blasfeo_dvec); // tmp_nux tmp_nbg
    size += 2*(N+1)*sizeof(struct blasfeo_dmat); // tmp_ny_ny tmp_nv_ny
	size += 3*(N+1)*sizeof(double *); // d_tmp_nv d_tmp_ny_nv_1 d_tmp_ny_nv_2

    for (ii = 0; ii < N+1; ii++)
    {
        size += blasfeo_memsize_dvec(nx[ii]+nu[ii]); // tmp_nux
        size += blasfeo_memsize_dvec(nb[ii]+ng[ii]); // tmp_nbg
        size += blasfeo_memsize_dmat(ny[ii], ny[ii]); // tmp_ny_ny
        size += blasfeo_memsize_dmat(nv[ii], ny[ii]); // tmp_nv_ny
		size += 1*nv[ii]*sizeof(double); // d_tmp_nv
		size += 2*(ny[ii]+ny[ii]*nv[ii])*sizeof(double); // d_tmp_ny_nv_1 d_tmp_ny_nv_2
    }

    size += 8; // blasfeo_struct align
    size += 64; // blasfeo_mem align

//    make_int_multiple_of(64, &size);

    return size;
}



void ocp_nlp_gn_sqp_cast_workspace(ocp_nlp_gn_sqp_work *work, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_args *args)
{
    char *c_ptr = (char *)work;
    c_ptr += sizeof(ocp_nlp_gn_sqp_work);

	// extract dims
    int N = mem->dims->N;
    int *nx = mem->dims->nx;
    int *nu = mem->dims->nu;
    int *nb = mem->dims->nb;
    int *ng = mem->dims->ng;
	ocp_nlp_cost_ls_dims *cost_dims = mem->dims->cost_dims;
	int *nv = cost_dims->nv;
	int *ny = cost_dims->ny;

    ocp_qp_dims qp_dims;
    cast_nlp_dims_to_qp_dims(&qp_dims, mem->dims);

	// blasfeo_struct align
    align_char_to(8, &c_ptr);

	assign_double_ptrs(N+1, &work->d_tmp_nv, &c_ptr);
	assign_double_ptrs(N+1, &work->d_tmp_ny_nv_1, &c_ptr);
	assign_double_ptrs(N+1, &work->d_tmp_ny_nv_2, &c_ptr);

	for (int ii=0; ii<=N; ii++)
		assign_double(nv[ii], work->d_tmp_nv+ii, &c_ptr);
	for (int ii=0; ii<=N; ii++)
		assign_double(ny[ii]+ny[ii]*nv[ii], work->d_tmp_ny_nv_1+ii, &c_ptr);
	for (int ii=0; ii<=N; ii++)
		assign_double(ny[ii]+ny[ii]*nv[ii], work->d_tmp_ny_nv_2+ii, &c_ptr);

    // set up local SQP data
    assign_blasfeo_dvec_structs(N+1, &work->tmp_nux, &c_ptr);
    assign_blasfeo_dvec_structs(N+1, &work->tmp_nbg, &c_ptr);
    assign_blasfeo_dmat_structs(N+1, &work->tmp_ny_ny, &c_ptr);
    assign_blasfeo_dmat_structs(N+1, &work->tmp_nv_ny, &c_ptr);

	// blasfeo_mem align
    align_char_to(64, &c_ptr);

	// tmp_ny_ny
    for (int ii = 0; ii <= N; ii++)
        assign_blasfeo_dmat_mem(ny[ii], ny[ii], work->tmp_ny_ny+ii, &c_ptr);
	// tmp_nv_ny
    for (int ii = 0; ii <= N; ii++)
        assign_blasfeo_dmat_mem(nv[ii], ny[ii], work->tmp_nv_ny+ii, &c_ptr);
	// tmp_nux
    for (int ii = 0; ii <= N; ii++)
        assign_blasfeo_dvec_mem(nx[ii]+nu[ii], work->tmp_nux+ii, &c_ptr);
	// tmp_nbg
    for (int ii = 0; ii <= N; ii++)
        assign_blasfeo_dvec_mem(nb[ii]+ng[ii], work->tmp_nbg+ii, &c_ptr);

    // set up QP solver
    work->qp_in = assign_ocp_qp_in(&qp_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(&qp_dims);
    work->qp_out = assign_ocp_qp_out(&qp_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(&qp_dims);

    work->qp_work = (void *)c_ptr;
    c_ptr += args->qp_solver->calculate_workspace_size(&qp_dims, args->qp_solver_args);

    // set up integrators
    sim_dims sim_dims;

    work->sim_in = (sim_in **) c_ptr;
    c_ptr += mem->dims->N*sizeof(sim_in *);
    work->sim_out = (sim_out **) c_ptr;
    c_ptr += mem->dims->N*sizeof(sim_out *);
    work->sim_solvers_work = (void **) c_ptr;
    c_ptr += mem->dims->N*sizeof(void *);

    int max_sim_work_size = get_max_sim_workspace_size(mem->dims, args);

    work->sim_solvers_work[0] = (void *)c_ptr;
    c_ptr += max_sim_work_size;

    for (int ii = 0; ii < mem->dims->N; ii++)
    {
        cast_nlp_dims_to_sim_dims(&sim_dims, mem->dims, ii);

        work->sim_in[ii] = assign_sim_in(&sim_dims, c_ptr);
        c_ptr += sim_in_calculate_size(&sim_dims);
        work->sim_out[ii] = assign_sim_out(&sim_dims, c_ptr);
        c_ptr += sim_out_calculate_size(&sim_dims);

        if (ii > 0) work->sim_solvers_work[ii] = work->sim_solvers_work[0];
    }

    assert((char *)work + ocp_nlp_gn_sqp_calculate_workspace_size(mem->dims, args) >= c_ptr);
}



/************************************************
* solver
************************************************/



static void initialize_objective(ocp_nlp_in *nlp_in, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *gn_sqp_mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i;

    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
	ocp_nlp_cost_ls_dims *cost_dims = nlp_in->dims->cost_dims;
	int *nv = cost_dims->nv;
	int *ny = cost_dims->ny;

    ocp_nlp_cost_ls *cost = (ocp_nlp_cost_ls*) nlp_in->cost;

	struct blasfeo_dmat *RSQrq = work->qp_in->RSQrq;

    for (i = 0; i <= N; i++)
	{

#if 0

		// identity Cyt
		blasfeo_dgecp(nu[i]+nx[i], nu[i]+nx[i], cost->W+i, 0, 0, RSQrq+i, 0, 0);

#else
		
		// general Cyt

		// TODO recompute factorization only if W are re-tuned ???
		blasfeo_dpotrf_l(ny[i], cost->W+i, 0, 0, work->tmp_ny_ny+i, 0, 0); // TODO move in the memory

		// linear ls
		// TODO avoid recomputing the Hessian if both W and Cyt do not change
		if (cost->nls_mask[i]==0)
		{
			blasfeo_dtrmm_rlnn(nv[i], ny[i], 1.0, work->tmp_ny_ny+i, 0, 0, cost->Cyt+i, 0, 0, work->tmp_nv_ny+i, 0, 0);
			blasfeo_dsyrk_ln(nv[i], ny[i], 1.0, work->tmp_nv_ny+i, 0, 0, work->tmp_nv_ny+i, 0, 0, 0.0, RSQrq+i, 0, 0, RSQrq+i, 0, 0);
		}

#endif

    }

	return;

}



static void initialize_constraints(ocp_nlp_in *nlp_in, ocp_nlp_gn_sqp_memory *gn_sqp_mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i, j;

    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    int *nb = nlp_in->dims->nb;
    int *ng = nlp_in->dims->ng;

	struct blasfeo_dmat *DCt = work->qp_in->DCt;
	int **idxb = work->qp_in->idxb;

	// initialize idxb
	for (i=0; i<=N; i++)
	{
		for (j=0; j<nb[i]; j++)
		{
			idxb[i][j] = nlp_in->idxb[i][j];
		}
	}

	// initialize general constraints matrix
    for (i=0; i<=N; i++)
	{
		blasfeo_dgecp(nu[i]+nx[i], ng[i], nlp_in->DCt+i, 0, 0, DCt+i, 0, 0);
	}

	return;

}



static void multiple_shooting(ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i;

	// extract dims
    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    int *nb = nlp_in->dims->nb;
    int *ng = nlp_in->dims->ng;
	ocp_nlp_cost_ls_dims *cost_dims = nlp_in->dims->cost_dims;
	int *nv = cost_dims->nv;
	int *ny = cost_dims->ny;

    struct blasfeo_dvec *tmp_nux = work->tmp_nux;
    struct blasfeo_dvec *tmp_nbg = work->tmp_nbg;

    ocp_nlp_cost_ls *cost = (ocp_nlp_cost_ls *) nlp_in->cost;
    struct blasfeo_dmat *Cyt = cost->Cyt;
    struct blasfeo_dvec *y_ref = cost->y_ref;

    struct blasfeo_dmat *BAbt = work->qp_in->BAbt;
    struct blasfeo_dmat *RSQrq = work->qp_in->RSQrq;
	struct blasfeo_dmat *DCt = work->qp_in->DCt;

	ocp_nlp_mem *nlp_mem = mem->nlp_mem;

	// initial stages
    for (i = 0; i <= N; i++)
    {

		// dynamics

		// XXX nlp mem: dyn_adj
		blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_mem->dyn_adj+i, 0);

		if (i<N)
		{
			// pass state and control to integrator
			blasfeo_unpack_dvec(nu[i], nlp_out->ux+i, 0, work->sim_in[i]->u);
			blasfeo_unpack_dvec(nx[i], nlp_out->ux+i, nu[i], work->sim_in[i]->x);

			// call integrator
			args->sim_solvers[i]->fun(work->sim_in[i], work->sim_out[i], args->sim_solvers_args[i],
				mem->sim_solvers_mem[i], work->sim_solvers_work[i]);

			// TODO(rien): transition functions for changing dimensions not yet implemented!

			// B
			blasfeo_pack_tran_dmat(nx[i+1], nu[i], &work->sim_out[i]->S_forw[nx[i+1]*nx[i]], nx[i+1], &BAbt[i], 0, 0);
			// A
			blasfeo_pack_tran_dmat(nx[i+1], nx[i], &work->sim_out[i]->S_forw[0], nx[i+1], &BAbt[i], nu[i], 0);

			// XXX nlp mem: dyn_fun
			blasfeo_pack_dvec(nx[i+1], work->sim_out[i]->xn, nlp_mem->dyn_fun+i, 0);
			blasfeo_daxpy(nx[i+1], -1.0, nlp_out->ux+i+1, nu[i+1], nlp_mem->dyn_fun+i, 0, nlp_mem->dyn_fun+i, 0);
			// XXX nlp mem: dyn_adj
			blasfeo_dgemv_n(nu[i]+nx[i], nx[i+1], -1.0, BAbt+i, 0, 0, nlp_out->pi+i, 0, 0.0, nlp_mem->dyn_adj+i, 0, nlp_mem->dyn_adj+i, 0);
		}

		// XXX nlp mem: dyn_adj
		if(i>0)
			blasfeo_daxpy(nx[i], 1.0, nlp_out->pi+i-1, 0, nlp_mem->dyn_adj+i, nu[i], nlp_mem->dyn_adj+i, nu[i]);



		// constraints
		// TODO merge dgemv_n and dgemv_t

		// XXX nlp_mem: ineq_fun
		blasfeo_dvecex_sp(nb[i], 1.0, nlp_in->idxb[i], nlp_out->ux+i, 0, tmp_nbg+i, 0);
		blasfeo_dgemv_t(nu[i]+nx[i], ng[i], 1.0, DCt+i, 0, 0, nlp_out->ux+i, 0, 0.0, tmp_nbg+i, nb[i], tmp_nbg+i, nb[i]);
		blasfeo_daxpy(nb[i]+ng[i], -1.0, tmp_nbg+i, 0, nlp_in->d+i, 0, nlp_mem->ineq_fun+i, 0);
		blasfeo_daxpy(nb[i]+ng[i], -1.0, nlp_in->d+i, nb[i]+ng[i], tmp_nbg+i, 0, nlp_mem->ineq_fun+i, nb[i]+ng[i]);

		// XXX nlp_mem: ineq_adj
		blasfeo_dvecse(nu[i]+nx[i], 0.0, nlp_mem->ineq_adj+i, 0);
		blasfeo_daxpy(nb[i]+ng[i], -1.0, nlp_out->lam+i, nb[i]+ng[i], nlp_out->lam+i, 0, tmp_nbg+i, 0);
		blasfeo_dvecad_sp(nb[i], 1.0, tmp_nbg+i, 0, nlp_in->idxb[i], nlp_mem->ineq_adj+i, 0);
		blasfeo_dgemv_n(nu[i]+nx[i], ng[i], 1.0, DCt+i, 0, 0, tmp_nbg+i, nb[i], 1.0, nlp_mem->ineq_adj+i, 0, nlp_mem->ineq_adj+i, 0);



        // cost
		// general Cyt

		// nonlinear ls
		if (cost->nls_mask[i]!=0)
		{
			// TODO evaluate nonlinear constraints into Cyt
			blasfeo_unpack_dvec(nu[i], nlp_out->ux+i, 0, work->d_tmp_nv[i]+nx[i]);
			blasfeo_unpack_dvec(nx[i], nlp_out->ux+i, nu[i], work->d_tmp_nv[i]);
			ls_cost_fun(nx[i], nu[i], work->d_tmp_nv[i], work->d_tmp_ny_nv_1[i], cost->nls_cost[i]);
			densify(work->d_tmp_ny_nv_1[i]+nv[i], work->d_tmp_ny_nv_2[i]+nv[i], cost->nls_cost_sparsity_jac[i]);
			blasfeo_pack_tran_dmat(ny[i], nv[i], work->d_tmp_ny_nv_2[i]+ny[i], ny[i], cost->Cyt+i, 0, 0);

			// TODO do something with the function evaluation !!!!!!!!!!!!!!!!!!!!!!!!!!!

			blasfeo_dtrmm_rlnn(nv[i], ny[i], 1.0, work->tmp_ny_ny+i, 0, 0, cost->Cyt+i, 0, 0, work->tmp_nv_ny+i, 0, 0);
			blasfeo_dsyrk_ln(nv[i], ny[i], 1.0, work->tmp_nv_ny+i, 0, 0, work->tmp_nv_ny+i, 0, 0, 0.0, RSQrq+i, 0, 0, RSQrq+i, 0, 0);
		}

		// TODO add nlp_mem cost_fun !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		// XXX nlp_mem: cost_grad

		blasfeo_dgemv_t(nv[i], ny[i], 1.0, Cyt+i, 0, 0, nlp_out->ux+i, 0, -1.0, y_ref+i, 0, tmp_nux+i, 0);

        blasfeo_dsymv_l(nu[i]+nx[i], nu[i]+nx[i], 1.0, RSQrq+i, 0, 0, tmp_nux+i, 0, 0.0, nlp_mem->cost_grad+i, 0, nlp_mem->cost_grad+i, 0);

		// XXX move to qp update ???
		if(i<N)
		{
			sim_rk_opts *opts = (sim_rk_opts*) args->sim_solvers_args[i];
			if (opts->scheme != NULL && opts->scheme->type != exact)
			{
				for (int_t j = 0; j < nx[i]; j++)
					DVECEL_LIBSTR(nlp_mem->cost_grad+i, nu[i]+j) += work->sim_out[i]->grad[j];
				for (int_t j = 0; j < nu[i]; j++)
					DVECEL_LIBSTR(nlp_mem->cost_grad+i, j) += work->sim_out[i]->grad[nx[i]+j];
			}
		}

    }

	// return
	return;

}



// update QP rhs for SQP (step prim var, abs dual var)
static void sqp_update_qp_vectors(ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i;

	// extract dims
    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    int *nb = nlp_in->dims->nb;
    int *ng = nlp_in->dims->ng;

    struct blasfeo_dmat *RSQrq = work->qp_in->RSQrq;
    struct blasfeo_dvec *rq = work->qp_in->rq;
    struct blasfeo_dmat *BAbt = work->qp_in->BAbt;
    struct blasfeo_dvec *b = work->qp_in->b;
    struct blasfeo_dvec *d = work->qp_in->d;

	ocp_nlp_mem *nlp_mem = mem->nlp_mem;

	// g
	for (i=0; i<=N; i++)
	{

		blasfeo_dveccp(nu[i]+nx[i], nlp_mem->cost_grad+i, 0, rq+i, 0);
        blasfeo_drowin(nu[i]+nx[i], 1.0, rq+i, 0, RSQrq+i, nu[i]+nx[i], 0); // XXX needed ???

	}

	// b
	for (i=0; i<N; i++)
	{
		blasfeo_dveccp(nx[i+1], nlp_mem->dyn_fun+i, 0, b+i, 0);
		blasfeo_drowin(nx[i+1], 1.0, b+i, 0, BAbt+i, nu[i]+nx[i], 0); // XXX needed ???
	}

	// d
	for (i=0; i<=N; i++)
	{

		blasfeo_dveccp(2*nb[i]+2*ng[i], nlp_mem->ineq_fun+i, 0, d+i, 0);

	}

	return;

}



static void update_variables(const ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i, j;

	// extract dims
    int N = nlp_out->dims->N;
    int *nx = nlp_out->dims->nx;
    int *nu = nlp_out->dims->nu;
    int *nb = nlp_out->dims->nb;
    int *ng = nlp_out->dims->ng;


    for (i = 0; i < N; i++)
    {
        for (j = 0; j < nx[i+1]; j++)
        {
            work->sim_in[i]->S_adj[j] = -DVECEL_LIBSTR(&work->qp_out->pi[i], j);
        }
    }

	// (full) step in primal variables
	for (i=0; i<=N; i++)
		blasfeo_daxpy(nu[i]+nx[i], 1.0, work->qp_out->ux+i, 0, nlp_out->ux+i, 0, nlp_out->ux+i, 0);

	// absolute in dual variables
	for (i=0; i<N; i++)
		blasfeo_dveccp(nx[i+1], work->qp_out->pi+i, 0, nlp_out->pi+i, 0);

	for (i=0; i<=N; i++)
		blasfeo_dveccp(2*nb[i]+2*ng[i], work->qp_out->lam+i, 0, nlp_out->lam+i, 0);

	for (i=0; i<=N; i++)
		blasfeo_dveccp(2*nb[i]+2*ng[i], work->qp_out->t+i, 0, nlp_out->t+i, 0);

	return;

}



// Simple fixed-step Gauss-Newton based SQP routine
int ocp_nlp_gn_sqp(ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_args *args, ocp_nlp_gn_sqp_memory *mem, void *work_)
{
    ocp_nlp_gn_sqp_work *work = (ocp_nlp_gn_sqp_work*) work_;
    ocp_nlp_gn_sqp_cast_workspace(work, mem, args);

    int N = nlp_in->dims->N;
    int *nx = nlp_in->dims->nx;
    int *nu = nlp_in->dims->nu;
    int *nb = nlp_in->dims->nb;
    int *ng = nlp_in->dims->ng;

    sim_rk_opts *sim_opts;

    // set up integrators
	ocp_nlp_model_expl *model = (ocp_nlp_model_expl *) nlp_in->model;
    for (int ii = 0; ii < N; ii++)
    {
        sim_opts = args->sim_solvers_args[ii];
        work->sim_in[ii]->step = sim_opts->interval/sim_opts->num_steps;

        work->sim_in[ii]->vde = model->vde[ii];
        work->sim_in[ii]->jac = model->jac[ii];
        work->sim_in[ii]->vde_adj = model->vde_adj[ii];
        work->sim_in[ii]->forward_vde_wrapper = &vde_fun;
        work->sim_in[ii]->jacobian_wrapper = &jac_fun;
        work->sim_in[ii]->adjoint_vde_wrapper = &vde_hess_fun;

        // TODO(dimitris): REVISE IF THIS IS CORRECT FOR VARYING DIMENSIONS!
        for (int jj = 0; jj < nx[ii+1] * (nx[ii] + nu[ii]); jj++)
            work->sim_in[ii]->S_forw[jj] = 0.0;
        for (int jj = 0; jj < nx[ii+1]; jj++)
            work->sim_in[ii]->S_forw[jj * (nx[ii] + 1)] = 1.0;
        for (int jj = 0; jj < nx[ii] + nu[ii]; jj++)
            work->sim_in[ii]->S_adj[jj] = 0.0;
        // for (int jj = 0; jj < nlp_in->dims->num_stages[ii] * nx[ii+1]; jj++)
            // work->sim_in[ii]->grad_K[jj] = 0.0;
    }

    initialize_objective(nlp_in, args, mem, work);

    initialize_constraints(nlp_in, mem, work);

	// start timer
    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);

	// main sqp loop
    int max_sqp_iterations =  args->maxIter;
	int sqp_iter = 0;
    for ( ; sqp_iter < max_sqp_iterations; sqp_iter++)
    {

        multiple_shooting(nlp_in, nlp_out, args, mem, work);

		// update QP rhs for SQP (step prim var, abs dual var)
        sqp_update_qp_vectors(nlp_in, nlp_out, args, mem, work);

		// compute nlp residuals
		ocp_nlp_res_compute(nlp_in, nlp_out, mem->nlp_res, mem->nlp_mem);

		// TODO exit conditions on residuals
		if( mem->nlp_res->inf_norm_res_g < args->min_res_g &
			mem->nlp_res->inf_norm_res_b < args->min_res_b &
			mem->nlp_res->inf_norm_res_d < args->min_res_d &
			mem->nlp_res->inf_norm_res_m < args->min_res_m )
		{

			// save sqp iterations number
			mem->sqp_iter = sqp_iter;

			// stop timer
			total_time += acados_toc(&timer);

			return 0;

		}

//print_ocp_qp_in(work->qp_in);
//exit(1);

        int_t qp_status = args->qp_solver->fun(work->qp_in, work->qp_out,
            args->qp_solver_args, mem->qp_solver_mem, work->qp_work);

//print_ocp_qp_out(work->qp_out);
//exit(1);

        if (qp_status != 0)
        {
            printf("QP solver returned error status %d\n", qp_status);
            return -1;
        }

        update_variables(nlp_out, args, mem, work);

//ocp_nlp_dims_print(nlp_out->dims);
//ocp_nlp_out_print(nlp_out);
//exit(1);

        for (int_t i = 0; i < N; i++)
        {
            sim_rk_opts *opts = (sim_rk_opts*) args->sim_solvers_args[i];
            if (opts->scheme == NULL)
                continue;
            opts->sens_adj = (opts->scheme->type != exact);
            if (nlp_in->freezeSens) {
                // freeze inexact sensitivities after first SQP iteration !!
                opts->scheme->freeze = true;
            }
        }

    }

	// stop timer
    total_time += acados_toc(&timer);

//	ocp_nlp_out_print(nlp_out);

	// save sqp iterations number
	mem->sqp_iter = sqp_iter;

	// maximum number of iterations reached
    return 1;

}
