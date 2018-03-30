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

// external
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// acados
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_erk_integrator.h>
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/external_function_generic.h"

// crane model
#include "examples/c/crane_model/crane_model.h"

// blasfeo
#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>



int main()
{

/************************************************
* bla bla bla
************************************************/

    int NREP = 500;
    acados_timer timer;

	/* double Time1, Time2, Time3; */

    int ii;
    int jj;

    int nx = 4;
    int nu = 1;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05;
    // int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;

/************************************************
* external functions (explicit model)
************************************************/

	// forward explicit VDE

	external_function_casadi expl_vde_for;
	expl_vde_for.casadi_fun = &vdeFun;
	expl_vde_for.casadi_work = &vdeFun_work;
	expl_vde_for.casadi_sparsity_in = &vdeFun_sparsity_in;
	expl_vde_for.casadi_sparsity_out = &vdeFun_sparsity_out;
	expl_vde_for.casadi_n_in = &vdeFun_n_in;
	expl_vde_for.casadi_n_out = &vdeFun_n_out;

	int forw_vde_size = external_function_casadi_calculate_size(&expl_vde_for);
	void *forw_vde_mem = malloc(forw_vde_size);
	external_function_casadi_assign(&expl_vde_for, forw_vde_mem);

	// adjoint explicit VDE

	external_function_casadi expl_vde_adj;
	expl_vde_adj.casadi_fun = &adjFun;
	expl_vde_adj.casadi_work = &adjFun_work;
	expl_vde_adj.casadi_sparsity_in = &adjFun_sparsity_in;
	expl_vde_adj.casadi_sparsity_out = &adjFun_sparsity_out;
	expl_vde_adj.casadi_n_in = &adjFun_n_in;
	expl_vde_adj.casadi_n_out = &adjFun_n_out;

	int adj_vde_size = external_function_casadi_calculate_size(&expl_vde_adj);
	void *adj_vde_mem = malloc(adj_vde_size);
	external_function_casadi_assign(&expl_vde_adj, adj_vde_mem);

	// jacobian explicit ODE

	external_function_casadi expl_ode_jac;
	expl_ode_jac.casadi_fun = &jacFun;
	expl_ode_jac.casadi_work = &jacFun_work;
	expl_ode_jac.casadi_sparsity_in = &jacFun_sparsity_in;
	expl_ode_jac.casadi_sparsity_out = &jacFun_sparsity_out;
	expl_ode_jac.casadi_n_in = &jacFun_n_in;
	expl_ode_jac.casadi_n_out = &jacFun_n_out;

	int jac_size = external_function_casadi_calculate_size(&expl_ode_jac);
	void *jac_mem = malloc(jac_size);
	external_function_casadi_assign(&expl_ode_jac, jac_mem);

	// hessian explicit ODE

	external_function_casadi expl_ode_hes;
	expl_ode_hes.casadi_fun = &hessFun;
	expl_ode_hes.casadi_work = &hessFun_work;
	expl_ode_hes.casadi_sparsity_in = &hessFun_sparsity_in;
	expl_ode_hes.casadi_sparsity_out = &hessFun_sparsity_out;
	expl_ode_hes.casadi_n_in = &hessFun_n_in;
	expl_ode_hes.casadi_n_out = &hessFun_n_out;

	int hess_ode_size = external_function_casadi_calculate_size(&expl_ode_hes);
	void *hess_ode_mem = malloc(hess_ode_size);
	external_function_casadi_assign(&expl_ode_hes, hess_ode_mem);

/************************************************
* external functions (implicit model)
************************************************/

	// implicit ODE

	external_function_casadi impl_ode_fun;
	impl_ode_fun.casadi_fun = &casadi_impl_ode_fun;
	impl_ode_fun.casadi_work = &casadi_impl_ode_fun_work;
	impl_ode_fun.casadi_sparsity_in = &casadi_impl_ode_fun_sparsity_in;
	impl_ode_fun.casadi_sparsity_out = &casadi_impl_ode_fun_sparsity_out;
	impl_ode_fun.casadi_n_in = &casadi_impl_ode_fun_n_in;
	impl_ode_fun.casadi_n_out = &casadi_impl_ode_fun_n_out;

	int impl_ode_fun_size = external_function_casadi_calculate_size(&impl_ode_fun);
	void *impl_ode_fun_mem = malloc(impl_ode_fun_size);
	external_function_casadi_assign(&impl_ode_fun, impl_ode_fun_mem);

	//

	external_function_casadi impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_fun = &casadi_impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_work = &casadi_impl_ode_fun_jac_x_xdot_work;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_sparsity_in;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_sparsity_out;
	impl_ode_fun_jac_x_xdot.casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_n_in;
	impl_ode_fun_jac_x_xdot.casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_n_out;

	int impl_ode_fun_jac_x_xdot_size = external_function_casadi_calculate_size(&impl_ode_fun_jac_x_xdot);
	void *impl_ode_fun_jac_x_xdot_mem = malloc(impl_ode_fun_jac_x_xdot_size);
	external_function_casadi_assign(&impl_ode_fun_jac_x_xdot, impl_ode_fun_jac_x_xdot_mem);

	//

	external_function_casadi impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_fun = &casadi_impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_work = &casadi_impl_ode_jac_x_xdot_u_work;
	impl_ode_jac_x_xdot_u.casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_sparsity_in;
	impl_ode_jac_x_xdot_u.casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_sparsity_out;
	impl_ode_jac_x_xdot_u.casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_n_in;
	impl_ode_jac_x_xdot_u.casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_n_out;

	int impl_ode_jac_x_xdot_u_size = external_function_casadi_calculate_size(&impl_ode_jac_x_xdot_u);
	void *impl_ode_jac_x_xdot_u_mem = malloc(impl_ode_jac_x_xdot_u_size);
	external_function_casadi_assign(&impl_ode_jac_x_xdot_u, impl_ode_jac_x_xdot_u_mem);

	//

	external_function_casadi impl_ode_jac_x_u;
	impl_ode_jac_x_u.casadi_fun = &casadi_impl_ode_jac_x_u;
	impl_ode_jac_x_u.casadi_work = &casadi_impl_ode_jac_x_u_work;
	impl_ode_jac_x_u.casadi_sparsity_in = &casadi_impl_ode_jac_x_u_sparsity_in;
	impl_ode_jac_x_u.casadi_sparsity_out = &casadi_impl_ode_jac_x_u_sparsity_out;
	impl_ode_jac_x_u.casadi_n_in = &casadi_impl_ode_jac_x_u_n_in;
	impl_ode_jac_x_u.casadi_n_out = &casadi_impl_ode_jac_x_u_n_out;

	int impl_ode_jac_x_u_size = external_function_casadi_calculate_size(&impl_ode_jac_x_u);
	void *impl_ode_jac_x_u_mem = malloc(impl_ode_jac_x_u_size);
	external_function_casadi_assign(&impl_ode_jac_x_u, impl_ode_jac_x_u_mem);



	int number_sim_solvers = 3;
	int nss;
	for (nss=0; nss<number_sim_solvers; nss++)
	{

/************************************************
* sim config
************************************************/

		int config_size = sim_solver_config_calculate_size();
		void *config_mem = malloc(config_size);
		sim_solver_config *config = sim_solver_config_assign(config_mem);

		switch (nss)
		{
			case 0: // erk
				printf("\n\nsim solver: ERK\n");
				sim_erk_config_initialize_default(config);
				break;

			case 1: // irk
				printf("\n\nsim solver: IRK\n");
				sim_irk_config_initialize_default(config);
				break;

			case 2: // lifted_irk
				printf("\n\nsim solver: Lifted_IRK\n");
				sim_lifted_irk_config_initialize_default(config);
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);
		}

/************************************************
* sim dims
************************************************/

		int dims_size = sim_dims_calculate_size();
		void *dims_mem = malloc(dims_size);
		sim_dims *dims = sim_dims_assign(dims_mem);

		dims->nx = nx;
		dims->nu = nu;

/************************************************
* sim opts
************************************************/

		int opts_size = config->opts_calculate_size(config, dims);
		void *opts_mem = malloc(opts_size);
		sim_rk_opts *opts = config->opts_assign(config, dims, opts_mem);
		config->opts_initialize_default(config, dims, opts);

		switch (nss)
		{
			case 0: // erk
				opts->ns = 4;
				opts->sens_adj = true;
				break;

			case 1: // irk
				opts->ns = 2;
				opts->sens_adj = true;
				break;

			case 2: // lifted_irk
				opts->ns = 2;
				opts->sens_adj = true;
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);
		}
		// recompute Butcher tableau after selecting ns
		config->opts_update(config, dims, opts);


/************************************************
* sim memory
************************************************/

		int mem_size = config->memory_calculate_size(config, dims, opts);
		void *mem_mem = malloc(mem_size);
		void *mem = config->memory_assign(config, dims, opts, mem_mem);

/************************************************
* sim workspace
************************************************/

		int work_size = config->workspace_calculate_size(config, dims, opts);
		void *work = malloc(work_size);

/************************************************
* sim in
************************************************/

		int in_size = sim_in_calculate_size(config, dims);
		void *in_mem = malloc(in_size);
		sim_in *in = sim_in_assign(config, dims, in_mem);

		in->T = T;

		// external functions
		switch (nss)
		{
			case 0: // erk
			{
				config->model_set_function(in->model, EXPL_VDE_FOR, &expl_vde_for);
				config->model_set_function(in->model, EXPL_VDE_ADJ, &expl_vde_adj);
				config->model_set_function(in->model, EXPL_ODE_HES, &expl_ode_hes);
				break;
			}
			case 1: // irk
			{
				config->model_set_function(in->model, IMPL_ODE_FUN, &impl_ode_fun);
				config->model_set_function(in->model, IMPL_ODE_FUN_JAC_X_XDOT, &impl_ode_fun_jac_x_xdot);
				config->model_set_function(in->model, IMPL_ODE_JAC_X_XDOT_U, &impl_ode_jac_x_xdot_u);
				config->model_set_function(in->model, IMPL_ODE_JAC_X_U, &impl_ode_jac_x_u);
				break;
			}
			case 2: // lifted_irk
			{
				config->model_set_function(in->model, EXPL_ODE_JAC, &expl_ode_jac);
				config->model_set_function(in->model, EXPL_VDE_FOR, &expl_vde_for);
				break;
			}
			default :
			{
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);
			}
		}

		// x
		for (ii = 0; ii < nx; ii++) {
			in->x[ii] = xref[ii];
		}

		// p
		for (ii = 0;ii < nu; ii++){
			in->u[ii] = 1.0;
		}

		// seeds forw
		for (ii = 0; ii < nx * NF; ii++)
			in->S_forw[ii] = 0.0;
		for (ii = 0; ii < nx; ii++)
			in->S_forw[ii * (nx + 1)] = 1.0;

		// seeds adj
		for (ii = 0; ii < nx; ii++)
			in->S_adj[ii] = 1.0;

/************************************************
* sim out
************************************************/

		int out_size = sim_out_calculate_size(config, dims);
		void *out_mem = malloc(out_size);
		sim_out *out = sim_out_assign(config, dims, out_mem);

/************************************************
* sim solver
************************************************/

		acados_tic(&timer);

		for (ii=0;ii<NREP;ii++)
			config->evaluate(config, in, out, opts, mem, work);

		/* Time1 = acados_toc(&timer)/NREP; */

		double *xn = out->xn;

/************************************************
* printing
************************************************/

		printf("\nxn: \n");
		for (ii=0;ii<nx;ii++)
			printf("%8.5f ",xn[ii]);
		printf("\n");

		double *S_forw_out;
		S_forw_out = NULL;
		if(opts->sens_forw){
			S_forw_out = out->S_forw;
			printf("\nS_forw_out: \n");
			for (ii=0;ii<nx;ii++){
				for (jj=0;jj<NF;jj++)
					printf("%8.5f ", S_forw_out[jj*nx+ii]);
				printf("\n");
			}
		}

		double *S_adj_out;
		if(opts->sens_adj){
			S_adj_out = out->S_adj;
			printf("\nS_adj_out: \n");
			for (ii=0;ii<nx+nu;ii++){
				printf("%8.5f ", S_adj_out[ii]);
			}
			printf("\n");
		}

		double *S_hess_out;
		if(opts->sens_hess){
			double zero = 0.0;
			S_hess_out = out->S_hess;
			printf("\nS_hess_out: \n");
			for (ii=0;ii<NF;ii++){
				for (jj=0;jj<NF;jj++){
					if (jj>ii){
						printf("%8.5f ", zero);
					}else{
						printf("%8.5f ", S_hess_out[jj*NF+ii]);
					}
				}
				printf("\n");
			}
		}


		printf("\n");
		printf("cpt: %8.4f [ms]\n", 1000*out->info->CPUtime);
		printf("AD cpt: %8.4f [ms]\n", 1000*out->info->ADtime);

		if(opts->sens_adj){
			struct blasfeo_dmat sA;
			blasfeo_allocate_dmat(nx, nx+nu, &sA);
			blasfeo_pack_dmat(nx, nx+nu, S_forw_out, nx, &sA, 0, 0);

			struct blasfeo_dvec sx;
			blasfeo_allocate_dvec(nx, &sx);
			blasfeo_pack_dvec(nx, in->S_adj, &sx, 0);

			struct blasfeo_dvec sz;
			blasfeo_allocate_dvec(nx+nu, &sz);
//			blasfeo_print_dmat(nx, nx+nu, &sA, 0, 0);
//			blasfeo_print_tran_dvec(nx, &sx, 0);
			blasfeo_dgemv_t(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);

			printf("\nJac times lambdaX:\n");
			blasfeo_print_tran_dvec(nx+nu, &sz, 0);

			blasfeo_free_dmat(&sA);
			blasfeo_free_dvec(&sx);
			blasfeo_free_dvec(&sz);
		}

/************************************************
* free
************************************************/

		free(config_mem);
		free(dims_mem);
		free(opts_mem);
		free(mem_mem);
		free(work);
		free(in_mem);
		free(out_mem);

	}

	// explicit model
	free(forw_vde_mem);
	free(adj_vde_mem);
	free(hess_ode_mem);
	// implicit model
	free(impl_ode_fun_mem);
	free(impl_ode_fun_jac_x_xdot_mem);
	free(impl_ode_jac_x_xdot_u_mem);
	free(impl_ode_jac_x_u_mem);

/************************************************
* return
************************************************/

	printf("\nsuccess!\n\n");

    return 0;
}
