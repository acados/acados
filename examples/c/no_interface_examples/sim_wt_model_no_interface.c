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
#include "acados/utils/external_function_generic.h"

// wt model
#include "examples/c/wt_model_nx3/wt_model.h"

// blasfeo
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>



int main()
{

/************************************************
* bla bla bla
************************************************/

    int NREP = 1;//500;
    acados_timer timer;

	/* double Time1, Time2, Time3; */

    int ii;
    int jj;

    int nx = 3;
    int nu = 4;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05; // simulation time
    // int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;

/************************************************
* external functions (explicit model)
************************************************/

	// forward explicit VDE

	external_function_casadi exfun_forw_vde;
	exfun_forw_vde.casadi_fun = &vde_energy_balanced_model;
	exfun_forw_vde.casadi_work = &vde_energy_balanced_model_work;
	exfun_forw_vde.casadi_sparsity_in = &vde_energy_balanced_model_sparsity_in;
	exfun_forw_vde.casadi_sparsity_out = &vde_energy_balanced_model_sparsity_out;
	exfun_forw_vde.casadi_n_in = &vde_energy_balanced_model_n_in;
	exfun_forw_vde.casadi_n_out = &vde_energy_balanced_model_n_out;

	int forw_vde_size = external_function_casadi_calculate_size(&exfun_forw_vde);
	void *forw_vde_mem = malloc(forw_vde_size);
	external_function_casadi_assign(&exfun_forw_vde, forw_vde_mem);

	// adjoint explicit VDE

	external_function_casadi exfun_adj_vde;
	exfun_adj_vde.casadi_fun = &vde_adj_energy_balanced_model;
	exfun_adj_vde.casadi_work = &vde_adj_energy_balanced_model_work;
	exfun_adj_vde.casadi_sparsity_in = &vde_adj_energy_balanced_model_sparsity_in;
	exfun_adj_vde.casadi_sparsity_out = &vde_adj_energy_balanced_model_sparsity_out;
	exfun_adj_vde.casadi_n_in = &vde_adj_energy_balanced_model_n_in;
	exfun_adj_vde.casadi_n_out = &vde_adj_energy_balanced_model_n_out;

	int adj_vde_size = external_function_casadi_calculate_size(&exfun_adj_vde);
	void *adj_vde_mem = malloc(adj_vde_size);
	external_function_casadi_assign(&exfun_adj_vde, adj_vde_mem);

	// jacobian explicit ODE

	external_function_casadi exfun_jac;
	exfun_jac.casadi_fun = &jac_energy_balanced_model;
	exfun_jac.casadi_work = &jac_energy_balanced_model_work;
	exfun_jac.casadi_sparsity_in = &jac_energy_balanced_model_sparsity_in;
	exfun_jac.casadi_sparsity_out = &jac_energy_balanced_model_sparsity_out;
	exfun_jac.casadi_n_in = &jac_energy_balanced_model_n_in;
	exfun_jac.casadi_n_out = &jac_energy_balanced_model_n_out;

	int jac_size = external_function_casadi_calculate_size(&exfun_jac);
	void *jac_mem = malloc(jac_size);
	external_function_casadi_assign(&exfun_jac, jac_mem);

	// hessian explicit ODE

	external_function_casadi exfun_hess_ode;
	exfun_hess_ode.casadi_fun = &vde_hess_energy_balanced_model;
	exfun_hess_ode.casadi_work = &vde_hess_energy_balanced_model_work;
	exfun_hess_ode.casadi_sparsity_in = &vde_hess_energy_balanced_model_sparsity_in;
	exfun_hess_ode.casadi_sparsity_out = &vde_hess_energy_balanced_model_sparsity_out;
	exfun_hess_ode.casadi_n_in = &vde_hess_energy_balanced_model_n_in;
	exfun_hess_ode.casadi_n_out = &vde_hess_energy_balanced_model_n_out;

	int hess_ode_size = external_function_casadi_calculate_size(&exfun_hess_ode);
	void *hess_ode_mem = malloc(hess_ode_size);
	external_function_casadi_assign(&exfun_hess_ode, hess_ode_mem);

/************************************************
* external functions (implicit model)
************************************************/

	// implicit ODE

	external_function_casadi exfun_ode;
	exfun_ode.casadi_fun = &impl_odeFun_energy_balanced_model;
	exfun_ode.casadi_work = &impl_odeFun_energy_balanced_model_work;
	exfun_ode.casadi_sparsity_in = &impl_odeFun_energy_balanced_model_sparsity_in;
	exfun_ode.casadi_sparsity_out = &impl_odeFun_energy_balanced_model_sparsity_out;
	exfun_ode.casadi_n_in = &impl_odeFun_energy_balanced_model_n_in;
	exfun_ode.casadi_n_out = &impl_odeFun_energy_balanced_model_n_out;

	int ode_size = external_function_casadi_calculate_size(&exfun_ode);
	void *ode_mem = malloc(ode_size);
	external_function_casadi_assign(&exfun_ode, ode_mem);

	// jac_x implicit ODE

	external_function_casadi exfun_jac_x_ode;
	exfun_jac_x_ode.casadi_fun = &impl_jacFun_x_energy_balanced_model;
	exfun_jac_x_ode.casadi_work = &impl_jacFun_x_energy_balanced_model_work;
	exfun_jac_x_ode.casadi_sparsity_in = &impl_jacFun_x_energy_balanced_model_sparsity_in;
	exfun_jac_x_ode.casadi_sparsity_out = &impl_jacFun_x_energy_balanced_model_sparsity_out;
	exfun_jac_x_ode.casadi_n_in = &impl_jacFun_x_energy_balanced_model_n_in;
	exfun_jac_x_ode.casadi_n_out = &impl_jacFun_x_energy_balanced_model_n_out;

	int jac_x_ode_size = external_function_casadi_calculate_size(&exfun_jac_x_ode);
	void *jac_x_ode_mem = malloc(jac_x_ode_size);
	external_function_casadi_assign(&exfun_jac_x_ode, jac_x_ode_mem);

	// jac_xdot implicit ODE

	external_function_casadi exfun_jac_xdot_ode;
	exfun_jac_xdot_ode.casadi_fun = &impl_jacFun_xdot_energy_balanced_model;
	exfun_jac_xdot_ode.casadi_work = &impl_jacFun_xdot_energy_balanced_model_work;
	exfun_jac_xdot_ode.casadi_sparsity_in = &impl_jacFun_xdot_energy_balanced_model_sparsity_in;
	exfun_jac_xdot_ode.casadi_sparsity_out = &impl_jacFun_xdot_energy_balanced_model_sparsity_out;
	exfun_jac_xdot_ode.casadi_n_in = &impl_jacFun_xdot_energy_balanced_model_n_in;
	exfun_jac_xdot_ode.casadi_n_out = &impl_jacFun_xdot_energy_balanced_model_n_out;

	int jac_xdot_ode_size = external_function_casadi_calculate_size(&exfun_jac_xdot_ode);
	void *jac_xdot_ode_mem = malloc(jac_xdot_ode_size);
	external_function_casadi_assign(&exfun_jac_xdot_ode, jac_xdot_ode_mem);

	// jac_u implicit ODE

	external_function_casadi exfun_jac_u_ode;
	exfun_jac_u_ode.casadi_fun = &impl_jacFun_u_energy_balanced_model;
	exfun_jac_u_ode.casadi_work = &impl_jacFun_u_energy_balanced_model_work;
	exfun_jac_u_ode.casadi_sparsity_in = &impl_jacFun_u_energy_balanced_model_sparsity_in;
	exfun_jac_u_ode.casadi_sparsity_out = &impl_jacFun_u_energy_balanced_model_sparsity_out;
	exfun_jac_u_ode.casadi_n_in = &impl_jacFun_u_energy_balanced_model_n_in;
	exfun_jac_u_ode.casadi_n_out = &impl_jacFun_u_energy_balanced_model_n_out;

	int jac_u_ode_size = external_function_casadi_calculate_size(&exfun_jac_u_ode);
	void *jac_u_ode_mem = malloc(jac_u_ode_size);
	external_function_casadi_assign(&exfun_jac_u_ode, jac_u_ode_mem);



	int number_sim_solvers = 3;
	int nss;
	for (nss=0; nss<number_sim_solvers; nss++)
// for (nss=1; nss<2; nss++)
	{

/************************************************
* sim config
************************************************/

		int config_size = sim_config_calculate_size();
		void *config_mem = malloc(config_size);
		sim_config *config = sim_config_assign(config_mem);

		switch (nss)
		{

			case 0: // erk
				printf("\n\nsim solver: ERK\n");
				sim_erk_config_initialize_default(config);
				config->ns = 4; // number of integration stages
				break;

			case 1: // irk
				printf("\n\nsim solver: IRK\n");
				sim_irk_config_initialize_default(config);
				config->ns = 2; // number of integration stages
				break;

			case 2: // old_lifted_irk
				printf("\n\nsim solver: Old_Lifted_IRK\n");
				sim_lifted_irk_config_initialize_default(config);
				config->ns = 2; // number of integration stages
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
		sim_opts *opts = config->opts_assign(config, dims, opts_mem);
		config->opts_initialize_default(config, dims, opts);

		opts->num_steps = 10; // integration steps
		opts->sens_adj = true;

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

		in->T = T; // simulation time

		// external functions
		switch (nss)
		{
			case 0: // erk
			{
				erk_model *model = in->model;
				model->forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
				model->adj_vde_expl = (external_function_generic *) &exfun_adj_vde;
				model->hess_ode_expl = (external_function_generic *) &exfun_hess_ode;
				break;
			}
			case 1: // irk
			{
				irk_model *model = in->model;
				model->ode_impl = (external_function_generic *) &exfun_ode;
				model->jac_x_ode_impl = (external_function_generic *) &exfun_jac_x_ode;
				model->jac_xdot_ode_impl = (external_function_generic *) &exfun_jac_xdot_ode;
				model->jac_u_ode_impl = (external_function_generic *) &exfun_jac_u_ode;
				break;
			}
			case 2: // old_lifted_irk
			{
				lifted_irk_model *model = in->model;
				model->forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
				model->jac_ode_expl = (external_function_generic *) &exfun_jac;
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
			in->x[ii] = 1.0;//xref[ii];
		}

		// p
		for (ii = 0;ii < nu; ii++){
			in->u[ii] = 1.0;//1.0;
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

// 	acados_tic(&timer);

		for (ii=0;ii<NREP;ii++)
			config->evaluate(config, in, out, opts, mem, work);



// 	Time1 = acados_toc(&timer)/NREP;

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
// 		blasfeo_print_dmat(nx, nx+nu, &sA, 0, 0);
// 		blasfeo_print_tran_dvec(nx, &sx, 0);
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
	free(ode_mem);
	free(jac_x_ode_mem);
	free(jac_xdot_ode_mem);
	free(jac_u_ode_mem);

/************************************************
* return
************************************************/

	printf("\nsuccess!\n\n");

    return 0;
}
