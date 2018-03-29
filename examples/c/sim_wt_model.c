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
#define M_PI 3.14159265358979323846
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// acados
// TODO(dimitris): remove most includes
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// wt model
#include "examples/c/wt_model_nx3/wt_model.h"

// blasfeo
#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>

// x0 and u for simulation
#include "examples/c/wt_model_nx3/u_x0.c"



int main()
{

/************************************************
* bla bla bla
************************************************/

    int ii;
    int jj;

    int nx = 3;
    int nu = 4;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05; // simulation time

	double x_sim[nx*(nsim+1)];
	for (ii=0; ii<nx; ii++)
		x_sim[ii] = x0[ii];

/************************************************
* external functions (explicit model)
************************************************/

	// explicit ODE

	external_function_casadi expl_ode;
	expl_ode.casadi_fun = &ode_energy_balanced_model;
	expl_ode.casadi_work = &ode_energy_balanced_model_work;
	expl_ode.casadi_sparsity_in = &ode_energy_balanced_model_sparsity_in;
	expl_ode.casadi_sparsity_out = &ode_energy_balanced_model_sparsity_out;
	expl_ode.casadi_n_in = &ode_energy_balanced_model_n_in;
	expl_ode.casadi_n_out = &ode_energy_balanced_model_n_out;
	external_function_casadi_create(&expl_ode);

	// forward explicit VDE

	external_function_casadi expl_forw_vde;
	expl_forw_vde.casadi_fun = &vde_energy_balanced_model;
	expl_forw_vde.casadi_work = &vde_energy_balanced_model_work;
	expl_forw_vde.casadi_sparsity_in = &vde_energy_balanced_model_sparsity_in;
	expl_forw_vde.casadi_sparsity_out = &vde_energy_balanced_model_sparsity_out;
	expl_forw_vde.casadi_n_in = &vde_energy_balanced_model_n_in;
	expl_forw_vde.casadi_n_out = &vde_energy_balanced_model_n_out;
	external_function_casadi_create(&expl_forw_vde);

	// adjoint explicit VDE

	external_function_casadi expl_adj_vde;
	expl_adj_vde.casadi_fun = &vde_adj_energy_balanced_model;
	expl_adj_vde.casadi_work = &vde_adj_energy_balanced_model_work;
	expl_adj_vde.casadi_sparsity_in = &vde_adj_energy_balanced_model_sparsity_in;
	expl_adj_vde.casadi_sparsity_out = &vde_adj_energy_balanced_model_sparsity_out;
	expl_adj_vde.casadi_n_in = &vde_adj_energy_balanced_model_n_in;
	expl_adj_vde.casadi_n_out = &vde_adj_energy_balanced_model_n_out;
	external_function_casadi_create(&expl_adj_vde);

	// jacobian explicit ODE

	external_function_casadi expl_jac;
	expl_jac.casadi_fun = &jac_energy_balanced_model;
	expl_jac.casadi_work = &jac_energy_balanced_model_work;
	expl_jac.casadi_sparsity_in = &jac_energy_balanced_model_sparsity_in;
	expl_jac.casadi_sparsity_out = &jac_energy_balanced_model_sparsity_out;
	expl_jac.casadi_n_in = &jac_energy_balanced_model_n_in;
	expl_jac.casadi_n_out = &jac_energy_balanced_model_n_out;
	external_function_casadi_create(&expl_jac);

	// hessian explicit ODE

	external_function_casadi expl_hess_ode;
	expl_hess_ode.casadi_fun = &vde_hess_energy_balanced_model;
	expl_hess_ode.casadi_work = &vde_hess_energy_balanced_model_work;
	expl_hess_ode.casadi_sparsity_in = &vde_hess_energy_balanced_model_sparsity_in;
	expl_hess_ode.casadi_sparsity_out = &vde_hess_energy_balanced_model_sparsity_out;
	expl_hess_ode.casadi_n_in = &vde_hess_energy_balanced_model_n_in;
	expl_hess_ode.casadi_n_out = &vde_hess_energy_balanced_model_n_out;
	external_function_casadi_create(&expl_hess_ode);

/************************************************
* external functions (implicit model)
************************************************/

	// implicit ODE

	external_function_casadi impl_ode;
	impl_ode.casadi_fun = &impl_odeFun_energy_balanced_model;
	impl_ode.casadi_work = &impl_odeFun_energy_balanced_model_work;
	impl_ode.casadi_sparsity_in = &impl_odeFun_energy_balanced_model_sparsity_in;
	impl_ode.casadi_sparsity_out = &impl_odeFun_energy_balanced_model_sparsity_out;
	impl_ode.casadi_n_in = &impl_odeFun_energy_balanced_model_n_in;
	impl_ode.casadi_n_out = &impl_odeFun_energy_balanced_model_n_out;
	external_function_casadi_create(&impl_ode);

	// jac_x implicit ODE

	external_function_casadi impl_jac_x_ode;
	impl_jac_x_ode.casadi_fun = &impl_jacFun_x_energy_balanced_model;
	impl_jac_x_ode.casadi_work = &impl_jacFun_x_energy_balanced_model_work;
	impl_jac_x_ode.casadi_sparsity_in = &impl_jacFun_x_energy_balanced_model_sparsity_in;
	impl_jac_x_ode.casadi_sparsity_out = &impl_jacFun_x_energy_balanced_model_sparsity_out;
	impl_jac_x_ode.casadi_n_in = &impl_jacFun_x_energy_balanced_model_n_in;
	impl_jac_x_ode.casadi_n_out = &impl_jacFun_x_energy_balanced_model_n_out;
	external_function_casadi_create(&impl_jac_x_ode);

	// jac_xdot implicit ODE

	external_function_casadi impl_jac_xdot_ode;
	impl_jac_xdot_ode.casadi_fun = &impl_jacFun_xdot_energy_balanced_model;
	impl_jac_xdot_ode.casadi_work = &impl_jacFun_xdot_energy_balanced_model_work;
	impl_jac_xdot_ode.casadi_sparsity_in = &impl_jacFun_xdot_energy_balanced_model_sparsity_in;
	impl_jac_xdot_ode.casadi_sparsity_out = &impl_jacFun_xdot_energy_balanced_model_sparsity_out;
	impl_jac_xdot_ode.casadi_n_in = &impl_jacFun_xdot_energy_balanced_model_n_in;
	impl_jac_xdot_ode.casadi_n_out = &impl_jacFun_xdot_energy_balanced_model_n_out;
	external_function_casadi_create(&impl_jac_xdot_ode);

	// jac_u implicit ODE

	external_function_casadi impl_jac_u_ode;
	impl_jac_u_ode.casadi_fun = &impl_jacFun_u_energy_balanced_model;
	impl_jac_u_ode.casadi_work = &impl_jacFun_u_energy_balanced_model_work;
	impl_jac_u_ode.casadi_sparsity_in = &impl_jacFun_u_energy_balanced_model_sparsity_in;
	impl_jac_u_ode.casadi_sparsity_out = &impl_jacFun_u_energy_balanced_model_sparsity_out;
	impl_jac_u_ode.casadi_n_in = &impl_jacFun_u_energy_balanced_model_n_in;
	impl_jac_u_ode.casadi_n_out = &impl_jacFun_u_energy_balanced_model_n_out;
	external_function_casadi_create(&impl_jac_u_ode);

/************************************************
* gathered functions
************************************************/	

	// implicit ODE including jacobian w.r.t. x, xdot

	external_function_casadi impl_ode_inc_J_xxdot;
	impl_ode_inc_J_xxdot.casadi_fun 		 = &impl_ode_inc_J_xxdot_energy_balanced_model;
	impl_ode_inc_J_xxdot.casadi_work 		 = &impl_ode_inc_J_xxdot_energy_balanced_model_work;
	impl_ode_inc_J_xxdot.casadi_sparsity_in  = &impl_ode_inc_J_xxdot_energy_balanced_model_sparsity_in;
	impl_ode_inc_J_xxdot.casadi_sparsity_out = &impl_ode_inc_J_xxdot_energy_balanced_model_sparsity_out;
	impl_ode_inc_J_xxdot.casadi_n_in 		 = &impl_ode_inc_J_xxdot_energy_balanced_model_n_in;
	impl_ode_inc_J_xxdot.casadi_n_out 		 = &impl_ode_inc_J_xxdot_energy_balanced_model_n_out;

	external_function_casadi_create(&impl_ode_inc_J_xxdot);

	// jacobian of implicit ODE w.r.t. x, u

	external_function_casadi impl_ode_J_xu;
	impl_ode_J_xu.casadi_fun 		  = &impl_ode_J_xu_energy_balanced_model;
	impl_ode_J_xu.casadi_work 		  = &impl_ode_J_xu_energy_balanced_model_work;
	impl_ode_J_xu.casadi_sparsity_in  = &impl_ode_J_xu_energy_balanced_model_sparsity_in;
	impl_ode_J_xu.casadi_sparsity_out = &impl_ode_J_xu_energy_balanced_model_sparsity_out;
	impl_ode_J_xu.casadi_n_in 		  = &impl_ode_J_xu_energy_balanced_model_n_in;
	impl_ode_J_xu.casadi_n_out 		  = &impl_ode_J_xu_energy_balanced_model_n_out;

	external_function_casadi_create(&impl_ode_J_xu);


	// jacobian of implicit ODE w.r.t. x, xdot, u

	external_function_casadi impl_ode_J_xxdotu;
	impl_ode_J_xxdotu.casadi_fun 		  = &impl_ode_J_xxdotu_energy_balanced_model;
	impl_ode_J_xxdotu.casadi_work 		  = &impl_ode_J_xxdotu_energy_balanced_model_work;
	impl_ode_J_xxdotu.casadi_sparsity_in  = &impl_ode_J_xxdotu_energy_balanced_model_sparsity_in;
	impl_ode_J_xxdotu.casadi_sparsity_out = &impl_ode_J_xxdotu_energy_balanced_model_sparsity_out;
	impl_ode_J_xxdotu.casadi_n_in 		  = &impl_ode_J_xxdotu_energy_balanced_model_n_in;
	impl_ode_J_xxdotu.casadi_n_out 		  = &impl_ode_J_xxdotu_energy_balanced_model_n_out;

	external_function_casadi_create(&impl_ode_J_xxdotu);


	int number_sim_solvers = 2;
	int nss;
	for (nss = 0; nss < number_sim_solvers; nss++)
	{
		/************************************************
		* sim plan & config
		************************************************/
		// printf("using solver no. %d\n",nss);
		// choose plan
		sim_solver_plan plan;

		switch (nss)
		{

			case 0:
				printf("\n\nsim solver: ERK\n");
				plan.sim_solver = ERK;
				break;

			case 1:
				printf("\n\nsim solver: IRK\n");
				plan.sim_solver = IRK;
				break;

			case 2:
				printf("\n\nsim solver: Lifted_IRK\n");
				plan.sim_solver = LIFTED_IRK;
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);

		}

		// create correct config based on plan
		sim_solver_config *config = sim_config_create(plan);

		/************************************************
		* sim dims
		************************************************/

		sim_dims *dims = sim_dims_create();

		dims->nx = nx;
		dims->nu = nu;

		/************************************************
		* sim opts
		************************************************/

		sim_rk_opts *opts = sim_opts_create(config, dims);

//		opts->ns = 4; // number of stages in rk integrator
//		opts->num_steps = 5; // number of integration steps
		opts->sens_adj = true;
		opts->sens_forw = true;

		switch (nss)
		{

			case 0:
				// ERK
				opts->ns = 4; // number of stages in rk integrator
				break;

			case 1:
				// IRK
				opts->ns = 2; // number of stages in rk integrator
				break;

			case 2:
				// lifted IRK
				opts->ns = 2; // number of stages in rk integrator
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);

		}


		/************************************************
		* sim in / out
		************************************************/

		sim_in *in = sim_in_create(config, dims);
		sim_out *out = sim_out_create(config, dims);

		in->T = T;

		// external functions
		switch (nss)
		{
			case 0:
			{
				erk_model *model = in->model;
				model->ode_expl = (external_function_generic *) &expl_ode;
				sim_set_model(config, in, "forward_vde", &expl_forw_vde);
				sim_set_model(config, in, "adjoint_vde", &expl_adj_vde);
				// model->hess_ode_expl = (external_function_generic *) &expl_hess_ode;
				break;
			}
			case 1:
			{
				irk_model *model = in->model;
				model->ode_impl = (external_function_generic *) &impl_ode;
				model->jac_x_ode_impl = (external_function_generic *) &impl_jac_x_ode;
				model->jac_xdot_ode_impl = (external_function_generic *) &impl_jac_xdot_ode;
				model->jac_u_ode_impl = (external_function_generic *) &impl_jac_u_ode;

				model->impl_ode_inc_J_xxdot = (external_function_generic *) &impl_ode_inc_J_xxdot;
				model->impl_ode_J_xu = (external_function_generic *) &impl_ode_J_xu;
				model->impl_ode_J_xxdotu = (external_function_generic *) &impl_ode_J_xxdotu;
				break;
			}
			case 2:
			{
				lifted_irk_model *model = in->model;
				model->forw_vde_expl = (external_function_generic *) &expl_forw_vde;
				model->jac_ode_expl = (external_function_generic *) &expl_jac;
				break;
			}
			default :
			{
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);
			}
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
		* sim solver
		************************************************/

		sim_solver *sim_solver = sim_create(config, dims, opts);

		int acados_return;

    	acados_timer timer;
		acados_tic(&timer);

		int nsim0 = nsim;

		double cpu_time = 0.0;
		double la_time = 0.0;
		double ad_time = 0.0;

//		for (ii=0; ii<nsim; ii++)
		for (ii=0; ii<nsim0; ii++)
		{
			// x
			for (jj = 0; jj < nx; jj++)
				in->x[jj] = x_sim[ii*nx+jj];

			// p
			for (jj = 0; jj < 2; jj++)
				in->u[jj] = u_sim[ii*2+jj];
			for (jj = 0; jj < nu; jj++)
				in->u[2+jj] = 0.1;

//			d_print_mat(1, nx, in->x, 1);
//			d_print_mat(1, nu, in->u, 1);

		    acados_return = sim_solve(sim_solver, in, out);
			if (acados_return != 0)
            	printf("error in sim solver\n");

			cpu_time += out->info->CPUtime;
			la_time += out->info->LAtime;
			ad_time += out->info->ADtime;

//			d_print_mat(1, nx, out->xn, 1);

			// x_out
			for (jj = 0; jj < nx; jj++)
				x_sim[(ii+1)*nx+jj] = out->xn[jj];

		}
		double total_cpu_time = acados_toc(&timer);

		/************************************************
		* printing
		************************************************/

		printf("\nxn: \n");
		for (ii=0; ii<nx; ii++)
			printf("%8.5f ", x_sim[nsim0*nx+ii]);
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
		if(opts->sens_adj)
		{
			S_adj_out = out->S_adj;
			printf("\nS_adj_out: \n");
			for (ii=0;ii<nx+nu;ii++){
				printf("%8.5f ", S_adj_out[ii]);
			}
			printf("\n");
		}

#if 0

		double *S_hess_out;
		if(opts->sens_hess)
		{
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

		if(opts->sens_adj)
		{
			struct blasfeo_dmat sA;
			blasfeo_allocate_dmat(nx, nx+nu, &sA);
			blasfeo_pack_dmat(nx, nx+nu, S_forw_out, nx, &sA, 0, 0);

			struct blasfeo_dvec sx;
			blasfeo_allocate_dvec(nx, &sx);
			blasfeo_pack_dvec(nx, in->S_adj, &sx, 0);

			struct blasfeo_dvec sz;
			blasfeo_allocate_dvec(nx+nu, &sz);
			// blasfeo_print_dmat(nx, nx+nu, &sA, 0, 0);
			// blasfeo_print_tran_dvec(nx, &sx, 0);
			blasfeo_dgemv_t(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);

			printf("\nJac times lambdaX:\n");
			blasfeo_print_tran_dvec(nx+nu, &sz, 0);

			blasfeo_free_dmat(&sA);
			blasfeo_free_dvec(&sx);
			blasfeo_free_dvec(&sz);
		}
#endif

//		printf("time split: %f ms CPU, %f ms LA, %f ms AD\n\n", cpu_time, la_time, ad_time);
		printf("time for %d simulation steps: %f ms (AD time: %f ms (%5.2f%%))\n\n", nsim, 1e3*total_cpu_time, 1e3*ad_time, 1e2*ad_time/cpu_time);

		free(sim_solver);
		free(in);
		free(out);

		free(opts);
		free(config);
	}

	// TODO(dimitris): free all external functions (or write a free_model)
	// explicit model
	external_function_casadi_free(&expl_forw_vde);
	external_function_casadi_free(&expl_adj_vde);
	external_function_casadi_free(&expl_jac);
	external_function_casadi_free(&expl_hess_ode);
	// implicit model
	external_function_casadi_free(&impl_ode);
	external_function_casadi_free(&impl_jac_x_ode);
	external_function_casadi_free(&impl_jac_xdot_ode);
	external_function_casadi_free(&impl_jac_u_ode);
	// gathered functions
	external_function_casadi_free(&impl_ode_inc_J_xxdot);
	external_function_casadi_free(&impl_ode_J_xu);
	external_function_casadi_free(&impl_ode_J_xxdotu);
	

	/************************************************
	* return
	************************************************/

	printf("\nsuccess! (RESULT NOT CHECKED) \n\n");

    return 0;
}
