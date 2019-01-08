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
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// wt model
#include "examples/c/wt_model_nx3/wt_model.h"

// x0 and u for simulation
#include "examples/c/wt_model_nx3/u_x0.c"

// blasfeo
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"


int main()
{

/************************************************
* initialize common stuff (for all integrators)
************************************************/

    int nx = 3;
    int nu = 4;
	int nz = 0;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05; // simulation time

	double x_sim[nx*(nsim+1)];
	for (int ii = 0; ii < nx; ii++)
		x_sim[ii] = x0[ii];

	/************************************************
	* external functions (explicit model)
	************************************************/

	// expl_ode_fun
	external_function_casadi expl_ode_fun;
	expl_ode_fun.casadi_fun = &casadi_expl_ode_fun;
	expl_ode_fun.casadi_work = &casadi_expl_ode_fun_work;
	expl_ode_fun.casadi_sparsity_in = &casadi_expl_ode_fun_sparsity_in;
	expl_ode_fun.casadi_sparsity_out = &casadi_expl_ode_fun_sparsity_out;
	expl_ode_fun.casadi_n_in = &casadi_expl_ode_fun_n_in;
	expl_ode_fun.casadi_n_out = &casadi_expl_ode_fun_n_out;
	external_function_casadi_create(&expl_ode_fun);

	// expl_vde_for
	external_function_casadi expl_vde_for;
	expl_vde_for.casadi_fun = &casadi_expl_vde_for;
	expl_vde_for.casadi_work = &casadi_expl_vde_for_work;
	expl_vde_for.casadi_sparsity_in = &casadi_expl_vde_for_sparsity_in;
	expl_vde_for.casadi_sparsity_out = &casadi_expl_vde_for_sparsity_out;
	expl_vde_for.casadi_n_in = &casadi_expl_vde_for_n_in;
	expl_vde_for.casadi_n_out = &casadi_expl_vde_for_n_out;
	external_function_casadi_create(&expl_vde_for);

	// expl_vde_adj
	external_function_casadi expl_vde_adj;
	expl_vde_adj.casadi_fun = &casadi_expl_vde_adj;
	expl_vde_adj.casadi_work = &casadi_expl_vde_adj_work;
	expl_vde_adj.casadi_sparsity_in = &casadi_expl_vde_adj_sparsity_in;
	expl_vde_adj.casadi_sparsity_out = &casadi_expl_vde_adj_sparsity_out;
	expl_vde_adj.casadi_n_in = &casadi_expl_vde_adj_n_in;
	expl_vde_adj.casadi_n_out = &casadi_expl_vde_adj_n_out;
	external_function_casadi_create(&expl_vde_adj);

	/************************************************
	* external functions (implicit model)
	************************************************/

	// impl_ode_fun
	external_function_casadi impl_ode_fun;
	impl_ode_fun.casadi_fun = &casadi_impl_ode_fun;
	impl_ode_fun.casadi_work = &casadi_impl_ode_fun_work;
	impl_ode_fun.casadi_sparsity_in = &casadi_impl_ode_fun_sparsity_in;
	impl_ode_fun.casadi_sparsity_out = &casadi_impl_ode_fun_sparsity_out;
	impl_ode_fun.casadi_n_in = &casadi_impl_ode_fun_n_in;
	impl_ode_fun.casadi_n_out = &casadi_impl_ode_fun_n_out;
	external_function_casadi_create(&impl_ode_fun);

	// impl_ode_fun_jac_x_xdot
	external_function_casadi impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_fun = &casadi_impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_work = &casadi_impl_ode_fun_jac_x_xdot_work;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_sparsity_in;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_sparsity_out;
	impl_ode_fun_jac_x_xdot.casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_n_in;
	impl_ode_fun_jac_x_xdot.casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_n_out;
	external_function_casadi_create(&impl_ode_fun_jac_x_xdot);

	// impl_ode_jac_x_xdot_u
	external_function_casadi impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_fun = &casadi_impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_work = &casadi_impl_ode_jac_x_xdot_u_work;
	impl_ode_jac_x_xdot_u.casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_sparsity_in;
	impl_ode_jac_x_xdot_u.casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_sparsity_out;
	impl_ode_jac_x_xdot_u.casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_n_in;
	impl_ode_jac_x_xdot_u.casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_n_out;
	external_function_casadi_create(&impl_ode_jac_x_xdot_u);

	// impl_ode_jac_x_xdot_u
	external_function_casadi impl_ode_fun_jac_x_xdot_u;
	impl_ode_fun_jac_x_xdot_u.casadi_fun = &casadi_impl_ode_fun_jac_x_xdot_u;
	impl_ode_fun_jac_x_xdot_u.casadi_work = &casadi_impl_ode_fun_jac_x_xdot_u_work;
	impl_ode_fun_jac_x_xdot_u.casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_u_sparsity_in;
	impl_ode_fun_jac_x_xdot_u.casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_u_sparsity_out;
	impl_ode_fun_jac_x_xdot_u.casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_u_n_in;
	impl_ode_fun_jac_x_xdot_u.casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_u_n_out;
	external_function_casadi_create(&impl_ode_fun_jac_x_xdot_u);

	/************************************************
	* external functions (Generalized Nonlinear Static Feedback (GNSF) model)
	************************************************/
    // phi_fun
    external_function_casadi phi_fun;
    phi_fun.casadi_fun            = &casadi_phi_fun;
    phi_fun.casadi_work           = &casadi_phi_fun_work;
    phi_fun.casadi_sparsity_in    = &casadi_phi_fun_sparsity_in;
    phi_fun.casadi_sparsity_out   = &casadi_phi_fun_sparsity_out;
    phi_fun.casadi_n_in           = &casadi_phi_fun_n_in;
    phi_fun.casadi_n_out          = &casadi_phi_fun_n_out;
	external_function_casadi_create(&phi_fun);

    // phi_fun_jac_y
    external_function_casadi phi_fun_jac_y;
    phi_fun_jac_y.casadi_fun            = &casadi_phi_fun_jac_y;
    phi_fun_jac_y.casadi_work           = &casadi_phi_fun_jac_y_work;
    phi_fun_jac_y.casadi_sparsity_in    = &casadi_phi_fun_jac_y_sparsity_in;
    phi_fun_jac_y.casadi_sparsity_out   = &casadi_phi_fun_jac_y_sparsity_out;
    phi_fun_jac_y.casadi_n_in           = &casadi_phi_fun_jac_y_n_in;
    phi_fun_jac_y.casadi_n_out          = &casadi_phi_fun_jac_y_n_out;
	external_function_casadi_create(&phi_fun_jac_y);

    // phi_jac_y_uhat
    external_function_casadi phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_fun                = &casadi_phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_work               = &casadi_phi_jac_y_uhat_work;
    phi_jac_y_uhat.casadi_sparsity_in        = &casadi_phi_jac_y_uhat_sparsity_in;
    phi_jac_y_uhat.casadi_sparsity_out       = &casadi_phi_jac_y_uhat_sparsity_out;
    phi_jac_y_uhat.casadi_n_in               = &casadi_phi_jac_y_uhat_n_in;
    phi_jac_y_uhat.casadi_n_out              = &casadi_phi_jac_y_uhat_n_out;
	external_function_casadi_create(&phi_jac_y_uhat);

    // f_lo_fun_jac_x1k1uz
    external_function_casadi f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_fun            = &casadi_f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_work           = &casadi_f_lo_fun_jac_x1k1uz_work;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_in    = &casadi_f_lo_fun_jac_x1k1uz_sparsity_in;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_out   = &casadi_f_lo_fun_jac_x1k1uz_sparsity_out;
    f_lo_fun_jac_x1k1uz.casadi_n_in           = &casadi_f_lo_fun_jac_x1k1uz_n_in;
    f_lo_fun_jac_x1k1uz.casadi_n_out          = &casadi_f_lo_fun_jac_x1k1uz_n_out;
	external_function_casadi_create(&f_lo_fun_jac_x1k1uz);

    // get_matrices_fun
    external_function_casadi get_matrices_fun;
    get_matrices_fun.casadi_fun            = &casadi_get_matrices_fun;
    get_matrices_fun.casadi_work           = &casadi_get_matrices_fun_work;
    get_matrices_fun.casadi_sparsity_in    = &casadi_get_matrices_fun_sparsity_in;
    get_matrices_fun.casadi_sparsity_out   = &casadi_get_matrices_fun_sparsity_out;
    get_matrices_fun.casadi_n_in           = &casadi_get_matrices_fun_n_in;
    get_matrices_fun.casadi_n_out          = &casadi_get_matrices_fun_n_out;
	external_function_casadi_create(&get_matrices_fun);



	int number_sim_solvers = 4;
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
				printf("\n\nsim solver: LIFTED_IRK\n");
				plan.sim_solver = LIFTED_IRK;
				break;


			case 3:
				printf("\n\nsim solver: GNSF\n");
				plan.sim_solver = GNSF;
				break;
			

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);

		}

		// create correct config based on plan
		sim_config *config = sim_config_create(plan);

		/************************************************
		* sim dims
		************************************************/

		void *dims = sim_dims_create(config);
		
		sim_dims_set(config, dims, "nx", &nx);
		sim_dims_set(config, dims, "nu", &nu);
		sim_dims_set(config, dims, "nz", &nz);

		/************************************************
		* sim opts
		************************************************/

		sim_opts *opts = sim_opts_create(config, dims);

		// opts->ns = 4; // number of stages in rk integrator
		// opts->num_steps = 3; // number of integration steps

		opts->sens_forw = true;
		opts->sens_adj = true;

		switch (nss)
		{

			case 0:
				// ERK
				opts->ns = 4; // number of stages in rk integrator
				opts->sens_adj = false;

				break;

			case 1:
				// IRK
				opts->ns = 2; // number of stages in rk integrator
				break;

			case 2:
				// new lifted IRK
				opts->ns = 2; // number of stages in rk integrator
				opts->sens_adj = false; // not implemented yet
				break;

			case 3:
				// GNSF
				opts->ns = 2; // number of stages in rk integrator
				opts->jac_reuse = true; // jacobian reuse
				opts->newton_iter = 3; // number of newton iterations per integration step

				// set additional dimensions
				int nx1 = nx;
				int nz1 = 0;
				int nout = 1;
				int ny = nx;
				int nuhat = nu;

				sim_dims_set(config, dims, "nx1", &nx1);
				sim_dims_set(config, dims, "nz1", &nz1);
				sim_dims_set(config, dims, "nout", &nout);
				sim_dims_set(config, dims, "ny", &ny);
				sim_dims_set(config, dims, "nuhat", &nuhat);

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
				config->model_set(in->model, "expl_ode_fun", &expl_ode_fun);
				config->model_set(in->model, "expl_vde_for", &expl_vde_for);
				config->model_set(in->model, "expl_vde_adj", &expl_vde_adj);
				break;
			}
			case 1:
			{
				config->model_set(in->model, "impl_ode_fun", &impl_ode_fun);
				config->model_set(in->model, "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot);
				config->model_set(in->model, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
				break;
			}
			case 2: // lifted_irk
			{
				config->model_set(in->model, "impl_ode_fun", &impl_ode_fun);
				config->model_set(in->model, "impl_ode_fun_jac_x_xdot_u", &impl_ode_fun_jac_x_xdot_u);
				break;
			}
			case 3: // GNSF
			{
				// set model funtions
				config->model_set(in->model, "phi_fun", &phi_fun);
				config->model_set(in->model, "phi_fun_jac_y", &phi_fun_jac_y);
				config->model_set(in->model, "phi_jac_y_uhat", &phi_jac_y_uhat);
				config->model_set(in->model, "f_lo_jac_x1_x1dot_u_z", &f_lo_fun_jac_x1k1uz);
				config->model_set(in->model, "get_gnsf_matrices", &get_matrices_fun);
				break;
			}
			default :
			{
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);
			}
		}

		// seeds forw
		for (int ii = 0; ii < nx * NF; ii++)
			in->S_forw[ii] = 0.0;
		for (int ii = 0; ii < nx; ii++)
			in->S_forw[ii * (nx + 1)] = 1.0;

		// seeds adj
		for (int ii = 0; ii < nx; ii++)
			in->S_adj[ii] = 1.0;

		/************************************************
		* sim solver
		************************************************/

		sim_solver *sim_solver = sim_solver_create(config, dims, opts);

		int acados_return;

		sim_precompute(sim_solver, in, out);

    	acados_timer timer;
		acados_tic(&timer);

		int nsim0 = 1;//nsim;

		double cpu_time = 0.0;
		double la_time = 0.0;
		double ad_time = 0.0;

		printf("\n---> testing integrator %d (num_steps = %d, num_stages = %d, jac_reuse = %d, newton_iter = %d )\n",
					nss, opts->num_steps, opts->ns, opts->jac_reuse, opts->newton_iter);

		for (int ii = 0; ii < nsim0; ii++)
		{
			// x
			for (int jj = 0; jj < nx; jj++)
				in->x[jj] = x_sim[ii*nx+jj];

			// u
			for (int jj = 0; jj < nu; jj++)
				in->u[jj] = u_sim[ii*nu+jj];

		    acados_return = sim_solve(sim_solver, in, out);
			if (acados_return != 0)
            	printf("error in sim solver\n");

			cpu_time += out->info->CPUtime;
			la_time += out->info->LAtime;
			ad_time += out->info->ADtime;

			// x_out
			for (int jj = 0; jj < nx; jj++)
				x_sim[(ii+1)*nx+jj] = out->xn[jj];

		}
		double total_cpu_time = acados_toc(&timer);

		/************************************************
		* printing
		************************************************/

		printf("\nxn: \n");
		for (int ii = 0; ii < nx; ii++)
			printf("%8.5f ", x_sim[nsim0*nx+ii]);
		printf("\n");

		double *S_forw_out;
		S_forw_out = NULL;
		if(opts->sens_forw){
			S_forw_out = out->S_forw;
			printf("\nS_forw_out: \n");
			for (int ii = 0; ii < nx; ii++){
				for (int jj = 0; jj < NF; jj++)
					printf("%8.5f ", S_forw_out[jj*nx+ii]);
				printf("\n");
			}
		}

		double *S_adj_out;
		if(opts->sens_adj)
		{
			S_adj_out = out->S_adj;
			printf("\nS_adj_out: \n");
			for (int ii = 0; ii < nx + nu; ii++){
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

// 	printf("time split: %f ms CPU, %f ms LA, %f ms AD\n\n", cpu_time, la_time, ad_time);
		printf("time for %d simulation steps: %f ms (AD time: %f ms (%5.2f%%))\n\n", nsim0, 1e3*total_cpu_time, 1e3*ad_time, 1e2*ad_time/cpu_time);

		free(sim_solver);
		free(in);
		free(out);

		free(opts);
		free(config);
		free(dims);
	}

	// TODO(dimitris): free all external functions (or write a free_model)
	// explicit model
    external_function_casadi_free(&expl_ode_fun);
	external_function_casadi_free(&expl_vde_for);
	external_function_casadi_free(&expl_vde_adj);
	// implicit model
	external_function_casadi_free(&impl_ode_fun);
	external_function_casadi_free(&impl_ode_fun_jac_x_xdot);
	external_function_casadi_free(&impl_ode_jac_x_xdot_u);
	external_function_casadi_free(&impl_ode_fun_jac_x_xdot_u);
	// gnsf functions:
	external_function_casadi_free(&f_lo_fun_jac_x1k1uz);
	external_function_casadi_free(&phi_fun);
	external_function_casadi_free(&phi_fun_jac_y);
	external_function_casadi_free(&phi_jac_y_uhat);
	external_function_casadi_free(&get_matrices_fun);	

	/************************************************
	* return
	************************************************/

	printf("\nsuccess! (RESULT NOT CHECKED) \n\n");

    return 0;
}
