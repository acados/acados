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
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

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

	external_function_casadi exfun_forw_vde;
	exfun_forw_vde.casadi_fun = &vdeFun;
	exfun_forw_vde.casadi_work = &vdeFun_work;
	exfun_forw_vde.casadi_sparsity_in = &vdeFun_sparsity_in;
	exfun_forw_vde.casadi_sparsity_out = &vdeFun_sparsity_out;
	exfun_forw_vde.casadi_n_in = &vdeFun_n_in;
	exfun_forw_vde.casadi_n_out = &vdeFun_n_out;
	external_function_casadi_create(&exfun_forw_vde);


	// adjoint explicit VDE

	external_function_casadi exfun_adj_vde;
	exfun_adj_vde.casadi_fun = &adjFun;
	exfun_adj_vde.casadi_work = &adjFun_work;
	exfun_adj_vde.casadi_sparsity_in = &adjFun_sparsity_in;
	exfun_adj_vde.casadi_sparsity_out = &adjFun_sparsity_out;
	exfun_adj_vde.casadi_n_in = &adjFun_n_in;
	exfun_adj_vde.casadi_n_out = &adjFun_n_out;
	external_function_casadi_create(&exfun_adj_vde);


	// jacobian explicit ODE

	external_function_casadi exfun_jac;
	exfun_jac.casadi_fun = &jacFun;
	exfun_jac.casadi_work = &jacFun_work;
	exfun_jac.casadi_sparsity_in = &jacFun_sparsity_in;
	exfun_jac.casadi_sparsity_out = &jacFun_sparsity_out;
	exfun_jac.casadi_n_in = &jacFun_n_in;
	exfun_jac.casadi_n_out = &jacFun_n_out;
	external_function_casadi_create(&exfun_jac);


	// hessian explicit ODE

	external_function_casadi exfun_hess_ode;
	exfun_hess_ode.casadi_fun = &hessFun;
	exfun_hess_ode.casadi_work = &hessFun_work;
	exfun_hess_ode.casadi_sparsity_in = &hessFun_sparsity_in;
	exfun_hess_ode.casadi_sparsity_out = &hessFun_sparsity_out;
	exfun_hess_ode.casadi_n_in = &hessFun_n_in;
	exfun_hess_ode.casadi_n_out = &hessFun_n_out;
	external_function_casadi_create(&exfun_hess_ode);


	/************************************************
	* external functions (implicit model)
	************************************************/

	// implicit ODE

	external_function_casadi exfun_ode;
	exfun_ode.casadi_fun = &impl_odeFun;
	exfun_ode.casadi_work = &impl_odeFun_work;
	exfun_ode.casadi_sparsity_in = &impl_odeFun_sparsity_in;
	exfun_ode.casadi_sparsity_out = &impl_odeFun_sparsity_out;
	exfun_ode.casadi_n_in = &impl_odeFun_n_in;
	exfun_ode.casadi_n_out = &impl_odeFun_n_out;
	external_function_casadi_create(&exfun_ode);

	// jac_x implicit ODE

	external_function_casadi exfun_jac_x_ode;
	exfun_jac_x_ode.casadi_fun = &impl_jacFun_x;
	exfun_jac_x_ode.casadi_work = &impl_jacFun_x_work;
	exfun_jac_x_ode.casadi_sparsity_in = &impl_jacFun_x_sparsity_in;
	exfun_jac_x_ode.casadi_sparsity_out = &impl_jacFun_x_sparsity_out;
	exfun_jac_x_ode.casadi_n_in = &impl_jacFun_x_n_in;
	exfun_jac_x_ode.casadi_n_out = &impl_jacFun_x_n_out;
	external_function_casadi_create(&exfun_jac_x_ode);

	// jac_xdot implicit ODE

	external_function_casadi exfun_jac_xdot_ode;
	exfun_jac_xdot_ode.casadi_fun = &impl_jacFun_xdot;
	exfun_jac_xdot_ode.casadi_work = &impl_jacFun_xdot_work;
	exfun_jac_xdot_ode.casadi_sparsity_in = &impl_jacFun_xdot_sparsity_in;
	exfun_jac_xdot_ode.casadi_sparsity_out = &impl_jacFun_xdot_sparsity_out;
	exfun_jac_xdot_ode.casadi_n_in = &impl_jacFun_xdot_n_in;
	exfun_jac_xdot_ode.casadi_n_out = &impl_jacFun_xdot_n_out;
	external_function_casadi_create(&exfun_jac_xdot_ode);

	// jac_u implicit ODE

	external_function_casadi exfun_jac_u_ode;
	exfun_jac_u_ode.casadi_fun = &impl_jacFun_u;
	exfun_jac_u_ode.casadi_work = &impl_jacFun_u_work;
	exfun_jac_u_ode.casadi_sparsity_in = &impl_jacFun_u_sparsity_in;
	exfun_jac_u_ode.casadi_sparsity_out = &impl_jacFun_u_sparsity_out;
	exfun_jac_u_ode.casadi_n_in = &impl_jacFun_u_n_in;
	exfun_jac_u_ode.casadi_n_out = &impl_jacFun_u_n_out;
	external_function_casadi_create(&exfun_jac_u_ode);


	int number_sim_solvers = 3;
	int nss;
	for (nss = 0; nss < number_sim_solvers; nss++)
	{
		/************************************************
		* sim config
		************************************************/

		sim_solver_plan plan;
		sim_solver_config *config;

		switch (nss)
		{

			case 0:
				printf("\n\nsim solver: ERK\n");
				plan.sim_solver = ERK;
				config = sim_config_create(&plan);
				// TODO(dimitris): move back to opts and then only create config outside switch once
				config->ns = 4;
				break;

			case 1:
				printf("\n\nsim solver: IRK\n");
				plan.sim_solver = IRK;
				config = sim_config_create(&plan);
				config->ns = 2;
				break;

			case 2:
				printf("\n\nsim solver: Lifted_IRK\n");
				plan.sim_solver = LIFTED_IRK;
				config = sim_config_create(&plan);
				config->ns = 2;
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);

		}

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

		opts->sens_adj = true;

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
				sim_set_model(config, in, "forward_vde", &exfun_forw_vde);
				sim_set_model(config, in, "adjoint_vde", &exfun_adj_vde);
				// model->hess_ode_expl = (external_function_generic *) &exfun_hess_ode;
				break;
			}
			case 1:
			{
				irk_model *model = in->model;
				model->ode_impl = (external_function_generic *) &exfun_ode;
				model->jac_x_ode_impl = (external_function_generic *) &exfun_jac_x_ode;
				model->jac_xdot_ode_impl = (external_function_generic *) &exfun_jac_xdot_ode;
				model->jac_u_ode_impl = (external_function_generic *) &exfun_jac_u_ode;
				break;
			}
			case 2:
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
		for (ii = 0; ii < nx; ii++)
			in->x[ii] = xref[ii];

		// p
		for (ii = 0;ii < nu; ii++)
			in->u[ii] = 1.0;

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

    	// acados_timer timer;
		// acados_tic(&timer);

		for (ii=0;ii<NREP;ii++)
		    acados_return = sim_solve(sim_solver, in, out);

		// double cpu_time = acados_toc(&timer)/NREP;

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
		if(opts->sens_adj)
		{
			S_adj_out = out->S_adj;
			printf("\nS_adj_out: \n");
			for (ii=0;ii<nx+nu;ii++){
				printf("%8.5f ", S_adj_out[ii]);
			}
			printf("\n");
		}

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
		free(sim_solver);
		free(in);
		free(out);

		free(opts);
		free(config);
	}

	// TODO(dimitris): free all external functions (or write a free_model)
	external_function_casadi_free(&exfun_adj_vde);

	/************************************************
	* return
	************************************************/

	printf("\nsuccess! (RESULT NOT CHECKED) \n\n");

    return 0;
}
