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

// blasfeo
#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>

// wt model
#include "examples/c/wt_model_nx6/wt_model.h"

// x0 and u for simulation
#include "examples/c/wt_model_nx6/u_x0.c"



int main()
{

	/************************************************
	* initialization
	************************************************/

    int ii, jj;

    int nx = 6;
    int nu = 2;
	int np = 1;

    int NF = nx + nu; // columns of forward seed

    double Ts = 0.2; // simulation time

	double *x_sim = malloc(sizeof(double)*nx*(nsim+1));

	for (ii=0; ii<nx; ii++)
		x_sim[ii] = x_ref[ii];

	/************************************************
	* external functions (explicit model)
	************************************************/

	// explicit ODE
	external_function_param_casadi expl_ode_fun;
	expl_ode_fun.casadi_fun = &casadi_expl_ode_fun;
	expl_ode_fun.casadi_work = &casadi_expl_ode_fun_work;
	expl_ode_fun.casadi_sparsity_in = &casadi_expl_ode_fun_sparsity_in;
	expl_ode_fun.casadi_sparsity_out = &casadi_expl_ode_fun_sparsity_out;
	expl_ode_fun.casadi_n_in = &casadi_expl_ode_fun_n_in;
	expl_ode_fun.casadi_n_out = &casadi_expl_ode_fun_n_out;
	external_function_param_casadi_create(&expl_ode_fun, np);

	// forward explicit VDE
	external_function_param_casadi expl_vde_for;
	expl_vde_for.casadi_fun = &casadi_expl_vde_for;
	expl_vde_for.casadi_work = &casadi_expl_vde_for_work;
	expl_vde_for.casadi_sparsity_in = &casadi_expl_vde_for_sparsity_in;
	expl_vde_for.casadi_sparsity_out = &casadi_expl_vde_for_sparsity_out;
	expl_vde_for.casadi_n_in = &casadi_expl_vde_for_n_in;
	expl_vde_for.casadi_n_out = &casadi_expl_vde_for_n_out;
	external_function_param_casadi_create(&expl_vde_for, np);

	/************************************************
	* sim plan & config
	************************************************/

	// choose plan
	sim_solver_plan plan;
	plan.sim_solver = ERK;

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
	opts->sens_adj = false;
	opts->sens_forw = true;

	opts->ns = 4; // number of stages in rk integrator
	opts->num_steps = 10; // number of integrationsteps

	/************************************************
	* sim in / out
	************************************************/

	sim_in *in = sim_in_create(config, dims);
	sim_out *out = sim_out_create(config, dims);

	in->T = Ts;

	// external functions
	sim_set_model(config, in, "expl_ode_fun", &expl_ode_fun);
	sim_set_model(config, in, "expl_vde_for", &expl_vde_for);

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

	// to avoid unstable behavior introduce a small pi-controller for rotor speed tracking
	double uctrl = 0.0;
	double uctrlI = 0.0;
	double kI = 1e-1;
	double kP = 10;
	double tmp, ctrlErr;

	for (ii=0; ii<nsim; ii++)
	{
		// update initial state
		for (jj = 0; jj < nx; jj++)
			in->x[jj] = x_sim[ii*nx+jj];

		// compute inputs
		for (jj = 0; jj < nu; jj++)
			in->u[jj] = u_sim[ii*nu+jj];
		tmp = in->u[1] - uctrl;
		in->u[1] = tmp>0.0 ? tmp : 0.0;

		// update parameters
		expl_ode_fun.set_param(&expl_ode_fun, p_sim+ii*np);
		expl_vde_for.set_param(&expl_vde_for, p_sim+ii*np);


		// d_print_mat(1, nx, in->x, 1);
		// d_print_mat(1, nu, in->u, 1);

		// execute simulation step with current input and state
		acados_return = sim_solve(sim_solver, in, out);
		if (acados_return != 0)
		{
			printf("error in sim solver\n");
			return ACADOS_FAILURE;
		}

		cpu_time += out->info->CPUtime;
		la_time += out->info->LAtime;
		ad_time += out->info->ADtime;

		// d_print_mat(1, nx, out->xn, 1);
		// d_print_mat(1, nx, x_ref+ii*nx, 1);

		// extract state at next time step
		for (jj = 0; jj < nx; jj++)
			x_sim[(ii+1)*nx+jj] = out->xn[jj];

		// update PI-controller
		ctrlErr = x_ref[nx*(ii+1)] - x_sim[nx*(ii+1)];
		uctrlI = uctrlI + kI*ctrlErr*Ts;
		uctrl = kP*ctrlErr + uctrlI;

		// if (ii < nsim-1)
		// 	printf("\nii = %d, sim error = %e\n", ii, ctrlErr);
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

#if 0
	printf("\n");
	printf("cpt: %8.4f [ms]\n", 1000*out->info->CPUtime);
	printf("AD cpt: %8.4f [ms]\n", 1000*out->info->ADtime);

#endif

	// printf("time split: %f ms CPU, %f ms LA, %f ms AD\n\n", cpu_time, la_time, ad_time);
	printf("\n\ntime for %d simulation steps: %f ms (AD time: %f ms (%5.2f%%))\n\n", nsim, 1e3*total_cpu_time, 1e3*ad_time, 1e2*ad_time/cpu_time);

	/************************************************
	* free memory
	************************************************/

	free(sim_solver);
	free(in);
	free(out);

	free(opts);
	free(config);

	external_function_param_casadi_free(&expl_ode_fun);
	external_function_param_casadi_free(&expl_vde_for);

	free(x_sim);

	printf("\nsuccess!\n\n");

    return 0;
}
