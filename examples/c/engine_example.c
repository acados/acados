/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include <stdlib.h>
// #include <xmmintrin.h>

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// acados
#include "acados/utils/print.h"
#include "acados_c/external_function_interface.h"
#include "acados_c/ocp_nlp_interface.h"

// example specific includes
#include "examples/c/engine_model/engine_impl_dae_fun.h"
#include "examples/c/engine_model/engine_impl_dae_fun_jac_x_xdot_z.h"
#include "examples/c/engine_model/engine_impl_dae_jac_x_xdot_u_z.h"
#include "examples/c/engine_model/engine_impl_dae_fun_jac_x_xdot_u_z.h"
#include "examples/c/engine_model/engine_ls_cost.h"
#include "examples/c/engine_model/engine_ls_cost_N.h"
#include "examples/c/engine_model/reference.c"

int main()
{
    // _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

    FILE *out_file = fopen("engine_ux.txt", "w");
    if(ferror(out_file))
        exit(1);

    int jj, sqp_iter;
    double time_tot;
    // double *stat;

    const int n_sim = 601;
    const int N = 20;

	int nx_ = 4;
	int nu_ = 2;
	int nz_ = 2;

    int nx[N+1], nu[N+1], nz[N+1], ny[N+1], nbx[N+1], nbu[N+1];
    for (int i = 0; i < N; i++)
    {
        nx[i] = nx_;
        nu[i] = nu_;
        nz[i] = nz_;
        ny[i] = 1 + nx_ + nu_;
        nbx[i] = nx_;
        nbu[i] = nu_;
    }
	nx[N] = nx_;
	nu[N] = 0;
	nz[N] = 0;
	ny[N] = 1 + nx_;
	nbx[N] = nx_;
	nbu[N] = 0;

    // sampling time (s)
    double T = 0.05;

    // x: u1, u2, xD1, xD2
    double x0[] = {50, 50, 1.14275, 1.53787};
    // z: xA1, xA2
    double z0[] = {1.28976, 1.78264};
    // u: u1_r, u2_r
    double u[] = {0, 0};

    // reference
    double y_ref[] = {1500, 50, 50, 1.14275, 1.53787, 0, 0};

    // weighting matrices
    double W[ny[0]*ny[0]];
    for (int i = 0; i < ny[0]*ny[0]; ++i)
        W[i] = 0;
    W[0*(ny[0]+1)] = 100;
    W[1*(ny[0]+1)] = 1;
    W[2*(ny[0]+1)] = 1;
    W[3*(ny[0]+1)] = 1;
    W[4*(ny[0]+1)] = 1;
    W[5*(ny[0]+1)] = 0.1;
    W[6*(ny[0]+1)] = 0.1;

    double W_N[ny[N]*ny[N]];
    for (int i = 0; i < ny[N]*ny[N]; ++i)
        W_N[i] = 0;
    W_N[0*(ny[N]+1)] = 100;
    W_N[1*(ny[N]+1)] = 1;
    W_N[2*(ny[N]+1)] = 1;
    W_N[3*(ny[N]+1)] = 1;
    W_N[4*(ny[N]+1)] = 1;

    // bounds
    double lbu[] = {-10000, -10000};
    double ubu[] = {10000, 10000};
	int idxbu[] = {0, 1};

    double lbx[] = {0, 0, 0.5, 0.5};
    double ubx[] = {100, 100, 1.757, 2.125};
	int idxbx[] = {0, 1, 2, 3};

    // implicit dae
    external_function_casadi impl_dae_fun;
    impl_dae_fun.casadi_fun = &engine_impl_dae_fun;
    impl_dae_fun.casadi_work = &engine_impl_dae_fun_work;
    impl_dae_fun.casadi_sparsity_in = &engine_impl_dae_fun_sparsity_in;
    impl_dae_fun.casadi_sparsity_out = &engine_impl_dae_fun_sparsity_out;
    impl_dae_fun.casadi_n_in = &engine_impl_dae_fun_n_in;
    impl_dae_fun.casadi_n_out = &engine_impl_dae_fun_n_out;
    external_function_casadi_create(&impl_dae_fun);

    external_function_casadi impl_dae_fun_jac_x_xdot_z;
    impl_dae_fun_jac_x_xdot_z.casadi_fun = &engine_impl_dae_fun_jac_x_xdot_z;
    impl_dae_fun_jac_x_xdot_z.casadi_work = &engine_impl_dae_fun_jac_x_xdot_z_work;
    impl_dae_fun_jac_x_xdot_z.casadi_sparsity_in = &engine_impl_dae_fun_jac_x_xdot_z_sparsity_in;
    impl_dae_fun_jac_x_xdot_z.casadi_sparsity_out = &engine_impl_dae_fun_jac_x_xdot_z_sparsity_out;
    impl_dae_fun_jac_x_xdot_z.casadi_n_in = &engine_impl_dae_fun_jac_x_xdot_z_n_in;
    impl_dae_fun_jac_x_xdot_z.casadi_n_out = &engine_impl_dae_fun_jac_x_xdot_z_n_out;
    external_function_casadi_create(&impl_dae_fun_jac_x_xdot_z);

    external_function_casadi impl_dae_jac_x_xdot_u_z;
    impl_dae_jac_x_xdot_u_z.casadi_fun = &engine_impl_dae_jac_x_xdot_u_z;
    impl_dae_jac_x_xdot_u_z.casadi_work = &engine_impl_dae_jac_x_xdot_u_z_work;
    impl_dae_jac_x_xdot_u_z.casadi_sparsity_in = &engine_impl_dae_jac_x_xdot_u_z_sparsity_in;
    impl_dae_jac_x_xdot_u_z.casadi_sparsity_out = &engine_impl_dae_jac_x_xdot_u_z_sparsity_out;
    impl_dae_jac_x_xdot_u_z.casadi_n_in = &engine_impl_dae_jac_x_xdot_u_z_n_in;
    impl_dae_jac_x_xdot_u_z.casadi_n_out = &engine_impl_dae_jac_x_xdot_u_z_n_out;
    external_function_casadi_create(&impl_dae_jac_x_xdot_u_z);

    // Only needed for lifted IRK:

    // external_function_casadi engine_impl_dae_fun_jac_x_xdot_u_z;
    // engine_impl_dae_fun_jac_x_xdot_u_z.casadi_fun = &engine_impl_dae_fun_jac_x_xdot_u_z;
    // engine_impl_dae_fun_jac_x_xdot_u_z.casadi_work = &engine_impl_dae_fun_jac_x_xdot_u_z_work;
    // engine_impl_dae_fun_jac_x_xdot_u_z.casadi_sparsity_in = &engine_impl_dae_fun_jac_x_xdot_u_z_sparsity_in;
    // engine_impl_dae_fun_jac_x_xdot_u_z.casadi_sparsity_out = &engine_impl_dae_fun_jac_x_xdot_u_z_sparsity_out;
    // engine_impl_dae_fun_jac_x_xdot_u_z.casadi_n_in = &engine_impl_dae_fun_jac_x_xdot_u_z_n_in;
    // engine_impl_dae_fun_jac_x_xdot_u_z.casadi_n_out = &engine_impl_dae_fun_jac_x_xdot_u_z_n_out;
    // external_function_casadi_create(&engine_impl_dae_fun_jac_x_xdot_u_z);

    external_function_casadi nls_cost_residual;
    nls_cost_residual.casadi_fun = &engine_ls_cost;
    nls_cost_residual.casadi_work = &engine_ls_cost_work;
    nls_cost_residual.casadi_sparsity_in = &engine_ls_cost_sparsity_in;
    nls_cost_residual.casadi_sparsity_out = &engine_ls_cost_sparsity_out;
    nls_cost_residual.casadi_n_in = &engine_ls_cost_n_in;
    nls_cost_residual.casadi_n_out = &engine_ls_cost_n_out;
    external_function_casadi_create(&nls_cost_residual);

    external_function_casadi nls_cost_N_residual;
    nls_cost_N_residual.casadi_fun = &engine_ls_cost_N;
    nls_cost_N_residual.casadi_work = &engine_ls_cost_N_work;
    nls_cost_N_residual.casadi_sparsity_in = &engine_ls_cost_N_sparsity_in;
    nls_cost_N_residual.casadi_sparsity_out = &engine_ls_cost_N_sparsity_out;
    nls_cost_N_residual.casadi_n_in = &engine_ls_cost_N_n_in;
    nls_cost_N_residual.casadi_n_out = &engine_ls_cost_N_n_out;
    external_function_casadi_create(&nls_cost_N_residual);

	ocp_nlp_plan_t *plan = ocp_nlp_plan_create(N);

	plan->nlp_solver = SQP;
//	plan->nlp_solver = SQP_RTI;

	for (int i = 0; i <= N; i++)
		plan->nlp_cost[i] = NONLINEAR_LS;

	plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;
//	plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_HPIPM;
//	plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_QPOASES;

	for (int i = 0; i < N; i++)
    {
		plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
		plan->sim_solver_plan[i].sim_solver = IRK;
	}

	for (int i = 0; i <= N; i++)
		plan->nlp_constraints[i] = BGH;

	ocp_nlp_config *config = ocp_nlp_config_create(*plan);

    // dimensions
    ocp_nlp_dims *dims = ocp_nlp_dims_create(config);

    ocp_nlp_dims_set_opt_vars(config, dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(config, dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(config, dims, "nz", nz);

	for (int i = 0; i <= N; i++)
    {
        ocp_nlp_dims_set_cost(config, dims, i, "ny", &ny[i]);

        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", &nbu[i]);
    }


    // in
	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);

	for (int i = 0; i < N; ++i)
		ocp_nlp_in_set(config, dims, nlp_in, i, "Ts", &T);

    // cost
    // ocp_nlp_cost_nls_model **cost = (ocp_nlp_cost_nls_model **) nlp_in->cost;
    int status = ACADOS_SUCCESS;

	for (int i = 0; i < N; ++i) {
        if(ocp_nlp_cost_model_set(config, dims, nlp_in, i, "nls_res_jac", &nls_cost_residual)) exit(1);
        if(ocp_nlp_cost_model_set(config, dims, nlp_in, i, "y_ref", y_ref)) exit(1);
        if(ocp_nlp_cost_model_set(config, dims, nlp_in, i, "W", W)) exit(1);
    }

    if(ocp_nlp_cost_model_set(config, dims, nlp_in, N, "nls_res_jac", &nls_cost_N_residual)) exit(1);
    if(ocp_nlp_cost_model_set(config, dims, nlp_in, N, "y_ref", y_ref)) exit(1);
    if(ocp_nlp_cost_model_set(config, dims, nlp_in, N, "W", W_N)) exit(1);

    // dynamics
    for (int i = 0; i < N; ++i)
    {
        if(ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun", &impl_dae_fun)) exit(1);
        if(ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_fun_jac_x_xdot", &impl_dae_fun_jac_x_xdot_z)) exit(1);
        if(ocp_nlp_dynamics_model_set(config, dims, nlp_in, i, "impl_ode_jac_x_xdot_u", &impl_dae_jac_x_xdot_u_z)) exit(1);
    }

    // constraints
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lbx", x0);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ubx", x0);
	ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "idxbx", idxbx);

    for (int i = 1; i <= N; i++)
    {
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lbx", lbx);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "ubx", ubx);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "idxbx", idxbx);
    }

    for (int i = 0; i < N; i++)
    {
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lbu", lbu);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "ubu", ubu);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "idxbu", idxbu);
    }

    // options
    void *nlp_opts = ocp_nlp_solver_opts_create(config, dims);
    int max_iter = 20;
    ocp_nlp_solver_opts_set(config, nlp_opts, "max_iter", &max_iter);
	int ext_qp_res = 1;
    ocp_nlp_solver_opts_set(config, nlp_opts, "ext_qp_res", &ext_qp_res);
	int qp_warm_start = 2;
    ocp_nlp_solver_opts_set(config, nlp_opts, "qp_warm_start", &qp_warm_start);
	double tol = 1e-6;
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_stat", &tol);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_eq", &tol);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_ineq", &tol);
    ocp_nlp_solver_opts_set(config, nlp_opts, "tol_comp", &tol);

    // out
    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

    // solver
	ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, nlp_opts);

    // initialize
    for (int i = 0; i < N; ++i)
    {
		ocp_nlp_out_set(config, dims, nlp_out, i, "u", u);
		ocp_nlp_out_set(config, dims, nlp_out, i, "x", x0);
		ocp_nlp_out_set(config, dims, nlp_out, i, "z", z0);
    }
	ocp_nlp_out_set(config, dims, nlp_out, N, "x", x0);

    status = ocp_nlp_precompute(solver, nlp_in, nlp_out);

	double *x_sol = malloc(nx_*(N+1)*sizeof(double));
	double *u_sol = malloc(nu_*(N+1)*sizeof(double));
	double *z_sol = malloc(nz_*(N+1)*sizeof(double));

    for (int i = 0; i < n_sim; ++i)
    {
        printf("\n-----\n");
        y_ref[0] = reference[i];

        for (int j = 0; j <= N; ++j)
            status = ocp_nlp_cost_model_set(config, dims, nlp_in, j, "y_ref", y_ref);

        status = ocp_nlp_solve(solver, nlp_in, nlp_out);

        blasfeo_print_to_file_dvec(out_file, nu[0]+nx[0], nlp_out->ux, 0);

		// get solution
		for(jj=0; jj<=N; jj++)
			ocp_nlp_out_get(config, dims, nlp_out, jj, "x", x_sol+jj*nx_);
		for(jj=0; jj<N; jj++)
			ocp_nlp_out_get(config, dims, nlp_out, jj, "u", u_sol+jj*nu_);
		for(jj=0; jj<N; jj++)
			ocp_nlp_out_get(config, dims, nlp_out, jj, "z", z_sol+jj*nz_);

		// set x0
        ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lbx", x_sol+1*nx_);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ubx", x_sol+1*nx_);

		// shift guess
		for(jj=0; jj<N; jj++)
			ocp_nlp_out_set(config, dims, nlp_out, jj, "x", x_sol+(jj+1)*nx_);
		for(jj=0; jj<N-1; jj++)
			ocp_nlp_out_set(config, dims, nlp_out, jj, "u", u_sol+(jj+1)*nu_);
		for(jj=0; jj<N-1; jj++)
			ocp_nlp_out_set(config, dims, nlp_out, jj, "z", z_sol+(jj+1)*nz_);

		config->get(config, dims, solver->mem, "sqp_iter", &sqp_iter);
		config->get(config, dims, solver->mem, "time_tot", &time_tot);
#if 0
		config->get(config, dims, solver->mem, "stat_m", &stat_m);
		config->get(config, dims, solver->mem, "stat_n", &stat_n);
		config->get(config, dims, solver->mem, "stat", &stat);
		printf("\niter\tres_g\t\tres_b\t\tres_d\t\tres_m\t\tqp_stat\tqp_iter\tqp_res_g\tqp_res_b\tqp_res_d\tqp_res_m\t");
		for(jj=0; jj<sqp_iter+1; jj++)
		{
			if(stat_n==10)
			{
				printf("\n%d\t%e\t%e\t%e\t%e\t%d\t%d\t%e\t%e\t%e\t%e\n", jj, stat[jj*stat_n+0], stat[jj*stat_n+1], stat[jj*stat_n+2], stat[jj*stat_n+3], (int) stat[jj*stat_n+4], (int) stat[jj*stat_n+5], stat[jj*stat_n+6], stat[jj*stat_n+7], stat[jj*stat_n+8], stat[jj*stat_n+9]);
			}
		}
#endif

        printf("n_sim %d, time: %e [s], status: %d, iter: %d, res = %g\n", i, time_tot, status, sqp_iter, nlp_out->inf_norm_res);

#if 0
	if(i==5)
		exit(1);
#endif
    }

    fclose(out_file);

    // free memory
	free(x_sol);
	free(u_sol);
	free(z_sol);

    ocp_nlp_solver_destroy(solver);
    ocp_nlp_out_destroy(nlp_out);
    ocp_nlp_solver_opts_destroy(nlp_opts);
    ocp_nlp_in_destroy(nlp_in);
    ocp_nlp_dims_destroy(dims);
    ocp_nlp_config_destroy(config);
    ocp_nlp_plan_destroy(plan);

    /* free external function */
    // implicit model
    external_function_casadi_free(&impl_dae_fun);
    external_function_casadi_free(&impl_dae_fun_jac_x_xdot_z);
    external_function_casadi_free(&impl_dae_jac_x_xdot_u_z);
    // external_function_casadi_free(&impl_ode_jac_x_xdot_u);
    external_function_casadi_free(&nls_cost_residual);
    external_function_casadi_free(&nls_cost_N_residual);

    return 0;
}
