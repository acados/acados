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
 *    Author: Jonathan Frey, Robin Verschueren
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

    const int n_sim = 601;
    const int N = 20;

    int nx[N+1], nu[N+1], nz[N+1], ny[N+1], nb[N+1], nbx[N+1], nbu[N+1], ng[N+1], nh[N+1], ns[N+1];
    for (int i = 0; i <= N; ++i)
    {
        nx[i] = 4;
        nu[i] = 2;
        nz[i] = 2;
        ny[i] = 1 + nx[i] + nu[i];
        nbx[i] = nx[i];
        nbu[i] = nu[i];
        nb[i] = nbx[i] + nbu[i];
        ng[i] = 0;
        nh[i] = 0;
        ns[i] = 0;
    }

    nu[N] = 0;
    nz[N] = 0;
    ny[N] = 1 + nx[N];
    nbu[N] = 0;
    nb[N] = nbx[N] + nbu[N];

    int idxb[nbu[0]+nbx[0]];
    for (int i = 0; i < nbu[0]+nbx[0]; ++i)
        idxb[i] = i;

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
    double lb_0[] = {-10000, -10000, 50, 50, 1.14275, 1.53787};
    double ub_0[] = {+10000, +10000, 50, 50, 1.14275, 1.53787};

    double lb[] = {-10000, -10000, 0, 0, 0.5, 0.5};
    double ub[] = {+10000, +10000, 100, 100, 1.757, 2.125};

    double lb_N[] = {0, 0, 0.5, 0.5};
    double ub_N[] = {100, 100, 1.757, 2.125};

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

	ocp_nlp_plan *plan = ocp_nlp_plan_create(N);

	plan->nlp_solver = SQP;

	for (int i = 0; i <= N; i++)
		plan->nlp_cost[i] = NONLINEAR_LS;

	plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_QPOASES;

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
    ocp_nlp_dims_set_opt_vars(config, dims, "ns", ns);

	for (int i = 0; i <= N; i++)
    {
        ocp_nlp_dims_set_cost(config, dims, i, "ny", &ny[i]);

        ocp_nlp_dims_set_constraints(config, dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(config, dims, i, "nh", &nh[i]);
    }


    // in
	ocp_nlp_in *nlp_in = ocp_nlp_in_create(config, dims);
	for (int i = 0; i < N; ++i)
    	nlp_in->Ts[i] = T;

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

    // bounds
	ocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) nlp_in->constraints;

    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lb", lb_0);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ub", ub_0);

    for (int i = 0; i < nb[0]; ++i)
        constraints[0]->idxb[i] = idxb[i];
    
    for (int i = 1; i < N; ++i)
    {
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "lb", lb);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, i, "ub", ub);

        for (int j = 0; j < nb[i]; ++j)
            constraints[i]->idxb[j] = idxb[j];
    }

    ocp_nlp_constraints_model_set(config, dims, nlp_in, N, "lb", lb_N);
    ocp_nlp_constraints_model_set(config, dims, nlp_in, N, "ub", ub_N);
    for (int i = 0; i < nb[N]; ++i)
        constraints[N]->idxb[i] = idxb[i];

    // options
    void *nlp_opts = ocp_nlp_opts_create(config, dims);
    int maxIter = 1;
    ocp_nlp_opts_set(config, nlp_opts, "maxIter", &maxIter);

    // out
    ocp_nlp_out *nlp_out = ocp_nlp_out_create(config, dims);

    // solver
	ocp_nlp_solver *solver = ocp_nlp_solver_create(config, dims, nlp_opts);

    // initialize
    for (int i = 0; i < N; ++i)
    {
        blasfeo_pack_dvec(nu[i], u, nlp_out->ux+i, 0);
        blasfeo_pack_dvec(nx[i], x0, nlp_out->ux+i, nu[i]);
        blasfeo_pack_dvec(nz[i], z0, nlp_out->z+i, 0);
    }
    blasfeo_pack_dvec(nx[N], x0, nlp_out->ux+N, nu[N]);

    status = ocp_nlp_precompute(solver, nlp_in, nlp_out);

    for (int i = 0; i < n_sim; ++i)
    {
        printf("\n-----\n");
        y_ref[0] = reference[i];

        for (int j = 0; j <= N; ++j)
            status = ocp_nlp_cost_model_set(config, dims, nlp_in, j, "y_ref", y_ref);

        status = ocp_nlp_solve(solver, nlp_in, nlp_out);
        
        blasfeo_print_to_file_dvec(out_file, nu[0]+nx[0], nlp_out->ux, 0);

        blasfeo_unpack_dvec(nx[1], nlp_out->ux+1, nu[1], &lb_0[nu[1]]);
        blasfeo_unpack_dvec(nx[1], nlp_out->ux+1, nu[1], &ub_0[nu[1]]);

        ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "lb", lb_0);
        ocp_nlp_constraints_model_set(config, dims, nlp_in, 0, "ub", ub_0);


        printf("iter %d, status: %d, res = %g\n", i, status, nlp_out->inf_norm_res);
    }

    fclose(out_file);

    // free memory
    ocp_nlp_solver_destroy(solver);
    ocp_nlp_out_destroy(nlp_out);
    ocp_nlp_opts_destroy(nlp_opts);
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
