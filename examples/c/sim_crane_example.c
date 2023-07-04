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


// external
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// acados
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// crane model
#include "examples/c/crane_model/crane_model.h"

// blasfeo
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"



int main()
{

    /************************************************
    * Initialize
    ************************************************/

    int NREP = 500;

    int ii, jj;

    int nx = 4;
    int nu = 1;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05;

    double *xref = calloc(nx, sizeof(double));
    xref[1] = M_PI;

	double *uref = calloc(nu, sizeof(double));
	for (ii=0; ii<nu; ii++) uref[ii] = 1.0;

	double *Sx = calloc(nx*nx, sizeof(double));
	for (ii = 0; ii < nx; ii++)
		Sx[ii*(nx+1)] = 1.0;

	double *Su = calloc(nx*nu, sizeof(double));

	double *seed_adj = calloc(nx, sizeof(double));
	for (ii = 0; ii < nx; ii++)
		seed_adj[ii] = 1.0;

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
    external_function_casadi_create(&expl_vde_for);


    // adjoint explicit VDE

    external_function_casadi expl_vde_adj;
    expl_vde_adj.casadi_fun = &adjFun;
    expl_vde_adj.casadi_work = &adjFun_work;
    expl_vde_adj.casadi_sparsity_in = &adjFun_sparsity_in;
    expl_vde_adj.casadi_sparsity_out = &adjFun_sparsity_out;
    expl_vde_adj.casadi_n_in = &adjFun_n_in;
    expl_vde_adj.casadi_n_out = &adjFun_n_out;
    external_function_casadi_create(&expl_vde_adj);


    // hessian explicit ODE

    external_function_casadi expl_hess_ode;
    expl_hess_ode.casadi_fun = &hessFun;
    expl_hess_ode.casadi_work = &hessFun_work;
    expl_hess_ode.casadi_sparsity_in = &hessFun_sparsity_in;
    expl_hess_ode.casadi_sparsity_out = &hessFun_sparsity_out;
    expl_hess_ode.casadi_n_in = &hessFun_n_in;
    expl_hess_ode.casadi_n_out = &hessFun_n_out;
    external_function_casadi_create(&expl_hess_ode);


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
    external_function_casadi_create(&impl_ode_fun);

    //

    external_function_casadi impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_fun = &casadi_impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_work = &casadi_impl_ode_fun_jac_x_xdot_work;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_sparsity_in;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_sparsity_out;
    impl_ode_fun_jac_x_xdot.casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_n_in;
    impl_ode_fun_jac_x_xdot.casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot);

    //

    external_function_casadi impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_fun = &casadi_impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_work = &casadi_impl_ode_jac_x_xdot_u_work;
    impl_ode_jac_x_xdot_u.casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_sparsity_in;
    impl_ode_jac_x_xdot_u.casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_sparsity_out;
    impl_ode_jac_x_xdot_u.casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_n_in;
    impl_ode_jac_x_xdot_u.casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_n_out;
    external_function_casadi_create(&impl_ode_jac_x_xdot_u);

    int number_sim_solvers = 2;
    int nss;
    for (nss = 0; nss < number_sim_solvers; nss++)
    {
        /************************************************
        * sim plan & config
        ************************************************/

        // choose plan
        sim_solver_plan_t plan;

        switch (nss)
        {

            case 0:
                printf("\nsim solver: ERK\n");
                plan.sim_solver = ERK;
                break;

            case 1:
                printf("\nim solver: IRK\n");
                plan.sim_solver = IRK;
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

        /************************************************
        * sim opts
        ************************************************/

        sim_opts *opts = sim_opts_create(config, dims);

        int ns = 4; // number of stages in rk integrator
        int num_steps = 5; // number of integration steps
        bool sens_adj = true;
        sim_opts_set(config, opts, "ns", &ns);
        sim_opts_set(config, opts, "num_steps", &num_steps);
        sim_opts_set(config, opts, "sens_adj", &sens_adj);

        /************************************************
        * sim in
        ************************************************/

        sim_in *in = sim_in_create(config, dims);

		// matrices
		sim_in_set(config, dims, in, "T", &T);
		sim_in_set(config, dims, in, "x", xref);
		sim_in_set(config, dims, in, "u", uref);
		sim_in_set(config, dims, in, "Sx", Sx);
		sim_in_set(config, dims, in, "Su", Su);
		sim_in_set(config, dims, in, "seed_adj", seed_adj);

        // external functions
        switch (nss)
        {
            case 0:
            {
                sim_in_set(config, dims, in, "expl_vde_for", &expl_vde_for);
                sim_in_set(config, dims, in, "expl_vde_adj", &expl_vde_adj);
                break;
            }
            case 1:
            {
                sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
                sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot);
                sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
                break;
            }
            default :
            {
                printf("\nnot enough sim solvers implemented!\n");
                exit(1);
            }
        }

        /************************************************
        * sim out
        ************************************************/

        sim_out *out = sim_out_create(config, dims);

        /************************************************
        * sim solver
        ************************************************/

        sim_solver *sim_solver = sim_solver_create(config, dims, opts);

        int acados_return;

        // acados_timer timer;
        // acados_tic(&timer);

        for (ii=0;ii<NREP;ii++)
        {
            acados_return = sim_solve(sim_solver, in, out);
            if (acados_return != 0)
                printf("error in sim solver\n");
        }
        // double cpu_time = acados_toc(&timer)/NREP;

        /************************************************
        * printing
        ************************************************/

        double *xn_out = calloc(nx, sizeof(double));
		sim_out_get(config, dims, out, "xn", xn_out);
        printf("\nxn: \n");
        d_print_exp_mat(1, nx, xn_out, 1);

	    double *S_forw_out = calloc(nx*(nx+nu), sizeof(double));
        if (opts->sens_forw){
		    sim_out_get(config, dims, out, "S_forw", S_forw_out);
            printf("\nS_forw_out: \n");
            d_print_exp_mat(nx, NF, S_forw_out, nx);
        }

	    double *S_adj_out = calloc(nx+nu, sizeof(double));
        if (opts->sens_adj){
		    sim_out_get(config, dims, out, "S_adj", S_adj_out);
            printf("\nS_adj_out: \n");
            d_print_exp_mat(1, nx+nu, S_adj_out, 1);
        }

        double *S_hess_out = calloc((nu+nx)*(nu+nx), sizeof(double));
        if (opts->sens_hess)
        {
		    sim_out_get(config, dims, out, "S_hess", S_hess_out);
            printf("\nS_hess_out: \n");
            for (ii=0;ii<NF;ii++){
                for (jj=0;jj<NF;jj++){
                    if (jj>ii){
                        printf("%8.5f ", 0.0);
                    }else{
                        printf("%8.5f ", S_hess_out[jj*NF+ii]);
                    }
                }
                printf("\n");
            }
        }

        if(opts->sens_adj)
        {
            struct blasfeo_dmat sA;
            blasfeo_allocate_dmat(nx, nx+nu, &sA);
            blasfeo_pack_dmat(nx, nx+nu, S_forw_out, nx, &sA, 0, 0);

            struct blasfeo_dvec sx;
            blasfeo_allocate_dvec(nx, &sx);
            blasfeo_pack_dvec(nx, in->S_adj, 1, &sx, 0);

            struct blasfeo_dvec sz;
            blasfeo_allocate_dvec(nx+nu, &sz);
            // blasfeo_print_dmat(nx, nx+nu, &sA, 0, 0);
            // blasfeo_print_tran_dvec(nx, &sx, 0);
            blasfeo_dgemv_t(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);

            printf("\nJac times lambdaX:\n");
            blasfeo_print_exp_tran_dvec(nx+nu, &sz, 0);

            blasfeo_free_dmat(&sA);
            blasfeo_free_dvec(&sx);
            blasfeo_free_dvec(&sz);
        }

        printf("\n");
        printf("cpt: %8.4f [ms]\n", 1000*out->info->CPUtime);
        printf("AD cpt: %8.4f [ms]\n", 1000*out->info->ADtime);
        printf("========================\n");

		free(xn_out);
		free(S_forw_out);
		free(S_adj_out);
		free(S_hess_out);

        sim_solver_destroy(sim_solver);
        sim_in_destroy(in);
        sim_out_destroy(out);
        sim_opts_destroy(opts);
        sim_dims_destroy(dims);
        sim_config_destroy(config);
    }

	free(xref);
	free(uref);
	free(Sx);
	free(Su);
	free(seed_adj);

    // explicit model
    external_function_casadi_free(&expl_vde_for);
    external_function_casadi_free(&expl_vde_adj);
    external_function_casadi_free(&expl_hess_ode);
    // implicit model
    external_function_casadi_free(&impl_ode_fun);
    external_function_casadi_free(&impl_ode_fun_jac_x_xdot);
    external_function_casadi_free(&impl_ode_jac_x_xdot_u);

    /************************************************
    * return
    ************************************************/

    printf("\nsuccess! (RESULT NOT CHECKED) \n\n");

    return 0;
}
