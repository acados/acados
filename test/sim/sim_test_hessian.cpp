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
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"

// acados
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/math.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// pendulum_model
#include "examples/c/pendulum_model/pendulum_model.h"


extern "C"
{
    // initialize dimensions
    const int nx = 4;
    const int nu = 1;
    const int nz = 0;

    const int NF = nx + nu;  // columns of forward seed

    const bool PRINT_HESS_RESULTS = false;  // can be used for debugging

    int acados_return;

    // generate x0, u_sim
    double x0_pendulum[nx];
    double u_sim_pendulum[nu];

    double x_ref_sol[nx];
    double S_forw_ref_sol[nx*NF];
    double S_adj_ref_sol[NF];
    double S_hess_ref_sol[NF * NF];

    double error[nx];
    double error_z[nz];
    double error_S_forw[nx*NF];
    double error_S_adj[NF];
    double error_S_hess[NF * NF];

    double norm_x_ref, norm_S_forw_ref, norm_S_adj_ref, norm_S_hess_ref;
    double norm_error, norm_error_forw, norm_error_adj, norm_error_hess;
    double rel_error_forw, rel_error_adj, rel_error_hess;
}

using std::vector;

sim_solver_t hashitsim_hess(std::string const& inString)
{
    if (inString == "ERK") return ERK;
    if (inString == "IRK") return IRK;

    return (sim_solver_t) -1;
}

double sim_solver_tolerance_sim(std::string const& inString)
{
    if (inString == "IRK")  return 1e-9;
    if (inString == "ERK") return 1e-5;

    return -1;
}

double sim_solver_tolerance_hess(std::string const& inString)
{
    if (inString == "IRK")  return 1e-8;
    if (inString == "ERK") return 1e-5;

    return -1;
}

TEST_CASE("pendulum_hessians", "[integrators]")
{
    vector<std::string> solvers = {"IRK", "ERK"};

    for (int ii = 0; ii < nx; ii++)
        x0_pendulum[ii] = 0.0;

    u_sim_pendulum[0] = 0.1;


    int nsim0 = 1;  // nsim;

    double T = 0.1;  // simulation time
    // former value 0.5

    double x_sim[nx*(nsim0+2)];

    for (int ii = 0; ii < nx; ii++)
        x_sim[ii] = x0_pendulum[ii];

/************************************************
* external functions
************************************************/
    /* IMPLICIT MODEL */
    // impl_ode_fun
    external_function_casadi impl_ode_fun;
    impl_ode_fun.casadi_fun = &pendulum_ode_impl_ode_fun;
    impl_ode_fun.casadi_work = &pendulum_ode_impl_ode_fun_work;
    impl_ode_fun.casadi_sparsity_in = &pendulum_ode_impl_ode_fun_sparsity_in;
    impl_ode_fun.casadi_sparsity_out = &pendulum_ode_impl_ode_fun_sparsity_out;
    impl_ode_fun.casadi_n_in = &pendulum_ode_impl_ode_fun_n_in;
    impl_ode_fun.casadi_n_out = &pendulum_ode_impl_ode_fun_n_out;
    external_function_casadi_create(&impl_ode_fun);

    // impl_ode_fun_jac_x_xdot
    external_function_casadi impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_fun = &pendulum_ode_impl_ode_fun_jac_x_xdot_z;
    impl_ode_fun_jac_x_xdot.casadi_work = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_work;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_sparsity_in;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_out =
                         &pendulum_ode_impl_ode_fun_jac_x_xdot_z_sparsity_out;
    impl_ode_fun_jac_x_xdot.casadi_n_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_n_in;
    impl_ode_fun_jac_x_xdot.casadi_n_out = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_fun = &pendulum_ode_impl_ode_jac_x_xdot_u_z;
    impl_ode_jac_x_xdot_u.casadi_work = &pendulum_ode_impl_ode_jac_x_xdot_u_z_work;
    impl_ode_jac_x_xdot_u.casadi_sparsity_in = &pendulum_ode_impl_ode_jac_x_xdot_u_z_sparsity_in;
    impl_ode_jac_x_xdot_u.casadi_sparsity_out = &pendulum_ode_impl_ode_jac_x_xdot_u_z_sparsity_out;
    impl_ode_jac_x_xdot_u.casadi_n_in = &pendulum_ode_impl_ode_jac_x_xdot_u_z_n_in;
    impl_ode_jac_x_xdot_u.casadi_n_out = &pendulum_ode_impl_ode_jac_x_xdot_u_z_n_out;
    external_function_casadi_create(&impl_ode_jac_x_xdot_u);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_fun = &pendulum_ode_impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_work = &pendulum_ode_impl_ode_fun_jac_x_xdot_u_work;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_in =
                            &pendulum_ode_impl_ode_fun_jac_x_xdot_u_sparsity_in;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_out =
                            &pendulum_ode_impl_ode_fun_jac_x_xdot_u_sparsity_out;
    impl_ode_fun_jac_x_xdot_u.casadi_n_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_u_n_in;
    impl_ode_fun_jac_x_xdot_u.casadi_n_out = &pendulum_ode_impl_ode_fun_jac_x_xdot_u_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot_u);

    // impl_ode_hess
    external_function_casadi impl_ode_hess;
    impl_ode_hess.casadi_fun = &pendulum_ode_impl_ode_hess;
    impl_ode_hess.casadi_work = &pendulum_ode_impl_ode_hess_work;
    impl_ode_hess.casadi_sparsity_in = &pendulum_ode_impl_ode_hess_sparsity_in;
    impl_ode_hess.casadi_sparsity_out = &pendulum_ode_impl_ode_hess_sparsity_out;
    impl_ode_hess.casadi_n_in = &pendulum_ode_impl_ode_hess_n_in;
    impl_ode_hess.casadi_n_out = &pendulum_ode_impl_ode_hess_n_out;
    external_function_casadi_create(&impl_ode_hess);

    /* EXPLICIT MODEL */
    // expl_ode_fun
    external_function_casadi expl_ode_fun;
    expl_ode_fun.casadi_fun = &pendulum_ode_expl_ode_fun;
    expl_ode_fun.casadi_work = &pendulum_ode_expl_ode_fun_work;
    expl_ode_fun.casadi_sparsity_in = &pendulum_ode_expl_ode_fun_sparsity_in;
    expl_ode_fun.casadi_sparsity_out = &pendulum_ode_expl_ode_fun_sparsity_out;
    expl_ode_fun.casadi_n_in = &pendulum_ode_expl_ode_fun_n_in;
    expl_ode_fun.casadi_n_out = &pendulum_ode_expl_ode_fun_n_out;
    external_function_casadi_create(&expl_ode_fun);

    // expl_vde_for
    external_function_casadi expl_vde_for;
    expl_vde_for.casadi_fun = &pendulum_ode_expl_vde_forw;
    expl_vde_for.casadi_work = &pendulum_ode_expl_vde_forw_work;
    expl_vde_for.casadi_sparsity_in = &pendulum_ode_expl_vde_forw_sparsity_in;
    expl_vde_for.casadi_sparsity_out = &pendulum_ode_expl_vde_forw_sparsity_out;
    expl_vde_for.casadi_n_in = &pendulum_ode_expl_vde_forw_n_in;
    expl_vde_for.casadi_n_out = &pendulum_ode_expl_vde_forw_n_out;
    external_function_casadi_create(&expl_vde_for);

    // expl_vde_adj
    external_function_casadi expl_vde_adj;
    expl_vde_adj.casadi_fun = &pendulum_ode_expl_vde_adj;
    expl_vde_adj.casadi_work = &pendulum_ode_expl_vde_adj_work;
    expl_vde_adj.casadi_sparsity_in = &pendulum_ode_expl_vde_adj_sparsity_in;
    expl_vde_adj.casadi_sparsity_out = &pendulum_ode_expl_vde_adj_sparsity_out;
    expl_vde_adj.casadi_n_in = &pendulum_ode_expl_vde_adj_n_in;
    expl_vde_adj.casadi_n_out = &pendulum_ode_expl_vde_adj_n_out;
    external_function_casadi_create(&expl_vde_adj);

    // expl_ode_hess
    external_function_casadi expl_ode_hess;
    expl_ode_hess.casadi_fun = &pendulum_ode_expl_ode_hess;
    expl_ode_hess.casadi_work = &pendulum_ode_expl_ode_hess_work;
    expl_ode_hess.casadi_sparsity_in = &pendulum_ode_expl_ode_hess_sparsity_in;
    expl_ode_hess.casadi_sparsity_out = &pendulum_ode_expl_ode_hess_sparsity_out;
    expl_ode_hess.casadi_n_in = &pendulum_ode_expl_ode_hess_n_in;
    expl_ode_hess.casadi_n_out = &pendulum_ode_expl_ode_hess_n_out;
    external_function_casadi_create(&expl_ode_hess);

/************************************************
* Create Reference Solution
************************************************/

    sim_solver_plan_t plan;
    plan.sim_solver = IRK;  // IRK

    sim_config *config = sim_config_create(plan);

    void *dims = sim_dims_create(config);

    /* set dimensions */
    sim_dims_set(config, dims, "nx", &nx);
    sim_dims_set(config, dims, "nu", &nu);
    sim_dims_set(config, dims, "nz", &nz);

    // set opts
    void *opts_ = sim_opts_create(config, dims);
    sim_opts *opts = (sim_opts *) opts_;
    config->opts_initialize_default(config, dims, opts);

    // opts reference solution
    opts->sens_forw         = true;
    opts->sens_adj          = true;
    opts->sens_algebraic    = false;
    opts->sens_hess         = true;
    opts->output_z          = false;
    opts->jac_reuse         = false;  // jacobian reuse
    opts->newton_iter       = 8;    // number of newton iterations per integration step
    opts->num_steps         = 200;  // number of steps
    opts->ns                = 4;    // number of stages in rk integrator

    sim_in *in = sim_in_create(config, dims);
    sim_out *out = sim_out_create(config, dims);

    in->T = T;

    // set model
    switch (plan.sim_solver)
    {
        case ERK:  // ERK
        {
            sim_in_set(config, dims, in, "expl_ode_fun", &expl_ode_fun);
            sim_in_set(config, dims, in, "expl_vde_for", &expl_vde_for);
            sim_in_set(config, dims, in, "expl_vde_adj", &expl_vde_adj);
            sim_in_set(config, dims, in, "expl_ode_hess", &expl_ode_hess);
            break;
        }
        case IRK:  // IRK
        {
            sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
            sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot",
                    &impl_ode_fun_jac_x_xdot);
            sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
            sim_in_set(config, dims, in, "impl_ode_hess", &impl_ode_hess);
            break;
        }
        default :
        {
            printf("\nnot plan.sim_solver not supported!\n");
            exit(1);
        }
    }

    // seeds forw
    for (int ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (int ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;
    in->identity_seed = true;


    // seeds adj
    for (int ii = 0; ii < nx; ii++)
        in->S_adj[ii] = 1.0;
    for (int ii = nx; ii < nx + nu; ii++)
        in->S_adj[ii] = 0.0;

    /************************************************
    * sim solver
    ************************************************/

    sim_solver *sim_solver = sim_solver_create(config, dims, opts);

    for (int ii = 0; ii < nsim0; ii++)
    {
        // x
        for (int jj = 0; jj < nx; jj++)
            in->x[jj] = x_sim[ii*nx+jj];

        // u
        for (int jj = 0; jj < nu; jj++)
            in->u[jj] = u_sim_pendulum[ii*nu+jj];

        acados_return = sim_solve(sim_solver, in, out);
        REQUIRE(acados_return == 0);

        for (int jj = 0; jj < nx; jj++)
            x_sim[(ii+1)*nx+jj] = out->xn[jj];
    }

    // store reference solution
    for (int jj = 0; jj < nx; jj++)
        x_ref_sol[jj] = out->xn[jj];

    for (int jj = 0; jj < nx*NF; jj++)
        S_forw_ref_sol[jj] = out->S_forw[jj];

    for (int jj = 0; jj < NF; jj++)
        S_adj_ref_sol[jj] = out->S_adj[jj];

    for (int jj = 0; jj < NF * NF; jj++)
        S_hess_ref_sol[jj] = out->S_hess[jj];

    norm_x_ref = onenorm(nx, 1, x_ref_sol);
    norm_S_forw_ref = onenorm(nx, nx + nu, S_forw_ref_sol);
    norm_S_adj_ref = onenorm(1, nx + nu, S_adj_ref_sol);
    norm_S_hess_ref = onenorm(nx + nu, nx + nu, S_hess_ref_sol);

    // printf("Reference xn \n");
    // d_print_e_mat(1, nx, &x_ref_sol[0], 1);

    // printf("Reference forward sensitivities \n");
    // d_print_e_mat(nx, NF, &S_forw_ref_sol[0], nx);

    // printf("reference adjoint sensitivities \n");
    // d_print_e_mat(1, nx + nu, &S_adj_ref_sol[0], 1);
    if ( PRINT_HESS_RESULTS ){
        printf("Reference Hessian \n");
        d_print_exp_mat(nx + nu, nx + nu, out->S_hess, nx + nu);
    }

    /* free */
    free(config);
    free(dims);
    free(opts);

    free(in);
    free(out);
    free(sim_solver);

/************************************************
* test solver loop
************************************************/



    for (int sens_forw = 0; sens_forw < 2; sens_forw++)
    {
    SECTION("sens_forw = " + std::to_string((bool)sens_forw))
    {
        for (int sens_adj = 0; sens_adj < 2; sens_adj++)
        {
        SECTION("sens_adj = " + std::to_string((bool)sens_adj))
        {
            for (std::string solver : solvers)
            {
            SECTION(solver)
            {
            for (int num_stages = 4; num_stages < 5; num_stages += 1)
            {
            SECTION("num_stages = " + std::to_string(num_stages))
            {
            for (int num_steps = 1; num_steps < 20; num_steps += 2)
            {
            SECTION("num_steps = " + std::to_string(num_steps))
            {


                double tol_hess = sim_solver_tolerance_hess(solver);
                double tol = sim_solver_tolerance_sim(solver);

                plan.sim_solver = hashitsim_hess(solver);

                // create correct config based on plan
                sim_config *config = sim_config_create(plan);

            /* sim dims */
                void *dims = sim_dims_create(config);
                /* set dimensions */
                sim_dims_set(config, dims, "nx", &nx);
                sim_dims_set(config, dims, "nu", &nu);
                sim_dims_set(config, dims, "nz", &nz);

            /* sim options */

                void *opts_ = sim_opts_create(config, dims);
                sim_opts *opts = (sim_opts *) opts_;
                config->opts_initialize_default(config, dims, opts);

                opts->jac_reuse = false;        // jacobian reuse
                opts->newton_iter = 4;          // number of newton iterations per integration step

                opts->ns                = num_stages;          // number of stages in rk integrator
                opts->num_steps         = num_steps;    // number of steps
                opts->sens_forw         = (bool) sens_forw;
                opts->sens_adj          = (bool) sens_adj;
                opts->output_z          = false;
                opts->sens_algebraic    = false;
                opts->sens_hess         = true;

                if (plan.sim_solver == ERK){
                    if (!(num_stages == 4 || num_stages == 2 || num_stages == 1 )){
                        std::cout << "\n\n ATTEMPTED to use ERK with num_stages not in {1,2,4}";
                        std::cout << "\n --->> NOT SUPPORTED -- corresponding test skipped \n";
                        break;
                    }
                    if ( !sens_forw )
                    {
                        std::cout << "\n ERK hessians only work with   ";
                        std::cout << "sens_forw = true; [known issue]";
                        std::cout << "\n --->> corresponding test skipped \n";
                        break;
                    }
                    else
                    {
                        opts->num_steps *= 2;  // use more steps as explict RK has lower order
                    }
                }


            /* sim in / out */

                sim_in *in = sim_in_create(config, dims);
                sim_out *out = sim_out_create(config, dims);

                in->T = T;

            /* set model */
                switch (plan.sim_solver)
                {
                    case ERK:  // ERK
                    {
                        sim_in_set(config, dims, in, "expl_ode_fun", &expl_ode_fun);
                        sim_in_set(config, dims, in, "expl_vde_for", &expl_vde_for);
                        sim_in_set(config, dims, in, "expl_vde_adj", &expl_vde_adj);
                        sim_in_set(config, dims, in, "expl_ode_hes", &expl_ode_hess);
                        break;
                    }
                    case IRK:  // IRK
                    {
                        sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
                        sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot",
                                &impl_ode_fun_jac_x_xdot);
                        sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
                        sim_in_set(config, dims, in, "impl_ode_hes", &impl_ode_hess);
                        break;
                    }
                    default :
                    {
                        printf("\nnot enough sim solvers implemented!\n");
                        exit(1);
                    }
                }

            /* seeds */
                for (int ii = 0; ii < nx * NF; ii++)
                    in->S_forw[ii] = 0.0;
                for (int ii = 0; ii < nx; ii++)
                    in->S_forw[ii * (nx + 1)] = 1.0;

                // seeds adj
                for (int ii = 0; ii < nx; ii++)
                    in->S_adj[ii] = 1.0;
                for (int ii = nx; ii < nx + nu; ii++)
                    in->S_adj[ii] = 0.0;

            /* sim solver  */
                sim_solver = sim_solver_create(config, dims, opts);

            /* print */
                std::cout << "\n---> testing integrator " << solver;
                std::cout << " OPTS: num_steps = " << opts->num_steps;
                std::cout << ", num_stages = " << opts->ns;
                std::cout << ", jac_reuse = " << opts->jac_reuse;
                std::cout << ", newton_iter = " << opts->newton_iter << ")\n";

                for (int ii = 0; ii < nsim0; ii++)
                {
                    // x
                    for (int jj = 0; jj < nx; jj++)
                        in->x[jj] = x_sim[ii*nx+jj];

                    // u
                    for (int jj = 0; jj < nu; jj++)
                        in->u[jj] = u_sim_pendulum[ii*nu+jj];

                    acados_return = sim_solve(sim_solver, in, out);
                    REQUIRE(acados_return == 0);

                    for (int jj = 0; jj < nx; jj++){
                        x_sim[(ii+1)*nx+jj] = out->xn[jj];
                    }

                }

            /************************************************
            * compute error w.r.t. reference solution
            ************************************************/

                // error sim
                for (int jj = 0; jj < nx; jj++){
                    error[jj] = fabs(out->xn[jj] - x_ref_sol[jj]);
                    REQUIRE(std::isnan(out->xn[jj]) == false);
                }
                norm_error = onenorm(nx, 1, error);
                double rel_error_x = norm_error / norm_x_ref;

                if ( opts->sens_forw ){     // error_S_forw
                    norm_error_forw = 0.0;
                    for (int jj = 0; jj < nx*NF; jj++){
                        REQUIRE(std::isnan(out->S_forw[jj]) == 0);
                        error_S_forw[jj] = fabs(S_forw_ref_sol[jj] - out->S_forw[jj]);
                    }
                    norm_error_forw = onenorm(nx, nx + nu, error_S_forw);
                    rel_error_forw = norm_error_forw / norm_S_forw_ref;
                }


                if ( opts->sens_adj ){               // error_S_adj
                    for (int jj = 0; jj < nx + nu; jj++){
                        REQUIRE(std::isnan(out->S_adj[jj]) == 0);
                        error_S_adj[jj] = S_adj_ref_sol[jj] - out->S_adj[jj];
                    }
                    norm_error_adj = onenorm(1, nx +nu, error_S_adj);
                    rel_error_adj = norm_error_adj / norm_S_adj_ref;
                }

                if ( opts->sens_hess ){               // error_S_hess
                    for (int jj = 0; jj < (nx + nu) * (nx + nu); jj++){
                        REQUIRE(std::isnan(out->S_hess[jj]) == 0);
                        error_S_hess[jj] = S_hess_ref_sol[jj] - out->S_hess[jj];
                    }
                    norm_error_hess = onenorm(nx + nu, nx + nu, error_S_hess);
                    rel_error_hess = norm_error_hess / norm_S_hess_ref;
                }



            /************************************************
            * printing
            ************************************************/

                std::cout  << "rel_error_sim    = " << rel_error_x <<  "\n";
                if ( opts->sens_forw )
                std::cout  << "rel_error_forw   = " << rel_error_forw << "\n";
                if ( opts->sens_adj )
                std::cout  << "rel_error_adj    = " << rel_error_adj  << "\n";
                if ( opts->sens_hess ){
                    std::cout  << "rel_error_hess   = " << rel_error_hess  << "\n";
                    if ( PRINT_HESS_RESULTS ){
                        printf("Tested Hessian \n");
                        d_print_exp_mat(nx + nu, nx + nu, out->S_hess, nx + nu);
                    }
                }
                // printf("Tested xn \n");
                // d_print_exp_mat(1, nx, out->xn, 1);
            /************************************************
            * asserts on erors
            ************************************************/
                REQUIRE(rel_error_x <= tol);

                if ( opts->sens_forw )
                    REQUIRE(rel_error_forw <= tol);

                if ( opts->sens_adj )
                    REQUIRE(rel_error_adj <= tol);

                if ( opts->sens_hess )
                    REQUIRE(rel_error_hess <= tol_hess);

            /************************************************
            * free tested solver
            ************************************************/
                free(config);
                free(dims);
                free(opts);

                free(in);
                free(out);
                free(sim_solver);
            }  // end SECTION
            }  // end for
            }  // end SECTION
            }  // end for
            }  // end SECTION
            }  // end for num_steps
            }  // end SECTION
            }  // end for num_stages
        }  // end section solver
    }  // END FOR SOLVERS

    // implicit model
    external_function_casadi_free(&impl_ode_fun);
    external_function_casadi_free(&impl_ode_fun_jac_x_xdot);
    external_function_casadi_free(&impl_ode_fun_jac_x_xdot_u);
    external_function_casadi_free(&impl_ode_jac_x_xdot_u);
    external_function_casadi_free(&impl_ode_hess);
    // explicit model
    external_function_casadi_free(&expl_ode_fun);
    external_function_casadi_free(&expl_vde_for);
    external_function_casadi_free(&expl_vde_adj);
    external_function_casadi_free(&expl_ode_hess);

}  // END_TEST_CASE


TEST_CASE("pendulum model hessians - Finite Differences", "compare against finite differences"){


    for (int ii = 0; ii < nx; ii++)
        x0_pendulum[ii] = 0.0;

    u_sim_pendulum[0] = 0.1;

    int NF = nx + nu;  // columns of forward seed

    double T = 0.1;  // simulation time

    double FD_EPS = 1e-8;
    double hess_FD[(nx+nu)*(nx+nu)];

/************************************************
* external functions
************************************************/
    /* IMPLICIT MODEL */
    // impl_ode_fun
    external_function_casadi impl_ode_fun;
    impl_ode_fun.casadi_fun = &pendulum_ode_impl_ode_fun;
    impl_ode_fun.casadi_work = &pendulum_ode_impl_ode_fun_work;
    impl_ode_fun.casadi_sparsity_in = &pendulum_ode_impl_ode_fun_sparsity_in;
    impl_ode_fun.casadi_sparsity_out = &pendulum_ode_impl_ode_fun_sparsity_out;
    impl_ode_fun.casadi_n_in = &pendulum_ode_impl_ode_fun_n_in;
    impl_ode_fun.casadi_n_out = &pendulum_ode_impl_ode_fun_n_out;
    external_function_casadi_create(&impl_ode_fun);

    // impl_ode_fun_jac_x_xdot
    external_function_casadi impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_fun = &pendulum_ode_impl_ode_fun_jac_x_xdot_z;
    impl_ode_fun_jac_x_xdot.casadi_work = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_work;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_sparsity_in;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_out =
                         &pendulum_ode_impl_ode_fun_jac_x_xdot_z_sparsity_out;
    impl_ode_fun_jac_x_xdot.casadi_n_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_n_in;
    impl_ode_fun_jac_x_xdot.casadi_n_out = &pendulum_ode_impl_ode_fun_jac_x_xdot_z_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_fun = &pendulum_ode_impl_ode_jac_x_xdot_u_z;
    impl_ode_jac_x_xdot_u.casadi_work = &pendulum_ode_impl_ode_jac_x_xdot_u_z_work;
    impl_ode_jac_x_xdot_u.casadi_sparsity_in = &pendulum_ode_impl_ode_jac_x_xdot_u_z_sparsity_in;
    impl_ode_jac_x_xdot_u.casadi_sparsity_out = &pendulum_ode_impl_ode_jac_x_xdot_u_z_sparsity_out;
    impl_ode_jac_x_xdot_u.casadi_n_in = &pendulum_ode_impl_ode_jac_x_xdot_u_z_n_in;
    impl_ode_jac_x_xdot_u.casadi_n_out = &pendulum_ode_impl_ode_jac_x_xdot_u_z_n_out;
    external_function_casadi_create(&impl_ode_jac_x_xdot_u);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_fun = &pendulum_ode_impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_work = &pendulum_ode_impl_ode_fun_jac_x_xdot_u_work;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_in =
                            &pendulum_ode_impl_ode_fun_jac_x_xdot_u_sparsity_in;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_out =
                            &pendulum_ode_impl_ode_fun_jac_x_xdot_u_sparsity_out;
    impl_ode_fun_jac_x_xdot_u.casadi_n_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_u_n_in;
    impl_ode_fun_jac_x_xdot_u.casadi_n_out = &pendulum_ode_impl_ode_fun_jac_x_xdot_u_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot_u);

    // impl_ode_hess
    external_function_casadi impl_ode_hess;
    impl_ode_hess.casadi_fun = &pendulum_ode_impl_ode_hess;
    impl_ode_hess.casadi_work = &pendulum_ode_impl_ode_hess_work;
    impl_ode_hess.casadi_sparsity_in = &pendulum_ode_impl_ode_hess_sparsity_in;
    impl_ode_hess.casadi_sparsity_out = &pendulum_ode_impl_ode_hess_sparsity_out;
    impl_ode_hess.casadi_n_in = &pendulum_ode_impl_ode_hess_n_in;
    impl_ode_hess.casadi_n_out = &pendulum_ode_impl_ode_hess_n_out;
    external_function_casadi_create(&impl_ode_hess);

/* generate adjoint sensitivities */
    sim_solver_plan_t plan;
    plan.sim_solver = IRK;  // IRK

    sim_config *config = sim_config_create(plan);

    void *dims = sim_dims_create(config);

    /* set dimensions */
    sim_dims_set(config, dims, "nx", &nx);
    sim_dims_set(config, dims, "nu", &nu);
    sim_dims_set(config, dims, "nz", &nz);

    // set opts
    void *opts_ = sim_opts_create(config, dims);
    sim_opts *opts = (sim_opts *) opts_;
    config->opts_initialize_default(config, dims, opts);

    // opts reference solution
    opts->sens_adj  = true;
    opts->jac_reuse = false;  // jacobian reuse
    opts->newton_iter = 8;  // number of newton iterations per integration step
    opts->num_steps = 100;  // number of steps
    opts->ns = 7;  // number of stages in rk integrator

    sim_in *in = sim_in_create(config, dims);
    sim_out *out = sim_out_create(config, dims);

    in->T = T;

    // set model
    switch (plan.sim_solver)
    {
        case IRK:  // IRK
        {
            sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
            sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot",
                    &impl_ode_fun_jac_x_xdot);
            sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
            sim_in_set(config, dims, in, "impl_ode_hess", &impl_ode_hess);
            break;
        }
        default :
        {
            printf("\nnot plan.sim_solver not supported!\n");
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
    for (int ii = nx; ii < nx + nu; ii++)
        in->S_adj[ii] = 0.0;

    /************************************************
    * sim solver
    ************************************************/

    sim_solver *sim_solver = sim_solver_create(config, dims, opts);

    // x
    for (int jj = 0; jj < nx; jj++)
        in->x[jj] = x0_pendulum[jj];

    // u
    for (int jj = 0; jj < nu; jj++)
        in->u[jj] = u_sim_pendulum[jj];

    acados_return = sim_solve(sim_solver, in, out);
    REQUIRE(acados_return == 0);

    // store reference solution
    for (int jj = 0; jj < NF; jj++)
        S_adj_ref_sol[jj] = out->S_adj[jj];

    /************************************************
    * finite differences loop
    ************************************************/

    // hessian FD:
    for (int s = 0; s < nx + nu; s++)
    {
        // set input
        for (int jj = 0; jj < nx; jj++)
            in->x[jj] = x0_pendulum[jj];
        for (int jj = 0; jj < nu; jj++)
            in->u[jj] = u_sim_pendulum[jj];

        if (s < nx)
        {
            in->x[s] = in->x[s] + FD_EPS;
        }
        else
        {
            in->u[s - nx] = in->u[s - nx] + FD_EPS;
        }

        acados_return = sim_solve(sim_solver, in, out);

        for (int i = 0; i < nx + nu; i++)
            hess_FD[s * (nx + nu) + i] = (out->S_adj[i] - S_adj_ref_sol[i]) / FD_EPS;
    }

    // if ( PRINT_HESS_RESULTS )
    printf("\n================================================================\n");
    printf("Finite differences Hessian result =\n");
    d_print_exp_mat(nx+nu, nx+nu, hess_FD, nx+nu);

    printf("IRK IND Hessian result =\n");
    d_print_exp_mat(nx+nu, nx+nu, S_hess_ref_sol, nx+nu);

    printf("ERROR Hessian (Finite differences vs IRK internal) =\n");
    for (int jj = 0; jj < (nx + nu) * (nx + nu); jj++){
        REQUIRE(std::isnan(out->S_hess[jj]) == 0);
        error_S_hess[jj] = S_hess_ref_sol[jj] - hess_FD[jj];
    }
    d_print_exp_mat(nx+nu, nx+nu, error_S_hess, nx+nu);

    norm_error_hess = onenorm(nx + nu, nx + nu, error_S_hess);
    rel_error_hess = norm_error_hess / norm_S_hess_ref;

    printf("relative errror hessian (Finite differences vs IRK internal) = %e\n", rel_error_hess);

    REQUIRE(rel_error_hess < 1e-4);

    // implicit model
    external_function_casadi_free(&impl_ode_fun);
    external_function_casadi_free(&impl_ode_fun_jac_x_xdot);
    external_function_casadi_free(&impl_ode_fun_jac_x_xdot_u);
    external_function_casadi_free(&impl_ode_jac_x_xdot_u);
    external_function_casadi_free(&impl_ode_hess);

    // free solver
    free(config);
    free(dims);
    free(opts);

    free(in);
    free(out);
    free(sim_solver);
}
