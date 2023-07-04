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

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// acados
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/math.h"

#include "acados_c/external_function_interface.h"
#include "interfaces/acados_c/sim_interface.h"

// crane dae model
#include "examples/c/crane_dae_model/crane_dae_model.h"


extern "C"
{

}

using std::vector;

sim_solver_t hashitsim_dae(std::string const& inString)
{
    if (inString == "IRK") return IRK;
    if (inString == "GNSF") return GNSF;

    return (sim_solver_t) -1;
}

double sim_solver_tolerance_dae(std::string const& inString)
{
    if (inString == "IRK")  return 1e-7;
    if (inString == "GNSF") return 1e-7;

    return -1;
}

double sim_solver_tolerance_algebraic_dae(std::string const& inString)
{
    if (inString == "IRK")  return 1e-3;
    if (inString == "GNSF") return 1e-3;

    return -1;
}



TEST_CASE("crane_dae_example", "[integrators]")
{
    vector<std::string> solvers = {"IRK", "GNSF"};
    // initialize dimensions

    int nx = 9;
    int nu = 2;
    int nz = 2;
    int nout = 1;    // gnsf split
    int ny = 4;
    int nuhat = 1;

    int nx1 = 5;
    int nz1 = 0;

    // generate x0, u_sim
    double x0[nx];
    double u_sim[nu];

    for (int ii = 0; ii < nx; ii++) {
        x0[ii] = 0.0;
    }
    // x0[2] = 0.8; // old (manual) gnsf split
    x0[0] = 0.8;  // value for xL

    u_sim[0] = 40.108149413030752;
    u_sim[1] = -50.446662212534974;

    int NF = nx + nu;  // columns of forward seed

    int nsim0 = 1;  // nsim;

    double T = 0.01;  // simulation time // former values .1, .05
    // reduced for faster test

    double x_sim[nx*(nsim0+2)];

    double x_ref_sol[nx];
    double S_forw_ref_sol[nx*NF];
    double S_adj_ref_sol[NF];
    double z_ref_sol[nz];
    double S_alg_ref_sol[NF*nz];

    double error[nx];
    double error_z[nz];
    double error_S_forw[nx*NF];
    double error_S_adj[NF];
    double error_S_alg[NF*nz];

    double norm_error, norm_error_forw, norm_error_adj, norm_error_z, norm_error_sens_alg;

    for (int ii = 0; ii < nx; ii++)
        x_sim[ii] = x0[ii];

/************************************************
* external functions
************************************************/

    // impl_ode_fun
    external_function_casadi impl_ode_fun;
    impl_ode_fun.casadi_fun = &crane_dae_impl_ode_fun;
    impl_ode_fun.casadi_work = &crane_dae_impl_ode_fun_work;
    impl_ode_fun.casadi_sparsity_in = &crane_dae_impl_ode_fun_sparsity_in;
    impl_ode_fun.casadi_sparsity_out = &crane_dae_impl_ode_fun_sparsity_out;
    impl_ode_fun.casadi_n_in = &crane_dae_impl_ode_fun_n_in;
    impl_ode_fun.casadi_n_out = &crane_dae_impl_ode_fun_n_out;
    external_function_casadi_create(&impl_ode_fun);

    // impl_ode_fun_jac_x_xdot
    external_function_casadi impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_fun = &crane_dae_impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_work = &crane_dae_impl_ode_fun_jac_x_xdot_work;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &crane_dae_impl_ode_fun_jac_x_xdot_sparsity_in;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &crane_dae_impl_ode_fun_jac_x_xdot_sparsity_out;
    impl_ode_fun_jac_x_xdot.casadi_n_in = &crane_dae_impl_ode_fun_jac_x_xdot_n_in;
    impl_ode_fun_jac_x_xdot.casadi_n_out = &crane_dae_impl_ode_fun_jac_x_xdot_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_fun = &crane_dae_impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_work = &crane_dae_impl_ode_jac_x_xdot_u_work;
    impl_ode_jac_x_xdot_u.casadi_sparsity_in = &crane_dae_impl_ode_jac_x_xdot_u_sparsity_in;
    impl_ode_jac_x_xdot_u.casadi_sparsity_out = &crane_dae_impl_ode_jac_x_xdot_u_sparsity_out;
    impl_ode_jac_x_xdot_u.casadi_n_in = &crane_dae_impl_ode_jac_x_xdot_u_n_in;
    impl_ode_jac_x_xdot_u.casadi_n_out = &crane_dae_impl_ode_jac_x_xdot_u_n_out;
    external_function_casadi_create(&impl_ode_jac_x_xdot_u);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_fun = &crane_dae_impl_ode_fun_jac_x_xdot_u;
    impl_ode_fun_jac_x_xdot_u.casadi_work = &crane_dae_impl_ode_fun_jac_x_xdot_u_work;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_in =
                            &crane_dae_impl_ode_fun_jac_x_xdot_u_sparsity_in;
    impl_ode_fun_jac_x_xdot_u.casadi_sparsity_out =
                            &crane_dae_impl_ode_fun_jac_x_xdot_u_sparsity_out;
    impl_ode_fun_jac_x_xdot_u.casadi_n_in = &crane_dae_impl_ode_fun_jac_x_xdot_u_n_in;
    impl_ode_fun_jac_x_xdot_u.casadi_n_out = &crane_dae_impl_ode_fun_jac_x_xdot_u_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot_u);

    /************************************************
    * external functions (Generalized Nonlinear Static Feedback (GNSF) model)
    ************************************************/
    // phi_fun
    external_function_casadi phi_fun;
    phi_fun.casadi_fun            = &crane_dae_phi_fun;
    phi_fun.casadi_work           = &crane_dae_phi_fun_work;
    phi_fun.casadi_sparsity_in    = &crane_dae_phi_fun_sparsity_in;
    phi_fun.casadi_sparsity_out   = &crane_dae_phi_fun_sparsity_out;
    phi_fun.casadi_n_in           = &crane_dae_phi_fun_n_in;
    phi_fun.casadi_n_out          = &crane_dae_phi_fun_n_out;
    external_function_casadi_create(&phi_fun);

    // phi_fun_jac_y
    external_function_casadi phi_fun_jac_y;
    phi_fun_jac_y.casadi_fun            = &crane_dae_phi_fun_jac_y;
    phi_fun_jac_y.casadi_work           = &crane_dae_phi_fun_jac_y_work;
    phi_fun_jac_y.casadi_sparsity_in    = &crane_dae_phi_fun_jac_y_sparsity_in;
    phi_fun_jac_y.casadi_sparsity_out   = &crane_dae_phi_fun_jac_y_sparsity_out;
    phi_fun_jac_y.casadi_n_in           = &crane_dae_phi_fun_jac_y_n_in;
    phi_fun_jac_y.casadi_n_out          = &crane_dae_phi_fun_jac_y_n_out;
    external_function_casadi_create(&phi_fun_jac_y);

    // phi_jac_y_uhat
    external_function_casadi phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_fun                = &crane_dae_phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_work               = &crane_dae_phi_jac_y_uhat_work;
    phi_jac_y_uhat.casadi_sparsity_in        = &crane_dae_phi_jac_y_uhat_sparsity_in;
    phi_jac_y_uhat.casadi_sparsity_out       = &crane_dae_phi_jac_y_uhat_sparsity_out;
    phi_jac_y_uhat.casadi_n_in               = &crane_dae_phi_jac_y_uhat_n_in;
    phi_jac_y_uhat.casadi_n_out              = &crane_dae_phi_jac_y_uhat_n_out;
    external_function_casadi_create(&phi_jac_y_uhat);

    // f_lo_fun_jac_x1k1uz
    external_function_casadi f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_fun            = &crane_dae_f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_work           = &crane_dae_f_lo_fun_jac_x1k1uz_work;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_in    = &crane_dae_f_lo_fun_jac_x1k1uz_sparsity_in;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_out   = &crane_dae_f_lo_fun_jac_x1k1uz_sparsity_out;
    f_lo_fun_jac_x1k1uz.casadi_n_in           = &crane_dae_f_lo_fun_jac_x1k1uz_n_in;
    f_lo_fun_jac_x1k1uz.casadi_n_out          = &crane_dae_f_lo_fun_jac_x1k1uz_n_out;
    external_function_casadi_create(&f_lo_fun_jac_x1k1uz);

    // get_matrices_fun
    external_function_casadi get_matrices_fun;
    get_matrices_fun.casadi_fun            = &crane_dae_get_matrices_fun;
    get_matrices_fun.casadi_work           = &crane_dae_get_matrices_fun_work;
    get_matrices_fun.casadi_sparsity_in    = &crane_dae_get_matrices_fun_sparsity_in;
    get_matrices_fun.casadi_sparsity_out   = &crane_dae_get_matrices_fun_sparsity_out;
    get_matrices_fun.casadi_n_in           = &crane_dae_get_matrices_fun_n_in;
    get_matrices_fun.casadi_n_out          = &crane_dae_get_matrices_fun_n_out;
    external_function_casadi_create(&get_matrices_fun);


/************************************************
* Create Reference Solution
************************************************/

    sim_solver_plan_t plan;
    plan.sim_solver = GNSF;  // IRK; -- works but super slow

    sim_config *config = sim_config_create(plan);

    void *dims = sim_dims_create(config);

    /* set dimensions */
    sim_dims_set(config, dims, "nx", &nx);
    sim_dims_set(config, dims, "nu", &nu);
    sim_dims_set(config, dims, "nz", &nz);

    // GNSF -- set additional dimensions
    if (plan.sim_solver == GNSF)
    {
        sim_dims_set(config, dims, "nx1", &nx1);
        sim_dims_set(config, dims, "nz1", &nz1);
        sim_dims_set(config, dims, "nout", &nout);
        sim_dims_set(config, dims, "ny", &ny);
        sim_dims_set(config, dims, "nuhat", &nuhat);
    }

    // set opts
    void *opts_ = sim_opts_create(config, dims);
    sim_opts *opts = (sim_opts *) opts_;
    config->opts_initialize_default(config, dims, opts);

    // opts reference solution
    opts->sens_forw = true;
    opts->sens_adj = true;
    opts->sens_algebraic = true;
    opts->sens_hess = false;
    opts->output_z = true;
    opts->jac_reuse = false;  // jacobian reuse
    opts->newton_iter = 8;  // number of newton iterations per integration step
    opts->num_steps = 100;  // number of steps
    opts->ns = 4;  // number of stages in rk integrator

    sim_in *in = sim_in_create(config, dims);
    sim_out *out = sim_out_create(config, dims);

    sim_in_set(config, dims, in, "T", &T);

    // set model
    switch (plan.sim_solver)
    {
        case IRK:  // IRK
        {
            sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
            sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot",
                    &impl_ode_fun_jac_x_xdot);
            sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
            break;
        }
        case GNSF:  // GNSF
        {
            // set model funtions
            sim_in_set(config, dims, in, "phi_fun", &phi_fun);
            sim_in_set(config, dims, in, "phi_fun_jac_y", &phi_fun_jac_y);
            sim_in_set(config, dims, in, "phi_jac_y_uhat", &phi_jac_y_uhat);
            sim_in_set(config, dims, in, "f_lo_jac_x1_x1dot_u_z", &f_lo_fun_jac_x1k1uz);
            sim_in_set(config, dims, in, "get_gnsf_matrices", &get_matrices_fun);
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
    sim_precompute(sim_solver, in, out);

    int acados_return;

    for (int ii = 0; ii < nsim0; ii++)
    {
        // x
        for (int jj = 0; jj < nx; jj++)
            in->x[jj] = x_sim[ii*nx+jj];

        // u
        for (int jj = 0; jj < nu; jj++)
            in->u[jj] = u_sim[ii*nu+jj];

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

    for (int jj = 0; jj < nz; jj++)
        z_ref_sol[jj] = out->zn[jj];

    for (int jj = 0; jj < nz*NF; jj++)
        S_alg_ref_sol[jj] = out->S_algebraic[jj];

    // compute one norms
    double norm_x_ref, norm_S_forw_ref, norm_S_adj_ref, norm_z_ref, norm_S_alg_ref = 0;

    norm_x_ref = onenorm(nx, 1, x_ref_sol);
    norm_S_forw_ref = onenorm(nx, nx + nu, S_forw_ref_sol);
    norm_S_adj_ref = onenorm(1, nx + nu, S_adj_ref_sol);
    norm_z_ref = onenorm(nz, 1, z_ref_sol);
    norm_S_alg_ref = onenorm(nz, nx + nu, S_alg_ref_sol);

    // printf("Reference xn \n");
    // d_print_exp_mat(1, nx, &x_ref_sol[0], 1);

    // printf("Reference zn \n");
    // d_print_exp_mat(1, nz, &z_ref_sol[0], 1);

    // printf("Reference forward sensitivities \n");
    // d_print_exp_mat(nx, NF, &S_forw_ref_sol[0], nx);

    // printf("reference algebraic sensitivities \n");
    // d_print_exp_mat(nz, nx + nu, &S_alg_ref_sol[0], nz);

    /* free */
    sim_config_destroy(config);
    sim_dims_destroy(dims);
    sim_opts_destroy(opts);

    sim_in_destroy(in);
    sim_out_destroy(out);
    sim_solver_destroy(sim_solver);

/************************************************
* test solver loop
************************************************/

    for (std::string solver : solvers)
    {
    SECTION(solver)
    {
        for (int sens_adj = 0; sens_adj < 2; sens_adj++)
        {
        SECTION("sens_adj = " + std::to_string((bool)sens_adj))
        {
            for (int output_z = 0; output_z < 2; output_z++)
            {
            SECTION("output_z = " + std::to_string((bool)output_z))
            {
            for (int sens_alg = 0; sens_alg < 2; sens_alg++)
            {
            SECTION("sens_alg = " + std::to_string((bool)sens_alg))
            {
            for (int num_stages = 3; num_stages < 4; num_stages++)
            {
            SECTION("num_stages = " + std::to_string(num_stages))
            {
            for (int num_steps = 3; num_steps < 5; num_steps += 2)
            {
            SECTION("num_steps = " + std::to_string(num_steps))
            {
            // for (int sens_forw = 0; sens_forw < 2; sens_forw++)
            // {
            // SECTION("sens_forw = " + std::to_string((bool)sens_forw))
            // {
                int sens_forw = 1;



                double tol = sim_solver_tolerance_dae(solver);
                double tol_algebraic = sim_solver_tolerance_algebraic_dae(solver);

                plan.sim_solver = hashitsim_dae(solver);

                // create correct config based on plan
                sim_config *config = sim_config_create(plan);
                void *dims = sim_dims_create(config);

            /* sim dims */
                sim_dims_set(config, dims, "nx", &nx);
                sim_dims_set(config, dims, "nu", &nu);
                sim_dims_set(config, dims, "nz", &nz);

                // GNSF -- set additional dimensions
                if (plan.sim_solver == GNSF)
                {
                    sim_dims_set(config, dims, "nx1", &nx1);
                    sim_dims_set(config, dims, "nz1", &nz1);
                    sim_dims_set(config, dims, "nout", &nout);
                    sim_dims_set(config, dims, "ny", &ny);
                    sim_dims_set(config, dims, "nuhat", &nuhat);
                }

            /* sim options */

                void *opts_ = sim_opts_create(config, dims);
                sim_opts *opts = (sim_opts *) opts_;
                config->opts_initialize_default(config, dims, opts);

                opts->jac_reuse = false;        // jacobian reuse
                opts->newton_iter = 2;          // number of newton iterations per integration step

                opts->ns                = num_stages;          // number of stages in rk integrator
                opts->num_steps         = num_steps;    // number of steps
                opts->sens_forw         = (bool) sens_forw;
                opts->sens_adj          = (bool) sens_adj;
                opts->output_z          = (bool) output_z;
                opts->sens_algebraic    = (bool) sens_alg;
                opts->sens_hess         = false;


            /* sim in / out */

                sim_in *in = sim_in_create(config, dims);
                sim_out *out = sim_out_create(config, dims);

                in->T = T;

            /* set model */
                switch (plan.sim_solver)
                {
                    case IRK:  // IRK
                    {
                        sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
                        sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot);
                        sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
                        break;
                    }
                    case GNSF:  // GNSF
                    {
                        // set model funtions
                        sim_in_set(config, dims, in, "phi_fun", &phi_fun);
                        sim_in_set(config, dims, in, "phi_fun_jac_y", &phi_fun_jac_y);
                        sim_in_set(config, dims, in, "phi_jac_y_uhat", &phi_jac_y_uhat);
                        sim_in_set(config, dims, in, "f_lo_jac_x1_x1dot_u_z", &f_lo_fun_jac_x1k1uz);
                        sim_in_set(config, dims, in, "get_gnsf_matrices", &get_matrices_fun);
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
                int acados_return;
                sim_precompute(sim_solver, in, out);

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
                        in->u[jj] = u_sim[ii*nu+jj];

                    acados_return = sim_solve(sim_solver, in, out);
                    REQUIRE(acados_return == 0);

                    for (int jj = 0; jj < nx; jj++){
                        x_sim[(ii+1)*nx+jj] = out->xn[jj];
                    }

                }

            /************************************************
            * compute error w.r.t. reference solution
            ************************************************/
                double rel_error_forw, rel_error_adj, rel_error_z, rel_error_sens_alg;

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

                if ( opts->output_z ){      // error_z
                    for (int jj = 0; jj < nz; jj++){
                        error_z[jj] = fabs(out->zn[jj] - z_ref_sol[jj]);
                        REQUIRE(std::isnan(out->zn[jj]) == 0);
                    }
                    norm_error_z = onenorm(nz, 1, error_z);
                    rel_error_z = norm_error_z / norm_z_ref;
                }

                if ( opts->sens_algebraic ){        // error_S_alg
                    for (int jj = 0; jj < nz * (nx + nu); jj++){
                        // REQUIRE(std::isnan(out->S_algebraic[jj]) == 0);
                        error_S_alg[jj] = fabs(out->S_algebraic[jj] - S_alg_ref_sol[jj]);
                    }
                    norm_error_sens_alg = onenorm(nz, nx + nu, error_S_alg);
                    rel_error_sens_alg = norm_error_sens_alg / norm_S_alg_ref;
                }




            /************************************************
            * printing
            ************************************************/

                std::cout  << "rel_error_sim      = " << rel_error_x <<  "\n";
                if ( opts->sens_forw )
                std::cout  << "rel_error_forw     = " << rel_error_forw << "\n";
                if ( opts->sens_adj )
                std::cout  << "rel_error_adj      = " << rel_error_adj  << "\n";
                if ( opts->output_z )
                std::cout  << "rel_error_z        = " << rel_error_z <<"\n";
                if ( opts->sens_algebraic )
                std::cout  << "rel_error_sens_alg = " << rel_error_sens_alg <<"\n";

                // printf("tested algebraic sensitivities \n");
                // d_print_exp_mat(nz, nx + nu, &S_alg_ref_sol[0], nz);
                // printf("tested xn \n");
                // d_print_exp_mat(1, nx, &out->xn[0], 1);
            /************************************************
            * asserts on erors
            ************************************************/
                REQUIRE(rel_error_x <= tol);

                if ( opts->sens_forw ){
                    // printf("tested forward sensitivities \n");
                    // d_print_exp_mat(nx, nx + nu, &out->S_forw[0], nx);
                    // printf("REF forward sensitivities \n");
                    // d_print_exp_mat(nx, nx + nu, S_forw_ref_sol, nx);
                    // printf("ERROR forward sensitivities \n");
                    // d_print_exp_mat(nx, nx + nu, error_S_forw, nx);
                    REQUIRE(rel_error_forw <= tol);
                }

                if ( opts->sens_adj ){
                    // printf("tested adjoint sensitivities \n");
                    // d_print_exp_mat(1, nx + nu, &out->S_adj[0], 1);
                    // printf("REF adjoint sensitivities \n");
                    // d_print_exp_mat(1, nx + nu, S_adj_ref_sol, 1);
                    REQUIRE(rel_error_adj <= tol);
                }

                if ( opts->output_z ){
                    // printf("tested z output \n");
                    // d_print_exp_mat(1, nz, &out->zn[0], 1);
                    // printf("REF z output \n");
                    // d_print_exp_mat(1, nz, z_ref_sol, 1);
                    REQUIRE(rel_error_z <= tol_algebraic);
                }

                if ( opts->sens_algebraic ){
                    // printf("tested algebraic sensitivities \n");
                    // d_print_exp_mat(nz, nx + nu, &out->S_algebraic[0], nz);
                    // printf("reference algebraic sensitivities \n");
                    // d_print_exp_mat(nz, nx + nu, &S_alg_ref_sol[0], nz);
                    REQUIRE(rel_error_sens_alg <= tol_algebraic);
                }

            /************************************************
            * free tested solver
            ************************************************/
                sim_config_destroy(config);
                sim_dims_destroy(dims);
                sim_opts_destroy(opts);

                sim_in_destroy(in);
                sim_out_destroy(out);
                sim_solver_destroy(sim_solver);
            // }  // end SECTION
            // }  // end for
            }  // end SECTION
            }  // end for
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
    // gnsf functions:
    external_function_casadi_free(&phi_fun);
    external_function_casadi_free(&phi_fun_jac_y);
    external_function_casadi_free(&phi_jac_y_uhat);
    external_function_casadi_free(&f_lo_fun_jac_x1k1uz);
    external_function_casadi_free(&get_matrices_fun);

}  // END_TEST_CASE
