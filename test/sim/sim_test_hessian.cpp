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
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"

// acados
#include "acados/sim/sim_common.h"
// #include "acados/sim/sim_gnsf.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/math.h"

#include "acados_c/external_function_interface.h"
#include "interfaces/acados_c/sim_interface.h"

#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// pendulum_model
#include "examples/c/pendulum_model/pendulum_model.h"


extern "C"
{

}

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

sim_solver_t hashitsim_hess(std::string const& inString)
{
    if (inString == "ERK") return ERK;
    if (inString == "IRK") return IRK;
    if (inString == "LIFTED_IRK") return LIFTED_IRK;
    if (inString == "GNSF") return GNSF;
    if (inString == "NEW_LIFTED_IRK") return NEW_LIFTED_IRK;

    return (sim_solver_t) -1;
}

double sim_solver_tolerance_hess(std::string const& inString)
{
    if (inString == "IRK")  return 1e-4;
    // if (inString == "LIFTED_IRK") return 1e-3;
    if (inString == "ERK") return 1e-4;
    // if (inString == "NEW_LIFTED_IRK") return 1e-3;

    return -1;
}


TEST_CASE("pendulum_hessians", "[integrators]")
{
    vector<std::string> solvers = {"IRK", "ERK"};
    // {"ERK", "IRK", "LIFTED_IRK", "GNSF", "NEW_LIFTED_IRK"};
    // initialize dimensions

    const int nx = 4;
    const int nu = 1;
    const int nz = 0;
    // const int nx1 = 5;  // gnsf split
    // const int nx2 = 4;
    // const int n_out = 3;
    // const int ny = 5;
    // const int nuhat = 1;

    // generate x0, u_sim
    double x0[nx];
    double u_sim[nu];

    for (int ii = 0; ii < nx; ii++)
        x0[ii] = 0.0;

    u_sim[0] = 0.1;

    int NF = nx + nu;  // columns of forward seed

    int nsim0 = 1;  // nsim;

    double T = 0.5;  // simulation time
    // reduced for faster test

    double x_sim[nx*(nsim0+2)];

    double x_ref_sol[nx];
    double S_forw_ref_sol[nx*NF];
    double S_adj_ref_sol[NF];

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
    impl_ode_fun_jac_x_xdot.casadi_fun = &pendulum_ode_impl_ode_fun_jac_x_xdot;
    impl_ode_fun_jac_x_xdot.casadi_work = &pendulum_ode_impl_ode_fun_jac_x_xdot_work;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_sparsity_in;
    impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &pendulum_ode_impl_ode_fun_jac_x_xdot_sparsity_out;
    impl_ode_fun_jac_x_xdot.casadi_n_in = &pendulum_ode_impl_ode_fun_jac_x_xdot_n_in;
    impl_ode_fun_jac_x_xdot.casadi_n_out = &pendulum_ode_impl_ode_fun_jac_x_xdot_n_out;
    external_function_casadi_create(&impl_ode_fun_jac_x_xdot);

    // impl_ode_jac_x_xdot_u
    external_function_casadi impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_fun = &pendulum_ode_impl_ode_jac_x_xdot_u;
    impl_ode_jac_x_xdot_u.casadi_work = &pendulum_ode_impl_ode_jac_x_xdot_u_work;
    impl_ode_jac_x_xdot_u.casadi_sparsity_in = &pendulum_ode_impl_ode_jac_x_xdot_u_sparsity_in;
    impl_ode_jac_x_xdot_u.casadi_sparsity_out = &pendulum_ode_impl_ode_jac_x_xdot_u_sparsity_out;
    impl_ode_jac_x_xdot_u.casadi_n_in = &pendulum_ode_impl_ode_jac_x_xdot_u_n_in;
    impl_ode_jac_x_xdot_u.casadi_n_out = &pendulum_ode_impl_ode_jac_x_xdot_u_n_out;
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

    // expl_ode_jac
    external_function_casadi expl_ode_jac;
    expl_ode_jac.casadi_fun = &pendulum_ode_expl_ode_jac;
    expl_ode_jac.casadi_work = &pendulum_ode_expl_ode_jac_work;
    expl_ode_jac.casadi_sparsity_in = &pendulum_ode_expl_ode_jac_sparsity_in;
    expl_ode_jac.casadi_sparsity_out = &pendulum_ode_expl_ode_jac_sparsity_out;
    expl_ode_jac.casadi_n_in = &pendulum_ode_expl_ode_jac_n_in;
    expl_ode_jac.casadi_n_out = &pendulum_ode_expl_ode_jac_n_out;
    external_function_casadi_create(&expl_ode_jac);

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

    sim_solver_plan plan;
    plan.sim_solver = ERK;  // IRK

    sim_solver_config *config = sim_config_create(plan);

    void *dims = sim_dims_create(config);

    /* set dimensions */
    config->set_nx(dims, nx);
    config->set_nu(dims, nu);
    config->set_nz(dims, nz);

    // set opts
    void *opts_ = sim_opts_create(config, dims);
    sim_rk_opts *opts = (sim_rk_opts *) opts_;
    config->opts_initialize_default(config, dims, opts);

    // opts reference solution
    opts->sens_forw = true;
    opts->sens_adj  = true;
    opts->sens_algebraic = false;
    opts->sens_hess = true;
    opts->output_z = false;
    opts->jac_reuse = false;  // jacobian reuse
    opts->newton_iter = 8;  // number of newton iterations per integration step
    opts->num_steps = 40;  // number of steps
    opts->ns = 4;  // number of stages in rk integrator

    sim_in *in = sim_in_create(config, dims);
    sim_out *out = sim_out_create(config, dims);

    in->T = T;

    // // import model matrices
    // external_function_generic *get_model_matrices =
    //         (external_function_generic *) &get_matrices_fun;
    // gnsf_model *model = (gnsf_model *) in->model;

    // set model
    switch (plan.sim_solver)
    {
        case ERK:  // ERK
        {
            sim_set_model(config, in, "expl_ode_fun", &expl_ode_fun);
            sim_set_model(config, in, "expl_vde_for", &expl_vde_for);
            sim_set_model(config, in, "expl_vde_adj", &expl_vde_adj);
            sim_set_model(config, in, "expl_ode_hess", &expl_ode_hess);
            break;
        }
        case IRK:  // IRK
        {
            sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
            sim_set_model(config, in, "impl_ode_fun_jac_x_xdot",
                    &impl_ode_fun_jac_x_xdot);
            sim_set_model(config, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
            sim_set_model(config, in, "impl_ode_hess", &impl_ode_hess);
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

    sim_solver *sim_solver = sim_create(config, dims, opts);

    // if (plan.sim_solver == GNSF){  // for gnsf: perform precomputation
    //     gnsf_model *model = (gnsf_model *) in->model;
    //     sim_gnsf_precompute(config, gnsf_dim, model, opts,
    //                 sim_solver->mem, sim_solver->work, in->T);
    // }

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

    // compute one norms
    double norm_x_ref, norm_S_forw_ref, norm_S_adj_ref = 0;

    norm_x_ref = onenorm(nx, 1, x_ref_sol);
    norm_S_forw_ref = onenorm(nx, nx + nu, S_forw_ref_sol);
    norm_S_adj_ref = onenorm(1, nx + nu, S_adj_ref_sol);

    // printf("Reference xn \n");
    // d_print_e_mat(1, nx, &x_ref_sol[0], 1);

    // printf("Reference forward sensitivities \n");
    // d_print_e_mat(nx, NF, &S_forw_ref_sol[0], nx);

    printf("reference adjoint sensitivities \n");
    d_print_e_mat(1, nx + nu, &S_adj_ref_sol[0], 1);

    printf("Reference Hessian \n");
    double *S_hess_out;
    if(opts->sens_hess)
    {
        double zero = 0.0;
        S_hess_out = out->S_hess;
        printf("\nS_hess_out: \n");
        for (int ii = 0;ii < NF; ii++){
            for (int jj = 0; jj < NF; jj++){
                if (jj > ii) {
                    printf("%8.5e \t", zero);
                } else {
                    printf("%8.5e \t", S_hess_out[jj*NF+ii]);
                }
            }
            printf("\n");
        }
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



    for (int sens_forw = 1; sens_forw < 2; sens_forw++)
    {
    SECTION("sens_forw = " + std::to_string((bool)sens_forw))
    {
        for (int sens_adj = 1; sens_adj < 2; sens_adj++)
        {
        SECTION("sens_adj = " + std::to_string((bool)sens_adj))
        {
            for (int num_stages = 2; num_stages < 5; num_stages += 2)
            {
            SECTION("num_stages = " + std::to_string(num_stages))
            {
            for (int num_steps = 3; num_steps < 5; num_steps += 2)
            {
            SECTION("num_steps = " + std::to_string(num_steps))
            {

            for (std::string solver : solvers)
            {
            SECTION(solver)
            {


                double tol = sim_solver_tolerance_hess(solver);

                plan.sim_solver = hashitsim_hess(solver);

                // create correct config based on plan
                sim_solver_config *config = sim_config_create(plan);

            /* sim dims */
                void *dims = sim_dims_create(config);
                config->set_nx(dims, nx);
                config->set_nu(dims, nu);
                config->set_nz(dims, nz);

            /* sim options */

                void *opts_ = sim_opts_create(config, dims);
                sim_rk_opts *opts = (sim_rk_opts *) opts_;
                config->opts_initialize_default(config, dims, opts);

                opts->jac_reuse = false;        // jacobian reuse
                opts->newton_iter = 2;          // number of newton iterations per integration step

                opts->ns                = num_stages;          // number of stages in rk integrator
                opts->num_steps         = num_steps;    // number of steps
                opts->sens_forw         = (bool) sens_forw;
                opts->sens_adj          = (bool) sens_adj;
                opts->output_z          = false;
                opts->sens_algebraic    = false;
                opts->sens_hess         = false;


            /* sim in / out */

                sim_in *in = sim_in_create(config, dims);
                sim_out *out = sim_out_create(config, dims);

                in->T = T;

            /* set model */
                switch (plan.sim_solver)
                {
                    case ERK:  // ERK
                    {
                        sim_set_model(config, in, "expl_ode_fun", &expl_ode_fun);
                        sim_set_model(config, in, "expl_vde_for", &expl_vde_for);
                        sim_set_model(config, in, "expl_vde_adj", &expl_vde_adj);
                        sim_set_model(config, in, "expl_ode_hess", &expl_ode_hess);
                        break;
                    }
                    case IRK:  // IRK
                    {
                        sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
                        sim_set_model(config, in, "impl_ode_fun_jac_x_xdot",
                                &impl_ode_fun_jac_x_xdot);
                        sim_set_model(config, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
                        sim_set_model(config, in, "impl_ode_hess", &impl_ode_hess);
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
                sim_solver = sim_create(config, dims, opts);
                int acados_return;

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
                double rel_error_forw, rel_error_adj, rel_error_z, rel_error_alg;

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




            /************************************************
            * printing
            ************************************************/

                std::cout  << "rel_error_sim    = " << rel_error_x <<  "\n";
                if ( opts->sens_forw )
                std::cout  << "rel_error_forw   = " << rel_error_forw << "\n";
                if ( opts->sens_adj )
                std::cout  << "rel_error_adj    = " << rel_error_adj  << "\n";

            /************************************************
            * asserts on erors
            ************************************************/
                REQUIRE(rel_error_x <= tol);

                if ( opts->sens_forw )
                    REQUIRE(rel_error_forw <= tol);

                if ( opts->sens_adj )
                    REQUIRE(rel_error_adj <= tol);

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
    // explicit model
    external_function_casadi_free(&expl_ode_fun);
    external_function_casadi_free(&expl_ode_jac);
    external_function_casadi_free(&expl_vde_for);
    external_function_casadi_free(&expl_vde_adj);
    external_function_casadi_free(&expl_ode_hess);

}  // END_TEST_CASE
