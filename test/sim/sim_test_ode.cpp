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

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "test/test_utils/eigen.h"
#include "catch/include/catch.hpp"

// acados
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// wt model
#include "examples/c/wt_model_nx3/wt_model.h"

// x0 and u for simulation
#include "examples/c/wt_model_nx3/u_x0.c"

extern "C"
{

}

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

sim_solver_t hashitsim(std::string const& inString)
{
    if (inString == "ERK") return ERK;
    if (inString == "IRK") return IRK;
    if (inString == "GNSF") return GNSF;
    if (inString == "LIFTED_IRK") return LIFTED_IRK;

    return (sim_solver_t) -1;
}

double sim_solver_tolerance(std::string const& inString)
{
    if (inString == "ERK") return 1e-7;
    if (inString == "IRK") return 1e-7;
    if (inString == "GNSF") return 1e-7;
    if (inString == "LIFTED_IRK") return 1e-5;

    return -1;
}




TEST_CASE("wt_nx3_example", "[integrators]")
{
    vector<std::string> solvers = {"ERK", "IRK", "GNSF", "LIFTED_IRK"};
    // initialize dimensions
    int ii, jj;

    const int nx = 3;
    const int nu = 4;
    int NF = nx + nu;  // columns of forward seed

    int nsim0 = 1;  // nsim;

    double T = 0.05;  // simulation time

    double x_sim[nx*(nsim0+1)];
    double x_ref_sol[nx];
    double S_forw_ref_sol[nx*NF];
    double S_adj_ref_sol[NF];

    double error[nx];
    double error_S_forw[nx*NF];
    double error_S_adj[NF];

    double max_error, max_error_forw, max_error_adj;

    for (ii=0; ii < nx; ii++)
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


    /************************************************
    * Create Reference Solution
    ************************************************/

    sim_solver_plan plan;
    plan.sim_solver = IRK;

    sim_config *config = sim_config_create(plan);

    void *dims = sim_dims_create(config);

    /* set dimensions */
    sim_dims_set(config, dims, "nx", &nx);
    sim_dims_set(config, dims, "nu", &nu);

    void *opts_ = sim_opts_create(config, dims);
    sim_opts *opts = (sim_opts *) opts_;

    opts->sens_forw = true;
    opts->sens_adj = true;

    opts->jac_reuse = false;  // jacobian reuse
    opts->newton_iter = 5;  // number of newton iterations per integration step
    opts->num_steps = 10;  // number of steps
    opts->ns = 5;  // number of stages in rk integrator

    sim_in *in = sim_in_create(config, dims);
    sim_out *out = sim_out_create(config, dims);

    in->T = T;

    sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
    sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot);
    sim_in_set(config, dims, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);

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

    sim_solver *sim_solver = sim_solver_create(config, dims, opts);

    int acados_return;

    for (ii=0; ii < nsim0; ii++)
    {
        // x
        for (jj = 0; jj < nx; jj++)
            in->x[jj] = x_sim[ii*nx+jj];

        // u
        for (jj = 0; jj < nu; jj++)
            in->u[jj] = u_sim[ii*nu+jj];

        acados_return = sim_solve(sim_solver, in, out);
        REQUIRE(acados_return == 0);

        for (jj = 0; jj < nx; jj++)
            x_sim[(ii+1)*nx+jj] = out->xn[jj];
    }

    for (jj = 0; jj < nx; jj++)
        x_ref_sol[jj] = out->xn[jj];

    for (jj = 0; jj < nx*NF; jj++)
        S_forw_ref_sol[jj] = out->S_forw[jj];

    for (jj = 0; jj < NF; jj++)
        S_adj_ref_sol[jj] = out->S_adj[jj];

    // printf("Reference forward sensitivities \n");
    // d_print_exp_mat(nx, NF, &S_forw_ref_sol[0], 1);


    free(config);
    free(dims);
    free(opts);

    free(in);
    free(out);
    free(sim_solver);

    for (std::string solver : solvers)
    {
        SECTION(solver)
        {
            for (int num_steps = 1; num_steps < 4; num_steps++)
            {
                double tol = sim_solver_tolerance(solver);

                plan.sim_solver = hashitsim(solver);

                // create correct config based on plan
                sim_config *config = sim_config_create(plan);

                /************************************************
                * sim dims
                ************************************************/

                void *dims = sim_dims_create(config);
                /* set dimensions */
                sim_dims_set(config, dims, "nx", &nx);
                sim_dims_set(config, dims, "nu", &nu);

                /************************************************
                * sim opts
                ************************************************/

                void *opts_ = sim_opts_create(config, dims);
                sim_opts *opts = (sim_opts *) opts_;

                if (plan.sim_solver != LIFTED_IRK)
                    opts->sens_adj = true;
                else
                    opts->sens_adj = false;

                opts->jac_reuse = true;  // jacobian reuse
                opts->newton_iter = 1;  // number of newton iterations per integration step
                opts->num_steps = num_steps;  // number of steps

                // gnsf dimension
                int nx1 = nx;
                int nz1 = 0;
                int ny = nx;
                int nuhat = nu;
                int nout = 1;
                int nz = 0;
                switch (plan.sim_solver)
                {

                    case ERK:
                         // ERK
                        opts->ns = 4;  // number of stages in rk integrator
                        break;

                    case IRK:
                         // IRK
                        opts->ns = 2;  // number of stages in rk integrator
                        break;

                    case GNSF:
                        // GNSF
                        opts->ns = 2;  // number of stages in rk integrator

                        // set additional dimensions
                        sim_dims_set(config, dims, "nx1", &nx1);
                        sim_dims_set(config, dims, "nz", &nz);
                        sim_dims_set(config, dims, "nz1", &nz1);
                        sim_dims_set(config, dims, "nout", &nout);
                        sim_dims_set(config, dims, "ny", &ny);
                        sim_dims_set(config, dims, "nuhat", &nuhat);

                        break;

                    case LIFTED_IRK:
                        // new lifted IRK
                        opts->ns = 2;  // number of stages in rk integrator
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
                switch (plan.sim_solver)
                {
                    case ERK:  // ERK
                    {
                        sim_in_set(config, dims, in, "expl_ode_fun", &expl_ode_fun);
                        sim_in_set(config, dims, in, "expl_vde_for", &expl_vde_for);
                        sim_in_set(config, dims, in, "expl_vde_adj", &expl_vde_adj);
                        break;
                    }
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
                    case LIFTED_IRK:  // lifted_irk
                    {
                        sim_in_set(config, dims, in, "impl_ode_fun", &impl_ode_fun);
                        sim_in_set(config, dims, in, "impl_ode_fun_jac_x_xdot_u",
                                 &impl_ode_fun_jac_x_xdot_u);
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
                std::cout << "\n---> testing integrator " << solver <<
                        " (num_steps = " << opts->num_steps << ", num_stages = " << opts->ns
                        << ", jac_reuse = " << opts->jac_reuse << ", newton_iter = "
                        << opts->newton_iter << ")\n";

                sim_solver = sim_solver_create(config, dims, opts);

                int acados_return;

                sim_precompute(sim_solver, in, out);

                for (ii=0; ii < nsim0; ii++)
                {
                    // x
                    for (jj = 0; jj < nx; jj++)
                        in->x[jj] = x_sim[ii*nx+jj];

                    // u
                    for (jj = 0; jj < nu; jj++)
                        in->u[jj] = u_sim[ii*nu+jj];

                    acados_return = sim_solve(sim_solver, in, out);
                    REQUIRE(acados_return == 0);

                    for (jj = 0; jj < nx; jj++){
                        x_sim[(ii+1)*nx+jj] = out->xn[jj];
                        REQUIRE(std::isnan(out->xn[jj]) == 0);
                    }

                }

                // error sim
                for (jj = 0; jj < nx; jj++)
                    error[jj] = fabs(out->xn[jj] - x_ref_sol[jj]);

                max_error = 0.0;
                for (int ii = 0; ii < nx; ii++)
                    max_error = (error[ii] >= max_error) ? error[ii] : max_error;

                // error_S_forw
                for (jj = 0; jj < nx*NF; jj++){
                    REQUIRE(std::isnan(out->S_forw[jj]) == 0);
                    error_S_forw[jj] = fabs(S_forw_ref_sol[jj] - out->S_forw[jj]);
                }

                max_error_forw = 0.0;
                for (jj = 0; jj < nx*NF; jj++)
                    max_error_forw = (error_S_forw[jj] >= max_error_forw)
                            ? error_S_forw[jj] : max_error_forw;

                // error_S_adj
                for (jj = 0; jj < NF; jj++){
                    REQUIRE(std::isnan(out->S_forw[jj]) == 0);
                    error_S_adj[jj] = S_adj_ref_sol[jj] - out->S_adj[jj];
                }
                max_error_adj = 0.0;
                for (jj = 0; jj < NF; jj++)
                    max_error_adj = (error_S_adj[jj] >= max_error_adj)
                            ? error_S_adj[jj] : max_error_adj;

                /************************************************
                * printing
                ************************************************/

                std::cout  << "error_sim   = " << max_error << "\n";
                std::cout  << "error_forw  = " << max_error_forw << "\n";

                REQUIRE(max_error <= tol);
                REQUIRE(max_error_forw <= tol);

                // TODO(FreyJo): implement adjoint sensitivites for these integrators!!!
                if ((plan.sim_solver != LIFTED_IRK)){
                    std::cout  << "error_adj   = " << max_error_adj  << "\n";
                    REQUIRE(max_error_adj <= tol);
                }


                free(config);
                free(dims);
                free(opts);

                free(in);
                free(out);
                free(sim_solver);
            }  // end for num_steps
        }  // end section
    }  // END FOR SOLVERS

    // explicit model
    external_function_casadi_free(&expl_ode_fun);
    external_function_casadi_free(&expl_vde_for);
    external_function_casadi_free(&expl_vde_adj);
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
