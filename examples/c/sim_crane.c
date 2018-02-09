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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// acados
#include <acados_c/sim.h>
#include <acados_c/options.h>
#include <acados_c/utils/casadi_wrapper.h>
// NOTE(nielsvd): required to cast memory etc. should go.
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_erk_integrator.h>
#include <acados/sim/sim_irk_integrator.h>
#include <acados/sim/sim_lifted_irk_integrator.h>
// #include <acados/sim/sim_casadi_wrapper.h>

#include "examples/c/crane_model/crane_model.h"

// blasfeo
#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>

// #define M_PI 3.14159265358979323846

int main() {

    int kk;
    int jj;
    int flag;

    int nx = 4;
    int nu = 1;
    int np = 1;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05;
    int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;

    /************************************************
     * sim solvers
     ************************************************/

    sim_solver_t sim_solvers[] =
    {
        ERK,
        IRK,
        LIFTED_IRK
    };

    int ii_max = 3;

    for (int ii = 0; ii < ii_max; ii++) {
        sim_solver_config config;
        config.sim_solver = sim_solvers[ii];
        config.extfun.type = CASADI_WRAPPER;

        if (config.sim_solver == IRK){
            num_stages = 3;  // TODO(nielsvd): IRK integrator only supports 3 stages atm
        } else {
            num_stages = 4;
        }

        sim_dims dims;
        dims.num_stages = num_stages;
        dims.nx = nx;
        dims.nu = nu;
        dims.np = np;

        sim_solver_fcn_ptrs *fcn_ptrs =
        create_sim_solver_fcn_ptrs(&config, &dims);

        void *args = sim_create_args(fcn_ptrs, &dims);
        int num_steps = 4;
        bool sens_forw = true;
        bool sens_adj = false;
        bool sens_hess = false;

        sim_erk_integrator_args *erk_args = (sim_erk_integrator_args *)args;
        sim_irk_integrator_args *irk_args = (sim_irk_integrator_args *)args;
        sim_lifted_irk_integrator_args *lifted_irk_args =
            (sim_lifted_irk_integrator_args *)args;
        switch (config.sim_solver) {
            case ERK:
                printf("\nCrane simulation using ERK:\n\n");
                erk_args->num_steps = num_steps;
                erk_args->sens_forw = sens_forw;
                erk_args->sens_adj = sens_adj;
                erk_args->sens_hess = sens_hess;
                // TODO(dimitris): SET IN DEFAULT ARGS
                erk_args->num_forw_sens = NF;
                // Forward VDE
                ((casadi_wrapper_args *)(erk_args->forward_vde_args))->fun =
                    &vdeFun;
                ((casadi_wrapper_args *)(erk_args->forward_vde_args))->dims =
                    &vdeFun_work;
                ((casadi_wrapper_args *)(erk_args->forward_vde_args))
                    ->sparsity = &vdeFun_sparsity_out;
                // Adjoint VDE
                ((casadi_wrapper_args *)(erk_args->adjoint_vde_args))->fun =
                    &adjFun;
                ((casadi_wrapper_args *)(erk_args->adjoint_vde_args))->dims =
                    &adjFun_work;
                ((casadi_wrapper_args *)(erk_args->adjoint_vde_args))
                    ->sparsity = &adjFun_sparsity_out;
                // Hessian VDE
                ((casadi_wrapper_args *)(erk_args->hess_vde_args))->fun =
                    &hessFun;
                ((casadi_wrapper_args *)(erk_args->hess_vde_args))->dims =
                    &hessFun_work;
                ((casadi_wrapper_args *)(erk_args->hess_vde_args))->sparsity =
                    &hessFun_sparsity_out;
                break;
            case IRK:
                printf("\nCrane simulation using IRK:\n\n");
                irk_args->num_steps = num_steps;
                irk_args->sens_forw = sens_forw;
                irk_args->sens_adj = sens_adj;
                irk_args->sens_hess = sens_hess;
                irk_args->jac_reuse = false;
                // TODO(dimitris): SET IN DEFAULT ARGS
                irk_args->num_forw_sens = NF;
                // Implicit ODE
                ((casadi_wrapper_args *)(irk_args->impl_res_args))->fun =
                    &impl_odeFun;
                ((casadi_wrapper_args *)(irk_args->impl_res_args))->dims =
                    &impl_odeFun_work;
                ((casadi_wrapper_args *)(irk_args->impl_res_args))->sparsity =
                    &impl_odeFun_sparsity_out;
                // Jacobians of implicit ODE
                ((casadi_wrapper_args *)(irk_args->impl_jac_args))->fun =
                    &impl_jacFun;
                ((casadi_wrapper_args *)(irk_args->impl_jac_args))->dims =
                    &impl_jacFun_work;
                ((casadi_wrapper_args *)(irk_args->impl_jac_args))->sparsity =
                    &impl_jacFun_sparsity_out;
                break;
            case LIFTED_IRK:
                printf("\nCrane simulation using Lifted IRK:\n\n");
                lifted_irk_args->num_steps = num_steps;
                lifted_irk_args->sens_forw = sens_forw;
                lifted_irk_args->sens_adj = sens_adj;
                lifted_irk_args->sens_hess = sens_hess;
                // TODO(dimitris): SET IN DEFAULT ARGS
                lifted_irk_args->num_forw_sens = NF;
                // Forward VDE
                ((casadi_wrapper_args *)(lifted_irk_args->forward_vde_args))
                    ->fun = &vdeFun;
                ((casadi_wrapper_args *)(lifted_irk_args->forward_vde_args))
                    ->dims = &vdeFun_work;
                ((casadi_wrapper_args *)(lifted_irk_args->forward_vde_args))
                    ->sparsity = &vdeFun_sparsity_out;
                // Adjoint VDE
                ((casadi_wrapper_args *)(lifted_irk_args->jacobian_ode_args))
                    ->fun = &jacFun;
                ((casadi_wrapper_args *)(lifted_irk_args->jacobian_ode_args))
                    ->dims = &jacFun_work;
                ((casadi_wrapper_args *)(lifted_irk_args->jacobian_ode_args))
                    ->sparsity = &jacFun_sparsity_out;
                break;
        }

        sim_in *in = create_sim_in(&dims);

        in->step = T / num_steps;

        for (kk = 0; kk < nx; kk++) {
            in->x[kk] = xref[kk];
        }
        for (kk = 0; kk < nu; kk++) {
            in->u[kk] = 1.0;
        }
        for (kk = 0; kk < np; kk++) {
            in->p[kk] = 1.0;
        }

        for (kk = 0; kk < nx * NF; kk++) in->S_forw[kk] = 0.0;
        for (kk = 0; kk < nx; kk++) in->S_forw[kk * (nx + 1)] = 1.0;

        for (kk = 0; kk < nx; kk++) in->S_adj[kk] = 1.0;

        sim_solver *solver = sim_create(fcn_ptrs, &dims, args);

        free(fcn_ptrs);

        sim_out *out = create_sim_out(&dims);

        flag = sim_solve(solver, in, out);

        double *xn = out->xn;

        printf("\nxn: \n");
        for (kk = 0; kk < nx; kk++) printf("%8.5f ", xn[kk]);
        printf("\n");

        double *S_forw_out;
        if (sens_forw) {
            S_forw_out = out->S_forw;
            printf("\nS_forw_out: \n");
            for (kk = 0; kk < nx; kk++) {
                for (jj = 0; jj < NF; jj++)
                    printf("%8.5f ", S_forw_out[jj * nx + kk]);
                printf("\n");
            }
        }

        double *S_adj_out;
        if (sens_adj) {
            S_adj_out = out->S_adj;
            printf("\nS_adj_out: \n");
            for (kk = 0; kk < nx + nu; kk++) {
                printf("%8.5f ", S_adj_out[kk]);
            }
            printf("\n");
        }

        double *S_hess_out;
        if (sens_hess) {
            double zero = 0.0;
            S_hess_out = out->S_hess;
            printf("\nS_hess_out: \n");
            for (kk = 0; kk < NF; kk++) {
                for (jj = 0; jj < NF; jj++) {
                    if (jj > kk) {
                        printf("%8.5f ", zero);
                    } else {
                        printf("%8.5f ", S_hess_out[jj * NF + kk]);
                    }
                }
                printf("\n");
            }
        }

        printf("\n");
        printf("cpt: %8.4f [ms]\n", out->info->CPUtime);
        printf("AD cpt: %8.4f [ms]\n", out->info->ADtime);

        if (sens_adj) {
            struct blasfeo_dmat sA;
            blasfeo_create_dmat(nx, nx + nu, &sA, S_forw_out);

            struct blasfeo_dvec sx;
            blasfeo_create_dvec(nx, &sx, in->S_adj);

            struct blasfeo_dvec sz;
            void *mz;
            v_zeros_align(&mz, blasfeo_memsize_dvec(nx + nu));
            blasfeo_create_dvec(nx + nu, &sz, mz);
            blasfeo_dgemv_t(nx, nx + nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0,
                            &sz, 0);

            printf("\nJac times lambdaX:\n");
            blasfeo_print_tran_dvec(nx + nu, &sz, 0);

            v_free_align(mz);
        }

        free(in);
        free(out);
        free(solver);       
        free(args);
    }

    free(xref);
    
    return flag;
}