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
// NOTE(nielsvd): required to cast memory etc. should go.
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_erk_integrator.h>
#include <acados/sim/sim_casadi_wrapper.h>

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

    int ii;
    int jj;

    int nx = 4;
    int nu = 1;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05;
    int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;

    sim_solver_plan plan;
    plan.sim_solver = ERK;

    sim_dims dims;
    dims.num_stages = num_stages;
    dims.nx = nx;
    dims.nu = nu;

    void *args = sim_create_args(&plan, &dims);
    
    sim_rk_opts *erk_opts = (sim_rk_opts *) args;
    erk_opts->num_steps = 4;
    erk_opts->sens_forw = true;
    erk_opts->sens_adj = false;
    erk_opts->sens_hess = false;
    // TODO(dimitris): SET IN DEFAULT ARGS
    erk_opts->num_forw_sens = NF;

    sim_in *in = create_sim_in(&dims);
    in->step = T / erk_opts->num_steps;
    in->vde = &vdeFun;
    in->vde_adj = &adjFun;
    in->hess = &hessFun;
    in->forward_vde_wrapper = &vde_fun;
    in->adjoint_vde_wrapper = &vde_adj_fun;
    in->Hess_fun = &vde_hess_fun;

    for (ii = 0; ii < nx; ii++) {
        in->x[ii] = xref[ii];
    }
    for (ii = 0;ii < nu; ii++){
        in->u[ii] = 1.0;
    }

    for (ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;

    for (ii = 0; ii < nx; ii++)
        in->S_adj[ii] = 1.0;

    sim_solver *solver = sim_create(&plan, &dims, args);

    sim_out *out = create_sim_out(&dims);

    int flag = sim_solve(solver, in, out);

    double *xn = out->xn;

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    double *S_forw_out;
    if(erk_opts->sens_forw){
        S_forw_out = out->S_forw;
        printf("\nS_forw_out: \n");
        for (ii=0;ii<nx;ii++){
            for (jj=0;jj<NF;jj++)
                printf("%8.5f ",S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }

    double *S_adj_out;
    if(erk_opts->sens_adj){
        S_adj_out = out->S_adj;
        printf("\nS_adj_out: \n");
        for (ii=0;ii<nx+nu;ii++){
            printf("%8.5f ",S_adj_out[ii]);
        }
        printf("\n");
    }

    double *S_hess_out;
    if(erk_opts->sens_hess){
        double zero = 0.0;
        S_hess_out = out->S_hess;
        printf("\nS_hess_out: \n");
        for (ii=0;ii<NF;ii++){
            for (jj=0;jj<NF;jj++){
                if (jj>ii){
                    printf("%8.5f ",zero);
                }else{
                    printf("%8.5f ",S_hess_out[jj*NF+ii]);
                }
            }
            printf("\n");
        }
    }


    printf("\n");
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", out->info->ADtime);

    if(erk_opts->sens_adj){
        struct blasfeo_dmat sA;
        blasfeo_create_dmat(nx, nx+nu, &sA, S_forw_out);

        struct blasfeo_dvec sx;
        blasfeo_create_dvec(nx, &sx, in->S_adj);

        struct blasfeo_dvec sz;
        void *mz;
        v_zeros_align(&mz, blasfeo_memsize_dvec(nx+nu));
        blasfeo_create_dvec(nx+nu, &sz, mz);
        blasfeo_dgemv_t(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);

        printf("\nJac times lambdaX:\n");
        blasfeo_print_tran_dvec(nx+nu, &sz, 0);

        v_free_align(mz);
    }

    free(xref);
    free(in);
    free(solver);
    free(out);

    return flag;
}