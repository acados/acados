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
//#include <acados/sim/sim_casadi_wrapper.h>
#include "acados/utils/external_function_generic.h"
// crane model
#include "examples/c/crane_model/crane_model.h"

// blasfeo
#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>

// #define M_PI 3.14159265358979323846



int main()
{

/************************************************
* external functions
************************************************/

	// forward explicit VDE

	external_function_casadi exfun_forw_vde;
	exfun_forw_vde.casadi_fun = &vdeFun;
	exfun_forw_vde.casadi_work = &vdeFun_work;
	exfun_forw_vde.casadi_sparsity_in = &vdeFun_sparsity_in;
	exfun_forw_vde.casadi_sparsity_out = &vdeFun_sparsity_out;

	int exfun_forw_vde_size = external_function_casadi_calculate_size(&exfun_forw_vde);
	void *exfun_forw_vde_mem = malloc(exfun_forw_vde_size);
	external_function_casadi_assign(&exfun_forw_vde, exfun_forw_vde_mem);

	// adjoint explicit VDE

	external_function_casadi exfun_adj_vde;
	exfun_adj_vde.casadi_fun = &adjFun;
	exfun_adj_vde.casadi_work = &adjFun_work;
	exfun_adj_vde.casadi_sparsity_in = &adjFun_sparsity_in;
	exfun_adj_vde.casadi_sparsity_out = &adjFun_sparsity_out;

	int exfun_adj_vde_size = external_function_casadi_calculate_size(&exfun_adj_vde);
	void *exfun_adj_vde_mem = malloc(exfun_adj_vde_size);
	external_function_casadi_assign(&exfun_adj_vde, exfun_adj_vde_mem);

	// hessian explicit ODE

	external_function_casadi exfun_hess_ode;
	exfun_hess_ode.casadi_fun = &hessFun;
	exfun_hess_ode.casadi_work = &hessFun_work;
	exfun_hess_ode.casadi_sparsity_in = &hessFun_sparsity_in;
	exfun_hess_ode.casadi_sparsity_out = &hessFun_sparsity_out;

	int exfun_hess_vde_size = external_function_casadi_calculate_size(&exfun_hess_ode);
	void *exfun_hess_vde_mem = malloc(exfun_hess_vde_size);
	external_function_casadi_assign(&exfun_hess_ode, exfun_hess_vde_mem);

/************************************************
* bla bla bla
************************************************/

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
    erk_opts->sens_adj = true;
    erk_opts->sens_hess = false;
    // TODO(dimitris): SET IN DEFAULT ARGS
    erk_opts->num_forw_sens = NF;

	// sim in
    sim_in *in = create_sim_in(&dims);
    in->step = T / erk_opts->num_steps;
	// casadi functins & wrappers
//    in->vde = &vdeFun;
//    in->vde_adj = &adjFun;
//    in->hess = &hessFun;
//    in->forward_vde_wrapper = &vde_fun;
//    in->adjoint_vde_wrapper = &vde_adj_fun;
//    in->Hess_fun = &vde_hess_fun;
	// external functions
	in->exfun_forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
	in->exfun_adj_vde_expl = (external_function_generic *) &exfun_adj_vde;
	in->exfun_hess_ode_expl = (external_function_generic *) &exfun_hess_ode;



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
                printf("%8.5f ", S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }

    double *S_adj_out;
    if(erk_opts->sens_adj){
        S_adj_out = out->S_adj;
        printf("\nS_adj_out: \n");
        for (ii=0;ii<nx+nu;ii++){
            printf("%8.5f ", S_adj_out[ii]);
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
                    printf("%8.5f ", zero);
                }else{
                    printf("%8.5f ", S_hess_out[jj*NF+ii]);
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
		blasfeo_allocate_dmat(nx, nx+nu, &sA);
		blasfeo_pack_dmat(nx, nx+nu, S_forw_out, nx, &sA, 0, 0);

        struct blasfeo_dvec sx;
		blasfeo_allocate_dvec(nx, &sx);
		blasfeo_pack_dvec(nx, in->S_adj, &sx, 0);

        struct blasfeo_dvec sz;
		blasfeo_allocate_dvec(nx+nu, &sz);
//		blasfeo_print_dmat(nx, nx+nu, &sA, 0, 0);
//		blasfeo_print_tran_dvec(nx, &sx, 0);
        blasfeo_dgemv_t(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);

        printf("\nJac times lambdaX:\n");
        blasfeo_print_tran_dvec(nx+nu, &sz, 0);

		blasfeo_free_dmat(&sA);
		blasfeo_free_dvec(&sx);
		blasfeo_free_dvec(&sz);
    }

	free(exfun_forw_vde_mem);
	free(exfun_adj_vde_mem);
	free(exfun_hess_vde_mem);

    free(xref);
    free(in);
    free(solver);
    free(out);

    return flag;
}
