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
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_erk_integrator.h>
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

// c interface
#ifdef ACADOS_WITH_C_INTERFACE
#include <acados_c/external_function_generic.h>
#include <acados_c/sim.h>
#include <acados_c/options.h>
#endif



int main()
{

/************************************************
* bla bla bla
************************************************/

    int NREP = 500;
    acados_timer timer;
    double Time1, Time2, Time3;

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



#ifdef ACADOS_WITH_C_INTERFACE



/************************************************
* external functions
************************************************/

	// forward explicit VDE

	external_function_casadi exfun_forw_vde;
	exfun_forw_vde.casadi_fun = &vdeFun;
	exfun_forw_vde.casadi_work = &vdeFun_work;
	exfun_forw_vde.casadi_sparsity_in = &vdeFun_sparsity_in;
	exfun_forw_vde.casadi_sparsity_out = &vdeFun_sparsity_out;

	create_external_function_casadi(&exfun_forw_vde);

	// adjoint explicit VDE

	external_function_casadi exfun_adj_vde;
	exfun_adj_vde.casadi_fun = &adjFun;
	exfun_adj_vde.casadi_work = &adjFun_work;
	exfun_adj_vde.casadi_sparsity_in = &adjFun_sparsity_in;
	exfun_adj_vde.casadi_sparsity_out = &adjFun_sparsity_out;

	create_external_function_casadi(&exfun_adj_vde);

	// hessian explicit ODE

	external_function_casadi exfun_hess_ode;
	exfun_hess_ode.casadi_fun = &hessFun;
	exfun_hess_ode.casadi_work = &hessFun_work;
	exfun_hess_ode.casadi_sparsity_in = &hessFun_sparsity_in;
	exfun_hess_ode.casadi_sparsity_out = &hessFun_sparsity_out;

	create_external_function_casadi(&exfun_hess_ode);

/************************************************
* sim dims
************************************************/

    sim_dims dims;
    dims.num_stages = num_stages;
    dims.nx = nx;
    dims.nu = nu;

/************************************************
* sim plan
************************************************/

    sim_solver_plan plan;
    plan.sim_solver = ERK;

/************************************************
* sim opts
************************************************/

    void *args = sim_create_args(&plan, &dims);
    
    sim_rk_opts *erk_opts = (sim_rk_opts *) args;
    erk_opts->num_steps = 4;
    erk_opts->sens_forw = true;
    erk_opts->sens_adj = true;
    erk_opts->sens_hess = false;
    // TODO(dimitris): SET IN DEFAULT ARGS
    erk_opts->num_forw_sens = NF;

/************************************************
* sim in
************************************************/

    sim_in *in = create_sim_in(&dims);

    in->T = T;

	// external functions
	in->forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
	in->adj_vde_expl = (external_function_generic *) &exfun_adj_vde;
	in->hess_ode_expl = (external_function_generic *) &exfun_hess_ode;



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

/************************************************
* sim_out
************************************************/

    sim_out *out = create_sim_out(&dims);

/************************************************
* sim solver
************************************************/

    sim_solver *solver = sim_create(&plan, &dims, args);

/************************************************
* call integrator
************************************************/

    acados_tic(&timer);
    for (ii=0;ii<NREP;ii++)
		sim_solve(solver, in, out);
    Time1 = acados_toc(&timer)/NREP;

    double *xn = out->xn;

/************************************************
* printing
************************************************/

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
    printf("cpt: %8.4f [ms]\n", 1000*out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", 1000*out->info->ADtime);

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

/************************************************
* free
************************************************/

	free_external_function_casadi(&exfun_forw_vde);
	free_external_function_casadi(&exfun_adj_vde);
	free_external_function_casadi(&exfun_hess_ode);

    free(xref);
    free(in);
    free(solver);
    free(out);



#else // ! ACADOS_WITH_C_INTERFACE



/************************************************
* external functions
************************************************/

	// forward explicit VDE

	external_function_casadi exfun_forw_vde;
	exfun_forw_vde.casadi_fun = &vdeFun;
	exfun_forw_vde.casadi_work = &vdeFun_work;
	exfun_forw_vde.casadi_sparsity_in = &vdeFun_sparsity_in;
	exfun_forw_vde.casadi_sparsity_out = &vdeFun_sparsity_out;

	int forw_vde_size = external_function_casadi_calculate_size(&exfun_forw_vde);
	void *forw_vde_mem = malloc(forw_vde_size);
	external_function_casadi_assign(&exfun_forw_vde, forw_vde_mem);

	// adjoint explicit VDE

	external_function_casadi exfun_adj_vde;
	exfun_adj_vde.casadi_fun = &adjFun;
	exfun_adj_vde.casadi_work = &adjFun_work;
	exfun_adj_vde.casadi_sparsity_in = &adjFun_sparsity_in;
	exfun_adj_vde.casadi_sparsity_out = &adjFun_sparsity_out;

	int adj_vde_size = external_function_casadi_calculate_size(&exfun_adj_vde);
	void *adj_vde_mem = malloc(adj_vde_size);
	external_function_casadi_assign(&exfun_adj_vde, adj_vde_mem);

	// hessian explicit ODE

	external_function_casadi exfun_hess_ode;
	exfun_hess_ode.casadi_fun = &hessFun;
	exfun_hess_ode.casadi_work = &hessFun_work;
	exfun_hess_ode.casadi_sparsity_in = &hessFun_sparsity_in;
	exfun_hess_ode.casadi_sparsity_out = &hessFun_sparsity_out;

	int hess_ode_size = external_function_casadi_calculate_size(&exfun_hess_ode);
	void *hess_ode_mem = malloc(hess_ode_size);
	external_function_casadi_assign(&exfun_hess_ode, hess_ode_mem);

/************************************************
* sim config
************************************************/

	sim_solver_config config;

	sim_erk_config_initialize_default(&config);

/************************************************
* sim dims
************************************************/

	int dims_size = sim_dims_calculate_size();
	void *dims_mem = malloc(dims_size);
	sim_dims *dims = sim_dims_assign(dims_mem);

    dims->num_stages = num_stages;
    dims->nx = nx;
    dims->nu = nu;

/************************************************
* sim opts
************************************************/

	int opts_size = config.opts_calculate_size(&config, dims);
	void *opts_mem = malloc(opts_size);
	sim_rk_opts *opts = config.opts_assign(&config, dims, opts_mem);
	config.opts_initialize_default(&config, dims, opts);

	opts->sens_adj = true;
	
/************************************************
* sim memory
************************************************/

	int mem_size = config.memory_calculate_size(&config, dims, opts);
	void *mem_mem = malloc(mem_size);
	void *mem = config.memory_assign(&config, dims, opts, mem_mem);

/************************************************
* sim workspace
************************************************/

	int work_size = config.workspace_calculate_size(&config, dims, opts);
	void *work = malloc(work_size);

/************************************************
* sim in
************************************************/

	int in_size = sim_in_calculate_size(dims, &config); // TODO move config as 1st
	void *in_mem = malloc(in_size);
	sim_in *in = sim_in_assign(dims, in_mem, &config); // TODO move config as 1st

    in->T = T;

	// external functions
	erk_model *model = in->model;
	model->forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
	model->adj_vde_expl = (external_function_generic *) &exfun_adj_vde;
	model->hess_ode_expl = (external_function_generic *) &exfun_hess_ode;

	// x
    for (ii = 0; ii < nx; ii++) {
        in->x[ii] = xref[ii];
    }

	// p
    for (ii = 0;ii < nu; ii++){
        in->u[ii] = 1.0;
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
* sim out
************************************************/

	int out_size = sim_out_calculate_size(dims);
	void *out_mem = malloc(out_size);
	sim_out *out = sim_out_assign(dims, out_mem);

/************************************************
* sim solver
************************************************/

    acados_tic(&timer);

    for (ii=0;ii<NREP;ii++)
		config.evaluate(&config, in, out, opts, mem, work);

    Time1 = acados_toc(&timer)/NREP;

    double *xn = out->xn;

/************************************************
* printing
************************************************/

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    double *S_forw_out;
    if(opts->sens_forw){
        S_forw_out = out->S_forw;
        printf("\nS_forw_out: \n");
        for (ii=0;ii<nx;ii++){
            for (jj=0;jj<NF;jj++)
                printf("%8.5f ", S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }

    double *S_adj_out;
    if(opts->sens_adj){
        S_adj_out = out->S_adj;
        printf("\nS_adj_out: \n");
        for (ii=0;ii<nx+nu;ii++){
            printf("%8.5f ", S_adj_out[ii]);
        }
        printf("\n");
    }

    double *S_hess_out;
    if(opts->sens_hess){
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
    printf("cpt: %8.4f [ms]\n", 1000*out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", 1000*out->info->ADtime);

    if(opts->sens_adj){
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

/************************************************
* free
************************************************/
	
	free(forw_vde_mem);
	free(adj_vde_mem);
	free(hess_ode_mem);

	free(dims_mem);
	free(opts_mem);
	free(mem_mem);
	free(work);
	free(in_mem);
	free(out_mem);



#endif // ACADOS_WITH_C_INTERFACE



	printf("\nsuccess!\n\n");

    return 0;
}
