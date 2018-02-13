#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// TODO(dimitris): add intergrator to c interface and clean this up
#include <acados_c/sim.h>
#include <acados_c/options.h>
#include "interfaces/acados_c/legacy_create.h"

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
//#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"
// crane model
#include "examples/c/crane_model/crane_model.h"

// blasfeo
#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"

// c interface
#include <acados_c/external_function_generic.h>



#define NO_C_INTERFACE



int main() {

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
    int num_stages = 3;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;



#ifdef NO_C_INTERFACE

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

	// jacobian explicit ODE

	external_function_casadi exfun_jac;
	exfun_jac.casadi_fun = &jacFun;
	exfun_jac.casadi_work = &jacFun_work;
	exfun_jac.casadi_sparsity_in = &jacFun_sparsity_in;
	exfun_jac.casadi_sparsity_out = &jacFun_sparsity_out;

	int jac_size = external_function_casadi_calculate_size(&exfun_jac);
	void *jac_mem = malloc(jac_size);
	external_function_casadi_assign(&exfun_jac, jac_mem);

/************************************************
* sim config
************************************************/

	sim_solver_config config;

	sim_lifted_irk_config_initialize_default(&config);

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

	int opts_size = config.opts_calculate_size(dims);
	void *opts_mem = malloc(opts_size);
	sim_rk_opts *opts = config.opts_assign(dims, opts_mem);
	config.opts_initialize_default(dims, opts);

	opts->sens_adj = true;
	
/************************************************
* sim memory
************************************************/

	int mem_size = config.memory_calculate_size(dims, opts);
	void *mem_mem = malloc(mem_size);
	void *mem = config.memory_assign(dims, opts, mem_mem);

/************************************************
* sim workspace
************************************************/

	int work_size = config.workspace_calculate_size(dims, opts);
	void *work = malloc(work_size);

/************************************************
* sim in
************************************************/

	int in_size = sim_in_calculate_size(dims, &config);
	void *in_mem = malloc(in_size);
	sim_in *in = sim_in_assign(dims, in_mem, &config);

    in->T = T;

	// external functions
	lifted_irk_model *model = in->model;
	model->forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
	model->jac_ode_expl = (external_function_generic *) &exfun_jac;

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
		config.fun(in, out, opts, mem, work);

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
	free(jac_mem);

	free(dims_mem);
	free(opts_mem);
	free(mem_mem);
	free(work);
	free(in_mem);
	free(out_mem);



#else



/************************************************
* create external functions
************************************************/

	// forward explicit VDE

	external_function_casadi exfun_forw_vde;
	exfun_forw_vde.casadi_fun = &vdeFun;
	exfun_forw_vde.casadi_work = &vdeFun_work;
	exfun_forw_vde.casadi_sparsity_in = &vdeFun_sparsity_in;
	exfun_forw_vde.casadi_sparsity_out = &vdeFun_sparsity_out;

	create_external_function_casadi(&exfun_forw_vde);

	// jacobian explicit ODE

	external_function_casadi exfun_jac;
	exfun_jac.casadi_fun = &jacFun;
	exfun_jac.casadi_work = &jacFun_work;
	exfun_jac.casadi_sparsity_in = &jacFun_sparsity_in;
	exfun_jac.casadi_sparsity_out = &jacFun_sparsity_out;

	create_external_function_casadi(&exfun_jac);

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
    int num_stages = 3;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;

/************************************************
* sim dims
************************************************/

    sim_dims dims;
    dims.num_stages = num_stages;
    dims.nx = nx;
    dims.nu = nu;

/************************************************
* sim opts
************************************************/

	int sim_opts_size = sim_lifted_irk_opts_calculate_size(&dims);
	void *sim_opts_mem = malloc(sim_opts_size);

	sim_rk_opts *sim_opts = sim_lifted_irk_assign_opts(&dims, sim_opts_mem);

//    sim_rk_opts *rk_opts = create_sim_lifted_irk_opts(&dims);

	sim_lifted_irk_initialize_default_args(&dims, sim_opts);

    sim_opts->num_steps = 2;
    sim_opts->sens_forw = true;
    sim_opts->sens_adj = true;
    sim_opts->sens_hess = false;

/************************************************
* sim in
************************************************/

    sim_in *in = create_sim_in(&dims);

    in->T = T;

	// external functions
	in->forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
	in->jac_ode_expl = (external_function_generic *) &exfun_jac;



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
* out
************************************************/

    sim_out *out = create_sim_out(&dims);

/************************************************
* workspace
************************************************/

    int workspace_size = sim_lifted_irk_calculate_workspace_size(&dims, sim_opts);
    void *workspace = malloc(workspace_size);

/************************************************
* memory
************************************************/

	int memory_size = sim_lifted_irk_calculate_memory_size(&dims, sim_opts);
	void *raw_memory = malloc(memory_size);

	sim_lifted_irk_memory *memory = sim_lifted_irk_assign_memory(&dims, sim_opts, raw_memory);

/************************************************
* call integrator
************************************************/

    acados_tic(&timer);
    for (ii=0;ii<NREP;ii++)
        sim_lifted_irk(in, out, sim_opts, memory, workspace);
    Time1 = acados_toc(&timer)/NREP;

/************************************************
* printing
************************************************/

    double *xn = out->xn;

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    double *S_forw_out;
    if(sim_opts->sens_forw){
        S_forw_out = out->S_forw;
        printf("\nS_forw_out: \n");
        for (ii=0;ii<nx;ii++){
            for (jj=0;jj<NF;jj++)
                printf("%8.5f ",S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }

    double *S_adj_out;
    if(sim_opts->sens_adj){
        S_adj_out = out->S_adj;
        printf("\nS_adj_out: \n");
        for (ii=0;ii<nx+nu;ii++){
            printf("%8.5f ",S_adj_out[ii]);
        }
        printf("\n");
    }

    double *S_hess_out;
    if(sim_opts->sens_hess){
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
    printf("cpt: %8.4f [ms]\n", 1000*out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", 1000*out->info->ADtime);

    if(sim_opts->sens_adj){
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
	free_external_function_casadi(&exfun_jac);

    free(xref);
	free(sim_opts_mem);
    free(in);
    free(workspace);
    free(raw_memory);
    free(out);

#endif

	printf("\nsuccess!\n\n");

    return 0;
}
