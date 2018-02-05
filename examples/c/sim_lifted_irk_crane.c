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
#include <acados/sim/sim_casadi_wrapper.h>
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

// #define M_PI 3.14159265358979323846

int main() {

/************************************************
* create external functions
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

	// jacobian explicit ODE

	external_function_casadi exfun_jac;
	exfun_jac.casadi_fun = &jacFun;
	exfun_jac.casadi_work = &jacFun_work;
	exfun_jac.casadi_sparsity_in = &jacFun_sparsity_in;
	exfun_jac.casadi_sparsity_out = &jacFun_sparsity_out;

	int exfun_jac_size = external_function_casadi_calculate_size(&exfun_jac);
	void *exfun_jac_mem = malloc(exfun_jac_size);
	external_function_casadi_assign(&exfun_jac, exfun_jac_mem);

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

    in->step = T / sim_opts->num_steps;

	// casadi functins & wrappers
    in->vde = &vdeFun;
    in->jac = &jacFun;
//    in->hess = &hessFun;
    in->forward_vde_wrapper = &vde_fun;
    in->jacobian_wrapper = &jac_fun;
//    in->Hess_fun = &vde_hess_fun;
	// external functions
	in->exfun_forw_vde_expl = (external_function_generic *) &exfun_forw_vde;
	in->exfun_jac_ode_expl = (external_function_generic *) &exfun_jac;



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

	free(exfun_forw_vde_mem);
	free(exfun_jac_mem);

    free(xref);
	free(sim_opts_mem);
    free(in);
    free(workspace);
    free(raw_memory);
    free(out);

    return 0;
}
