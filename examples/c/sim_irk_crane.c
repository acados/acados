#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// TODO(dimitris): add intergrator to c interface and clean this up
#include <acados_c/sim.h>
#include <acados_c/options.h>
#include "interfaces/acados_c/legacy_create.h"

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_irk_integrator.h"
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

// #define M_PI 3.14159265358979323846

int main() {

/************************************************
* external functions
************************************************/

	// implicit ODE

	external_function_casadi exfun_ode;
	exfun_ode.casadi_fun = &impl_odeFun;
	exfun_ode.casadi_work = &impl_odeFun_work;
	exfun_ode.casadi_sparsity_in = &impl_odeFun_sparsity_in;
	exfun_ode.casadi_sparsity_out = &impl_odeFun_sparsity_out;

	int exfun_ode_size = external_function_casadi_calculate_size(&exfun_ode);
	void *exfun_ode_mem = malloc(exfun_ode_size);
	external_function_casadi_assign(&exfun_ode, exfun_ode_mem);

	// jac_x implicit ODE

	external_function_casadi exfun_jac_x_ode;
	exfun_jac_x_ode.casadi_fun = &impl_jacFun_x;
	exfun_jac_x_ode.casadi_work = &impl_jacFun_x_work;
	exfun_jac_x_ode.casadi_sparsity_in = &impl_jacFun_x_sparsity_in;
	exfun_jac_x_ode.casadi_sparsity_out = &impl_jacFun_x_sparsity_out;

	int exfun_jac_x_ode_size = external_function_casadi_calculate_size(&exfun_jac_x_ode);
	void *exfun_jac_x_ode_mem = malloc(exfun_jac_x_ode_size);
	external_function_casadi_assign(&exfun_jac_x_ode, exfun_jac_x_ode_mem);

	// jac_xdot implicit ODE

	external_function_casadi exfun_jac_xdot_ode;
	exfun_jac_xdot_ode.casadi_fun = &impl_jacFun_xdot;
	exfun_jac_xdot_ode.casadi_work = &impl_jacFun_xdot_work;
	exfun_jac_xdot_ode.casadi_sparsity_in = &impl_jacFun_xdot_sparsity_in;
	exfun_jac_xdot_ode.casadi_sparsity_out = &impl_jacFun_xdot_sparsity_out;

	int exfun_jac_xdot_ode_size = external_function_casadi_calculate_size(&exfun_jac_xdot_ode);
	void *exfun_jac_xdot_ode_mem = malloc(exfun_jac_xdot_ode_size);
	external_function_casadi_assign(&exfun_jac_xdot_ode, exfun_jac_xdot_ode_mem);

	// jac_u implicit ODE

	external_function_casadi exfun_jac_u_ode;
	exfun_jac_u_ode.casadi_fun = &impl_jacFun_u;
	exfun_jac_u_ode.casadi_work = &impl_jacFun_u_work;
	exfun_jac_u_ode.casadi_sparsity_in = &impl_jacFun_u_sparsity_in;
	exfun_jac_u_ode.casadi_sparsity_out = &impl_jacFun_u_sparsity_out;

	int exfun_jac_u_ode_size = external_function_casadi_calculate_size(&exfun_jac_u_ode);
	void *exfun_jac_u_ode_mem = malloc(exfun_jac_u_ode_size);
	external_function_casadi_assign(&exfun_jac_u_ode, exfun_jac_u_ode_mem);

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

    sim_dims dims;
    dims.num_stages = num_stages;
    dims.nx = nx;
    dims.nu = nu;

    // the function to create opts should not be named by erk
    sim_rk_opts *irk_opts = create_sim_irk_opts(&dims);

    sim_in *in = create_sim_in(&dims);

    irk_opts->num_steps = 2;
    in->step = T / irk_opts->num_steps;
    irk_opts->sens_forw = true;
    irk_opts->sens_adj = true;
    irk_opts->sens_hess = false;
    irk_opts->jac_reuse = false;

	// casadi functins & wrappers
//    in->impl_ode = &impl_odeFun;
//    in->eval_impl_res = &impl_ode_fun;
//    in->impl_jac_x = &impl_jacFun_x;
//    in->eval_impl_jac_x = &impl_jac_x_fun;
//    in->impl_jac_xdot = &impl_jacFun_xdot;
//    in->eval_impl_jac_xdot = &impl_jac_xdot_fun;
//    in->impl_jac_u = &impl_jacFun_u;
//    in->eval_impl_jac_u = &impl_jac_u_fun;
	// external functions
	in->exfun_ode_impl = (external_function_generic *) &exfun_ode;
	in->exfun_jac_x_ode_impl = (external_function_generic *) &exfun_jac_x_ode;
	in->exfun_jac_xdot_ode_impl = (external_function_generic *) &exfun_jac_xdot_ode;
	in->exfun_jac_u_ode_impl = (external_function_generic *) &exfun_jac_u_ode;


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

    // TODO(dimitris): SET IN DEFAULT ARGS
    irk_opts->num_forw_sens = NF;

    int workspace_size = sim_irk_calculate_workspace_size(&dims, irk_opts);
    void *workspace = malloc(workspace_size);

    sim_out *out = create_sim_out(&dims);

    acados_tic(&timer);
    for (ii=0;ii<NREP;ii++)
        sim_irk(in, out, irk_opts, NULL, workspace);
    Time1 = acados_toc(&timer)/NREP;

    // printf("\nS_forw_out (blasfeo): \n");
    // blasfeo_print_dmat(nx, NF, ((sim_irk_workspace *)workspace)->S_forw,0 ,0);
    // exit(1);

#if 0
    irk_opts->jac_reuse = true;
    acados_tic(&timer);
    for (ii=0;ii<NREP;ii++)
        sim_irk(in, out, irk_opts, NULL, workspace);
    Time2 = acados_toc(&timer)/NREP;

    irk_opts->jac_reuse = false;
    irk_opts->sens_forw = false;
    irk_opts->sens_adj = true;
    acados_tic(&timer);
    for (ii=0;ii<NREP;ii++)
        sim_irk(in, out, irk_opts, NULL, workspace);
    Time3 = acados_toc(&timer)/NREP;
#endif

    printf("IRK Example for Inverted Pendulum:\n");
    printf("No. of Integration Steps: %d", irk_opts->num_steps);

    printf("\n");

    printf("IRK with exact Jacobian:\n");
    printf("NO. of run: %d   Avg cpt: %8.5f[ms]\n", NREP, Time1*1000);

    printf("IRK with Jacobian re-use:\n");
    printf("NO. of run: %d   Avg cpt: %8.5f[ms]\n", NREP, Time2*1000);

    printf("Speed Up Factor: %8.5f\n", Time1/Time2);

    printf("IRK with Adjoint Sensitivity:\n");
    printf("NO. of run: %d   Avg cpt: %8.5f[ms]\n", NREP, Time3*1000);

    double *xn = out->xn;

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    double *S_forw_out;
    if(irk_opts->sens_forw){
        S_forw_out = out->S_forw;
        printf("\nS_forw_out: \n");
        for (ii=0;ii<nx;ii++){
            for (jj=0;jj<NF;jj++)
                printf("%8.5f ",S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }

    double *S_adj_out;
    if(irk_opts->sens_adj){
        S_adj_out = out->S_adj;
        printf("\nS_adj_out: \n");
        for (ii=0;ii<nx+nu;ii++){
            printf("%8.5f ",S_adj_out[ii]);
        }
        printf("\n");
    }

    double *S_hess_out;
    if(irk_opts->sens_hess){
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


    // printf("\n");
    // printf("cpt: %8.5f [ms]\n", out->info->CPUtime*1000);
    // printf("AD cpt: %8.5f [ms]\n", out->info->ADtime*1000);

    // if(irk_opts->sens_adj && irk_opts->sens_forw){
    //     struct d_strmat sA;
    //     d_create_strmat(nx, nx+nu, &sA, S_forw_out);

    //     struct d_strvec sx;
    //     d_create_strvec(nx, &sx, in->S_adj);

    //     struct d_strvec sz;
    //     void *mz;
    //     v_zeros_align(&mz, d_size_strvec(nx+nu));
    //     d_create_strvec(nx+nu, &sz, mz);
    //     dgemv_t_libstr(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);

    //     printf("\nJac times lambdaX:\n");
    //     d_print_tran_strvec(nx+nu, &sz, 0);

    //     v_free_align(mz);
    // }

    free(xref);
    free(irk_opts);
    free(in);
    free(workspace);
    free(out);

    return 0;
}
