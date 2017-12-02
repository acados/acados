#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/create.h"

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

    sim_dims dims;
    dims.num_stages = num_stages;
    dims.nx = nx;
    dims.nu = nu;

    sim_rk_opts *erk_opts = create_sim_erk_opts(&dims);

    sim_in *in = create_sim_in(&dims);

    erk_opts->num_steps = 4;
    in->step = T / erk_opts->num_steps;
    erk_opts->sens_forw = true;
    erk_opts->sens_adj = false;
    erk_opts->sens_hess = false;

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

    // TODO(dimitris): SET IN DEFAULT ARGS
    erk_opts->num_forw_sens = NF;

    int workspace_size = sim_erk_calculate_workspace_size(&dims, erk_opts);
    void *workspace = malloc(workspace_size);

    sim_out *out = create_sim_out(&dims);

    int flag = sim_erk(in, out, erk_opts, NULL, workspace);

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
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime*1000);
    printf("AD cpt: %8.4f [ms]\n", out->info->ADtime*1000);

    if(erk_opts->sens_adj){
        struct d_strmat sA;
        d_create_strmat(nx, nx+nu, &sA, S_forw_out);

        struct d_strvec sx;
        d_create_strvec(nx, &sx, in->S_adj);

        struct d_strvec sz;
        void *mz;
        v_zeros_align(&mz, d_size_strvec(nx+nu));
        d_create_strvec(nx+nu, &sz, mz);
        dgemv_t_libstr(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);

        printf("\nJac times lambdaX:\n");
        d_print_tran_strvec(nx+nu, &sz, 0);

        v_free_align(mz);
    }

    free(xref);
    free(erk_opts);
    free(in);
    free(workspace);
    free(out);

    return flag;
}