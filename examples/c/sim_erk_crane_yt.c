#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_erk_integrator_yt.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/yutao_model/yutao_model.h"

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

    sim_RK_opts *erk_opts = create_sim_RK_opts(&dims);

    sim_in *in = create_sim_in(nx, nu, NF);

    in->num_steps = 4;
    in->step = T / in->num_steps;
    in->sens_forw = true;
    in->sens_adj = true;
    in->sens_hess = false;

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

    sim_erk_memory *erk_mem = sim_erk_create_memory(in, erk_opts);

    sim_out *out = create_sim_out(nx, nu, NF);

    int flag = sim_erk_yt(in, out, erk_opts, erk_mem, NULL);

    double *xn = out->xn;

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    double *S_forw_out = out->S_forw;
    printf("\nS_forw_out: \n");
    for (ii=0;ii<nx;ii++){
        for (jj=0;jj<NF;jj++)
            printf("%8.5f ",S_forw_out[jj*nx+ii]);
        printf("\n");
    }

    double *S_adj_out = out->S_adj;
    printf("\nS_adj_out: \n");
    for (ii=0;ii<nx+nu;ii++){
        printf("%8.5f ",S_adj_out[ii]);
    }
    printf("\n");

    // double zero = 0.0;
    // if(in->sens_hess){
    //     double *S_hess_out = out->S_hess;
    //     printf("\nS_hess_out: \n");
    //     for (ii=0;ii<NF;ii++){
    //         for (jj=0;jj<NF;jj++){
    //             if (jj>ii){
    //                 printf("%8.5f ",zero);
    //             }else{
    //                 printf("%8.5f ",S_hess_out[jj*NF+ii]);
    //             }
    //         }
    //         printf("\n");
    //     }
    // }


    printf("\n");
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", out->info->ADtime);

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

    free(xref);
    free(erk_opts);
    free(in);
    free(erk_mem);
    free(out);
    v_free_align(mz);

    return flag;
}