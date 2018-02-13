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

// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_gnsf.h"
//#include "acados/sim/sim_erk_integrator.h"
//#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/sim/sim_gnsf_casadi_wrapper.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux.h"




void print_gnsf_dims( gnsf_dims *dims)
{
    printf("\n");
    printf("nx          = %d \n", dims->nx);
    printf("nu          = %d \n", dims->nu);
    printf("nz          = %d \n", dims->nz);
    printf("nx1         = %d \n", dims->nx1);
    printf("nx2         = %d \n", dims->nx2);
    printf("n_out       = %d \n", dims->n_out);
    printf("n_stages    = %d \n", dims->num_stages);
    printf("n_steps     = %d \n", dims->num_steps);
    printf("\n");
}

void print_gnsf_res_in( gnsf_dims *dims, double *res_in )
{
    int nff = dims->num_stages * dims->n_out;
    printf("res_in ff = \n");    
    for (int i=0; i<nff; i++)
        printf("\t%f\n", res_in[i]);
    printf("res_in x1 = \n");
    for (int i=nff; i< nff + dims->nx1; i++)
        printf("\t%f\n", res_in[i]);
    printf("res_in u = \n");
    for (int i= nff + dims->nx1; i< nff+ dims->nx1 + dims->nu; i++)
        printf("\t%f\n", res_in[i]);   
}


void print_gnsf_res_out( gnsf_dims *dims, double *res_out )
{
    int nff = dims->num_stages * dims->n_out;
    printf("res_out res_val = \n");    
    for (int i=0; i<nff; i++)
        printf("\t%5.5f\n", res_out[i]);
    printf("\nres_out J_res_ff = \n");
    for (int i=0; i<nff; i++){
        for (int j=1; j<nff+1; j++){
            printf("\t %5.5f", res_out[i+j*nff]);}
        printf("\n");
    }
}

void gnsf_get_dims( gnsf_dims *dims, casadi_function_t get_ints_fun)
{
    double *ints_out;
    ints_out = (double*) calloc(8,sizeof(double));
    export_from_ML_wrapped(ints_out, ints_out, get_ints_fun);

    dims->nx = (int) ints_out[0];
    dims->nu = (int) ints_out[1]; 
    dims->nz = (int) ints_out[2];
    dims->nx1 = (int) ints_out[3];
    dims->nx2 = (int) ints_out[4];
    dims->num_stages = (int) ints_out[5];
    dims->num_steps = (int) ints_out[6];
    dims->n_out = (int) ints_out[7];
}

void gnsf_get_KK_mat(gnsf_dims *dims, gnsf_fixed *fix, casadi_function_t KK_mat_fun)
{
    int nK1 = dims->nx1 * dims->num_stages;
    int nff = dims->num_stages * dims->n_out;
    double *out;
    out = (double*) calloc((nK1) * (nff + dims->nx1 + dims->nu),sizeof(double));
    // export_from_ML_wrapped(fix->KKf, fix->KKf, KK_mat_fun);
    export_from_ML_wrapped(out, out, KK_mat_fun);
    // for (int i=0; i<dims->nx1 *dims->num_stages*nff; i++)
    // fix->KKf[i] = out[i];
    // fix->KKf = &out[0];
    // printf("test\n");
    blasfeo_allocate_dmat(nK1, nff, &fix->KKf);
    blasfeo_allocate_dmat(nK1, dims->nx1, &fix->KKx);
    blasfeo_allocate_dmat(nK1, dims->nu , &fix->KKu);
    // fix->KKx = &out[dims->nx1 *dims->num_stages*nff];
    // fix->KKu = &out[(dims->nx1 *dims->num_stages)*(nff + dims->nx1)];
    blasfeo_pack_dmat(nK1, nff, out                                                        , dims->num_stages * dims->nx1, &fix->KKf, 0, 0);
    blasfeo_pack_dmat(nK1, dims->nx1, out + dims->nx1 *dims->num_stages*nff                , dims->num_stages * dims->nx1, &fix->KKx, 0, 0);
    blasfeo_pack_dmat(nK1, dims->nu, out + (dims->nx1 *dims->num_stages)*(nff + dims->nx1) , dims->num_stages * dims->nx1, &fix->KKu, 0, 0);
    printf("KKf Mat\n");
    blasfeo_print_dmat(nK1, nff, &fix->KKf, 0,0);
    printf("KKx Mat\n");
    blasfeo_print_dmat(nK1, dims->nx1, &fix->KKx, 0,0);
    printf("KKu Mat \n");
    blasfeo_print_dmat(nK1, dims->nu , &fix->KKu, 0,0);
    // printf("KKx Mat\n");
    // d_print_mat(dims->num_stages* dims->nx1, dims->nx1, fix->KKx, dims->nx1 *dims->num_stages);
}


// void gnsf_allocate_fixed( gnsf_dims *dims, gnsf_fixed *fix){ // Didnt really work..
//     fix->KKf = calloc(dims->nx1 * dims->num_stages * (nff + dims->nu + dims->nx1), sizeof(double));
//     fix->KKx = fix->KKf + (dims->nx1 * dims->num_stages *nff)* sizeof(double*);
//     fix->KKu = fix->KKx + dims->nx1 * dims->num_stages *dims->nx1;
// }

void gnsf_simulate( gnsf_dims *dims, gnsf_fixed *fix, gnsf_in *in, gnsf_out out)
{
    int nff = dims->n_out * dims->num_stages;
    int nK1 = dims->num_stages * dims->nx1;
    printf("GENERALIZED NONLINEAR STATIC FEEDBACK (GNSF) SIMULATION \n");
    print_gnsf_dims(dims);
    printf("KKu Mat \n");
    blasfeo_print_dmat(nK1, dims->nu , &fix->KKu, 0,0);

    double *res_in;
    double *K1_traj;
    double *ff_traj;
    int res_in_size = nff + dims->nx1 + dims->nu;
    int K1_traj_size = dims->num_stages * dims->nx1 * dims->num_steps;
    int ff_traj_size = nff * dims->num_steps;
    res_in  = (double*) calloc(res_in_size , sizeof(double));
    K1_traj = (double*) calloc(K1_traj_size, sizeof(double));
    ff_traj = (double*) calloc(ff_traj_size, sizeof(double));


    double *res_out;
    int res_out_size = nff * (1 + nff);
    res_out = (double*) calloc(res_out_size, sizeof(double));

    int *ipiv = (int*) calloc(nff, sizeof(int));
    
    // res_in[2+nff] = 0.8;

    // in->A_dt[0] = 1.0;
    // printf("%f",in->A_dt[0]);
    // printf("u = \t%f\n", in->u[0]);
    // printf("u = \t%f\n", in->u[1]);
    int newton_max = 1;
    struct blasfeo_dmat J_r_ff; // store the the jacobian of the residual w.r.t. ff
    struct blasfeo_dvec ff_val[dims->num_steps];
    struct blasfeo_dvec K1_val[dims->num_steps];
    struct blasfeo_dvec res_val;
    struct blasfeo_dvec u0;
    struct blasfeo_dvec x0_1;    
    blasfeo_allocate_dmat(nff, nff, &J_r_ff);
    blasfeo_allocate_dvec(nff, &res_val);
    blasfeo_allocate_dvec(dims->nu, &u0);
    blasfeo_pack_dvec(dims->nu, in->u, &u0, 0);
    blasfeo_allocate_dvec(dims->nx1, &x0_1);
    blasfeo_pack_dvec(dims->nx1, in->x, &x0_1, 0);  
    for (int ss = 0; ss < dims->num_steps; ss++) {
        blasfeo_allocate_dvec(nK1, &K1_val[ss]);
        blasfeo_allocate_dvec(nff, &ff_val[ss]);
    }

    for (int ss = 0; ss < 1; ss++) { // TODO: replace 1 with dim num_steps
        for (int iter = 0; iter < newton_max; iter++) { // NEWTON-ITERATION
            // set input for residual function
            blasfeo_unpack_dvec(nff, &ff_val[ss], 0, res_in);
            for (int i = 0; i<dims->nx1; i++) {
                res_in[i+nff] = in->x[i];
            }
            for (int i = 0; i<dims->nu; i++) {
                res_in[i+nff+dims->nx1] = in->u[i];
            }
            // evaluate residual and neccessary jacobians & pack into blasfeo mat/vec
            print_gnsf_res_in( dims, res_in );
            res_inc_Jff_wrapped(dims->nx1, dims->nu, dims->n_out, dims->num_stages, res_in, res_out, in->res_inc_Jff);
            // print_gnsf_res_out( *dims, res_out );
            blasfeo_pack_dvec(nff, res_out, &res_val, 0);
            // void blasfeo_pack_dmat(int m, int n, double *A, int lda, struct blasfeo_dmat *sB, int bi, int bj);
            blasfeo_pack_dmat(nff, nff, res_out+nff, nff, &J_r_ff, 0, 0); // pack residual result into blasfeo struct
            if (0); {
                printf("\nJ_r_ff = \n");
                blasfeo_print_dmat(nff, nff, &J_r_ff, 0,0);
                printf("\n residual value = \n");
                blasfeo_print_dvec(nff, &res_val, 0);
            }
            // // D <= lu( C ) ; no pivoting
            // void blasfeo_dgetrf_nopivot(int m, int n, struct blasfeo_dmat *sC, int ci, int cj, struct blasfeo_dmat *sD, int di, int dj);

            blasfeo_dgetrf_nopivot(nff, nff, &J_r_ff, 0,0, &J_r_ff, 0, 0); // invert J_r_ff

            // for (int i = 0; i < nff; i++) {
            //     printf("%f \t", ipiv[i]);
            // }
            // printf("\n After LU-fac; r_J_r_ffx1u = \n");
            // blasfeo_print_dmat(nff, 1+nff+dims->nx1+dims->nu, &r_J_r_ffx1u, 0,0);
            // z <= inv(A) * x, A (m)x(m) upper, not_transposed, not_unit
            blasfeo_dtrsv_unn(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            // z <= inv( A ) * x, A (m)x(m) lower, not_transposed, unit
            // void blasfeo_dtrsv_lnu(int m, struct blasfeo_dmat *sA, int ai, int aj, struct blasfeo_dvec *sx, int xi, struct blasfeo_dvec *sz, int zi);
            blasfeo_dtrsv_lnu(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            // z = y + alpha*x
            // void blasfeo_daxpy(int kmax, double alpha, struct blasfeo_dvec *sx, int xi, struct blasfeo_dvec *sy, int yi, struct blasfeo_dvec *sz, int zi);
            blasfeo_daxpy(nff, -1.0, &res_val, 0, &ff_val[ss], 0, &ff_val[ss], 0);
            // blasfeo_dtrsm_lunn(nff, 1, 1.0, &r_J_r_ffx1u, 0, 1, &r_J_r_ffx1u, 0, 0, &delta_ff, 0, 0);
            // D <= alpha * A^{-1} * B , with A lower triangular with unit diagonal
            // void blasfeo_dtrsm_llnu(int m, int n, double alpha, struct blasfeo_dmat *sA, int ai, int aj, struct blasfeo_dmat *sB, int bi, int bj, struct blasfeo_dmat *sD, int di, int dj);
            // blasfeo_dtrsm_llnu(nff, 1, 1.0, &r_J_r_ffx1u, 0, 1, &delta_ff, 0, 0, &delta_ff, 0, 0); // TODO put alpha = -1; but not implemented in blasfeo HP yet
            printf("\n ff =  \n");
            blasfeo_print_dvec(nff, &ff_val[ss], 0);
        }
        // K1_val = s.KKf * fftraj(:,ss) + s.KKu * u0 + s.KKx * x0_1;
        // z <= beta * y + alpha * A * x
// void blasfeo_dgemv_n(int m, int n, double alpha, struct blasfeo_dmat *sA, int ai, int aj, struct blasfeo_dvec *sx, int xi, double beta, struct blasfeo_dvec *sy, int yi, struct blasfeo_dvec *sz, int zi);
    blasfeo_dgemv_n(nK1, nff,       1.0, &fix->KKf, 0, 0, &ff_val[ss], 0, 0.0, &K1_val[ss], 0, &K1_val[ss], 0);

    printf("KKf Mat\n");
    blasfeo_print_dmat(nK1, nff, &fix->KKf, 0,0);
    printf("KKx Mat\n");
    blasfeo_print_dmat(nK1, dims->nx1, &fix->KKx, 0,0);
    printf("KKu Mat \n");
    blasfeo_print_dmat(nK1, dims->nu , &fix->KKu, 0,0);

    blasfeo_print_dvec(nK1, &K1_val[ss],0);
    blasfeo_print_dvec(dims->nu, &u0, 0);


    blasfeo_dgemv_n(nK1, dims->nu , 1.0, &fix->KKu, 0, 0, &u0        , 0, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
    blasfeo_dgemv_n(nK1, dims->nx1, 1.0, &fix->KKx, 0, 0, &x0_1, 0, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
    blasfeo_print_dvec(nK1, &K1_val[ss],0);



    }


// free memory
    for (int ss = 0; ss < dims->num_steps; ss++) {
        blasfeo_free_dvec(&ff_val[ss]);        
    }
    blasfeo_free_dmat(&J_r_ff);
    blasfeo_free_dvec(&res_val);
}