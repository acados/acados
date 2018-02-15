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

void gnsf_get_ZZ_mat(gnsf_dims *dims, gnsf_fixed *fix, casadi_function_t KK_mat_fun)
{
    int nZ = dims->nz * dims->num_stages;
    int nff = dims->num_stages * dims->n_out;
    double *out;
    out = (double*) calloc((nZ) * (nff + dims->nx1 + dims->nu),sizeof(double));
    // export_from_ML_wrapped(fix->KKf, fix->KKf, KK_mat_fun);
    export_from_ML_wrapped(out, out, KK_mat_fun);
    blasfeo_allocate_dmat(nZ, nff      , &fix->ZZf);
    blasfeo_allocate_dmat(nZ, dims->nx1, &fix->ZZx);
    blasfeo_allocate_dmat(nZ, dims->nu , &fix->ZZu);
    blasfeo_pack_dmat(nZ, nff, out                                                        , nZ, &fix->ZZf, 0, 0);
    blasfeo_pack_dmat(nZ, dims->nx1, out + dims->nx1 *dims->num_stages*nff                , nZ, &fix->ZZx, 0, 0);
    blasfeo_pack_dmat(nZ, dims->nu, out + (dims->nx1 *dims->num_stages)*(nff + dims->nx1) , nZ, &fix->ZZu, 0, 0);
    printf("ZZf Mat\n");
    blasfeo_print_dmat(nZ, nff, &fix->ZZf, 0,0);
    printf("ZZx Mat\n");
    blasfeo_print_dmat(nZ, dims->nx1, &fix->ZZx, 0,0);
    printf("ZZu Mat \n");
    blasfeo_print_dmat(nZ, dims->nu , &fix->ZZu, 0,0);
}

void gnsf_get_ALO_M2_dK2dx2(gnsf_dims *dims, gnsf_fixed *fix, casadi_function_t ALO_M2_dK2dx2_fun)
{
    int nK2 = dims->nx2 * dims->num_stages;
    int nff = dims->num_stages * dims->n_out;
    double *out;
    out = (double*) calloc(nK2 * (nK2 + dims->nx2) + dims->nx2 * dims->nx2, sizeof(double));
    export_from_ML_wrapped(out, out, ALO_M2_dK2dx2_fun);
    blasfeo_allocate_dmat(dims->nx2, dims->nx2    , &fix->ALO);
    blasfeo_allocate_dmat( nK2, nK2, &fix->M2inv);
    blasfeo_allocate_dmat( nK2, dims->nx2 , &fix->dK2_dx2);
    blasfeo_pack_dmat(dims->nx2, dims->nx2 , &out[0], dims->nx2, &fix->ALO, 0, 0);
    blasfeo_pack_dmat( nK2, nK2, &out[dims->nx2 *dims->nx2], nK2, &fix->M2inv, 0, 0);
    blasfeo_pack_dmat( nK2, dims->nx2, &out[dims->nx2 *dims->nx2 + nK2*nK2] , nK2, &fix->dK2_dx2, 0, 0);
    printf("ALO Mat\n");
    blasfeo_print_dmat(dims->nx2, dims->nx2, &fix->ALO, 0,0);
    printf("M2inv\n");
    blasfeo_print_dmat(nK2, nK2, &fix->M2inv, 0,0);
    printf("dK2_dx2 \n");
    blasfeo_print_dmat(nK2, dims->nx2 , &fix->dK2_dx2, 0,0);
}

void gnsf_get_butcher(gnsf_dims* dims, gnsf_fixed *fix, casadi_function_t Butcher_fun)
{
    double *out;
    out = (double*) calloc(dims->num_stages * (2+ dims->num_stages),sizeof(double));
    export_from_ML_wrapped(out, out, Butcher_fun);
    fix->A_dt = &out[0];
    fix->b_dt = &out[dims->num_stages * dims->num_stages];
    fix->c    = &out[dims->num_stages * (dims->num_stages+1)];
    d_print_e_mat(dims->num_stages, dims->num_stages, fix->A_dt, dims->num_stages);
}

void gnsf_simulate( gnsf_dims *dims, gnsf_fixed *fix, gnsf_in *in, gnsf_out out)
{
    printf("GENERALIZED NONLINEAR STATIC FEEDBACK (GNSF) SIMULATION \n");
    print_gnsf_dims(dims);

    // helpful integers
    int nff = dims->n_out * dims->num_stages;
    int nK1 = dims->num_stages * dims->nx1;
    int nK2 = dims->num_stages * dims->nx2;
    int nZ  = dims->num_stages * dims->nz;

    // allocate doubles 
    double *res_in;
    double *res_out;
    double *f_LO_in;
    double *f_LO_out;

    int res_in_size = nff + dims->nx1 + dims->nu;
    int res_out_size = nff * (nff + dims->nx1 + dims->nu); // size(out_res_inc_Jff) = nff* (1+nff), size(out_jac_res_ffx1u) = nff(nff+nx1+nu)
    int f_LO_in_size = 2*dims->nx1 + dims->nu + dims->nz;
    int f_LO_out_size = dims->nx2 * (1 + 2*dims->nx1 + dims->nu + dims->nz);

    res_in   = (double*) calloc(res_in_size  , sizeof(double));
    res_out  = (double*) calloc(res_out_size , sizeof(double));
    f_LO_in  = (double*) calloc(f_LO_in_size , sizeof(double));
    f_LO_out = (double*) calloc(f_LO_out_size, sizeof(double));

    int *ipiv = (int*) calloc(nff, sizeof(int));
    int newton_max = 10;

    // allocate blasfeo structures
    struct blasfeo_dmat J_r_ff; // store the the jacobian of the residual w.r.t. ff
    struct blasfeo_dmat J_r_x1u;
    struct blasfeo_dmat dK1_dx1;
    struct blasfeo_dmat dK1_du;
    struct blasfeo_dmat dZ_dx1;
    struct blasfeo_dmat dZ_du;
    struct blasfeo_dmat f_LO_jac[dims->num_steps];

    struct blasfeo_dvec ff_val[dims->num_steps];
    struct blasfeo_dvec K1_val[dims->num_steps];
    struct blasfeo_dvec x1_val[dims->num_steps];
    struct blasfeo_dvec Z_val[dims->num_steps];
    struct blasfeo_dvec f_LO_val[dims->num_steps];
    struct blasfeo_dvec K2_val;
    struct blasfeo_dvec x0_traj;
    struct blasfeo_dvec res_val;
    struct blasfeo_dvec u0;
    struct blasfeo_dvec x0_1;
    struct blasfeo_dvec x0_2;

    blasfeo_allocate_dvec((dims->num_steps +1) * dims->nx, &x0_traj);   // x0_traj
    blasfeo_allocate_dmat(nff, nff, &J_r_ff); //J_r_ff
    blasfeo_allocate_dmat(nff, dims->nx1 + dims->nu, &J_r_x1u); //    J_r_x1u
    
    blasfeo_allocate_dmat(nK1, dims->nx1, &dK1_dx1); //   dK1_dx1
    blasfeo_allocate_dmat(nK1, dims->nu,  &dK1_du); //    dK1_du
    blasfeo_allocate_dmat(nZ, dims->nx1, &dZ_dx1); //   dZ_dx1
    blasfeo_allocate_dmat(nZ, dims->nu,  &dZ_du); //    dZ_du

    blasfeo_allocate_dvec(nff, &res_val);   //res_val
    blasfeo_allocate_dvec(dims->nu, &u0); // u0
    blasfeo_pack_dvec(dims->nu, in->u, &u0, 0);
    blasfeo_allocate_dvec(dims->nx1, &x0_1); // x0_1
    blasfeo_allocate_dvec(dims->nx2, &x0_2); // x0_2
    blasfeo_pack_dvec(dims->nx1, &in->x[0], &x0_1, 0);
    blasfeo_pack_dvec(dims->nx2, &in->x[dims->nx1], &x0_2, 0);
    blasfeo_allocate_dvec(nK2, &K2_val);
    
    blasfeo_pack_dvec(dims->nx, &in->x[0], &x0_traj, 0);

    // printf("x0_1 = \n");
    // blasfeo_print_dvec(dims->nx1, &x0_1, 0);
    // printf("x0_2 = \n");
    // blasfeo_print_dvec(dims->nx2, &x0_2, 0);

    for (int ss = 0; ss < dims->num_steps; ss++) {
        blasfeo_allocate_dvec(nK1, &K1_val[ss]);
        blasfeo_allocate_dvec(nK1, &x1_val[ss]);
        blasfeo_allocate_dvec(nff, &ff_val[ss]);
        blasfeo_allocate_dvec(nZ,  &Z_val[ss]);
        blasfeo_allocate_dvec(nK2,  &f_LO_val[ss]);
        blasfeo_allocate_dmat(nK2, 2*dims->nx1 + dims->nu + dims->nz,  &f_LO_jac[ss]);
    }

    for (int ss = 0; ss < 1; ss++) { // TODO: replace 1 with dim num_steps
        // todo: usage of x0_1/ x0_2 can be avoided  % Initialization inside
        blasfeo_daxpy(dims->nx1, 0.0, &x0_traj, 0, &x0_traj, dims->nx * ss, &x0_1, 0);
        blasfeo_daxpy(dims->nx2, 0.0, &x0_traj, 0, &x0_traj, dims->nx * ss + dims->nx1, &x0_2, 0);
        for (int iter = 0; iter < newton_max; iter++) { // NEWTON-ITERATION
            // set input for residual function
            blasfeo_unpack_dvec(nff, &ff_val[ss], 0, &res_in[0]);
            blasfeo_unpack_dvec(dims->nx1, &x0_1, 0, &res_in[nff]);
            for (int i = 0; i<dims->nu; i++) {
                res_in[i+nff+dims->nx1] = in->u[i];
            }
            // evaluate residual and neccessary jacobians & pack into blasfeo mat/vec
            // print_gnsf_res_in( dims, res_in );
                printf("res_in = \n");
                d_print_e_mat(res_in_size, 1, res_in, res_in_size);
                
            res_inc_Jff_wrapped(dims->nx1, dims->nu, dims->n_out, dims->num_stages, res_in, res_out, in->res_inc_Jff);
            blasfeo_pack_dvec(nff, &res_out[0], &res_val, 0);
            blasfeo_pack_dmat(nff, nff, &res_out[nff], nff, &J_r_ff, 0, 0); // pack residual result into blasfeo struct
                // print_gnsf_res_out( *dims, res_out );
                printf("\nJ_r_ff = \n");
                blasfeo_print_exp_dmat(nff, nff, &J_r_ff, 0,0);
                // printf("\n residual value = \n");
                // blasfeo_print_dvec(nff, &res_val, 0);
            // // D <= lu( C ) ; no pivoting
            // void blasfeo_dgetrf_nopivot(int m, int n, struct blasfeo_dmat *sC, int ci, int cj, struct blasfeo_dmat *sD, int di, int dj);
            printf("Jacobian= \n");
            blasfeo_print_exp_dmat(nff, nff, &J_r_ff, 0, 0);
            blasfeo_dgetrf_nopivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0); // invert J_r_ff
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
        }
        printf("\n ff =  \n");
        blasfeo_print_exp_dvec(nff, &ff_val[ss], 0);
        // K1_val = s.KKf * fftraj(:,ss) + s.KKu * u0 + s.KKx * x0_1;
        // z <= beta * y + alpha * A * x
// void blasfeo_dgemv_n(int m, int n, double alpha, struct blasfeo_dmat *sA, int ai, int aj, struct blasfeo_dvec *sx, int xi, double beta, struct blasfeo_dvec *sy, int yi, struct blasfeo_dvec *sz, int zi);
        blasfeo_dgemv_n(nK1, nff,       1.0, &fix->KKf, 0, 0, &ff_val[ss], 0, 0.0, &K1_val[ss], 0, &K1_val[ss], 0);
        blasfeo_dgemv_n(nK1, dims->nu , 1.0, &fix->KKu, 0, 0, &u0        , 0, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
        blasfeo_dgemv_n(nK1, dims->nx1, 1.0, &fix->KKx, 0, 0, &x0_1,       0, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
        printf("\n K1_val =  \n");
        blasfeo_print_exp_dvec(nK1, &K1_val[ss],0);

        blasfeo_dgemv_n(nZ, nff,       1.0, &fix->ZZf, 0, 0, &ff_val[ss], 0, 0.0, &Z_val[ss], 0, &Z_val[ss], 0);
        blasfeo_dgemv_n(nZ, dims->nu , 1.0, &fix->ZZu, 0, 0, &u0        , 0, 1.0, &Z_val[ss], 0, &Z_val[ss], 0);
        blasfeo_dgemv_n(nZ, dims->nx1, 1.0, &fix->ZZx, 0, 0, &x0_1,       0, 1.0, &Z_val[ss], 0, &Z_val[ss], 0);
        printf("\n Z_val =  \n");
        blasfeo_print_exp_dvec(nZ, &Z_val[ss],0);
        // printf("\n adress \n%p",(void*)&Z_val[0]);
        // build x1 stage values
        for (int ii = 0; ii < dims->num_stages; ii++){
            blasfeo_daxpy(dims->nx1, 0.0, &x1_val[ss], 0, &x0_1, 0, &x1_val[ss], dims->nx1 * ii);
            for (int jj = 0; jj <dims->num_stages; jj++) {
                blasfeo_daxpy(dims->nx1, fix->A_dt[ii+dims->num_stages*jj], &K1_val[ss], dims->nx1*jj, &x1_val[ss], dims->nx1*ii, &x1_val[ss], dims->nx1*ii);
            }
        }
        printf("x1_val = \n");
        blasfeo_print_exp_dvec(nK1, &x1_val[ss], 0);

        // SIMULATE LINEAR OUTPUT SYSTEM
        // printf("\n adress \n%p",(void*)&Z_val[0]);
        // printf("\n nx1= %d",dims->nx1);
        // printf("\n adress %p\n",(void*)&f_LO_jac[ss]);
        for (int ii = 0; ii < dims->num_stages; ii++) {
            blasfeo_unpack_dvec(dims->nx1, &x1_val[ss], ii*dims->nx1, &f_LO_in[0]);
            blasfeo_unpack_dvec(dims->nx1, &K1_val[ss], ii*dims->nx1, &f_LO_in[dims->nx1]);
            blasfeo_unpack_dvec(dims->nz,  &Z_val[ss] , ii*dims->nz , &f_LO_in[2*dims->nx1]);
            blasfeo_unpack_dvec(dims->nu,  &u0        ,  0          , &f_LO_in[2*dims->nx1 +dims->nz]);
            // printf("f_LO_in = \n");
            // d_print_mat(f_LO_in_size, 1, &f_LO_in[0], f_LO_in_size);
            // d_print_mat_e(f_LO_in_size, 1, &f_LO_in[0], f_LO_in_size);
            // printf("f_LO_in = \n");
            // d_print_mat(f_LO_in_size,1,f_LO_in,f_LO_in_size);
            f_LO_inc_J_x1k1uz_wrapped(dims->nx1, dims->nz, f_LO_in, f_LO_out, in->f_LO_inc_J_x1k1uz);
            // printf("f_LO_out= \n");
            // d_print_mat(f_LO_out_size,1, &f_LO_out[0],f_LO_out_size);
            // printf("f_LO_out_size= %d \n", f_LO_out_size);
            blasfeo_pack_dvec(dims->nx2, &f_LO_out[0], &f_LO_val[ss], dims->nx2 * ii);
            // printf("\n adress %p\n",(void*)&f_LO_jac[ss]);
            blasfeo_pack_dmat(dims->nx2, 2*dims->nx1 + dims->nu + dims->nz, &f_LO_out[dims->nx2], dims->nx2, &f_LO_jac[ss], dims->nx2 * ii, 0);
            blasfeo_dgemv_n(dims->nx2, dims->nx2, 1.0, &fix->ALO, 0, 0, &x0_2, 0, -1.0, &f_LO_val[ss], dims->nx2 * ii, &f_LO_val[ss], dims->nx2 * ii); // todo: repmat( - s.ALO * x0_2, q, 1); could be translated more efficient
        }
        // blasfeo_print_exp_dvec(dims->num_stages * dims->nx2, &f_LO_val[ss], 0);
        // blasfeo_print_dmat(    dims->num_stages * dims->nx2, 2*dims->nx1 + dims->nz + dims->nu, &f_LO_jac[ss],0,0);

        // z <= beta * y + alpha * A * x
        blasfeo_dgemv_n( nK2, nK2, -1.0, &fix->M2inv, 0, 0, &f_LO_val[ss], 0, 0.0, &K2_val, 0, &K2_val, 0);
        // printf("K2_val = \n");
        // blasfeo_print_exp_dvec( nK2, &K2_val, 0);
        // Get simulation result
        blasfeo_daxpy(dims->nx, 0.0, &x0_traj, 0, &x0_traj, dims->nx * ss, &x0_traj, dims->nx * (ss+1));
        for (int ii = 0; ii < dims->num_stages; ii++) {
            blasfeo_daxpy(dims->nx1, fix->b_dt[ii], &K1_val[ss], ii*dims->nx1, &x0_traj, dims->nx * (ss+1) , &x0_traj, dims->nx * (ss+1));
            blasfeo_daxpy(dims->nx2, fix->b_dt[ii], &K2_val    , ii*dims->nx2, &x0_traj, dims->nx1 + dims->nx * (ss+1),  &x0_traj, dims->nx1 + dims->nx * (ss+1));
        }
        // blasfeo_print_exp_dvec(dims->nx, &x0_traj, (ss+1) * dims->nx);
        // set input for residual function
        blasfeo_unpack_dvec(nff, &ff_val[ss], 0, &res_in[0]);
        blasfeo_unpack_dvec(dims->nx1, &x0_1, 0, &res_in[nff]);
        for (int i = 0; i<dims->nu; i++) {
            res_in[i+nff+dims->nx1] = in->u[i];
        }

        printf("res_in = \n");
        d_print_e_mat(res_in_size, 1, res_in, res_in_size);

        jac_res_ffx1u_wrapped(dims->nx1, dims->nu, dims->n_out, dims->num_stages, res_in, res_out, in->jac_res_ffx1u);
        blasfeo_pack_dmat(nff, nff, &res_out[0], nff, &J_r_ff, 0, 0); // pack residual result into blasfeo struct
        blasfeo_pack_dmat(nff, dims->nx1+ dims->nu, &res_out[nff*nff], nff, &J_r_x1u, 0, 0); // pack residual result into blasfeo struct
        // blasfeo_print_exp_dmat(nff, dims->nx1+ dims->nu, &J_r_x1u, 0, 0);

        blasfeo_dgetrf_nopivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0); // invert J_r_ff
        blasfeo_dtrsm_lunn(nff, dims->nx1 + dims->nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
        blasfeo_dtrsm_llnu(nff, dims->nx1 + dims->nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
        // blasfeo_dgemm_nn(nff, dims->nx1 + dims->nu, 1, 0.0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0, -1.0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0); // J_r_x1u = - J_r_x1u, because alpha=-1.0 not supported above;
        // printf("-dff_dx1u= \n");
        // blasfeo_print_exp_dmat(nff, dims->nx1 + dims->nu, &J_r_x1u, 0, 0);

        // D <= beta * C + alpha * A * B
        blasfeo_dgemm_nn(nK1, dims->nx1, nff, -1.0, &fix->KKf, 0, 0, &J_r_x1u, 0,  0,        1.0, &fix->KKx, 0, 0, &dK1_dx1, 0, 0);
        blasfeo_dgemm_nn(nK1, dims->nu,  nff, -1.0, &fix->KKf, 0, 0, &J_r_x1u, 0, dims->nx1, 1.0, &fix->KKu, 0, 0, &dK1_du , 0, 0);

        blasfeo_dgemm_nn(nZ, dims->nx1, nff, -1.0, &fix->ZZf, 0, 0, &J_r_x1u, 0, 0,         1.0, &fix->ZZx, 0, 0, &dZ_dx1, 0, 0);
        blasfeo_dgemm_nn(nZ, dims->nu, nff, -1.0,  &fix->ZZf, 0, 0, &J_r_x1u, 0, dims->nx1, 1.0, &fix->ZZu, 0, 0, &dZ_du, 0, 0);

        // printf("dK1_dx1 = \n");
        // blasfeo_print_exp_dmat(nK1, dims->nx1, &dK1_dx1, 0, 0);
        // printf("dK1_du = \n");
        // blasfeo_print_exp_dmat(nK1, dims->nu, &dK1_du, 0, 0);

        // printf("dZ_dx1 = \n");
        // blasfeo_print_exp_dmat(nZ, dims->nx1, &dZ_dx1, 0, 0);
        // printf("dZ_du = \n"); // TODO: check this out, is it wrong?! error w.r.t. matlab solution is E-9, actually i double checked...



    }
    // printf("x0_traj 1st:\n");
    // blasfeo_print_exp_dvec(dims->nx , &x0_traj, dims->nx);
    // printf("x0_traj last:\n");
    // blasfeo_print_exp_dvec(dims->nx , &x0_traj, dims->num_steps*dims->nx);

// free memory
    for (int ss = 0; ss < dims->num_steps; ss++) {
        blasfeo_free_dvec(&ff_val[ss]);
        blasfeo_free_dvec(&x1_val[ss]);
        blasfeo_free_dvec(&K1_val[ss]);
        blasfeo_free_dvec(&Z_val[ss]);
        blasfeo_free_dvec(&f_LO_val[ss]);
        blasfeo_free_dmat(&f_LO_jac[ss]);
    }
    blasfeo_free_dmat(&J_r_ff);
    blasfeo_free_dmat(&J_r_x1u);
    blasfeo_free_dmat(&dK1_dx1);
    blasfeo_free_dmat(&dK1_du);
    blasfeo_free_dmat(&dZ_dx1);
    blasfeo_free_dmat(&dZ_du);


    blasfeo_free_dvec(&res_val);
    blasfeo_free_dvec(&u0);
    blasfeo_free_dvec(&x0_1);
    blasfeo_free_dvec(&x0_2);

}