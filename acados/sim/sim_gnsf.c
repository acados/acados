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
#include "acados/utils/timing.h"

#include "acados/sim/sim_common.h"
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


int gnsf_calculate_workspace_size(gnsf_dims *dims, gnsf_opts* opts)
{
    int nx  = dims->nx;
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int num_stages = dims->num_stages;
    int num_steps = dims->num_steps;
    int nff = n_out * num_stages;

    int size = sizeof(gnsf_workspace);

    int res_in_size = nff + nx1 + nu;
    int res_out_size = nff * (nff + nx1 + nu); // size(out_res_inc_Jff) = nff* (1+nff), size(out_jac_res_ffx1u) = nff(nff+nx1+nu)
    int f_LO_in_size = 2*nx1 + nu + nz;
    int f_LO_out_size = nx2 * (1 + 2*nx1 + nu + nz);

    size += (res_in_size + res_out_size + f_LO_in_size + f_LO_out_size) * sizeof(double); // rhs_forw_in

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}


void *gnsf_cast_workspace(gnsf_dims* dims, void *raw_memory)
{
    int nx  = dims->nx;
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int num_stages = dims->num_stages;
    int num_steps = dims->num_steps;
    int nff = n_out * num_stages;

    int res_in_size = nff + nx1 + nu;
    int res_out_size = nff * (nff + nx1 + nu); // size(out_res_inc_Jff) = nff* (1+nff), size(out_jac_res_ffx1u) = nff(nff+nx1+nu)
    int f_LO_in_size = 2*nx1 + nu + nz;
    int f_LO_out_size = nx2 * (1 + 2*nx1 + nu + nz);

    char *c_ptr = (char *)raw_memory;
    gnsf_workspace *workspace = (gnsf_workspace *) c_ptr;
    c_ptr += sizeof(gnsf_workspace);
    align_char_to(8, &c_ptr);

    assign_double(res_in_size, &workspace->res_in, &c_ptr);
    assign_double(res_out_size, &workspace->res_out, &c_ptr);
    assign_double(f_LO_in_size, &workspace->f_LO_in, &c_ptr);
    assign_double(f_LO_out_size, &workspace->f_LO_out, &c_ptr);



    return (void *)workspace;
}

void gnsf_simulate( gnsf_dims *dims, gnsf_fixed *fix, gnsf_in *in, sim_out *out, gnsf_opts *opts, void *work_)
{
    acados_timer tot_timer, casadi_timer;
    acados_tic(&tot_timer);
    printf("GENERALIZED NONLINEAR STATIC FEEDBACK (GNSF) SIMULATION \n");
    print_gnsf_dims(dims);

    gnsf_workspace *workspace = (gnsf_workspace *) gnsf_cast_workspace(dims, work_);
    // helpful integers
    int nx  = dims->nx;
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int num_stages = dims->num_stages;
    int num_steps = dims->num_steps;

    int nff = n_out * num_stages;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;

    // allocate doubles
    int res_in_size = nff + nx1 + nu;
    int res_out_size = nff * (nff + nx1 + nu); // size(out_res_inc_Jff) = nff* (1+nff), size(out_jac_res_ffx1u) = nff(nff+nx1+nu)
    int f_LO_in_size = 2*nx1 + nu + nz;
    int f_LO_out_size = nx2 * (1 + 2*nx1 + nu + nz);

    int *ipiv = (int*) calloc(nff, sizeof(int));
    int newton_max = 10;

    // allocate blasfeo structures
    struct blasfeo_dmat J_r_ff; // store the the jacobian of the residual w.r.t. ff

    struct blasfeo_dmat J_r_x1u;  // needed for sensitivity propagation

    struct blasfeo_dmat dK1_dx1;
    struct blasfeo_dmat dK1_du;
    struct blasfeo_dmat dZ_dx1;
    struct blasfeo_dmat dZ_du;
    struct blasfeo_dmat aux_G2_x1;
    struct blasfeo_dmat aux_G2_u;
    struct blasfeo_dmat J_G2_K1;
    struct blasfeo_dmat dK2_dx1;
    struct blasfeo_dmat dK2_du;
    struct blasfeo_dmat dK2_dff;
    struct blasfeo_dmat dxf_dwn;
    struct blasfeo_dmat S_forw_new; // used to avoid side effects

    struct blasfeo_dmat S_forw;

    struct blasfeo_dmat f_LO_jac[num_steps];

    struct blasfeo_dvec ff_val[num_steps];
    struct blasfeo_dvec K1_val[num_steps];
    struct blasfeo_dvec x1_val[num_steps];
    struct blasfeo_dvec Z_val[num_steps];
    struct blasfeo_dvec f_LO_val[num_steps];
    struct blasfeo_dvec K2_val;
    struct blasfeo_dvec x0_traj;
    struct blasfeo_dvec res_val;
    struct blasfeo_dvec u0;
    struct blasfeo_dvec x0_1;
    struct blasfeo_dvec x0_2;

    struct blasfeo_dmat aux_G2_ff;
    struct blasfeo_dmat dPsi_dff;
    struct blasfeo_dmat dPsi_dx;
    struct blasfeo_dmat dPsi_du;

    double *res_in = workspace->res_in;
    double *res_out = workspace->res_out;
    double *f_LO_in = workspace->f_LO_in;
    double *f_LO_out = workspace->f_LO_out;

    if (opts->sens_forw || opts->sens_adj) {
        blasfeo_allocate_dmat(nK1, nx1, &dK1_dx1); //   dK1_dx1
        blasfeo_allocate_dmat(nK1, nu,  &dK1_du); //    dK1_du
        blasfeo_allocate_dmat(nZ, nx1, &dZ_dx1); //   dZ_dx1
        blasfeo_allocate_dmat(nZ, nu,  &dZ_du); //    dZ_du
        blasfeo_allocate_dmat(nK2, nx1,  &aux_G2_x1); //    aux_G2_x1
        blasfeo_allocate_dmat(nK2, nu,  &aux_G2_u); //    aux_G2_u
        blasfeo_allocate_dmat(nK2, nK1,  &J_G2_K1); //    J_G2_K1
        blasfeo_allocate_dmat(nK2, nx1,  &dK2_dx1); //    dK2_dx1
        blasfeo_allocate_dmat(nK2, nu ,  &dK2_du); //    dK2_du
        blasfeo_allocate_dmat(nK2, nff ,  &dK2_dff); //    dK2_dff
        blasfeo_allocate_dmat(nx, nx + nu,  &dxf_dwn); //    dxf_dwn
        blasfeo_allocate_dmat(nx, nx + nu,  &S_forw_new); //    S_forw_new
    }

    blasfeo_allocate_dmat(nK2, nff,  &aux_G2_ff); //    aux_G2_ff

    blasfeo_allocate_dvec((num_steps +1) * nx, &x0_traj);   // x0_traj
    blasfeo_allocate_dmat(nff, nff, &J_r_ff); //J_r_ff
    blasfeo_allocate_dmat(nff, nx1 + nu, &J_r_x1u); //    J_r_x1u

    blasfeo_allocate_dmat(nx, nx + nu,  &S_forw); //    S_forw

    blasfeo_allocate_dmat(nx, nff, &dPsi_dff);
    blasfeo_allocate_dmat(nx, nx, &dPsi_dx);
    blasfeo_allocate_dmat(nx, nu, &dPsi_du);

    blasfeo_allocate_dvec(nff, &res_val);   //res_val
    blasfeo_allocate_dvec(nu, &u0); // u0
    blasfeo_pack_dvec(nu, in->u, &u0, 0);
    blasfeo_allocate_dvec(nx1, &x0_1); // x0_1
    blasfeo_allocate_dvec(nx2, &x0_2); // x0_2
    blasfeo_pack_dvec(nx1, &in->x[0], &x0_1, 0);
    blasfeo_pack_dvec(nx2, &in->x[nx1], &x0_2, 0);
    blasfeo_allocate_dvec(nK2, &K2_val);
    
    blasfeo_pack_dvec(nx, &in->x[0], &x0_traj, 0);
    blasfeo_pack_dmat(nx, nx + nu, &in->S_forw[0], nx, &S_forw, 0, 0);

    for (int ss = 0; ss < num_steps; ss++) {
        blasfeo_allocate_dvec(nK1, &K1_val[ss]);
        blasfeo_allocate_dvec(nK1, &x1_val[ss]);
        blasfeo_allocate_dvec(nff, &ff_val[ss]);
        blasfeo_allocate_dvec(nZ,  &Z_val[ss]);
        blasfeo_allocate_dvec(nK2,  &f_LO_val[ss]);
        blasfeo_allocate_dmat(nK2, 2*nx1 + nu + nz,  &f_LO_jac[ss]);
    }
    out->info->ADtime = 0;

    for (int ss = 0; ss < num_steps; ss++) { // TODO: replace 1 with dim num_steps
        // todo: usage of x0_1/ x0_2 can be avoided  % Initialization inside
        blasfeo_dveccp(nx1, &x0_traj, nx * ss, &x0_1, 0);
        blasfeo_dveccp(nx2, &x0_traj, nx * ss+ nx1, &x0_2, 0);
        for (int iter = 0; iter < newton_max; iter++) { // NEWTON-ITERATION
            // set input for residual function
            blasfeo_unpack_dvec(nff, &ff_val[ss], 0, &res_in[0]);
            blasfeo_unpack_dvec(nx1, &x0_1, 0, &res_in[nff]);
            for (int i = 0; i<nu; i++) {
                res_in[i+nff+nx1] = in->u[i];
            }
            // evaluate residual and neccessary jacobians & pack into blasfeo mat/vec         
            acados_tic(&casadi_timer);
            res_inc_Jff_wrapped(nx1, nu, n_out, num_stages, res_in, res_out, in->res_inc_Jff);
            out->info->ADtime += acados_toc(&casadi_timer);

            blasfeo_pack_dvec(nff, &res_out[0], &res_val, 0);
            blasfeo_pack_dmat(nff, nff, &res_out[nff], nff, &J_r_ff, 0, 0); // pack residual result into blasfeo struct
            blasfeo_dgetrf_nopivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0); // factorize J_r_ff
            blasfeo_dtrsv_unn(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            blasfeo_dtrsv_lnu(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            blasfeo_daxpy(nff, -1.0, &res_val, 0, &ff_val[ss], 0, &ff_val[ss], 0);
        }
        // printf("\n ff =  \n");
        // blasfeo_print_exp_dvec(nff, &ff_val[ss], 0);
        // K1_val = s.KKf * fftraj(:,ss) + s.KKu * u0 + s.KKx * x0_1;
        blasfeo_dgemv_n(nK1, nff,       1.0, &fix->KKf, 0, 0, &ff_val[ss], 0, 0.0, &K1_val[ss], 0, &K1_val[ss], 0);
        blasfeo_dgemv_n(nK1, nu , 1.0, &fix->KKu, 0, 0, &u0        , 0, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
        blasfeo_dgemv_n(nK1, nx1, 1.0, &fix->KKx, 0, 0, &x0_1,       0, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
        // printf("\n K1_val =  \n");
        // blasfeo_print_exp_dvec(nK1, &K1_val[ss],0);

        blasfeo_dgemv_n(nZ, nff,       1.0, &fix->ZZf, 0, 0, &ff_val[ss], 0, 0.0, &Z_val[ss], 0, &Z_val[ss], 0);
        blasfeo_dgemv_n(nZ, nu , 1.0, &fix->ZZu, 0, 0, &u0        , 0, 1.0, &Z_val[ss], 0, &Z_val[ss], 0);
        blasfeo_dgemv_n(nZ, nx1, 1.0, &fix->ZZx, 0, 0, &x0_1,       0, 1.0, &Z_val[ss], 0, &Z_val[ss], 0);
        // printf("\n Z_val =  \n");
        // blasfeo_print_exp_dvec(nZ, &Z_val[ss],0);
        // printf("\n adress \n%p",(void*)&Z_val[0]);
        // build x1 stage values
        for (int ii = 0; ii < num_stages; ii++){
            blasfeo_daxpy(nx1, 0.0, &x1_val[ss], 0, &x0_1, 0, &x1_val[ss], nx1 * ii);
            for (int jj = 0; jj <num_stages; jj++) {
                blasfeo_daxpy(nx1, fix->A_dt[ii+num_stages*jj], &K1_val[ss], nx1*jj, &x1_val[ss], nx1*ii, &x1_val[ss], nx1*ii);
            }
        }
        // printf("x1_val = \n");
        // blasfeo_print_exp_dvec(nK1, &x1_val[ss], 0);

        // SIMULATE LINEAR OUTPUT SYSTEM
        for (int ii = 0; ii < num_stages; ii++) {
            blasfeo_unpack_dvec(nx1, &x1_val[ss], ii*nx1, &f_LO_in[0]);
            blasfeo_unpack_dvec(nx1, &K1_val[ss], ii*nx1, &f_LO_in[nx1]);
            blasfeo_unpack_dvec(nz,  &Z_val[ss] , ii*nz , &f_LO_in[2*nx1]);
            blasfeo_unpack_dvec(nu,  &u0        ,  0          , &f_LO_in[2*nx1 +nz]);
            // printf("f_LO_in = \n");
            // d_print_mat(f_LO_in_size, 1, &f_LO_in[0], f_LO_in_size);
            // d_print_mat_e(f_LO_in_size, 1, &f_LO_in[0], f_LO_in_size);
            // printf("f_LO_in = \n");
            // d_print_mat(f_LO_in_size,1,f_LO_in,f_LO_in_size);
            acados_tic(&casadi_timer);
            f_LO_inc_J_x1k1uz_wrapped(nx1, nz, f_LO_in, f_LO_out, in->f_LO_inc_J_x1k1uz);
            out->info->ADtime += acados_toc(&casadi_timer);
            // printf("f_LO_out= \n");
            // d_print_mat(f_LO_out_size,1, &f_LO_out[0],f_LO_out_size);
            // printf("f_LO_out_size= %d \n", f_LO_out_size);
            blasfeo_pack_dvec(nx2, &f_LO_out[0], &f_LO_val[ss], nx2 * ii);
            // printf("\n adress %p\n",(void*)&f_LO_jac[ss]);
            blasfeo_pack_dmat(nx2, 2*nx1 + nu + nz, &f_LO_out[nx2], nx2, &f_LO_jac[ss], nx2 * ii, 0); // NOTE: f_LO_jac has different sign compared to Matlab prototype
            blasfeo_dgemv_n(nx2, nx2, 1.0, &fix->ALO, 0, 0, &x0_2, 0, -1.0, &f_LO_val[ss], nx2 * ii, &f_LO_val[ss], nx2 * ii); // todo: repmat( - s.ALO * x0_2, q, 1); could be translated more efficient
        }
        // printf("f_LO = \n");
        // blasfeo_print_exp_dvec(num_stages * nx2, &f_LO_val[ss], 0);
        // blasfeo_print_dmat(    num_stages * nx2, 2*nx1 + nz + nu, &f_LO_jac[ss],0,0);

        // z <= beta * y + alpha * A * x
        blasfeo_dgemv_n( nK2, nK2, -1.0, &fix->M2inv, 0, 0, &f_LO_val[ss], 0, 0.0, &K2_val, 0, &K2_val, 0);
        // printf("K2_val = \n");
        // blasfeo_print_exp_dvec( nK2, &K2_val, 0);
        // Get simulation result
        blasfeo_daxpy(nx, 0.0, &x0_traj, 0, &x0_traj, nx * ss, &x0_traj, nx * (ss+1));
        for (int ii = 0; ii < num_stages; ii++) {
            blasfeo_daxpy(nx1, fix->b_dt[ii], &K1_val[ss], ii*nx1, &x0_traj, nx * (ss+1) , &x0_traj, nx * (ss+1));
            blasfeo_daxpy(nx2, fix->b_dt[ii], &K2_val    , ii*nx2, &x0_traj, nx1 + nx * (ss+1),  &x0_traj, nx1 + nx * (ss+1));
        }
        // blasfeo_print_exp_dvec(nx, &x0_traj, (ss+1) * nx);
        // set input for residual function
        blasfeo_unpack_dvec(nff, &ff_val[ss], 0, &res_in[0]);
        blasfeo_unpack_dvec(nx1, &x0_1, 0, &res_in[nff]);
        for (int i = 0; i<nu; i++) {
            res_in[i+nff+nx1] = in->u[i];
        }

        // printf("res_in = \n");
        // d_print_e_mat(res_in_size, 1, res_in, res_in_size);
        if (opts->sens_forw) {
            acados_tic(&casadi_timer);
            jac_res_ffx1u_wrapped(nx1, nu, n_out, num_stages, res_in, res_out, in->jac_res_ffx1u);
            out->info->ADtime += acados_toc(&casadi_timer);
            blasfeo_pack_dmat(nff, nff, &res_out[0], nff, &J_r_ff, 0, 0); // pack residual result into blasfeo struct
            blasfeo_pack_dmat(nff, nx1+ nu, &res_out[nff*nff], nff, &J_r_x1u, 0, 0); // pack residual result into blasfeo struct
            // blasfeo_print_exp_dmat(nff, nx1+ nu, &J_r_x1u, 0, 0);

            blasfeo_dgetrf_nopivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0); // factorize J_r_ff
            blasfeo_dtrsm_lunn(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
            blasfeo_dtrsm_llnu(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
            // blasfeo_dgemm_nn(nff, nx1 + nu, 1, 0.0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0, -1.0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0); // J_r_x1u = - J_r_x1u, because alpha=-1.0 not supported above;
            // printf("-dff_dx1u= \n");
            // blasfeo_print_exp_dmat(nff, nx1 + nu, &J_r_x1u, 0, 0);

            // D <= beta * C + alpha * A * B
            blasfeo_dgemm_nn(nK1, nx1, nff, -1.0, &fix->KKf, 0, 0, &J_r_x1u, 0,  0,        1.0, &fix->KKx, 0, 0, &dK1_dx1, 0, 0);
            blasfeo_dgemm_nn(nK1, nu,  nff, -1.0, &fix->KKf, 0, 0, &J_r_x1u, 0, nx1, 1.0, &fix->KKu, 0, 0, &dK1_du , 0, 0); // Blasfeo HP & Reference differ here

            blasfeo_dgemm_nn(nZ, nx1, nff, -1.0, &fix->ZZf, 0, 0, &J_r_x1u, 0, 0,         1.0, &fix->ZZx, 0, 0, &dZ_dx1, 0, 0);
            blasfeo_dgemm_nn(nZ, nu , nff, -1.0, &fix->ZZf, 0, 0, &J_r_x1u, 0, nx1, 1.0, &fix->ZZu, 0, 0, &dZ_du, 0, 0);
            // blasfeo_print_exp_dmat(nZ, nu, &dZ_du, 0, 0);
            // BUILD J_G2_wn, J_G2_K1    
            for (int ii = 0; ii < num_stages; ii++) {
                for (int jj = 0; jj < num_stages; jj++) {
                    blasfeo_dgecpsc( nx2, nx1, -fix->A_dt[ii+ jj*num_stages], &f_LO_jac[ss], ii*nx2, 0, &J_G2_K1, ii*nx2, jj*nx1);
                }
                blasfeo_dgead(nx2, nx1, 1.0, &f_LO_jac[ss], ii*nx2, nx1, &J_G2_K1, ii*nx2, ii*nx1);
                blasfeo_dgemm_nn( nx2, nx1, nz, 1.0, &f_LO_jac[ss], ii* nx2, 2 * nx1, &dZ_dx1, ii*nz, 0, 0.0, &aux_G2_x1, ii*nx2, 0, &aux_G2_x1, ii*nx2, 0);
                blasfeo_dgemm_nn( nx2, nu , nz, 1.0, &f_LO_jac[ss], ii* nx2, 2 * nx1, &dZ_du , ii*nz, 0, 0.0, &aux_G2_u,  ii*nx2, 0, &aux_G2_u,  ii*nx2, 0);
            } // TODO: aux_G2_x1u seem correct but should be tested with nontrivial values i.e. not zeros..

            // BUILD dK2_dwn // dK2_dx1
            blasfeo_dgemm_nn(nK2, nx1, nK1, 1.0, &J_G2_K1, 0, 0, &dK1_dx1, 0, 0, 1.0, &aux_G2_x1, 0, 0, &aux_G2_x1, 0, 0);
            blasfeo_dgead(nK2, nx1, -1.0, &f_LO_jac[ss], 0, 0, &aux_G2_x1, 0, 0);
            blasfeo_dgemm_nn(nK2, nx1, nK2, -1.0, &fix->M2inv, 0, 0, &aux_G2_x1, 0, 0, 0.0, &dK2_dx1, 0, 0, &dK2_dx1, 0, 0);
            // dK2_du
            blasfeo_dgemm_nn(nK2, nu, nK1, 1.0, &J_G2_K1, 0, 0, &dK1_du, 0, 0, 1.0, &aux_G2_u, 0, 0, &aux_G2_u, 0, 0);
            blasfeo_dgead(nK2, nu, -1.0, &f_LO_jac[ss], 0, 2*nx1 + nz, &aux_G2_u, 0, 0);
            blasfeo_dgemm_nn(nK2, nu, nK2, -1.0, &fix->M2inv, 0, 0, &aux_G2_u, 0, 0, 0.0, &dK2_du, 0, 0, &dK2_du, 0, 0);
            // BUILD dxf_dwn
            blasfeo_dgese(nx, nx + nu, 0.0, &dxf_dwn, 0, 0); // Initialize as unit matrix
            for (int ii = 0; ii < nx; ii++) {
                blasfeo_dgein1(1.0, &dxf_dwn, ii,ii);            
            }
            for (int ii = 0; ii < num_stages; ii++) {
                blasfeo_dgead(nx1, nx1, fix->b_dt[ii], &dK1_dx1, ii * nx1, 0, &dxf_dwn, 0, 0);  // derivatives w.r.t. x1
                blasfeo_dgead(nx2, nx1, fix->b_dt[ii], &dK2_dx1, ii * nx2, 0, &dxf_dwn, nx1, 0);

                blasfeo_dgead(nx2, nx2, fix->b_dt[ii], &fix->dK2_dx2, ii * nx2, 0, &dxf_dwn, nx1, nx1);  // derivatives w.r.t. x2
                
                blasfeo_dgead(nx1, nu, fix->b_dt[ii], &dK1_du, ii * nx1, 0, &dxf_dwn, 0, nx);  // derivatives w.r.t. u
                blasfeo_dgead(nx2, nu, fix->b_dt[ii], &dK2_du, ii * nx2, 0, &dxf_dwn, nx1, nx);
            }
            // printf("dxf_dwn = \n");
            // blasfeo_print_exp_dmat(nx, nx + nu, &dxf_dwn, 0, 0);
            blasfeo_dgemm_nn(nx, nx, nx, 1.0, &dxf_dwn, 0, 0, &S_forw, 0, 0, 0.0, &S_forw_new, 0, 0, &S_forw_new, 0, 0);
            blasfeo_dgemm_nn(nx, nu, nx, 1.0, &dxf_dwn, 0, 0, &S_forw, 0, nx, 1.0, &dxf_dwn, 0, nx, &S_forw_new, 0, nx);
            blasfeo_dgecp(nx, nx +nu, &S_forw_new, 0, 0, &S_forw, 0, 0);
        }
    }
    if (opts->sens_adj) {
        // ADJOINT SENSITIVITY PROPAGATION:
        for (int ss = num_steps-1; ss >= num_steps-1; ss--) {
            blasfeo_dveccp(nx1, &x0_traj, nx * ss, &x0_1, 0);
            for (int ii = 0; ii < num_stages; ii++) {
                blasfeo_dgemm_nn(nx2, nff      , nz, -1.0, &f_LO_jac[ss], nx2 * ii, 2*nx1, &fix->ZZf, ii* nz, 0, 0.0, &fix->ZZf, 0, 0, &aux_G2_ff, ii * nx2, 0);
                blasfeo_dgemm_nn(nx2, nx1, nz, -1.0, &f_LO_jac[ss], nx2 * ii, 2*nx1, &fix->ZZx, ii* nz, 0, 0.0, &fix->ZZx, 0, 0, &aux_G2_x1, ii * nx2, 0);
                blasfeo_dgemm_nn(nx2, nu, nz, -1.0, &f_LO_jac[ss], nx2 * ii, 2*nx1, &fix->ZZu, ii* nz, 0, 0.0, &fix->ZZu, 0, 0, &aux_G2_u, ii * nx2, 0);
                for (int jj = 0; jj < num_stages; jj++) {
                    blasfeo_dgecpsc(nx2, nx1, -fix->A_dt[ii+jj*num_stages], &f_LO_jac[ss], nx2 * ii, 0, &J_G2_K1, ii*nx2, jj*nx1);
                }
                blasfeo_dgead(nx2, nx1, -1.0, &f_LO_jac[ss], nx2*ii, nx1, &J_G2_K1, nx2*ii, nx1*ii);
            }
            // printf("J_G2_K1 = \n");
            // blasfeo_print_exp_dmat(nK2,nK1,&J_G2_K1,0,0);
            blasfeo_dgemm_nn(nK2, nff, nK1, 1.0, &J_G2_K1, 0, 0, &fix->KKf, 0, 0, 1.0, &aux_G2_ff, 0, 0,&aux_G2_ff, 0, 0);
            blasfeo_dgemm_nn(nK2, nx1, nK1, 1.0, &J_G2_K1, 0, 0, &fix->KKx, 0, 0, 1.0, &aux_G2_x1, 0, 0,&aux_G2_x1, 0, 0);
            blasfeo_dgemm_nn(nK2, nu , nK1, 1.0, &J_G2_K1, 0, 0, &fix->KKu, 0, 0, 1.0, &aux_G2_u , 0, 0,&aux_G2_u , 0, 0);

            blasfeo_dgead(nK2, nx1, -1.0, &f_LO_jac[ss], 0, 0, &aux_G2_x1, 0, 0); //TODO: stattdessen vorher kopieren und dann addieren oben möglich
            blasfeo_dgead(nK2, nu , -1.0, &f_LO_jac[ss], 0, 2*nx1 + nz, &aux_G2_u, 0, 0);

            blasfeo_dgemm_nn(nK2, nff, nK2, -1.0, &fix->M2inv, 0, 0, &aux_G2_ff, 0, 0, 0.0, &dK2_dff, 0, 0, &dK2_dff, 0, 0);
            blasfeo_dgemm_nn(nK2, nx1, nK2, -1.0, &fix->M2inv, 0, 0, &aux_G2_x1, 0, 0, 0.0, &dK2_dx1, 0, 0, &dK2_dx1, 0, 0);
            blasfeo_dgemm_nn(nK2, nu , nK2, -1.0, &fix->M2inv, 0, 0, &aux_G2_u , 0, 0, 0.0, &dK2_du, 0, 0, &dK2_du, 0, 0);

            blasfeo_dgese(nx, nff, 0.0, &dPsi_dff, 0, 0); // initialize dPsi_d.. 
            blasfeo_dgese(nx, nx, 0.0, &dPsi_dx, 0, 0);
            blasfeo_dgese(nx, nu, 0.0, &dPsi_du, 0, 0);
            for (int ii = 0; ii < nx; ii++) {
                blasfeo_dgein1(1.0, &dPsi_dx, ii, ii);
            }
            // compute dPsi_d..
            for (int ii = 0; ii < num_stages; ii++) {
                blasfeo_dgead(nx1, nff, fix->b_dt[ii], &fix->KKf, ii*nx1, 0, &dPsi_dff, 0, 0);
                blasfeo_dgead(nx1, nx1, fix->b_dt[ii], &fix->KKx, ii*nx1, 0, &dPsi_dx, 0, 0);
                blasfeo_dgead(nx1, nu, fix->b_dt[ii], &fix->KKu, ii*nx1, 0, &dPsi_du, 0, 0);
                
                blasfeo_dgead(nx2, nff, fix->b_dt[ii], &dK2_dff, ii*nx2, 0, &dPsi_dff, nx1, 0);
                blasfeo_dgead(nx2, nx1, fix->b_dt[ii], &dK2_dx1, ii*nx2, 0, &dPsi_dx, nx1, 0);
                blasfeo_dgead(nx2, nx2, fix->b_dt[ii], &fix->dK2_dx2, ii*nx2, 0, &dPsi_dx, nx1, nx1);
                blasfeo_dgead(nx2, nu, fix->b_dt[ii], &dK2_du, ii*nx2, 0, &dPsi_du, nx1, 0);            
            }
            acados_tic(&casadi_timer);
            jac_res_ffx1u_wrapped(nx1, nu, n_out, num_stages, res_in, res_out, in->jac_res_ffx1u);
            out->info->ADtime += acados_toc(&casadi_timer);
            blasfeo_pack_dmat(nff, nff, &res_out[0], nff, &J_r_ff, 0, 0); // pack residual result into blasfeo struct
            blasfeo_pack_dmat(nff, nx1+ nu, &res_out[nff*nff], nff, &J_r_x1u, 0, 0); // pack residual result into blasfeo struct

            blasfeo_dgetrf_nopivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0); // factorize J_r_ff
            blasfeo_dtrsm_lunn(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
            // blasfeo_dtrsm_llnu(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
        }
    }
    printf("forw_Sensitivities = \n");
    blasfeo_print_exp_dmat(nx, nx + nu, &S_forw, 0, 0);

    printf("x0_traj last:\n");
    blasfeo_print_exp_dvec(nx , &x0_traj, num_steps*nx);

    out->info->CPUtime = acados_toc(&tot_timer);
    printf("tot_time = %f\n", out->info->CPUtime);

// free memory
    
    for (int ss = 0; ss < num_steps; ss++) {
        blasfeo_free_dvec(&ff_val[ss]);
        blasfeo_free_dvec(&x1_val[ss]);
        blasfeo_free_dvec(&K1_val[ss]);
        blasfeo_free_dvec(&Z_val[ss]);
        blasfeo_free_dvec(&f_LO_val[ss]);
        blasfeo_free_dmat(&f_LO_jac[ss]);
    }
    blasfeo_free_dmat(&J_r_ff);
    blasfeo_free_dmat(&J_r_x1u);

    if (opts->sens_forw || opts->sens_adj) {
        // blasfeo_free_dmat(&dK1_dx1);
        blasfeo_free_dmat(&dK1_du);
        blasfeo_free_dmat(&dZ_dx1);
        blasfeo_free_dmat(&dZ_du);
        blasfeo_free_dmat(&aux_G2_x1);
        blasfeo_free_dmat(&aux_G2_u);
        blasfeo_free_dmat(&J_G2_K1);
        blasfeo_free_dmat(&dK2_dx1);
        blasfeo_free_dmat(&dK2_du);
        blasfeo_free_dmat(&dK2_dff);
        blasfeo_free_dmat(&dxf_dwn);
        blasfeo_free_dmat(&S_forw_new);
    }
    blasfeo_free_dvec(&res_val);
    blasfeo_free_dvec(&u0);
    blasfeo_free_dvec(&x0_1);
    blasfeo_free_dvec(&x0_2);
}