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

void print_gnsf_res_in( gnsf_dims dims, double *res_in )
{
    printf("res_in ff = \n");    
    for (int i=0; i<dims.n_out *dims.num_stages; i++)
        printf("\t%f\n", res_in[i]);
    printf("res_in x1 = \n");
    for (int i=dims.n_out *dims.num_stages; i<dims.n_out *dims.num_stages + dims.nx1; i++)
        printf("\t%f\n", res_in[i]);
    printf("res_in u = \n");
    for (int i=dims.n_out *dims.num_stages + dims.nx1; i<dims.n_out *dims.num_stages + dims.nx1 + dims.nu; i++)
        printf("\t%f\n", res_in[i]);   
}


void print_gnsf_res_out( gnsf_dims dims, double *res_out )
{
    printf("res_out res_val = \n");    
    for (int i=0; i<dims.nff; i++)
        printf("\t%5.8f\n", res_out[i]);
    printf("\nres_out J_res_ff = \n");
    for (int i=0; i<dims.nff; i++){
        for (int j=1; j<dims.nff+1; j++){
            printf("\t %5.8f", res_out[i+j*dims.nff]);}
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
    dims->nff = dims->n_out * dims->num_stages;
}

void gnsf_get_KK_mat(gnsf_dims *dims, gnsf_fixed *fix, casadi_function_t KK_mat_fun)
{
    double *out;
    out = (double*) calloc((dims->nx1 * dims->num_stages) * (dims->nff + dims->nx1 + dims->nu),sizeof(double));
    // export_from_ML_wrapped(fix->KKf, fix->KKf, KK_mat_fun);
    export_from_ML_wrapped(out, out, KK_mat_fun);
    // for (int i=0; i<dims->nx1 *dims->num_stages*dims->nff; i++)
    // fix->KKf[i] = out[i];
    fix->KKf = &out[0];
    fix->KKx = &out[dims->nx1 *dims->num_stages*dims->nff];
    fix->KKu = &out[(dims->nx1 *dims->num_stages)*(dims->nff + dims->nx1)];
    printf("KKf Mat\n");
    d_print_mat(dims->num_stages* dims->nx1, dims->nff, fix->KKf, dims->nx1 *dims->num_stages);
    printf("KKx Mat\n");
    d_print_mat(dims->num_stages* dims->nx1, dims->nx1, fix->KKx, dims->nx1 *dims->num_stages);
        printf("KKu Mat\n");
    d_print_mat(dims->num_stages* dims->nx1, dims->nu, fix->KKu, dims->nx1 *dims->num_stages);

}


void gnsf_allocate_fixed( gnsf_dims *dims, gnsf_fixed *fix){ // Didnt really work..
    fix->KKf = calloc(dims->nx1 * dims->num_stages * (dims->nff + dims->nu + dims->nx1), sizeof(double));
    fix->KKx = fix->KKf + (dims->nx1 * dims->num_stages *dims->nff)* sizeof(double*);
    fix->KKu = fix->KKx + dims->nx1 * dims->num_stages *dims->nx1;
}

void gnsf_simulate( gnsf_dims *dims, gnsf_in *in, gnsf_out out)
{
    printf("GENERALIZED NONLINEAR STATIC FEEDBACK (GNSF) SIMULATION \n");
    print_gnsf_dims(dims);

    double *res_in;
    int res_in_size = dims->nff + dims->nx1 + dims->nu;
    res_in = (double*) calloc(res_in_size, sizeof(double));


    double *res_out;
    int res_out_size = dims->nff * (1 + dims->nff);
    res_out = (double*) calloc(res_out_size, sizeof(double));

    res_in[2+dims->nff] = 0.8;

    in->A_dt[0] = 1.0;
    printf("%f",in->A_dt[0]);

    print_gnsf_res_in( *dims, res_in );
    res_inc_Jff_wrapped(dims->nx1, dims->nu, dims->n_out, dims->num_stages, res_in, res_out, in->res_inc_Jff);
    print_gnsf_res_out( *dims, res_out );

    struct blasfeo_dmat res_inc_jac;
    blasfeo_allocate_dmat(dims->nff, dims->nff+1, &res_inc_jac);
    // void blasfeo_pack_dmat(int m, int n, double *A, int lda, struct blasfeo_dmat *sB, int bi, int bj);

    blasfeo_pack_dmat(1, 1, in->A_dt, 1, &res_inc_jac, 0, 0);

    blasfeo_print_dmat(dims->nff, dims->nff+1, &res_inc_jac, 0,0);

    blasfeo_free_dmat(&res_inc_jac);
}