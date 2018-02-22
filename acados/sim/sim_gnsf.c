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
#include "acados/sim/sim_gnsf_casadi_wrapper.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux.h"


void print_gnsf_dims( gnsf_dims *dims )
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

int gnsf_dims_calculate_size()
{
    int size = sizeof(gnsf_dims);
    return size;
}

gnsf_dims *assign_gnsf_dims(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    gnsf_dims *dims = (gnsf_dims *) c_ptr;
    c_ptr += sizeof(gnsf_dims);
    assert((char *) raw_memory + gnsf_dims_calculate_size() == c_ptr);
    return dims;
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
    free(ints_out);
}

int gnsf_in_calculate_size(gnsf_dims *dims)
{
    int size = sizeof(gnsf_in);

    int nx = dims->nx;
    int nu = dims->nu;

    size += nx * sizeof(double);  // x
    size += nu * sizeof(double);  // u
    size += nx * (nx+nu) * sizeof(double);  // S_forw (max dimension)
    // size += (nx + nu) * sizeof(double);  // S_adj

    make_int_multiple_of(8, &size);
    return size;
}


gnsf_in *assign_gnsf_in(gnsf_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    // printf("address c_ptr : %p\n",c_ptr);
    gnsf_in *in = (gnsf_in *) c_ptr;
    c_ptr += sizeof(gnsf_in);

    int nx = dims->nx;
    int nu = dims->nu;

    align_char_to(8, &c_ptr);

    assign_double(nx, &in->x, &c_ptr);
    // printf("address in_x : %p\n",in->x);
    assign_double(nu, &in->u, &c_ptr);
    assign_double(nx * (nx+nu), &in->S_forw, &c_ptr);
    assert((char*)raw_memory + gnsf_in_calculate_size(dims) == c_ptr);

    // printf("address c_ptr : %p\n",c_ptr);
    return in;
}


int gnsf_opts_calculate_size(gnsf_dims *dims)
{
    int size = sizeof(gnsf_opts);
    make_int_multiple_of(8, &size);
    return size;
}


gnsf_opts *assign_gnsf_opts(gnsf_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    gnsf_opts *in = (gnsf_opts *) c_ptr;
    c_ptr += sizeof(gnsf_opts);

    int nx = dims->nx;
    int nu = dims->nu;

    align_char_to(8, &c_ptr);
    assert((char*)raw_memory + gnsf_opts_calculate_size(dims) == c_ptr);
    return in;
}


void gnsf_import(gnsf_dims* dims, gnsf_fixed *fix, casadi_function_t But_KK_ZZ_LO_fun)
{
    acados_timer atimer;
    acados_tic(&atimer);
    int nx  = dims->nx;
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int num_stages = dims->num_stages;
    int num_steps = dims->num_steps;

    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;
    int nff = n_out * num_stages;

    // double *out;
    int exported_doubles = 0;
    exported_doubles += num_stages * (num_stages +2); // Butcher matrices
    exported_doubles += nK1 * (nff + nx1 + nu); //KK* matrices
    exported_doubles += nZ * (nff + nx1 + nu); //ZZ* matrices
    exported_doubles += nK2 * (nK2 + nx2) + nx2 * nx2; //LO matrices

    // printf("exported_%d \n",exported_doubles);
    double *export_in  = (double*) malloc(1*sizeof(double));
    double *exp_out = (double*) malloc(exported_doubles*sizeof(double));
    export_from_ML_wrapped(export_in, exp_out, But_KK_ZZ_LO_fun);

    // printf("\n adress \n%p",(void*)&exp_out);
    double *read_mem = exp_out;

    // IMPORT BUTCHER
    for (int ii = 0; ii < num_stages*num_stages; ii++) {
        fix->A_dt[ii] = read_mem[ii];
    }
    read_mem += num_stages*num_stages;

    for (int ii = 0; ii < num_stages; ii++) {
        fix->b_dt[ii] = read_mem[ii];
    }
    read_mem += num_stages;

    for (int ii = 0; ii < num_stages; ii++) {
        fix->c[ii] = read_mem[ii];
    }
    read_mem += num_stages;

    // IMPORT KKmat
    blasfeo_pack_dmat(nK1, nff, read_mem, nK1, &fix->KKf, 0, 0);
    read_mem += nK1 * nff;
    blasfeo_pack_dmat(nK1, nx1, read_mem, nK1, &fix->KKx, 0, 0);
    read_mem += nK1 * nx1;
    blasfeo_pack_dmat(nK1, nu,  read_mem, nK1, &fix->KKu, 0, 0);
    read_mem += nK1 * nu;

    // IMPORT ZZmat
    blasfeo_pack_dmat(nZ, nff, read_mem, nZ, &fix->ZZf, 0, 0);
    read_mem += nZ * nff;
    blasfeo_pack_dmat(nZ, nx1, read_mem, nZ, &fix->ZZx, 0, 0);
    read_mem += nZ * nx1;
    blasfeo_pack_dmat(nZ, nu,  read_mem, nZ, &fix->ZZu, 0, 0);
    read_mem += nZ * nu;

    // IMPORT LO matrices
    blasfeo_pack_dmat(nx2, nx2, read_mem, nx2, &fix->ALO, 0, 0);
    read_mem += nx2 * nx2;
    blasfeo_pack_dmat(nK2, nK2, read_mem, nK2, &fix->M2inv, 0, 0);
    read_mem += nK2 * nK2;
    blasfeo_pack_dmat(nK2, nx2, read_mem, nx2, &fix->dK2_dx2, 0, 0);
    read_mem += nK2 * nx2;

    free(exp_out);
    free(export_in);

    // d_print_e_mat(num_stages, num_stages, fix->A_dt, num_stages);
    // d_print_e_mat(num_stages, 1, fix->b_dt, num_stages);
    // d_print_e_mat(num_stages, 1, fix->c, num_stages);

    // printf("KKf Mat\n");
    // blasfeo_print_exp_dmat(nK1, nff, &fix->KKf, 0,0);
    // printf("KKx Mat\n");
    // blasfeo_print_exp_dmat(nK1, nx1, &fix->KKx, 0,0);
    // printf("KKu Mat \n");
    // blasfeo_print_exp_dmat(nK1, nu , &fix->KKu, 0,0);

    // printf("ZZf Mat\n");
    // blasfeo_print_exp_dmat(nZ, nff, &fix->ZZf, 0,0);
    // printf("ZZx Mat\n");
    // blasfeo_print_exp_dmat(nZ, nx1, &fix->ZZx, 0,0);
    // printf("ZZu Mat \n");
    // blasfeo_print_exp_dmat(nZ, nu , &fix->ZZu, 0,0);

    // printf("ALO Mat\n");
    // blasfeo_print_dmat(dims->nx2, dims->nx2, &fix->ALO, 0,0);
    // printf("M2inv\n");
    // blasfeo_print_dmat(nK2, nK2, &fix->M2inv, 0,0);
    // printf("dK2_dx2 \n");
    // blasfeo_print_dmat(nK2, dims->nx2 , &fix->dK2_dx2, 0,0);
}

int gnsf_fixed_calculate_size(gnsf_dims *dims, gnsf_opts* opts)
{
    int nx  = dims->nx;
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int num_stages = dims->num_stages;
    int num_steps = dims->num_steps;

    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;
    int nff = n_out * num_stages;

    int size = 8;
    size += sizeof(gnsf_fixed);

    size += num_stages * (num_stages + 2) * sizeof(double); // A_dt, b_dt, c_butcher;

	// size += 9*sizeof(struct blasfeo_dmat); // 3 * KK*, 3* ZZ*, 3* LO_mat // not needed because KK*,... are part of gnsf_fixed

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    size += blasfeo_memsize_dmat(nK1, nff); // KKf
    size += blasfeo_memsize_dmat(nK1, nx1); // KKx
    size += blasfeo_memsize_dmat(nK1, nu ); // KKu

    size += blasfeo_memsize_dmat(nZ, nff); // ZZf
    size += blasfeo_memsize_dmat(nZ, nx1); // ZZx
    size += blasfeo_memsize_dmat(nZ, nu ); // ZZu

    size += blasfeo_memsize_dmat(nx2, nx2); // ALO
    size += blasfeo_memsize_dmat(nK2, nK2); // M2inv
    size += blasfeo_memsize_dmat(nK2, nx2); // dK2_dx2
    // printf("gnsf_fixed_size = %d \n", size);

    return size;
}

gnsf_fixed *gnsf_fixed_assign(gnsf_dims *dims, void *raw_memory, int memsize)
{
    char *c_ptr = (char *) raw_memory;

	// extract sizes
    int nx  = dims->nx;
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int num_stages = dims->num_stages;
    int num_steps = dims->num_steps;

    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;
    int nff = n_out * num_stages;

	// initial align
	align_char_to(8, &c_ptr);
    // printf("\n adress of cptr fixed %p\n",(void*)c_ptr);

	// struct
    gnsf_fixed *fix = (gnsf_fixed *) c_ptr;
    c_ptr += sizeof(gnsf_fixed);

    // assign butcher
    assign_double(num_stages * num_stages, &fix->A_dt, &c_ptr);
    assign_double(num_stages, &fix->b_dt, &c_ptr);
    assign_double(num_stages, &fix->c,    &c_ptr);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

    // blasfeo_dmat_mem
    assign_blasfeo_dmat_mem(nK1, nff, &fix->KKf, &c_ptr);
    // printf("\n adress of KKf %p\n",(void*)&fix->KKf);
    assign_blasfeo_dmat_mem(nK1, nx1, &fix->KKx, &c_ptr);
    assign_blasfeo_dmat_mem(nK1, nu,  &fix->KKu, &c_ptr);

    assign_blasfeo_dmat_mem(nZ,  nff, &fix->ZZf, &c_ptr);
    assign_blasfeo_dmat_mem(nZ,  nx1, &fix->ZZx, &c_ptr);
    assign_blasfeo_dmat_mem(nZ,  nu,  &fix->ZZu, &c_ptr);

    assign_blasfeo_dmat_mem(nx2, nx2, &fix->ALO, &c_ptr);
    assign_blasfeo_dmat_mem(nK2, nK2, &fix->M2inv, &c_ptr);
    assign_blasfeo_dmat_mem(nK2, nx2, &fix->dK2_dx2, &c_ptr);

	// // assert
    // assert((char *) raw_memory + memsize == c_ptr); TODO recheck..
	return fix;
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
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;

    int size = sizeof(gnsf_workspace);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    int res_in_size = nff + nx1 + nu;
    int res_out_size = nff * (nff + nx1 + nu); // size(out_res_inc_Jff) = nff* (1+nff), size(out_jac_res_ffx1u) = nff(nff+nx1+nu)
    int f_LO_in_size = 2*nx1 + nu + nz;
    int f_LO_out_size = nx2 * (1 + 2*nx1 + nu + nz);

    size += (res_in_size + res_out_size + f_LO_in_size + f_LO_out_size) * sizeof(double); // input and outputs of residual and LO-fcn;

    size += 5 *num_steps * sizeof(struct blasfeo_dvec); // K1_val, x1_val, ff_val, Z_val, f_LO_val
	size += num_steps * sizeof(struct blasfeo_dmat); // f_LO_jac

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    size += num_steps * blasfeo_memsize_dmat(nK2, 2*nx1 +nu +nz); //f_LO_jac
    size += 2* num_steps * blasfeo_memsize_dvec(nK1); //K1_val, x1_val
    size += num_steps * blasfeo_memsize_dvec(nff); // ff_val
    size += num_steps * blasfeo_memsize_dvec(nZ); // Z_val
    size += num_steps * blasfeo_memsize_dvec(nK2); // f_LO_val

    size += blasfeo_memsize_dmat(nff, nff); // J_r_ff
    size += blasfeo_memsize_dmat(nff, nx1+nu); // J_r_x1u

    size += blasfeo_memsize_dmat(nK1, nx1); // dK1_dx1
    size += blasfeo_memsize_dmat(nK1, nu ); // dK1_du
    size += blasfeo_memsize_dmat(nZ, nx1); // dZ_dx1
    size += blasfeo_memsize_dmat(nZ, nu); // dZ_du
    size += blasfeo_memsize_dmat(nK2, nx1); // aux_G2_x1
    size += blasfeo_memsize_dmat(nK2, nu); // aux_G2_u
    size += blasfeo_memsize_dmat(nK2, nK1); // J_G2_K1
    size += blasfeo_memsize_dmat(nK2, nx1); // dK2_dx1
    size += blasfeo_memsize_dmat(nK2, nu); // dK2_du
    size += blasfeo_memsize_dmat(nK2, nff); // dK2_dff
    size += blasfeo_memsize_dmat(nx, nx + nu); // dxf_dwn
    size += 2*blasfeo_memsize_dmat(nx, nx + nu); // S_forw_new, S_forw

    size += blasfeo_memsize_dvec(nK2); // K2_val
    size += blasfeo_memsize_dvec(nx*(num_steps+1)); // x0_traj
    size += blasfeo_memsize_dvec(nff); // res_val
    size += blasfeo_memsize_dvec(nu); // u0

    size += blasfeo_memsize_dmat(nK2,nff); // aux_G2_ff
    size += blasfeo_memsize_dmat(nx, nff); // dPsi_dff
    size += blasfeo_memsize_dmat(nx, nx ); // dPsi_dx
    size += blasfeo_memsize_dmat(nx, nu ); // dPsi_du

    make_int_multiple_of(8, &size);
    size += 1 * 8;
    // printf("workspace size = %d \n", size);
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
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;

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

    assign_blasfeo_dmat_structs(num_steps, &workspace->f_LO_jac, &c_ptr);

    assign_blasfeo_dvec_structs(num_steps, &workspace->K1_val, &c_ptr);
    assign_blasfeo_dvec_structs(num_steps, &workspace->x1_val, &c_ptr);
    assign_blasfeo_dvec_structs(num_steps, &workspace->ff_val, &c_ptr);
    assign_blasfeo_dvec_structs(num_steps, &workspace->Z_val, &c_ptr);
    assign_blasfeo_dvec_structs(num_steps, &workspace->f_LO_val, &c_ptr);

    // blasfeo_mem align
	align_char_to(64, &c_ptr);
    for (int ii=0; ii<num_steps; ii++){
        assign_blasfeo_dmat_mem(nK2, 2*nx1+nu+nz, workspace->f_LO_jac+ii, &c_ptr);     // f_LO_jac

        assign_blasfeo_dvec_mem(nK1, workspace->K1_val+ii, &c_ptr);     // K1_val
        assign_blasfeo_dvec_mem(nK1, workspace->x1_val+ii, &c_ptr);     // x1_val
        assign_blasfeo_dvec_mem(nff, workspace->ff_val+ii, &c_ptr);     // ff_val
        assign_blasfeo_dvec_mem(nZ , workspace->Z_val+ii, &c_ptr);     // Z_val
        assign_blasfeo_dvec_mem(nK2, workspace->f_LO_val+ii, &c_ptr);     // Z_val
    }

    assign_blasfeo_dmat_mem(nff, nx1+nu, &workspace->J_r_x1u , &c_ptr);

    assign_blasfeo_dmat_mem(nff, nff, &workspace->J_r_ff , &c_ptr);
    assign_blasfeo_dmat_mem(nK1, nx1, &workspace->dK1_dx1 , &c_ptr);
    assign_blasfeo_dmat_mem(nK1, nu , &workspace->dK1_du  , &c_ptr);
    assign_blasfeo_dmat_mem(nZ, nx1, &workspace->dZ_dx1 , &c_ptr);
    assign_blasfeo_dmat_mem(nZ, nu , &workspace->dZ_du  , &c_ptr);

    assign_blasfeo_dmat_mem(nK2, nx1, &workspace->aux_G2_x1, &c_ptr);
    assign_blasfeo_dmat_mem(nK2, nu , &workspace->aux_G2_u , &c_ptr);
    assign_blasfeo_dmat_mem(nK2, nK1 , &workspace->J_G2_K1 , &c_ptr);

    assign_blasfeo_dmat_mem(nK2, nx1, &workspace->dK2_dx1 , &c_ptr);
    assign_blasfeo_dmat_mem(nK2, nu, &workspace->dK2_du , &c_ptr);
    assign_blasfeo_dmat_mem(nK2, nff, &workspace->dK2_dff, &c_ptr);
    assign_blasfeo_dmat_mem(nx, nx+nu, &workspace->dxf_dwn , &c_ptr);
    assign_blasfeo_dmat_mem(nx, nx+nu, &workspace->S_forw_new , &c_ptr);
    assign_blasfeo_dmat_mem(nx, nx+nu, &workspace->S_forw, &c_ptr);

    assign_blasfeo_dvec_mem(nK2, &workspace->K2_val, &c_ptr);
    assign_blasfeo_dvec_mem((num_steps+1)*nx, &workspace->x0_traj, &c_ptr);
    assign_blasfeo_dvec_mem(nff, &workspace->res_val, &c_ptr);
    assign_blasfeo_dvec_mem(nu, &workspace->u0, &c_ptr);

    assign_blasfeo_dmat_mem(nK2, nff, &workspace->aux_G2_ff, &c_ptr);
    assign_blasfeo_dmat_mem(nx, nff, &workspace->dPsi_dff , &c_ptr);
    assign_blasfeo_dmat_mem(nx, nx, &workspace->dPsi_dx , &c_ptr);
    assign_blasfeo_dmat_mem(nx, nu, &workspace->dPsi_du, &c_ptr);

    return (void *)workspace;
}

void gnsf_simulate( gnsf_dims *dims, gnsf_fixed *fix, gnsf_in *in, sim_out *out, gnsf_opts *opts, void *work_)
{
    acados_timer tot_timer, casadi_timer;
    acados_tic(&tot_timer);
    printf("GENERALIZED NONLINEAR STATIC FEEDBACK (GNSF) SIMULATION \n");
    // print_gnsf_dims(dims);

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

    int newton_max = opts->newton_max;

    double *res_in = workspace->res_in;
    double *res_out = workspace->res_out;
    double *f_LO_in = workspace->f_LO_in;
    double *f_LO_out = workspace->f_LO_out;

    struct blasfeo_dmat J_r_ff = workspace->J_r_ff; // store the the jacobian of the residual w.r.t. ff
    struct blasfeo_dmat J_r_x1u = workspace->J_r_x1u;  // needed for sensitivity propagation

    struct blasfeo_dmat dK1_dx1    = workspace->dK1_dx1;
    struct blasfeo_dmat dK1_du     = workspace->dK1_du;
    struct blasfeo_dmat dZ_dx1     = workspace->dZ_dx1;
    struct blasfeo_dmat dZ_du      = workspace->dZ_du;
    struct blasfeo_dmat aux_G2_x1  = workspace->aux_G2_x1;
    struct blasfeo_dmat aux_G2_u   = workspace->aux_G2_u;
    struct blasfeo_dmat J_G2_K1    = workspace->J_G2_K1;
    struct blasfeo_dmat dK2_dx1    = workspace->dK2_dx1;
    struct blasfeo_dmat dK2_du     = workspace->dK2_du;
    struct blasfeo_dmat dK2_dff    = workspace->dK2_dff;
    struct blasfeo_dmat dxf_dwn    = workspace->dxf_dwn;
    struct blasfeo_dmat S_forw_new = workspace->S_forw_new; // used to avoid side effects
    struct blasfeo_dmat S_forw     = workspace->S_forw;

    struct blasfeo_dmat *f_LO_jac = workspace->f_LO_jac;

    struct blasfeo_dvec *ff_val  = workspace->ff_val; 
    struct blasfeo_dvec *K1_val  = workspace->K1_val; 
    struct blasfeo_dvec *x1_val  = workspace->x1_val; 
    struct blasfeo_dvec *Z_val   = workspace->Z_val;
    struct blasfeo_dvec *f_LO_val= workspace->f_LO_val;

    struct blasfeo_dvec K2_val  = workspace->K2_val;
    struct blasfeo_dvec x0_traj = workspace->x0_traj;
    struct blasfeo_dvec res_val = workspace->res_val;
    struct blasfeo_dvec u0      = workspace->u0;

    struct blasfeo_dmat aux_G2_ff = workspace->aux_G2_ff;
    struct blasfeo_dmat dPsi_dff  = workspace->dPsi_dff;
    struct blasfeo_dmat dPsi_dx   = workspace->dPsi_dx;
    struct blasfeo_dmat dPsi_du   = workspace->dPsi_du;

    blasfeo_pack_dvec(nu, in->u, &u0, 0);    
    blasfeo_pack_dvec(nx, &in->x[0], &x0_traj, 0);
    blasfeo_pack_dmat(nx, nx + nu, &in->S_forw[0], nx, &S_forw, 0, 0);

    out->info->ADtime = 0;

    for (int ss = 0; ss < num_steps; ss++) {
        for (int iter = 0; iter < newton_max; iter++) { // NEWTON-ITERATION
            // set input for residual function
            blasfeo_unpack_dvec(nff, &ff_val[ss], 0, &res_in[0]);
            blasfeo_unpack_dvec(nx1, &x0_traj, ss*nx, &res_in[nff]);
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
        // K1_val = s.KKf * fftraj(:,ss) + s.KKu * u0 + s.KKx * x0_1;
        blasfeo_dgemv_n(nK1, nff,       1.0, &fix->KKf, 0, 0, &ff_val[ss], 0, 0.0, &K1_val[ss], 0, &K1_val[ss], 0);
        blasfeo_dgemv_n(nK1, nu , 1.0, &fix->KKu, 0, 0, &u0        , 0, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
        blasfeo_dgemv_n(nK1, nx1, 1.0, &fix->KKx, 0, 0, &x0_traj, ss*nx, 1.0, &K1_val[ss], 0, &K1_val[ss], 0);
        // printf("\n K1_val =  \n");
        // blasfeo_print_exp_dvec(nK1, &K1_val[ss],0);
        blasfeo_dgemv_n(nZ, nff,       1.0, &fix->ZZf, 0, 0, &ff_val[ss], 0, 0.0, &Z_val[ss], 0, &Z_val[ss], 0);
        blasfeo_dgemv_n(nZ, nu , 1.0, &fix->ZZu, 0, 0, &u0        , 0, 1.0, &Z_val[ss], 0, &Z_val[ss], 0);
        blasfeo_dgemv_n(nZ, nx1, 1.0, &fix->ZZx, 0, 0, &x0_traj, ss*nx, 1.0, &Z_val[ss], 0, &Z_val[ss], 0);
        // build x1 stage values
        for (int ii = 0; ii < num_stages; ii++){
            blasfeo_daxpy(nx1, 0.0, &x1_val[ss], 0, &x0_traj, ss*nx, &x1_val[ss], nx1 * ii);
            for (int jj = 0; jj <num_stages; jj++) {
                blasfeo_daxpy(nx1, fix->A_dt[ii+num_stages*jj], &K1_val[ss], nx1*jj, &x1_val[ss], nx1*ii, &x1_val[ss], nx1*ii);
            }
        }

        // SIMULATE LINEAR OUTPUT SYSTEM
        for (int ii = 0; ii < num_stages; ii++) {
            blasfeo_unpack_dvec(nx1, &x1_val[ss], ii*nx1, &f_LO_in[0]);
            blasfeo_unpack_dvec(nx1, &K1_val[ss], ii*nx1, &f_LO_in[nx1]);
            blasfeo_unpack_dvec(nz,  &Z_val[ss] , ii*nz , &f_LO_in[2*nx1]);
            blasfeo_unpack_dvec(nu,  &u0        ,  0          , &f_LO_in[2*nx1 +nz]);
            // printf("f_LO_in = \n");
            // d_print_mat(f_LO_in_size, 1, &f_LO_in[0], f_LO_in_size);
            acados_tic(&casadi_timer);
            f_LO_inc_J_x1k1uz_wrapped(nx1, nz, f_LO_in, f_LO_out, in->f_LO_inc_J_x1k1uz);
            out->info->ADtime += acados_toc(&casadi_timer);
            // printf("f_LO_out= \n");
            // d_print_mat(f_LO_out_size,1, &f_LO_out[0],f_LO_out_size);
            // printf("f_LO_out_size= %d \n", f_LO_out_size);
            blasfeo_pack_dvec(nx2, &f_LO_out[0], &f_LO_val[ss], nx2 * ii);
            // printf("\n adress %p\n",(void*)&f_LO_jac[ss]);
            blasfeo_pack_dmat(nx2, 2*nx1 + nu + nz, &f_LO_out[nx2], nx2, &f_LO_jac[ss], nx2 * ii, 0); // NOTE: f_LO_jac has different sign compared to Matlab prototype
            blasfeo_dgemv_n(nx2, nx2, 1.0, &fix->ALO, 0, 0, &x0_traj, ss*nx+nx1, -1.0, &f_LO_val[ss], nx2 * ii, &f_LO_val[ss], nx2 * ii); // todo: repmat( - s.ALO * x0_2, q, 1); could be translated more efficient
        }
        // printf("f_LO = \n");
        // blasfeo_print_exp_dvec(num_stages * nx2, &f_LO_val[ss], 0);
        // blasfeo_print_dmat(    num_stages * nx2, 2*nx1 + nz + nu, &f_LO_jac[ss],0,0);
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
        blasfeo_unpack_dvec(nx1, &x0_traj, nx * ss, &res_in[nff]);
        for (int i = 0; i<nu; i++) {
            res_in[i+nff+nx1] = in->u[i];
        }
        if (opts->sens_forw) {
            acados_tic(&casadi_timer);
            jac_res_ffx1u_wrapped(nx1, nu, n_out, num_stages, res_in, res_out, in->jac_res_ffx1u);
            out->info->ADtime += acados_toc(&casadi_timer);
            blasfeo_pack_dmat(nff, nff, &res_out[0], nff, &J_r_ff, 0, 0); // pack residual result into blasfeo struct
            blasfeo_pack_dmat(nff, nx1+ nu, &res_out[nff*nff], nff, &J_r_x1u, 0, 0); // pack residual result into blasfeo struct

            blasfeo_dgetrf_nopivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0); // factorize J_r_ff
            blasfeo_dtrsm_lunn(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
            blasfeo_dtrsm_llnu(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);

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
            blasfeo_dgemm_nn(nx, nx, nx, 1.0, &dxf_dwn, 0, 0, &S_forw, 0, 0, 0.0, &S_forw_new, 0, 0, &S_forw_new, 0, 0);
            blasfeo_dgemm_nn(nx, nu, nx, 1.0, &dxf_dwn, 0, 0, &S_forw, 0, nx, 1.0, &dxf_dwn, 0, nx, &S_forw_new, 0, nx);
            blasfeo_dgecp(nx, nx +nu, &S_forw_new, 0, 0, &S_forw, 0, 0);
        }
    }
    if (opts->sens_adj) {
        // ADJOINT SENSITIVITY PROPAGATION:
        for (int ss = num_steps-1; ss >= num_steps-1; ss--) {
            for (int ii = 0; ii < num_stages; ii++) {
                blasfeo_dgemm_nn(nx2, nff      , nz, -1.0, &f_LO_jac[ss], nx2 * ii, 2*nx1, &fix->ZZf, ii* nz, 0, 0.0, &fix->ZZf, 0, 0, &aux_G2_ff, ii * nx2, 0);
                blasfeo_dgemm_nn(nx2, nx1, nz, -1.0, &f_LO_jac[ss], nx2 * ii, 2*nx1, &fix->ZZx, ii* nz, 0, 0.0, &fix->ZZx, 0, 0, &aux_G2_x1, ii * nx2, 0);
                blasfeo_dgemm_nn(nx2, nu, nz, -1.0, &f_LO_jac[ss], nx2 * ii, 2*nx1, &fix->ZZu, ii* nz, 0, 0.0, &fix->ZZu, 0, 0, &aux_G2_u, ii * nx2, 0);
                for (int jj = 0; jj < num_stages; jj++) {
                    blasfeo_dgecpsc(nx2, nx1, -fix->A_dt[ii+jj*num_stages], &f_LO_jac[ss], nx2 * ii, 0, &J_G2_K1, ii*nx2, jj*nx1);
                }
                blasfeo_dgead(nx2, nx1, -1.0, &f_LO_jac[ss], nx2*ii, nx1, &J_G2_K1, nx2*ii, nx1*ii);
            }
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
        }
    }
    out->info->CPUtime = acados_toc(&tot_timer);
    printf("tot_time = %f\n", out->info->CPUtime);
    blasfeo_unpack_dvec(nx, &x0_traj, nx * num_steps, out->xn); //TODO pack everythin in out
    blasfeo_unpack_dmat(nx, nx + nu, &S_forw, 0, 0, out->S_forw, nx);
}