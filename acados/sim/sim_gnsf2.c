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
 *    Author: Jonathan Frey
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
#include "acados/sim/sim_gnsf2.h"
#include "acados/sim/sim_gnsf_casadi_wrapper.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux.h"


int gnsf2_dims_calculate_size()
{
    int size = sizeof(gnsf2_dims);
    return size;
}


gnsf2_dims *gnsf2_dims_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    gnsf2_dims *dims = (gnsf2_dims *) c_ptr;
    c_ptr += sizeof(gnsf2_dims);
    assert((char *) raw_memory + gnsf2_dims_calculate_size() == c_ptr);
    return dims;
}

void gnsf2_get_dims( gnsf2_dims *dims, casadi_function_t get_ints_fun)
{
    double *ints_out;
    ints_out = (double*) calloc(9,sizeof(double));
    export_from_ML_wrapped(ints_out, ints_out, get_ints_fun);

    dims->nx = (int) ints_out[0];
    dims->nu = (int) ints_out[1]; 
    dims->nz = (int) ints_out[2];
    dims->nx1 = (int) ints_out[3];
    dims->nx2 = (int) ints_out[4];
    dims->num_stages = (int) ints_out[5];
    dims->num_steps = (int) ints_out[6];
    dims->n_out = (int) ints_out[7];
    dims->n_in = (int) ints_out[8];
    free(ints_out);
}

int gnsf2_in_calculate_size(gnsf2_dims *dims)
{
    int size = sizeof(gnsf2_in);

    int nx = dims->nx;
    int nu = dims->nu;

    size += nx * sizeof(double);  // x
    size += nu * sizeof(double);  // u
    size += nx * (nx+nu) * sizeof(double);  // S_forw (max dimension)
    size += (nx + nu) * sizeof(double);  // S_adj

    make_int_multiple_of(8, &size);
    return size;
}

gnsf2_in *gnsf2_in_assign(gnsf2_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    // printf("address c_ptr : %p\n",c_ptr);
    gnsf2_in *in = (gnsf2_in *) c_ptr;
    c_ptr += sizeof(gnsf2_in);

    int nx = dims->nx;
    int nu = dims->nu;

    align_char_to(8, &c_ptr);

    assign_double(nx, &in->x, &c_ptr);
    // printf("address in_x : %p\n",in->x);
    assign_double(nu, &in->u, &c_ptr);
    assign_double(nx * (nx+nu), &in->S_forw, &c_ptr);
    assign_double(nx+nu, &in->S_adj, &c_ptr);
    assert((char*)raw_memory + gnsf2_in_calculate_size(dims) == c_ptr);

    return in;
}

int gnsf2_opts_calculate_size(gnsf2_dims *dims)
{
    int size = sizeof(gnsf2_opts);
    make_int_multiple_of(8, &size);
    return size;
}


gnsf2_opts *gnsf2_opts_assign(gnsf2_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    gnsf2_opts *opts = (gnsf2_opts *) c_ptr;
    c_ptr += sizeof(gnsf2_opts);

    align_char_to(8, &c_ptr);
    assert((char*)raw_memory + gnsf_opts_calculate_size(dims) == c_ptr);
    return opts;
}

int gnsf2_fixed_calculate_size(gnsf2_dims *dims, gnsf2_opts* opts)
{
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int num_stages = dims->num_stages;
    int n_in = dims->n_in;

    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;
    int nff = n_out * num_stages;
    int nyy = n_in  * num_stages;

    int size = 8;
    size += sizeof(gnsf2_fixed);

    size += num_stages * (num_stages + 2) * sizeof(double); // A_dt, b_dt, c_butcher;

	// size += 9*sizeof(struct blasfeo_dmat); // 3 * KK*, 3* ZZ*, 3* LO_mat // not needed because KK*,... are part of gnsf_fixed

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    size += blasfeo_memsize_dmat(nK1, nff); // KKf
    size += blasfeo_memsize_dmat(nK1, nx1); // KKx
    size += blasfeo_memsize_dmat(nK1, nu ); // KKu

    size += blasfeo_memsize_dmat(nyy, nff); // YYf
    size += blasfeo_memsize_dmat(nyy, nx1); // YYx
    size += blasfeo_memsize_dmat(nyy, nu ); // YYu

    size += blasfeo_memsize_dmat(nZ, nff); // ZZf
    size += blasfeo_memsize_dmat(nZ, nx1); // ZZx
    size += blasfeo_memsize_dmat(nZ, nu ); // ZZu

    size += blasfeo_memsize_dmat(nx2, nx2); // ALO
    size += blasfeo_memsize_dmat(nK2, nK2); // M2inv
    size += blasfeo_memsize_dmat(nK2, nx2); // dK2_dx2
    // printf("gnsf_fixed_size = %d \n", size);

    return size;
}

gnsf2_fixed *gnsf2_fixed_assign(gnsf2_dims *dims, void *raw_memory, int memsize)
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
    int n_in = dims->n_in;

    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;
    int nff = n_out * num_stages;
    int nyy = n_in  * num_stages;

	// initial align
	align_char_to(8, &c_ptr);
    // printf("\n adress of cptr fixed %p\n",(void*)c_ptr);

	// struct
    gnsf2_fixed *fix = (gnsf2_fixed *) c_ptr;
    c_ptr += sizeof(gnsf2_fixed);

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

    assign_blasfeo_dmat_mem(nyy,  nff, &fix->YYf, &c_ptr);
    assign_blasfeo_dmat_mem(nyy,  nx1, &fix->YYx, &c_ptr);
    assign_blasfeo_dmat_mem(nyy,  nu,  &fix->YYu, &c_ptr);

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

void gnsf2_import(gnsf2_dims* dims, gnsf2_fixed *fix, casadi_function_t But_KK_YY_ZZ_LO_fun)
{
    acados_timer atimer;
    acados_tic(&atimer);
    int nu  = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int n_in = dims->n_in;
    int num_stages = dims->num_stages;

    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;
    int nff = n_out * num_stages;
    int nyy = n_in  * num_stages;

    // double *out;
    int exported_doubles = 0;
    exported_doubles += num_stages * (num_stages +2); // Butcher matrices
    exported_doubles += nK1 * (nff + nx1 + nu); //KK* matrices
    exported_doubles += nyy * (nff + nx1 + nu); //YY* matrices
    exported_doubles += nZ * (nff + nx1 + nu); //ZZ* matrices
    exported_doubles += nK2 * (nK2 + nx2) + nx2 * nx2; //LO matrices

    // printf("exported_%d \n",exported_doubles);
    double *export_in  = (double*) malloc(1*sizeof(double));
    double *exp_out = (double*) malloc(exported_doubles*sizeof(double));
    export_from_ML_wrapped(export_in, exp_out, But_KK_YY_ZZ_LO_fun);

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

    // IMPORT YYmat
    blasfeo_pack_dmat(nyy, nff, read_mem, nyy, &fix->YYf, 0, 0);
    read_mem += nyy * nff;
    blasfeo_pack_dmat(nyy, nx1, read_mem, nyy, &fix->YYx, 0, 0);
    read_mem += nyy * nx1;
    blasfeo_pack_dmat(nyy, nu,  read_mem, nyy, &fix->YYu, 0, 0);
    read_mem += nyy * nu;
    printf("YYu = \n");
    blasfeo_print_exp_dmat(nyy, nu, &fix->YYu, 0, 0);

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
}


void sim_gnsf2_config_initialize_default(void *config_)
{
	sim_solver_config *config = config_;
// TODO:!!
	// config->evaluate = &gnsf2_simulate;
	// config->opts_calculate_size = &gnsf2_opts_calculate_size;
	// config->opts_assign = &gnsf2_opts_assign;
	// config->opts_initialize_default = &sim_irk_opts_initialize_default; TODO
	// config->memory_calculate_size = &sim_irk_memory_calculate_size; TODO
	// config->memory_assign = &sim_irk_memory_assign;  TODO
	// config->workspace_calculate_size = &gnsf_calculate_workspace_size;
	// config->model_calculate_size = &gnsf2_model_calculate_size;
	// config->model_assign = &gnsf2_model_assign;

	return;
}