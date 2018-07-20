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
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/math.h"

#include "acados/sim/sim_collocation_utils.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_gnsf.h"

// blasfeo
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_target.h"
// #include "blasfeo/include/blasfeo_d_aux_ext_dep.h" // can be included for printing while
// debugging

/************************************************
 * dims
 ************************************************/

int sim_gnsf_dims_calculate_size()
{
    int size = sizeof(sim_gnsf_dims);
    return size;
}

void *sim_gnsf_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    sim_gnsf_dims *dims = (sim_gnsf_dims *) c_ptr;
    c_ptr += sizeof(sim_gnsf_dims);

    dims->nx = 0;
    dims->nu = 0;
    dims->nz = 0;
    dims->nx1 = 0;
    dims->nx2 = 0;
    dims->n_out = 0;
    dims->ny = 0;
    dims->nuhat = 0;

    assert((char *) raw_memory + sim_gnsf_dims_calculate_size() == c_ptr);
    return dims;
}

/************************************************
 * get & set functions
 ************************************************/

void sim_gnsf_set_nx(void *dims_, int nx)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    dims->nx = nx;
}

void sim_gnsf_set_nu(void *dims_, int nu)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    dims->nu = nu;
}

void sim_gnsf_set_nz(void *dims_, int nz)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    dims->nz = nz;
}

void sim_gnsf_get_nx(void *dims_, int *nx)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    *nx = dims->nx;
}

void sim_gnsf_get_nu(void *dims_, int *nu)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    *nu = dims->nu;
}

void sim_gnsf_get_nz(void *dims_, int *nz)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    *nz = dims->nz;
}

/************************************************
 * import functions
 ************************************************/

void sim_gnsf_import_matrices(sim_gnsf_dims *dims, gnsf_model *model,
                              external_function_generic *get_matrices_fun)
{
    // calling the external function
    ext_fun_arg_t ext_fun_type_in[1];
    void *ext_fun_in[1];
    ext_fun_arg_t ext_fun_type_out[10];
    void *ext_fun_out[10];

    ext_fun_type_in[0] = COLMAJ;
    ext_fun_in[0] = model->A;  // just to have some input;

    for (int ii = 0; ii < 10; ii++) {
        ext_fun_type_out[ii] = COLMAJ;
    }

    ext_fun_out[0] = model->A;
    ext_fun_out[1] = model->B;
    ext_fun_out[2] = model->C;
    ext_fun_out[3] = model->E;
    ext_fun_out[4] = model->L_x;
    ext_fun_out[5] = model->L_xdot;
    ext_fun_out[6] = model->L_z;
    ext_fun_out[7] = model->L_u;
    ext_fun_out[8] = model->A_LO;
    ext_fun_out[9] = model->c;

    get_matrices_fun->evaluate(get_matrices_fun, ext_fun_type_in, ext_fun_in, ext_fun_type_out,
                               ext_fun_out);
}

/************************************************
 * OUTDATED FUNCTIONS: the following functions were implemented
 * before some major changes were applied, maybe similar
 * functionalities will be implemented at some point:
 ************************************************/
// void gnsf_get_dims( sim_gnsf_dims *dims, casadi_function_t get_ints_fun)

// void gnsf_import_precomputed(sim_gnsf_dims* dims, gnsf_model *model, casadi_function_t
// But_KK_YY_ZZ_LO_fun) this function can be used instead of import_matrices + precompute (not
// working now, ATTENTION!!!)


/************************************************
 * opts
 ************************************************/

int sim_gnsf_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = NS_MAX;

    int size = 0;

    size += sizeof(sim_rk_opts);

    size += ns_max * ns_max * sizeof(double);  // A_mat
    size += ns_max * sizeof(double);           // b_vec
    size += ns_max * sizeof(double);           // c_vec

    int tmp0 = gauss_nodes_work_calculate_size(ns_max);
    int tmp1 = butcher_table_work_calculate_size(ns_max);
    int work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    size += work_size;  // work

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}

void *sim_gnsf_opts_assign(void *config_, void *dims, void *raw_memory)
{
    int ns_max = NS_MAX;

    char *c_ptr = (char *) raw_memory;

    sim_rk_opts *opts = (sim_rk_opts *) c_ptr;
    c_ptr += sizeof(sim_rk_opts);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(ns_max * ns_max, &opts->A_mat, &c_ptr);
    assign_and_advance_double(ns_max, &opts->b_vec, &c_ptr);
    assign_and_advance_double(ns_max, &opts->c_vec, &c_ptr);

    // work
    int tmp0 = gauss_nodes_work_calculate_size(ns_max);
    int tmp1 = butcher_table_work_calculate_size(ns_max);
    int work_size = tmp0 > tmp1 ? tmp0 : tmp1;
    opts->work = c_ptr;
    c_ptr += work_size;

    assert((char *) raw_memory + sim_gnsf_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}

void sim_gnsf_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_rk_opts *opts = opts_;

    opts->ns = 3;  // GL 3
    int ns = opts->ns;

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    // gauss collocation nodes
    gauss_nodes(ns, opts->c_vec, opts->work);

    // butcher tableau
    butcher_table(ns, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    // default options
    opts->newton_iter = 3;
    opts->scheme = NULL;
    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
    opts->jac_reuse = true;

    opts->output_z = false;
    opts->sens_algebraic = false;

    return;
}

void sim_gnsf_opts_update(void *config_, void *dims, void *opts_)
{
    sim_rk_opts *opts = opts_;

    int ns = opts->ns;

    assert(ns <= NS_MAX && "ns > NS_MAX!");

    // set tableau size
    opts->tableau_size = opts->ns;

    // gauss collocation nodes
    gauss_nodes(ns, opts->c_vec, opts->work);

    // butcher tableau
    butcher_table(ns, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    return;
}

/************************************************
 * model
 ************************************************/

int sim_gnsf_model_calculate_size(void *config, void *dims_)
{
    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;

    // necessary integers
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;

    int size = 0;
    size += sizeof(gnsf_model);

    // model defining matrices
    size += (nx1 + nz) * (nx1 + nu + n_out + (nx1 + nz)) * sizeof(double);  // A,B,C,E

    size += ny * (2 * nx1 + nz) * sizeof(double);  // L_x, L_xdot, L_z
    size += nuhat * nu * sizeof(double);           // L_u
    size += nx2 * nx2 * sizeof(double);            // A_LO
    size += (nx1 + nz) * sizeof(double);           // c

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}

void *sim_gnsf_model_assign(void *config, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;

    // necessary integers
    // int nx      = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    gnsf_model *model = (gnsf_model *) c_ptr;
    c_ptr += sizeof(gnsf_model);

    // assign model matrices
    assign_and_advance_double((nx1 + nz) * nx1, &model->A, &c_ptr);
    assign_and_advance_double((nx1 + nz) * nu, &model->B, &c_ptr);
    assign_and_advance_double((nx1 + nz) * n_out, &model->C, &c_ptr);
    assign_and_advance_double((nx1 + nz) * (nx1 + nz), &model->E, &c_ptr);

    assign_and_advance_double(ny * nx1, &model->L_x, &c_ptr);
    assign_and_advance_double(ny * nx1, &model->L_xdot, &c_ptr);
    assign_and_advance_double(ny * nz,  &model->L_z, &c_ptr);

    assign_and_advance_double(nuhat * nu, &model->L_u, &c_ptr);

    assign_and_advance_double(nx2 * nx2, &model->A_LO, &c_ptr);
    assign_and_advance_double(nx1 + nz , &model->c, &c_ptr);

    assert((char *) raw_memory + sim_gnsf_model_calculate_size(config, dims_) >= c_ptr);
    return model;
}

int sim_gnsf_model_set_function(void *model_, sim_function_t fun_type, void *fun)
{
    gnsf_model *model = model_;

    switch (fun_type)
    {
        case PHI_FUN:
            model->phi_fun = (external_function_generic *) fun;
            break;
        case PHI_FUN_JAC_Y:
            model->phi_fun_jac_y = (external_function_generic *) fun;
            break;
        case PHI_JAC_Y_UHAT:
            model->phi_jac_y_uhat = (external_function_generic *) fun;
            break;
        case LO_FUN:
            model->f_lo_fun_jac_x1_x1dot_u_z = (external_function_generic *) fun;
            break;
        default:
            return ACADOS_FAILURE;
    }
    return ACADOS_SUCCESS;
}

/************************************************
 * GNSF PRECOMPUTATION
 ************************************************/

static void *gnsf_cast_pre_workspace(void *config_, sim_gnsf_dims *dims_, void *opts_,
                                     void *raw_memory)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_rk_opts *opts = (sim_rk_opts *) opts_;
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;

    int num_stages = opts->ns;

    int nff = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;

    char *c_ptr = (char *) raw_memory;
    align_char_to(8, &c_ptr);

    // struct
    gnsf_pre_workspace *work = (gnsf_pre_workspace *) c_ptr;
    c_ptr += sizeof(gnsf_pre_workspace);

    assign_and_advance_int(nK1, &work->ipivEE1, &c_ptr);
    assign_and_advance_int(nZ, &work->ipivEE2, &c_ptr);
    assign_and_advance_int(nZ, &work->ipivQQ1, &c_ptr);
    assign_and_advance_int(nK2, &work->ipivM2, &c_ptr);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx1, nx1, &work->E11, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx1, nz, &work->E12, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, nx1, &work->E21, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, nz, &work->E22, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx1, nx1, &work->A1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, nx1, &work->A2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx1, nu, &work->B1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, nu, &work->B2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx1, n_out, &work->C1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, n_out, &work->C2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK1, nx1, &work->AA1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nx1, &work->AA2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nu, &work->BB1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nu, &work->BB2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK1, nff, &work->CC1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nff, &work->CC2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nZ, &work->DD1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nK1, &work->DD2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK1, nK1, &work->EE1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nZ, &work->EE2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nZ, nZ, &work->QQ1, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nyy, nZ, &work->LLZ, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nx1, &work->LLx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nK1, &work->LLK, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK2, nK2, &work->M2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nx2, &work->dK2_dx2_work, &c_ptr);

    if (opts->sens_algebraic){
        // for algebraic sensitivities
        assign_and_advance_blasfeo_dmat_mem(nx1, nx1  , &work->K0x, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx1, nu   , &work->K0u, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nx1, n_out, &work->K0f, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nz,  nz,    &work->Q1, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem(nK1, &work->cc1, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nZ , &work->cc2, &c_ptr);

    assert((char *) raw_memory + sim_gnsf_workspace_calculate_size(config_, dims, opts) >= c_ptr);
    return (void *) work;
}  // cast pre_workspace

void sim_gnsf_precompute(void *config, sim_gnsf_dims *dims, gnsf_model *model, sim_rk_opts *opts,
                         void *mem_, void *work_, double T)
{
    acados_timer atimer;
    acados_tic(&atimer);
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;
    int num_stages = opts->ns;
    int num_steps = opts->num_steps;

    int nff = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;

    // set up precomputation workspace
    gnsf_pre_workspace *work =
        (gnsf_pre_workspace *) gnsf_cast_pre_workspace(config, dims, opts, work_);
    // set up memory
    sim_gnsf_memory *mem = (sim_gnsf_memory *) mem_;

    double dt = T / num_steps;
    mem->dt = dt;

    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    double *c_vec = opts->c_vec;

    double *c_butcher = mem->c_butcher;
    double *b_dt = mem->b_dt;
    double *A_dt = mem->A_dt;

    double *A_LO = model->A_LO;

    // workspace stuff
    int *ipivEE1 = work->ipivEE1;
    int *ipivEE2 = work->ipivEE2;
    int *ipivQQ1 = work->ipivQQ1;
    int *ipivM2 = work->ipivM2;

    struct blasfeo_dmat E11 = work->E11;
    struct blasfeo_dmat E12 = work->E12;
    struct blasfeo_dmat E21 = work->E21;
    struct blasfeo_dmat E22 = work->E22;

    struct blasfeo_dmat A1 = work->A1;
    struct blasfeo_dmat A2 = work->A2;
    struct blasfeo_dmat B1 = work->B1;
    struct blasfeo_dmat B2 = work->B2;
    struct blasfeo_dmat C1 = work->C1;
    struct blasfeo_dmat C2 = work->C2;

    struct blasfeo_dmat AA1 = work->AA1;
    struct blasfeo_dmat AA2 = work->AA2;
    struct blasfeo_dmat BB1 = work->BB1;
    struct blasfeo_dmat BB2 = work->BB2;
    struct blasfeo_dmat CC1 = work->CC1;
    struct blasfeo_dmat CC2 = work->CC2;

    struct blasfeo_dvec cc1 = work->cc1;
    struct blasfeo_dvec cc2 = work->cc2;

    struct blasfeo_dmat DD1 = work->DD1;
    struct blasfeo_dmat DD2 = work->DD2;
    struct blasfeo_dmat EE1 = work->EE1;
    struct blasfeo_dmat EE2 = work->EE2;

    struct blasfeo_dmat LLZ = work->LLZ;
    struct blasfeo_dmat LLx = work->LLx;
    struct blasfeo_dmat LLK = work->LLK;

    struct blasfeo_dmat QQ1 = work->QQ1;

    struct blasfeo_dmat M2 = work->M2;
    struct blasfeo_dmat dK2_dx2_work = work->dK2_dx2_work;

    // memory - precomputed matrices
    struct blasfeo_dmat KKf = mem->KKf;
    struct blasfeo_dmat KKx = mem->KKx;
    struct blasfeo_dmat KKu = mem->KKu;

    struct blasfeo_dmat YYf = mem->YYf;
    struct blasfeo_dmat YYx = mem->YYx;
    struct blasfeo_dmat YYu = mem->YYu;

    struct blasfeo_dmat ZZf = mem->ZZf;
    struct blasfeo_dmat ZZx = mem->ZZx;
    struct blasfeo_dmat ZZu = mem->ZZu;

    struct blasfeo_dmat ALO = mem->ALO;
    struct blasfeo_dmat M2inv = mem->M2inv;
    struct blasfeo_dmat dK2_dx2 = mem->dK2_dx2;

    struct blasfeo_dmat Lu = mem->Lu;

    // ONLY for algebraic sensitivity propagation
    struct blasfeo_dmat Lx = mem->Lx;
    struct blasfeo_dmat Lxdot = mem->Lxdot;
    struct blasfeo_dmat Lz = mem->Lz;

    struct blasfeo_dmat K0x = work->K0x;
    struct blasfeo_dmat K0u = work->K0u;
    struct blasfeo_dmat K0f = work->K0f;
    struct blasfeo_dmat Q1  = work->Q1;

    struct blasfeo_dmat Z0x = mem->Z0x;
    struct blasfeo_dmat Z0u = mem->Z0u;
    struct blasfeo_dmat Z0f = mem->Z0f;

    struct blasfeo_dmat Y0x = mem->Y0x;
    struct blasfeo_dmat Y0u = mem->Y0u;
    struct blasfeo_dmat Y0f = mem->Y0f;

    // precomputed vectors
    struct blasfeo_dvec KK0 = mem->KK0;
    struct blasfeo_dvec YY0 = mem->YY0;
    struct blasfeo_dvec ZZ0 = mem->ZZ0;

    // set memory to zeros
    blasfeo_dgese(nK1, nff, 0.0, &CC1, 0, 0);
    blasfeo_dgese(nZ, nff, 0.0, &CC2, 0, 0);
    blasfeo_dgese(nK1, nZ, 0.0, &DD1, 0, 0);
    blasfeo_dgese(nZ, nK1, 0.0, &DD2, 0, 0);
    blasfeo_dgese(nK1, nK1, 0.0, &EE1, 0, 0);
    blasfeo_dgese(nZ, nZ, 0.0, &EE2, 0, 0);
    blasfeo_dgese(nyy, nK1, 0.0, &LLK, 0, 0);
    blasfeo_dgese(nK2, nK2, 0.0, &M2, 0, 0);
    blasfeo_dgese(nK2, nK2, 0.0, &M2inv, 0, 0);

    blasfeo_pack_dmat(nx1, nx1, model->E, nx1 + nz, &E11, 0, 0);

    blasfeo_pack_dmat(nx1, nz, &model->E[(nx1 + nz) * nx1], nx1 + nz, &E12, 0, 0);
    blasfeo_pack_dmat(nz, nx1, &model->E[nx1], nx1 + nz, &E21, 0, 0);
    blasfeo_pack_dmat(nz, nz, &model->E[nx1 + (nx1 + nz) * nx1], nx1 + nz, &E22, 0, 0);

    blasfeo_pack_dmat(nx1, nx1, model->A, nx1 + nz, &A1, 0, 0);
    blasfeo_pack_dmat(nz, nx1, &model->A[nx1], nx1 + nz, &A2, 0, 0);

    blasfeo_pack_dmat(nx1, nu, model->B, nx1 + nz, &B1, 0, 0);
    blasfeo_pack_dmat(nz, nu, &model->B[nx1], nx1 + nz, &B2, 0, 0);

    blasfeo_pack_dmat(nx1, n_out, model->C, nx1 + nz, &C1, 0, 0);
    blasfeo_pack_dmat(nz, n_out, &model->C[nx1], nx1 + nz, &C2, 0, 0);

    blasfeo_pack_dmat(nx2, nx2, A_LO, nx2, &ALO, 0, 0);
    blasfeo_pack_dmat(nuhat, nu, model->L_u, nuhat, &Lu, 0, 0);

    for (int ii = 0; ii < num_stages * num_stages; ii++)
    {
        A_dt[ii] = A_mat[ii] * dt;
    }
    for (int ii = 0; ii < num_stages; ii++)
    {
        b_dt[ii] = b_vec[ii] * dt;
        c_butcher[ii] = c_vec[ii];
    }

    // Build fat matrices AA1, AA2, ... , EE2, LLx, ..., LLZ, vectors cc1, cc2
    for (int ii = 0; ii < num_stages; ii++)
    {  // num_stages
        blasfeo_dgecp(nx1, nx1, &A1, 0, 0, &AA1, ii * nx1, 0);
        blasfeo_dgecp(nz, nx1, &A2, 0, 0, &AA2, ii * nz, 0);
        blasfeo_dgecp(nx1, nu, &B1, 0, 0, &BB1, ii * nx1, 0);
        blasfeo_dgecp(nz, nu, &B2, 0, 0, &BB2, ii * nz, 0);

        blasfeo_dgecp(nx1, n_out, &C1, 0, 0, &CC1, ii * nx1, ii * n_out);
        blasfeo_dgecp(nz, n_out, &C2, 0, 0, &CC2, ii * nz, ii * n_out);

        blasfeo_dgecpsc(nx1, nz, -1.0, &E12, 0, 0, &DD1, ii * nx1, ii * nz);
        blasfeo_dgecpsc(nz, nx1, -1.0, &E21, 0, 0, &DD2, ii * nz, ii * nx1);

        blasfeo_dgecp(nx1, nx1, &E11, 0, 0, &EE1, ii * nx1, ii * nx1);
        blasfeo_dgecp(nz, nz, &E22, 0, 0, &EE2, ii * nz, ii * nz);

        blasfeo_pack_dmat(ny, nz, model->L_z, ny, &LLZ, ii * ny, ii * nz);
        blasfeo_pack_dmat(ny, nx1, model->L_x, ny, &LLx, ii * ny, 0);
        blasfeo_pack_dmat(ny, nx1, model->L_xdot, ny, &LLK, ii * ny, ii * nx1);

        blasfeo_pack_dmat(nx2, nx2, A_LO, nx2, &dK2_dx2_work, ii * nx2,
                          0);  // dK2_dx2 = repmat(s.ALO,opts.n_stages,1);

        blasfeo_pack_dvec(nx1, model->c, &cc1, ii * nx1);
        blasfeo_pack_dvec(nz , &model->c[nx1], &cc2, ii * nz);
    }

    for (int ii = 0; ii < num_stages; ii++)
    {  // perform kronecker product
        for (int jj = 0; jj < num_stages; jj++)
        {
            blasfeo_dgead(nz, nx1, A_dt[ii * num_stages + jj], &A2, 0, 0, &DD2, jj * nz, ii * nx1);
            blasfeo_dgead(nx1, nx1, -A_dt[ii * num_stages + jj], &A1, 0, 0, &EE1, jj * nx1,
                          ii * nx1);
            blasfeo_dgead(ny, nx1, A_dt[ii * num_stages + jj], &LLx, 0, 0, &LLK, jj * ny, ii * nx1);
        }
    }

    /************************************************
     *   Compute QQ1 KK*, ZZ* via QQ1
     ************************************************/

    // SOLVE EE1 \ DD1, ... EE1 \ AA1;
    blasfeo_dgetrf_rowpivot(nK1, nK1, &EE1, 0, 0, &EE1, 0, 0, ipivEE1);  // factorize EE1

    blasfeo_drowpe(nK1, ipivEE1, &AA1);  // permute also rhs
    blasfeo_dtrsm_llnu(nK1, nx1, 1.0, &EE1, 0, 0, &AA1, 0, 0, &AA1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nx1, 1.0, &EE1, 0, 0, &AA1, 0, 0, &AA1, 0,
                       0);              // AA1 now contains EE1\AA1
    blasfeo_drowpe(nK1, ipivEE1, &BB1);  // permute also rhs
    blasfeo_dtrsm_llnu(nK1, nu, 1.0, &EE1, 0, 0, &BB1, 0, 0, &BB1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nu, 1.0, &EE1, 0, 0, &BB1, 0, 0, &BB1, 0,
                       0);               // BB1 now contains EE1\BB1
    blasfeo_drowpe(nK1, ipivEE1, &CC1);  // permute also rhs

    blasfeo_dtrsm_llnu(nK1, nff, 1.0, &EE1, 0, 0, &CC1, 0, 0, &CC1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nff, 1.0, &EE1, 0, 0, &CC1, 0, 0, &CC1, 0,
                       0);               // CC1 now contains EE1\CC1

    blasfeo_drowpe(nK1, ipivEE1, &DD1);  // permute also rhs
    blasfeo_dtrsm_llnu(nK1, nZ, 1.0, &EE1, 0, 0, &DD1, 0, 0, &DD1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nZ, 1.0, &EE1, 0, 0, &DD1, 0, 0, &DD1, 0,
                       0);  // DD1 now contains EE1\DD1

    blasfeo_dvecpe(nK1, ipivEE1, &cc1, 0);  // permute rhs
    blasfeo_dtrsv_lnu(nK1, &EE1, 0, 0, &cc1, 0, &cc1, 0);
    blasfeo_dtrsv_unn(nK1, &EE1, 0, 0, &cc1, 0, &cc1, 0);
    // cc1 now contains EE1\cc1


    // SOLVE EE2 \ DD2, ... EE2 \ AA2;
    blasfeo_dgetrf_rowpivot(nZ, nZ, &EE2, 0, 0, &EE2, 0, 0, ipivEE2);  // factorize EE2

    blasfeo_drowpe(nZ, ipivEE2, &AA2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ, nx1, 1.0, &EE2, 0, 0, &AA2, 0, 0, &AA2, 0, 0);
    blasfeo_dtrsm_lunn(nZ, nx1, 1.0, &EE2, 0, 0, &AA2, 0, 0, &AA2, 0, 0);
                                    // AA2 now contains EE2\AA2

    blasfeo_drowpe(nZ, ipivEE2, &BB2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ, nu, 1.0, &EE2, 0, 0, &BB2, 0, 0, &BB2, 0, 0);
    blasfeo_dtrsm_lunn(nZ, nu, 1.0, &EE2, 0, 0, &BB2, 0, 0, &BB2, 0, 0);
                                    // BB2 now contains EE2\BB2

    blasfeo_drowpe(nZ, ipivEE2, &CC2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ, nff, 1.0, &EE2, 0, 0, &CC2, 0, 0, &CC2, 0, 0);
    blasfeo_dtrsm_lunn(nZ, nff, 1.0, &EE2, 0, 0, &CC2, 0, 0, &CC2, 0, 0);
                                    // CC2 now contains EE2\CC2

    blasfeo_drowpe(nZ, ipivEE2, &DD2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ, nK1, 1.0, &EE2, 0, 0, &DD2, 0, 0, &DD2, 0, 0);
    blasfeo_dtrsm_lunn(nZ, nK1, 1.0, &EE2, 0, 0, &DD2, 0, 0, &DD2, 0, 0);
                                    // DD2 now contains EE2\DD2

    blasfeo_dvecpe(nZ, ipivEE2, &cc2, 0);  // permute rhs
    blasfeo_dtrsv_lnu(nZ, &EE2, 0, 0, &cc2, 0, &cc2, 0);
    blasfeo_dtrsv_unn(nZ, &EE2, 0, 0, &cc2, 0, &cc2, 0);
    // cc2 now contains EE2\cc2

    /* Build and factorize QQ1 */
    blasfeo_dgemm_nn(nZ, nZ, nK1, -1.0, &DD2, 0, 0, &DD1, 0, 0, 0.0, &QQ1, 0, 0, &QQ1, 0,
                     0);                                               // QQ1 = -DD2*DD1
    blasfeo_ddiare(nZ, 1.0, &QQ1, 0, 0);                               // add eye(nZ) to QQ1
    blasfeo_dgetrf_rowpivot(nZ, nZ, &QQ1, 0, 0, &QQ1, 0, 0, ipivQQ1);  // factorize QQ1

    /* build ZZf */
    blasfeo_dgemm_nn(nZ, nff, nK1, 1.0, &DD2, 0, 0, &CC1, 0, 0, 0.0, &ZZf, 0, 0, &ZZf, 0,
                     0);                                  // ZZf = DD2 * CC1
    blasfeo_dgead(nZ, nff, 1.0, &CC2, 0, 0, &ZZf, 0, 0);  // ZZf = ZZf + CC2;
    // solve QQ1\ZZf and store result in ZZf;
    blasfeo_drowpe(nZ, ipivQQ1, &ZZf);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ, nff, 1.0, &QQ1, 0, 0, &ZZf, 0, 0, &ZZf, 0, 0);
    blasfeo_dtrsm_lunn(nZ, nff, 1.0, &QQ1, 0, 0, &ZZf, 0, 0, &ZZf, 0,
                       0);  // ZZf now contains QQ1\ZZf

    /* build ZZu */
    blasfeo_dgemm_nn(nZ, nu, nK1, 1.0, &DD2, 0, 0, &BB1, 0, 0, 0.0, &ZZu, 0, 0, &ZZu, 0,
                     0);                                 // ZZu = DD2 * BB1
    blasfeo_dgead(nZ, nu, 1.0, &BB2, 0, 0, &ZZu, 0, 0);  // ZZu = ZZu + BB2;
    // solve QQ1\ZZu and store result in ZZu;
    blasfeo_drowpe(nZ, ipivQQ1, &ZZu);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ, nu, 1.0, &QQ1, 0, 0, &ZZu, 0, 0, &ZZu, 0, 0);
    blasfeo_dtrsm_lunn(nZ, nu, 1.0, &QQ1, 0, 0, &ZZu, 0, 0, &ZZu, 0,
                       0);  // ZZu now contains QQ1\ZZu

    /*  build ZZx  */
    blasfeo_dgemm_nn(nZ, nx1, nK1, 1.0, &DD2, 0, 0, &AA1, 0, 0, 0.0, &ZZx, 0, 0, &ZZx, 0,
                     0);                                  // ZZx = DD2 * AA1;
    blasfeo_dgead(nZ, nx1, 1.0, &AA2, 0, 0, &ZZx, 0, 0);  // ZZx = ZZx + AA2;

    // solve QQ1\ZZx and store result in ZZx;
    blasfeo_drowpe(nZ, ipivQQ1, &ZZu);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ, nx1, 1.0, &QQ1, 0, 0, &ZZx, 0, 0, &ZZx, 0, 0);
    blasfeo_dtrsm_lunn(nZ, nx1, 1.0, &QQ1, 0, 0, &ZZx, 0, 0, &ZZx, 0,
                       0);  // ZZx now contains QQ1\ZZx

    /* build ZZ0  */
    blasfeo_dgemv_n(nZ, nK1, 1.0, &DD2, 0, 0, &cc1, 0, 1.0, &cc2, 0, &ZZ0, 0);
                                                // ZZ0 = DD2 * cc1 + cc2;
    // solve QQ1\ZZ0 and store result in ZZ0;
    blasfeo_dvecpe(nZ, ipivQQ1, &ZZ0, 0);  // permute rhs
    blasfeo_dtrsv_lnu(nZ, &QQ1, 0, 0, &ZZ0, 0, &ZZ0, 0);
    blasfeo_dtrsv_unn(nZ, &QQ1, 0, 0, &ZZ0, 0, &ZZ0, 0);

    /* build KKf, KKu, KKx */
    blasfeo_dgemm_nn(nK1, nff, nZ, 1.0, &DD1, 0, 0, &ZZf, 0, 0, 1.0, &CC1, 0, 0, &KKf, 0,
                     0);  // KKf = DD1 * ZZf + CC1
    blasfeo_dgemm_nn(nK1, nu, nZ, 1.0, &DD1, 0, 0, &ZZu, 0, 0, 1.0, &BB1, 0, 0, &KKu, 0,
                     0);  // KKu = DD1 * ZZu + BB1
    blasfeo_dgemm_nn(nK1, nx1, nZ, 1.0, &DD1, 0, 0, &ZZx, 0, 0, 1.0, &AA1, 0, 0, &KKx, 0,
                     0);  // KKx = DD1 * ZZx + AA1
    blasfeo_dgemv_n(nK1, nZ, 1.0, &DD1, 0, 0, &ZZ0, 0, 1.0, &cc1, 0, &KK0, 0);
                          // KK0 = DD1 * ZZ0 + cc1;

    /* build YYx */
    blasfeo_dgemm_nn(nyy, nx1, nK1, 1.0, &LLK, 0, 0, &KKx, 0, 0, 1.0, &LLx, 0, 0, &YYx, 0, 0);
    blasfeo_dgemm_nn(nyy, nx1, nZ, 1.0, &LLZ, 0, 0, &ZZx, 0, 0, 1.0, &YYx, 0, 0, &YYx, 0, 0);

    /* build YYu */
    blasfeo_dgemm_nn(nyy, nu, nK1, 1.0, &LLK, 0, 0, &KKu, 0, 0, 0.0, &YYu, 0, 0, &YYu, 0, 0);
    blasfeo_dgemm_nn(nyy, nu, nZ, 1.0, &LLZ, 0, 0, &ZZu, 0, 0, 1.0, &YYu, 0, 0, &YYu, 0, 0);
    // printf("YYu = (in precompute) \n");
    // blasfeo_print_exp_dmat(nyy, nu, &YYu, 0, 0);

    /* build YYf */
    blasfeo_dgemm_nn(nyy, nff, nK1, 1.0, &LLK, 0, 0, &KKf, 0, 0, 0.0, &YYf, 0, 0, &YYf, 0, 0);
    blasfeo_dgemm_nn(nyy, nff, nZ, 1.0, &LLZ, 0, 0, &ZZf, 0, 0, 1.0, &YYf, 0, 0, &YYf, 0, 0);

    /* build YY0 */
    blasfeo_dgemv_n(nyy, nK1, 1.0, &LLK, 0, 0, &KK0, 0, 0.0, &YY0, 0, &YY0, 0);
                                                // YY0 = LLK * KK0;
    blasfeo_dgemv_n(nyy, nZ, 1.0, &LLZ, 0, 0, &ZZ0, 0, 1.0, &YY0, 0, &YY0, 0);
                                                // YY0 = YY0 + LLZ * ZZ0

    // build M2_inv
    blasfeo_ddiare(nK2, 1.0, &M2inv, 0, 0);  // add eye(nK2) to M2inv
    blasfeo_ddiare(nK2, 1.0, &M2, 0, 0);     // add eye(nK2) to M2, both are eye now
    for (int ii = 0; ii < num_stages; ii++)
    {
        for (int jj = 0; jj < num_stages; jj++)
        {
            blasfeo_dgead(nx2, nx2, - A_dt[ii * num_stages + jj], &ALO, 0, 0, &M2, jj * nx2,
                          ii * nx2);
        }
    }
    blasfeo_dgetrf_rowpivot(nK2, nK2, &M2, 0, 0, &M2, 0, 0, ipivM2);  // factorize M2

    // solve M2\M2inv, result is M2inv
    blasfeo_drowpe(nK2, ipivM2, &M2inv);  // permute also rhs
    blasfeo_dtrsm_llnu(nK2, nK2, 1.0, &M2, 0, 0, &M2inv, 0, 0, &M2inv, 0, 0);
    blasfeo_dtrsm_lunn(nK2, nK2, 1.0, &M2, 0, 0, &M2inv, 0, 0, &M2inv, 0, 0);
                  // M2inv now contains inv(M2)

    blasfeo_dgemm_nn(nK2, nx2, nK2, 1.0, &M2inv, 0, 0, &dK2_dx2_work, 0, 0, 0.0, &dK2_dx2, 0, 0,
                     &dK2_dx2, 0, 0);

    if (opts->sens_algebraic){

        /* ONLY FOR ALGEBRAIC SENSITIVITIES */
        /* Compute K0*, Z0*, Y0*, * \in {f,x,u} */
        // pack matrices to memory for algebraic sensitivities
        blasfeo_pack_dmat(ny, nx1, model->L_x, ny, &Lx, 0, 0);
        blasfeo_pack_dmat(ny, nx1, model->L_xdot, ny, &Lxdot, 0, 0);
        blasfeo_pack_dmat(ny, nz , model->L_z, ny, &Lz, 0, 0);

        // SOLVE E11 \ A1,   E11 \ B1,      E11 \ C1,   E11\E12;
        blasfeo_dgetrf_rowpivot(nx1, nx1, &E11, 0, 0, &E11, 0, 0, ipivEE1);  // factorize E11

        blasfeo_drowpe(nx1, ipivEE1, &A1);  // permute also rhs
        blasfeo_dtrsm_llnu(nx1, nx1, 1.0, &E11, 0, 0, &A1, 0, 0, &A1, 0, 0);
        blasfeo_dtrsm_lunn(nx1, nx1, 1.0, &E11, 0, 0, &A1, 0, 0, &A1, 0,
                        0);              // A1 now contains E11\A1
        blasfeo_drowpe(nx1, ipivEE1, &B1);  // permute also rhs
        blasfeo_dtrsm_llnu(nx1, nu, 1.0, &E11, 0, 0, &B1, 0, 0, &B1, 0, 0);
        blasfeo_dtrsm_lunn(nx1, nu, 1.0, &E11, 0, 0, &B1, 0, 0, &B1, 0,
                        0);               // B1 now contains E11\B1

        blasfeo_drowpe(nx1, ipivEE1, &C1);  // permute also rhs
        blasfeo_dtrsm_llnu(nx1, n_out, 1.0, &E11, 0, 0, &C1, 0, 0, &C1, 0, 0);
        blasfeo_dtrsm_lunn(nx1, n_out, 1.0, &E11, 0, 0, &C1, 0, 0, &C1, 0,
                        0);               // C1 now contains E11\C1

        blasfeo_drowpe(nx1, ipivEE1, &E12);  // permute also rhs
        blasfeo_dtrsm_llnu(nx1, nz, 1.0, &E11, 0, 0, &E12, 0, 0, &E12, 0, 0);
        blasfeo_dtrsm_lunn(nx1, nz, 1.0, &E11, 0, 0, &E12, 0, 0, &E12, 0,
                        0);  // E12 now contains E11\E12


        // SOLVE E22 \ A2,   E22 \ B2,      E22 \ C2,   E22\E21;
        blasfeo_dgetrf_rowpivot(nz, nz, &E22, 0, 0, &E22, 0, 0, ipivEE2);  // factorize E22

        blasfeo_drowpe(nz, ipivEE2, &A2);  // permute also rhs
        blasfeo_dtrsm_llnu(nz, nx1, 1.0, &E22, 0, 0, &A2, 0, 0, &A2, 0, 0);
        blasfeo_dtrsm_lunn(nz, nx1, 1.0, &E22, 0, 0, &A2, 0, 0, &A2, 0, 0);
                                        // A2 now contains E22\A2

        blasfeo_drowpe(nz, ipivEE2, &B2);  // permute also rhs
        blasfeo_dtrsm_llnu(nz, nu, 1.0, &E22, 0, 0, &B2, 0, 0, &B2, 0, 0);
        blasfeo_dtrsm_lunn(nz, nu, 1.0, &E22, 0, 0, &B2, 0, 0, &B2, 0, 0);
                                        // B2 now contains E22\B2

        blasfeo_drowpe(nz, ipivEE2, &C2);  // permute also rhs
        blasfeo_dtrsm_llnu(nz, n_out, 1.0, &E22, 0, 0, &C2, 0, 0, &C2, 0, 0);
        blasfeo_dtrsm_lunn(nz, n_out, 1.0, &E22, 0, 0, &C2, 0, 0, &C2, 0, 0);
                                        // C2 now contains E22\C2

        blasfeo_drowpe(nz, ipivEE2, &E21);  // permute also rhs
        blasfeo_dtrsm_llnu(nz, nx1, 1.0, &E22, 0, 0, &E21, 0, 0, &E21, 0, 0);
        blasfeo_dtrsm_lunn(nz, nx1, 1.0, &E22, 0, 0, &E21, 0, 0, &E21, 0, 0);
                                        // E21 now contains E22\E21

        /* Build and factorize Q1 */
        blasfeo_dgemm_nn(nz, nz, nx1, -1.0, &E21, 0, 0, &E12, 0, 0, 0.0, &Q1, 0, 0, &Q1, 0,
                        0);                                               // Q1 = -E21*E12
        blasfeo_ddiare(nz, 1.0, &Q1, 0, 0);                               // add eye(nz) to Q1
        blasfeo_dgetrf_rowpivot(nz, nz, &Q1, 0, 0, &Q1, 0, 0, ipivQQ1);  // factorize Q1

        /* build Z0f */
        blasfeo_dgemm_nn(nz, n_out, nx1, 1.0, &E21, 0, 0, &C1, 0, 0, 1.0, &C2, 0, 0, &Z0f, 0,
                        0);                                  // Z0f = E21 * C1 + C2
        // blasfeo_dgead(nz, n_out, 1.0, &C2, 0, 0, &Z0f, 0, 0);  // Z0f = Z0f + C2;

        // solve Q1\Z0f and store result in Z0f;
        blasfeo_drowpe(nz, ipivQQ1, &Z0f);  // permute also rhs
        blasfeo_dtrsm_llnu(nz, n_out, 1.0, &Q1, 0, 0, &Z0f, 0, 0, &Z0f, 0, 0);
        blasfeo_dtrsm_lunn(nz, n_out, 1.0, &Q1, 0, 0, &Z0f, 0, 0, &Z0f, 0,
                        0);  // Z0f now contains Q1\Z0f

        /* build Z0u */
        blasfeo_dgemm_nn(nz, nu, nx1, 1.0, &E21, 0, 0, &B1, 0, 0, 0.0, &Z0u, 0, 0, &Z0u, 0,
                        0);                                 // Z0u = E21 * B1
        blasfeo_dgead(nz, nu, 1.0, &B2, 0, 0, &Z0u, 0, 0);  // Z0u = Z0u + B2;
        // solve Q1\Z0u and store result in Z0u;
        blasfeo_drowpe(nz, ipivQQ1, &Z0u);  // permute also rhs
        blasfeo_dtrsm_llnu(nz, nu, 1.0, &Q1, 0, 0, &Z0u, 0, 0, &Z0u, 0, 0);
        blasfeo_dtrsm_lunn(nz, nu, 1.0, &Q1, 0, 0, &Z0u, 0, 0, &Z0u, 0,
                        0);  // Z0u now contains Q1\Z0u

        /*  build Z0x  */
        blasfeo_dgemm_nn(nz, nx1, nx1, 1.0, &E21, 0, 0, &A1, 0, 0, 0.0, &Z0x, 0, 0, &Z0x, 0,
                        0);                                  // Z0x = E21 * A1;
        blasfeo_dgead(nz, nx1, 1.0, &A2, 0, 0, &Z0x, 0, 0);  // Z0x = Z0x + A2;

        // solve Q1\Z0x and store result in Z0x;
        blasfeo_drowpe(nz, ipivQQ1, &Z0u);  // permute also rhs
        blasfeo_dtrsm_llnu(nz, nx1, 1.0, &Q1, 0, 0, &Z0x, 0, 0, &Z0x, 0, 0);
        blasfeo_dtrsm_lunn(nz, nx1, 1.0, &Q1, 0, 0, &Z0x, 0, 0, &Z0x, 0,
                        0);  // Z0x now contains Q1\Z0x

        /* build K0f, K0u, K0x */
        blasfeo_dgemm_nn(nx1, n_out, nz, -1.0, &E12, 0, 0, &Z0f, 0, 0, 1.0, &C1, 0, 0, &K0f, 0,
                        0);  // K0f = -E12 * Z0f + C1
        blasfeo_dgemm_nn(nx1, nu, nz, -1.0, &E12, 0, 0, &Z0u, 0, 0, 1.0, &B1, 0, 0, &K0u, 0,
                        0);  // K0u = -E12 * Z0u + B1
        blasfeo_dgemm_nn(nx1, nx1, nz, -1.0, &E12, 0, 0, &Z0x, 0, 0, 1.0, &A1, 0, 0, &K0x, 0,
                        0);  // K0x = -E12 * Z0x + A1

        /* build Y0x */
        blasfeo_dgemm_nn(ny, nx1, nx1, 1.0, &Lxdot, 0, 0, &K0x, 0, 0, 1.0, &Lx, 0, 0, &Y0x, 0, 0);
        blasfeo_dgemm_nn(ny, nx1, nz, 1.0, &Lz, 0, 0, &Z0x, 0, 0, 1.0, &Y0x, 0, 0, &Y0x, 0, 0);

        /* build Y0u */
        blasfeo_dgemm_nn(ny, nu, nx1, 1.0, &Lxdot, 0, 0, &K0u, 0, 0, 0.0, &Y0u, 0, 0, &Y0u, 0, 0);
        blasfeo_dgemm_nn(ny, nu, nz, 1.0, &Lz, 0, 0, &Z0u, 0, 0, 1.0, &Y0u, 0, 0, &Y0u, 0, 0);

        /* build Y0f */
        blasfeo_dgemm_nn(ny, n_out, nx1, 1.0, &Lxdot, 0, 0, &K0f, 0, 0, 0.0, &Y0f,
                         0, 0, &Y0f, 0, 0);
        blasfeo_dgemm_nn(ny, n_out, nz, 1.0, &Lz, 0, 0, &Z0f, 0, 0, 1.0, &Y0f, 0, 0, &Y0f, 0, 0);

    }

    // double precomputation_time = acados_toc(&atimer) * 1000;
    // printf("time 2 precompute = %f [ms]\n", precomputation_time);
}

/************************************************
 * memory & workspace
 ************************************************/

int sim_gnsf_memory_calculate_size(void *config, void *dims_, void *opts_)
{
    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_rk_opts *opts = opts_;

    // necessary integers
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;

    int num_stages = opts->ns;

    int nff = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ = num_stages * nz;

    int size = sizeof(sim_gnsf_memory);

    // scaled butcher table
    size += num_stages * num_stages * sizeof(double);  // A_dt
    size += 2 * num_stages * sizeof(double);           // b_dt, c_butcher;

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    // precomputed matrices
    size += blasfeo_memsize_dmat(nK1, nff);  // KKf
    size += blasfeo_memsize_dmat(nK1, nx1);  // KKx
    size += blasfeo_memsize_dmat(nK1, nu);   // KKu

    size += blasfeo_memsize_dmat(nyy, nff);  // YYf
    size += blasfeo_memsize_dmat(nyy, nx1);  // YYx
    size += blasfeo_memsize_dmat(nyy, nu);   // YYu

    size += blasfeo_memsize_dmat(nZ, nff);  // ZZf
    size += blasfeo_memsize_dmat(nZ, nx1);  // ZZx
    size += blasfeo_memsize_dmat(nZ, nu);   // ZZu

    size += blasfeo_memsize_dmat(nx2, nx2);  // ALO
    size += blasfeo_memsize_dmat(nK2, nK2);  // M2inv
    size += blasfeo_memsize_dmat(nK2, nx2);  // dK2_dx2

    size += blasfeo_memsize_dmat(nuhat, nu);  // Lu

    if (opts->sens_algebraic){
        // matrices only needed for algebraic sensitivities
        size += blasfeo_memsize_dmat(nz, nx1);       // Z0x
        size += blasfeo_memsize_dmat(nz, nu);        // Z0u
        size += blasfeo_memsize_dmat(nz, n_out);     // Z0f

        size += blasfeo_memsize_dmat(ny, nx1);       // Y0x
        size += blasfeo_memsize_dmat(ny, nu);        // Y0u
        size += blasfeo_memsize_dmat(ny, n_out);     // Y0f

        size += 2 * blasfeo_memsize_dmat(ny, nx1);   // Lx, Lxdot
        size += blasfeo_memsize_dmat(ny, nz);        // Lz
    }

    size += blasfeo_memsize_dvec(nZ);   // ZZ0
    size += blasfeo_memsize_dvec(nK1);  // KK0
    size += blasfeo_memsize_dvec(nyy);  // YY0

    return size;
}

void *sim_gnsf_memory_assign(void *config, void *dims_, void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_rk_opts *opts = opts_;

    // necessary integers
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;

    int num_stages = opts->ns;

    int nff = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ = num_stages * nz;

    // struct
    sim_gnsf_memory *mem = (sim_gnsf_memory *) c_ptr;
    c_ptr += sizeof(sim_gnsf_memory);

    // assign scaled butcher table
    assign_and_advance_double(num_stages * num_stages, &mem->A_dt, &c_ptr);
    assign_and_advance_double(num_stages, &mem->b_dt, &c_ptr);
    assign_and_advance_double(num_stages, &mem->c_butcher, &c_ptr);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dmat_mem
    assign_and_advance_blasfeo_dmat_mem(nK1, nff, &mem->KKf, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nx1, &mem->KKx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nu, &mem->KKu, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nyy, nff, &mem->YYf, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nx1, &mem->YYx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nu, &mem->YYu, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nZ, nff, &mem->ZZf, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nx1, &mem->ZZx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nu, &mem->ZZu, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx2, nx2, &mem->ALO, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nK2, &mem->M2inv, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nx2, &mem->dK2_dx2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nuhat, nu, &mem->Lu, &c_ptr);

    if (opts->sens_algebraic){
        // for algebraic sensitivity propagation
        assign_and_advance_blasfeo_dmat_mem(ny, nx1, &mem->Lx, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(ny, nx1, &mem->Lxdot, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(ny, nz , &mem->Lz, &c_ptr);

        assign_and_advance_blasfeo_dmat_mem(nz,   nx1, &mem->Z0x, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nz,   nu , &mem->Z0u, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nz, n_out, &mem->Z0f, &c_ptr);

        assign_and_advance_blasfeo_dmat_mem(ny,   nx1, &mem->Y0x, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(ny,   nu , &mem->Y0u, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(ny, n_out, &mem->Y0f, &c_ptr);
    }

    assign_and_advance_blasfeo_dvec_mem(nZ , &mem->ZZ0, &c_ptr);  // ZZ0
    assign_and_advance_blasfeo_dvec_mem(nyy, &mem->YY0, &c_ptr);  // YY0
    assign_and_advance_blasfeo_dvec_mem(nK1, &mem->KK0, &c_ptr);  // KK0


    assert((char *) raw_memory + sim_gnsf_memory_calculate_size(config, dims_, opts_) >= c_ptr);
    return mem;
}

int sim_gnsf_workspace_calculate_size(void *config, void *dims_, void *opts_)
{
    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_rk_opts *opts = opts_;

    // necessary integers
    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;

    int num_stages = opts->ns;
    int num_steps = opts->num_steps;

    int nff = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ = num_stages * nz;

    /* Calculate workspace size for precompute function */
    int pre_size = sizeof(gnsf_pre_workspace);

    make_int_multiple_of(8, &pre_size);
    pre_size += 1 * 8;

    pre_size += (2 * nZ + nK1) * sizeof(int);  // ipivEE1, ipivEE2, ipivQQ1
    pre_size += nK2 * sizeof(int);             // ipivM2

    make_int_multiple_of(64, &pre_size);
    pre_size += 1 * 64;

    pre_size += blasfeo_memsize_dmat(nx1, nx1);  // E11
    pre_size += blasfeo_memsize_dmat(nx1, nz);   // E12
    pre_size += blasfeo_memsize_dmat(nz, nx1);   // E21
    pre_size += blasfeo_memsize_dmat(nz, nz);    // E22

    pre_size += blasfeo_memsize_dmat(nx1, nx1);    // A1
    pre_size += blasfeo_memsize_dmat(nz, nx1);     // A2
    pre_size += blasfeo_memsize_dmat(nx1, nu);     // B1
    pre_size += blasfeo_memsize_dmat(nz, nu);      // B2
    pre_size += blasfeo_memsize_dmat(nx1, n_out);  // C1
    pre_size += blasfeo_memsize_dmat(nz, n_out);   // C2

    pre_size += blasfeo_memsize_dmat(nK1, nx1);  // AA1
    pre_size += blasfeo_memsize_dmat(nZ, nx1);   // AA2
    pre_size += blasfeo_memsize_dmat(nK1, nu);   // BB1
    pre_size += blasfeo_memsize_dmat(nZ, nu);    // BB2

    pre_size += blasfeo_memsize_dmat(nK1, nff);  // CC1
    pre_size += blasfeo_memsize_dmat(nZ, nff);   // CC2
    pre_size += blasfeo_memsize_dmat(nK1, nZ);   // DD1
    pre_size += blasfeo_memsize_dmat(nZ, nK1);   // DD2

    pre_size += blasfeo_memsize_dmat(nK1, nK1);  // EE1
    pre_size += blasfeo_memsize_dmat(nZ, nZ);    // EE2

    pre_size += blasfeo_memsize_dmat(nZ, nZ);  // QQ1

    pre_size += blasfeo_memsize_dmat(nyy, nZ);   // LLZ
    pre_size += blasfeo_memsize_dmat(nyy, nx1);  // LLx
    pre_size += blasfeo_memsize_dmat(nyy, nK1);  // LLK

    pre_size += blasfeo_memsize_dmat(nK2, nK2);  // M2
    pre_size += blasfeo_memsize_dmat(nK2, nx2);  // dK2_dx2_work

    pre_size += blasfeo_memsize_dvec(nK1);  // cc1
    pre_size += blasfeo_memsize_dvec(nZ);  // cc2

    if (opts->sens_algebraic){
        pre_size += blasfeo_memsize_dmat(nx1, nx1);  // K0x
        pre_size += blasfeo_memsize_dmat(nx1, nu);  // K0u
        pre_size += blasfeo_memsize_dmat(nx1, n_out);  // K0f
        pre_size += blasfeo_memsize_dmat(nz, nz);  // Q1
    }

    make_int_multiple_of(8, &pre_size);
    pre_size += 1 * 8;

    /* Calculate workspace size for simulation function */
    int size = sizeof(gnsf_workspace);
    make_int_multiple_of(8, &size);
    size += 1 * 8;

    size += num_stages * sizeof(double);  // Z_work

    size += 2 * num_steps * sizeof(struct blasfeo_dvec);  // ff_traj, yy_traj
    size += num_steps * sizeof(struct blasfeo_dmat);  // f_LO_jac_traj

    size += nff * sizeof(int);  // ipiv

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    size += num_steps * blasfeo_memsize_dvec(nff);      // ff_traj
    size += num_steps * blasfeo_memsize_dvec(nyy);      // yy_traj

    size += blasfeo_memsize_dvec(nZ);                    // Z_val
    size += 2 * blasfeo_memsize_dvec(nK1);               // K1_val, x1_stage_val
    size += blasfeo_memsize_dvec(nK2);                   // f_LO_val
    size += blasfeo_memsize_dvec(nK2);                   // K2_val
    size += blasfeo_memsize_dvec(nx * (num_steps + 1));  // x0_traj
    size += blasfeo_memsize_dvec(nff);                   // res_val
    size += blasfeo_memsize_dvec(nu);                    // u0
    size += blasfeo_memsize_dvec(nx + nu);               // lambda
    size += blasfeo_memsize_dvec(nx + nu);               // lambda_old

    size += blasfeo_memsize_dvec(nyy);              // yyu
    size += blasfeo_memsize_dvec(nyy * num_steps);  // yyss
    size += blasfeo_memsize_dvec(ny);  // y_one_stage

    size += blasfeo_memsize_dvec(nK1);  // K1u
    size += blasfeo_memsize_dvec(nZ);   // Zu
    size += blasfeo_memsize_dvec(nx2);  // ALOtimesx02

    size += blasfeo_memsize_dvec(nuhat);  // uhat

    if (opts->sens_algebraic){
        size += blasfeo_memsize_dvec(nx1);  // x0dot_1;
        size += blasfeo_memsize_dvec(nz);  // z0;
        size += n_out * sizeof(int);  // ipiv_ff0
    }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    size += num_steps * blasfeo_memsize_dmat(nK2, 2 * nx1 + nu + nz);  // f_LO_jac_traj

    size += blasfeo_memsize_dmat(nff, nff);       // J_r_ff
    size += blasfeo_memsize_dmat(nff, nx1 + nu);  // J_r_x1u

    if (opts->sens_algebraic){
        size += blasfeo_memsize_dmat(n_out, n_out);  // dr0_dff0
        size += blasfeo_memsize_dmat(nz, nx1 + nu);  // dz0_dx1u
    }

    size += blasfeo_memsize_dmat(nK1, nx1);  // dK1_dx1
    size += blasfeo_memsize_dmat(nK1, nu);   // dK1_du
    size += blasfeo_memsize_dmat(nZ, nx1);   // dZ_dx1
    size += blasfeo_memsize_dmat(nZ, nu);    // dZ_du
    size += blasfeo_memsize_dmat(nK2, nx1);  // aux_G2_x1
    size += blasfeo_memsize_dmat(nK2, nu);   // aux_G2_u
    size += blasfeo_memsize_dmat(nK2, nK1);  // J_G2_K1

    size += blasfeo_memsize_dmat(nK2, nx1);         // dK2_dx1
    size += blasfeo_memsize_dmat(nK2, nu);          // dK2_du
    size += blasfeo_memsize_dmat(nK2, nff);         // dK2_dff
    size += blasfeo_memsize_dmat(nx, nx + nu);      // dxf_dwn
    size += 2 * blasfeo_memsize_dmat(nx, nx + nu);  // S_forw_new, S_forw

    size += blasfeo_memsize_dmat(nK2, nff);  // aux_G2_ff
    size += blasfeo_memsize_dmat(nx, nff);   // dPsi_dff
    size += blasfeo_memsize_dmat(nx, nx);    // dPsi_dx
    size += blasfeo_memsize_dmat(nx, nu);    // dPsi_du

    size += blasfeo_memsize_dmat(nff, ny + nuhat);  // dPHI_dyuhat

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    /* take maximum of both workspace sizes */
    size = (size > pre_size) ? size : pre_size;

    return size;
}

static void *sim_gnsf_cast_workspace(void *config, void *dims_, void *opts_, void *raw_memory)
{
    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_rk_opts *opts = opts_;

    // necessary integers
    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;

    int num_stages = opts->ns;
    int num_steps = opts->num_steps;

    int nff = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ = num_stages * nz;

    char *c_ptr = (char *) raw_memory;
    gnsf_workspace *workspace = (gnsf_workspace *) c_ptr;
    c_ptr += sizeof(gnsf_workspace);
    align_char_to(8, &c_ptr);

    assign_and_advance_double(num_stages, &workspace->Z_work, &c_ptr);

    assign_and_advance_int(nff, &workspace->ipiv, &c_ptr);

    align_char_to(8, &c_ptr);

    assign_and_advance_blasfeo_dmat_structs(num_steps, &workspace->f_LO_jac_traj, &c_ptr);

    assign_and_advance_blasfeo_dvec_structs(num_steps, &workspace->ff_traj, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(num_steps, &workspace->yy_traj, &c_ptr);

    for (int ii = 0; ii < num_steps; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nff, workspace->ff_traj + ii, &c_ptr);    // ff_traj
        assign_and_advance_blasfeo_dvec_mem(nyy, workspace->yy_traj + ii, &c_ptr);    // yy_traj
    }

    assign_and_advance_blasfeo_dvec_mem(nK2, &workspace->f_LO_val, &c_ptr);  // f_LO_val
    assign_and_advance_blasfeo_dvec_mem(nK1, &workspace->K1_val, &c_ptr);  // K1_val
    assign_and_advance_blasfeo_dvec_mem(nK1, &workspace->x1_stage_val, &c_ptr);  // x1_stage_val
    assign_and_advance_blasfeo_dvec_mem(nZ,  &workspace->Z_val, &c_ptr);      // Z_val

    assign_and_advance_blasfeo_dvec_mem(nK2, &workspace->K2_val, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem((num_steps + 1) * nx, &workspace->x0_traj, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nff, &workspace->res_val, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nu, &workspace->u0, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, &workspace->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, &workspace->lambda_old, &c_ptr);

    assign_and_advance_blasfeo_dvec_mem(nyy, &workspace->yyu, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nyy * num_steps, &workspace->yyss, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(ny, &workspace->y_one_stage, &c_ptr);

    assign_and_advance_blasfeo_dvec_mem(nK1, &workspace->K1u, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nZ, &workspace->Zu, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx2, &workspace->ALOtimesx02, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nuhat, &workspace->uhat, &c_ptr);

    if (opts->sens_algebraic){
        assign_and_advance_blasfeo_dvec_mem(nx1, &workspace->x0dot_1, &c_ptr);
        assign_and_advance_blasfeo_dvec_mem(nz , &workspace->z0, &c_ptr);

        assign_and_advance_int(n_out, &workspace->ipiv_ff0, &c_ptr);

    }

    // blasfeo_dmat_mem align
    align_char_to(64, &c_ptr);
    for (int ii = 0; ii < num_steps; ii++)
    {
        assign_and_advance_blasfeo_dmat_mem(nK2, 2 * nx1 + nu + nz, workspace->f_LO_jac_traj + ii,
                                            &c_ptr);  // f_LO_jac_traj
    }

    assign_and_advance_blasfeo_dmat_mem(nff, nx1 + nu, &workspace->J_r_x1u, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nff, nff, &workspace->J_r_ff, &c_ptr);

    if (opts->sens_algebraic){
        assign_and_advance_blasfeo_dmat_mem(nz, nx1 + nu, &workspace->dz0_dx1u, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(n_out, n_out, &workspace->dr0_dff0, &c_ptr);
    }

    assign_and_advance_blasfeo_dmat_mem(nK1, nx1, &workspace->dK1_dx1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nu, &workspace->dK1_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nx1, &workspace->dZ_dx1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ, nu, &workspace->dZ_du, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK2, nx1, &workspace->aux_G2_x1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nu, &workspace->aux_G2_u, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nK1, &workspace->J_G2_K1, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nff, ny + nuhat, &workspace->dPHI_dyuhat, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK2, nx1, &workspace->dK2_dx1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nu, &workspace->dK2_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nff, &workspace->dK2_dff, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->dxf_dwn, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->S_forw_new, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->S_forw, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK2, nff, &workspace->aux_G2_ff, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nff, &workspace->dPsi_dff, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &workspace->dPsi_dx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nu, &workspace->dPsi_du, &c_ptr);

    assert((char *) raw_memory + sim_gnsf_workspace_calculate_size(config, dims_, opts) >= c_ptr);

    return (void *) workspace;
}

int sim_gnsf(void *config, sim_in *in, sim_out *out, void *args, void *mem_, void *work_)
{
    acados_timer tot_timer, casadi_timer, la_timer;
    acados_tic(&tot_timer);

    // typecast
    sim_gnsf_memory *mem = (sim_gnsf_memory *) mem_;
    sim_rk_opts *opts = (sim_rk_opts *) args;
    sim_gnsf_dims *dims = (sim_gnsf_dims *) in->dims;
    gnsf_model *model = in->model;
    gnsf_workspace *workspace =
        (gnsf_workspace *) sim_gnsf_cast_workspace(config, dims, opts, work_);

    // necessary integers
    int nx = dims->nx;
    int nu = dims->nu;
    int nx1 = dims->nx1;
    int nx2 = dims->nx2;
    int nz = dims->nz;
    int n_out = dims->n_out;
    int ny = dims->ny;
    int nuhat = dims->nuhat;

    int num_stages = opts->ns;
    int num_steps  = opts->num_steps;
    int newton_max = opts->newton_iter;

    int nff = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nx2;
    int nZ  = num_stages * nz;

    // assert - only use supported features
    assert(mem->dt == in->T / opts->num_steps &&
           "model->dt not equal to im_in.T/opts->num_steps, check initialization");

    // assign variables from workspace
    double *Z_work = workspace->Z_work;

    struct blasfeo_dmat J_r_ff =
        workspace->J_r_ff;  // store the the jacobian of the residual w.r.t. ff
    int *ipiv = workspace->ipiv;
    struct blasfeo_dmat J_r_x1u = workspace->J_r_x1u;  // needed for sensitivity propagation

    struct blasfeo_dmat dK1_dx1 = workspace->dK1_dx1;
    struct blasfeo_dmat dK1_du = workspace->dK1_du;
    struct blasfeo_dmat dZ_dx1 = workspace->dZ_dx1;
    struct blasfeo_dmat dZ_du = workspace->dZ_du;
    struct blasfeo_dmat aux_G2_x1 = workspace->aux_G2_x1;
    struct blasfeo_dmat aux_G2_u = workspace->aux_G2_u;
    struct blasfeo_dmat J_G2_K1 = workspace->J_G2_K1;
    struct blasfeo_dmat dK2_dx1 = workspace->dK2_dx1;
    struct blasfeo_dmat dK2_du = workspace->dK2_du;
    struct blasfeo_dmat dK2_dff = workspace->dK2_dff;
    struct blasfeo_dmat dxf_dwn = workspace->dxf_dwn;
    struct blasfeo_dmat S_forw_new = workspace->S_forw_new;  // used to avoid side effects
    struct blasfeo_dmat S_forw = workspace->S_forw;

    struct blasfeo_dmat *f_LO_jac_traj = workspace->f_LO_jac_traj;

    struct blasfeo_dvec *ff_traj = workspace->ff_traj;
    struct blasfeo_dvec *yy_traj = workspace->yy_traj;

    struct blasfeo_dvec f_LO_val = workspace->f_LO_val;
    struct blasfeo_dvec K1_val = workspace->K1_val;
    struct blasfeo_dvec x1_stage_val = workspace->x1_stage_val;
    struct blasfeo_dvec Z_val = workspace->Z_val;

    struct blasfeo_dvec yyu = workspace->yyu;
    struct blasfeo_dvec yyss = workspace->yyss;
    struct blasfeo_dvec y_one_stage = workspace->y_one_stage;
    struct blasfeo_dmat dPHI_dyuhat = workspace->dPHI_dyuhat;

    struct blasfeo_dvec K2_val = workspace->K2_val;
    struct blasfeo_dvec x0_traj = workspace->x0_traj;
    struct blasfeo_dvec res_val = workspace->res_val;
    struct blasfeo_dvec u0 = workspace->u0;
    struct blasfeo_dvec lambda = workspace->lambda;
    struct blasfeo_dvec lambda_old = workspace->lambda_old;

    struct blasfeo_dmat aux_G2_ff = workspace->aux_G2_ff;
    struct blasfeo_dmat dPsi_dff = workspace->dPsi_dff;
    struct blasfeo_dmat dPsi_dx = workspace->dPsi_dx;
    struct blasfeo_dmat dPsi_du = workspace->dPsi_du;

    struct blasfeo_dvec K1u = workspace->K1u;
    struct blasfeo_dvec Zu = workspace->Zu;
    struct blasfeo_dvec ALOtimesx02 = workspace->ALOtimesx02;
    struct blasfeo_dvec uhat = workspace->uhat;

    // memory only available if (opts->sens_algebraic)
    struct blasfeo_dmat dz0_dx1u = workspace->dz0_dx1u;
    struct blasfeo_dvec x0dot_1 = workspace->x0dot_1;
    struct blasfeo_dvec z0 = workspace->z0;

    // memory - precomputed matrices
    double *A_dt = mem->A_dt;
    double *b_dt = mem->b_dt;

    struct blasfeo_dmat KKf = mem->KKf;
    struct blasfeo_dmat KKx = mem->KKx;
    struct blasfeo_dmat KKu = mem->KKu;

    struct blasfeo_dmat YYf = mem->YYf;
    struct blasfeo_dmat YYx = mem->YYx;
    struct blasfeo_dmat YYu = mem->YYu;

    struct blasfeo_dmat ZZf = mem->ZZf;
    struct blasfeo_dmat ZZx = mem->ZZx;
    struct blasfeo_dmat ZZu = mem->ZZu;

    struct blasfeo_dmat ALO = mem->ALO;
    struct blasfeo_dmat M2inv = mem->M2inv;
    struct blasfeo_dmat dK2_dx2 = mem->dK2_dx2;

    struct blasfeo_dmat Lu = mem->Lu;

    // precomputed vectors
    struct blasfeo_dvec KK0 = mem->KK0;
    struct blasfeo_dvec YY0 = mem->YY0;
    struct blasfeo_dvec ZZ0 = mem->ZZ0;

    // ONLY available for algebraic sensitivity propagation
    struct blasfeo_dmat Z0x = mem->Z0x;
    struct blasfeo_dmat Z0u = mem->Z0u;
    struct blasfeo_dmat Z0f = mem->Z0f;

    struct blasfeo_dmat Y0x = mem->Y0x;
    struct blasfeo_dmat Y0u = mem->Y0u;
    struct blasfeo_dmat Y0f = mem->Y0f;

    struct blasfeo_dmat Lx = mem->Lx;
    struct blasfeo_dmat Lxdot = mem->Lxdot;
    struct blasfeo_dmat Lz = mem->Lz;

    int *ipiv_ff0 = workspace->ipiv_ff0;

    struct blasfeo_dmat dr0_dff0 = workspace->dr0_dff0;

    // transform inputs to blasfeo
    blasfeo_pack_dvec(nu, in->u, &u0, 0);
    blasfeo_pack_dvec(nx, &in->x[0], &x0_traj, 0);
    blasfeo_pack_dvec(nx + nu, &in->S_adj[0], &lambda, 0);
    blasfeo_pack_dvec(nx + nu, &in->S_adj[0], &lambda_old, 0);
    blasfeo_pack_dmat(nx, nx + nu, &in->S_forw[0], nx, &S_forw, 0, 0);

    // initialize ff_traj
    for (int ss = 0; ss < num_steps; ss++)
    {
        blasfeo_dvecse(nff, 0, &ff_traj[ss], 0);
    }

    /************************************************
     * Set up function input & outputs
     ************************************************/

    /* PHI - NONLINEARITY FUNCTION */
    ext_fun_arg_t phi_type_in[2];
    void *phi_in[2];

    ext_fun_arg_t phi_fun_type_out[1];
    void *phi_fun_out[1];

    ext_fun_arg_t phi_fun_jac_y_type_out[2];
    void *phi_fun_jac_y_out[2];

    ext_fun_arg_t phi_jac_yuhat_type_out[2];
    void *phi_jac_yuhat_out[2];

    // set up external function argument structs
    struct blasfeo_dvec_args y_in;              // input for y of phi;
    struct blasfeo_dvec_args phi_fun_val_arg;   // output arg for function value of phi
    struct blasfeo_dmat_args phi_jac_y_arg;     // output arg for jacobian of phi w.r.t. y
    struct blasfeo_dmat_args phi_jac_uhat_arg;  // output arg for jacobian of phi w.r.t. uhat

    // set input for phi
    phi_type_in[0] = BLASFEO_DVEC_ARGS;
    phi_type_in[1] = BLASFEO_DVEC;
    phi_in[0] = &y_in;
    phi_in[1] = &uhat;

    // set output for phi_fun
    phi_fun_type_out[0] = BLASFEO_DVEC_ARGS;
    phi_fun_out[0] = &phi_fun_val_arg;
    phi_fun_val_arg.x = &res_val;

    // set output for phi_fun_jac_y
    phi_fun_jac_y_type_out[0] = BLASFEO_DVEC_ARGS;
    phi_fun_jac_y_out[0] = &phi_fun_val_arg;

    phi_fun_jac_y_type_out[1] = BLASFEO_DMAT_ARGS;  // jac_y
    phi_jac_y_arg.A = &dPHI_dyuhat;
    phi_jac_y_arg.aj = 0;  // never changes
    phi_fun_jac_y_out[1] = &phi_jac_y_arg;

    // set output for phi_jac_yuhat // jac_y
    phi_jac_yuhat_type_out[0] = BLASFEO_DMAT_ARGS;
    phi_jac_yuhat_out[0] = &phi_jac_y_arg;
    // jac_uhat
    phi_jac_yuhat_type_out[1] = BLASFEO_DMAT_ARGS;
    phi_jac_uhat_arg.A = &dPHI_dyuhat;
    phi_jac_uhat_arg.aj = ny;  // never changes
    phi_jac_yuhat_out[1] = &phi_jac_uhat_arg;

    /* f_lo - LINEAR OUTPUT FUNCTION */
    ext_fun_arg_t f_lo_fun_type_in[4];
    void *f_lo_fun_in[4];
    ext_fun_arg_t f_lo_fun_type_out[2];
    void *f_lo_fun_out[2];

    f_lo_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    f_lo_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    f_lo_fun_type_in[2] = BLASFEO_DVEC_ARGS;
    f_lo_fun_type_in[3] = BLASFEO_DVEC;

    // f_lo_in[0]: x1;
    struct blasfeo_dvec_args f_lo_in_x1;
    f_lo_fun_in[0] = &f_lo_in_x1;
    // f_lo_in[1]: k1;
    struct blasfeo_dvec_args f_lo_in_k1;
    f_lo_fun_in[1] = &f_lo_in_k1;
    // f_lo_in[2]: z;
    struct blasfeo_dvec_args f_lo_in_z;
    f_lo_fun_in[2] = &f_lo_in_z;
    // f_lo_in[3]: u;
    f_lo_fun_in[3] = &u0;

    // output
    f_lo_fun_type_out[0] = BLASFEO_DVEC_ARGS;
    f_lo_fun_type_out[1] = BLASFEO_DMAT_ARGS;

    struct blasfeo_dvec_args f_lo_val_out;
    f_lo_fun_out[0] = &f_lo_val_out;

    struct blasfeo_dmat_args f_lo_jac_out;
    f_lo_fun_out[1] = &f_lo_jac_out;
    f_lo_jac_out.aj = 0;

    /* TIMINGS */
    out->info->ADtime = 0;
    out->info->LAtime = 0;
    out->info->CPUtime = 0;

    // PRECOMPUTE YY0 + YYu * u, KK0 + KKu * u, ZZ0 + ZZu * u;
    blasfeo_dgemv_n(nyy, nu, 1.0, &YYu, 0, 0, &u0, 0, 1.0, &YY0, 0, &yyu, 0);
    blasfeo_dgemv_n(nK1, nu, 1.0, &KKu, 0, 0, &u0, 0, 1.0, &KK0, 0, &K1u, 0);
    blasfeo_dgemv_n(nZ , nu, 1.0, &ZZu, 0, 0, &u0, 0, 1.0, &ZZ0, 0, &Zu , 0);

    // compute uhat
    blasfeo_dgemv_n(nuhat, nu, 1.0, &Lu, 0, 0, &u0, 0, 0.0, &uhat, 0, &uhat, 0);

    /************************************************
     * FORWARD LOOP
     ************************************************/

    for (int ss = 0; ss < num_steps; ss++)
    {  // STEP LOOP
        blasfeo_dgemv_n(nyy, nx1, 1.0, &YYx, 0, 0, &x0_traj, ss * nx, 1.0, &yyu, 0, &yyss,
                        nyy * ss);  // yyss = YY0 + YYu * u + YYx * x0(at_stage)

        y_in.x = &yy_traj[ss];

        f_lo_in_x1.x = &x1_stage_val;
        f_lo_in_k1.x = &K1_val;
        f_lo_in_z.x = &Z_val;

        f_lo_val_out.x = &f_LO_val;
        f_lo_jac_out.A = &f_LO_jac_traj[ss];

        for (int iter = 0; iter < newton_max; iter++)
        {  // NEWTON-ITERATION
            /* EVALUATE RESIDUAL FUNCTION & JACOBIAN */

            blasfeo_dgemv_n(nyy, nff, 1.0, &YYf, 0, 0, &ff_traj[ss], 0, 1.0, &yyss, nyy * ss,
                            &yy_traj[ss], 0);

            if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
            {
                // set J_r_ff to unit matrix
                blasfeo_dgese(nff, nff, 0.0, &J_r_ff, 0, 0);
                for (int ii = 0; ii < nff; ii++)
                {
                    blasfeo_dgein1(1.0, &J_r_ff, ii, ii);
                }
            }

            for (int ii = 0; ii < num_stages; ii++)
            {  // eval phi, respectively phi_fun_jac_y
                y_in.xi = ii * ny;
                phi_fun_val_arg.xi = ii * n_out;
                phi_jac_y_arg.ai = ii * n_out;
                if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
                {
                    // evaluate
                    acados_tic(&casadi_timer);
                    model->phi_fun_jac_y->evaluate(model->phi_fun_jac_y, phi_type_in, phi_in,
                                                   phi_fun_jac_y_type_out, phi_fun_jac_y_out);
                    out->info->ADtime += acados_toc(&casadi_timer);

                    // build jacobian J_r_ff
                    blasfeo_dgemm_nn(n_out, nff, ny, -1.0, &dPHI_dyuhat, ii * n_out, 0, &YYf,
                                     ii * ny, 0, 1.0, &J_r_ff, ii * n_out, 0, &J_r_ff, ii * n_out,
                                     0);
                }
                else
                {
                    acados_tic(&casadi_timer);
                    model->phi_fun->evaluate(model->phi_fun, phi_type_in, phi_in, phi_fun_type_out,
                                             phi_fun_out);
                    out->info->ADtime += acados_toc(&casadi_timer);
                }
            }

            // set up res_val
            blasfeo_dveccpsc(nff, -1.0, &res_val, 0, &res_val,
                             0);  // set res_val = - res_val; NOW: res_val = -PHI_val;
            blasfeo_dvecad(nff, 1.0, &ff_traj[ss], 0, &res_val, 0);
                    // set res_val = res_val + ff_traj;
                    // this is the actual value of the residual function!

            acados_tic(&la_timer);
            // factorize J_r_ff
            if ((opts->jac_reuse && (ss == 0) & (iter == 0)) || (!opts->jac_reuse))
            {
                blasfeo_dgetrf_rowpivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0, ipiv);
            }

            /* Solve linear system and update ff */
            blasfeo_dvecpe(nff, ipiv, &res_val, 0);  // permute r.h.s.
            blasfeo_dtrsv_lnu(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            blasfeo_dtrsv_unn(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            out->info->LAtime += acados_toc(&la_timer);

            blasfeo_daxpy(nff, -1.0, &res_val, 0, &ff_traj[ss], 0, &ff_traj[ss], 0);

        }  // END NEWTON-ITERATION

        // compute K1 and Z values
        blasfeo_dgemv_n(nK1, nff, 1.0, &KKf, 0, 0, &ff_traj[ss], 0, 1.0, &K1u, 0, &K1_val,
                        0);  // K1u contains KKu * u0 + KK0;
        blasfeo_dgemv_n(nK1, nx1, 1.0, &KKx, 0, 0, &x0_traj, ss * nx, 1.0, &K1_val, 0,
                        &K1_val, 0);
        if (nz)
        {
            blasfeo_dgemv_n(nZ, nff, 1.0, &ZZf, 0, 0, &ff_traj[ss], 0, 1.0, &Zu, 0, &Z_val,
                            0);  // Zu contains ZZu * u0 + ZZ0;
            blasfeo_dgemv_n(nZ, nx1, 1.0, &ZZx, 0, 0, &x0_traj, ss * nx, 1.0, &Z_val, 0,
                            &Z_val, 0);
        }
        // build x1 stage values
        for (int ii = 0; ii < num_stages; ii++)
        {
            blasfeo_daxpy(nx1, 0.0, &x1_stage_val, 0, &x0_traj, ss * nx, &x1_stage_val, nx1 * ii);
            for (int jj = 0; jj < num_stages; jj++)
            {
                blasfeo_daxpy(nx1, A_dt[ii + num_stages * jj], &K1_val, nx1 * jj, &x1_stage_val,
                              nx1 * ii, &x1_stage_val, nx1 * ii);
            }
        }
        if (nx2)
        {
            /* SIMULATE LINEAR OUTPUT SYSTEM */
            blasfeo_dgemv_n(nx2, nx2, 1.0, &ALO, 0, 0, &x0_traj, ss * nx + nx1, 0.0, &f_LO_val,
                            0, &ALOtimesx02, 0);
                    // ALOtimesx02 = ALO * x2part(x0_traj(step));
            for (int ii = 0; ii < num_stages; ii++)
            {  // Evaluate f_LO + jacobian and pack to blasfeo structs
                f_lo_in_x1.xi = ii * nx1;
                f_lo_in_k1.xi = ii * nx1;
                f_lo_in_z.xi = ii * nz;

                f_lo_val_out.xi = ii * nx2;
                f_lo_jac_out.ai = ii * nx2;

                acados_tic(&casadi_timer);
                model->f_lo_fun_jac_x1_x1dot_u_z->evaluate(model->f_lo_fun_jac_x1_x1dot_u_z,
                                                           f_lo_fun_type_in, f_lo_fun_in,
                                                           f_lo_fun_type_out, f_lo_fun_out);
                out->info->ADtime += acados_toc(&casadi_timer);

                blasfeo_dvecsc(nx2, -1.0, &f_LO_val, nx2 * ii);  // f_LO_val = - f_LO_val
                blasfeo_dvecad(nx2, -1.0, &ALOtimesx02, 0, &f_LO_val, nx2 * ii);
                        // f_LO_val = f_LO_val - ALOtimesx02 (actual right hand side)
            }
            // solve for K2
            blasfeo_dgemv_n(nK2, nK2, -1.0, &M2inv, 0, 0, &f_LO_val, 0, 0.0, &K2_val, 0,
                            &K2_val, 0);
        }

        /* Get simulation result */
        blasfeo_daxpy(nx, 0.0, &x0_traj, 0, &x0_traj, nx * ss, &x0_traj, nx * (ss + 1));
        for (int ii = 0; ii < num_stages; ii++)
        {
            blasfeo_daxpy(nx1, b_dt[ii], &K1_val, ii * nx1, &x0_traj, nx * (ss + 1), &x0_traj,
                          nx * (ss + 1));
            blasfeo_daxpy(nx2, b_dt[ii], &K2_val, ii * nx2, &x0_traj, nx1 + nx * (ss + 1), &x0_traj,
                          nx1 + nx * (ss + 1));
        }

        // Forward Sensitivities (via IND)
        if (opts->sens_forw)
        {
            // evaluate jacobian of residual function
            // update yy
            blasfeo_dgemv_n(nyy, nff, 1.0, &YYf, 0, 0, &ff_traj[ss], 0, 1.0, &yyss, nyy * ss,
                            &yy_traj[ss], 0);
            // set J_r_ff to unit matrix
            blasfeo_dgese(nff, nff, 0.0, &J_r_ff, 0, 0);
            for (int ii = 0; ii < nff; ii++)
            {
                blasfeo_dgein1(1.0, &J_r_ff, ii, ii);
            }

            for (int ii = 0; ii < num_stages; ii++)
            {                                      //
                y_in.xi = ii * ny;                 // set input of phi
                phi_jac_uhat_arg.ai = ii * n_out;  // set output
                phi_jac_y_arg.ai = ii * n_out;

                acados_tic(&casadi_timer);
                model->phi_jac_y_uhat->evaluate(model->phi_jac_y_uhat, phi_type_in, phi_in,
                                                phi_jac_yuhat_type_out, phi_jac_yuhat_out);
                out->info->ADtime += acados_toc(&casadi_timer);

                // build J_r_ff
                blasfeo_dgemm_nn(n_out, nff, ny, -1.0, &dPHI_dyuhat, ii * n_out, 0, &YYf, ii * ny,
                                 0, 1.0, &J_r_ff, ii * n_out, 0, &J_r_ff, ii * n_out, 0);
                // build J_r_x1u
                blasfeo_dgemm_nn(n_out, nx1, ny, -1.0, &dPHI_dyuhat, ii * n_out, 0, &YYx, ii * ny,
                                 0, 0.0, &J_r_x1u, ii * n_out, 0, &J_r_x1u, ii * n_out,
                                 0);  // w.r.t. x1
                blasfeo_dgemm_nn(n_out, nu, ny, -1.0, &dPHI_dyuhat, ii * n_out, 0, &YYu, ii * ny, 0,
                                 0.0, &J_r_x1u, ii * n_out, nx1, &J_r_x1u, ii * n_out,
                                 nx1);  // w.r.t. u
                blasfeo_dgemm_nn(n_out, nu, nuhat, -1.0, &dPHI_dyuhat, ii * n_out, ny, &Lu, 0, 0,
                                 1.0, &J_r_x1u, ii * n_out, nx1, &J_r_x1u, ii * n_out,
                                 nx1);  // + dPhi_duhat * L_u;
            }
            acados_tic(&la_timer);
            blasfeo_dgetrf_rowpivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0,
                                    ipiv);        // factorize J_r_ff
            // printf("dPHI_dyuhat = (forward, ss = %d) \n", ss);
            // blasfeo_print_exp_dmat(nff, ny+nuhat, &dPHI_dyuhat, 0, 0);

            blasfeo_drowpe(nff, ipiv, &J_r_x1u);  // permute also rhs
            blasfeo_dtrsm_llnu(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
            blasfeo_dtrsm_lunn(nff, nx1 + nu, 1.0, &J_r_ff, 0, 0, &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);
            out->info->LAtime += acados_toc(&la_timer);

            blasfeo_dgemm_nn(nK1, nx1, nff, -1.0, &KKf, 0, 0, &J_r_x1u, 0, 0, 1.0, &KKx, 0, 0,
                             &dK1_dx1, 0, 0);
            blasfeo_dgemm_nn(nK1, nu, nff, -1.0, &KKf, 0, 0, &J_r_x1u, 0, nx1, 1.0, &KKu, 0, 0,
                             &dK1_du, 0, 0);  // Blasfeo HP & Reference differ here
            blasfeo_dgemm_nn(nZ, nx1, nff, -1.0, &ZZf, 0, 0, &J_r_x1u, 0, 0, 1.0, &ZZx, 0, 0,
                             &dZ_dx1, 0, 0);
            blasfeo_dgemm_nn(nZ, nu, nff, -1.0, &ZZf, 0, 0, &J_r_x1u, 0, nx1, 1.0, &ZZu, 0, 0,
                             &dZ_du, 0, 0);

            if (nx2)
            {
                // BUILD J_G2_wn, J_G2_K1
                for (int ii = 0; ii < num_stages; ii++)
                {
                    for (int jj = 0; jj < num_stages; jj++)
                    {
                        blasfeo_dgecpsc(nx2, nx1, -A_dt[ii + jj * num_stages], &f_LO_jac_traj[ss],
                                        ii * nx2, 0, &J_G2_K1, ii * nx2, jj * nx1);
                            //  copy - T * a_ij * df_dx1(k1_i, x1_i, z_i, u) into (i,j)th block
                    }
                    blasfeo_dgead(nx2, nx1, -1.0, &f_LO_jac_traj[ss], ii * nx2, nx1,
                                 &J_G2_K1, ii * nx2, ii * nx1);
                            //  add - df_dx1dot(k1_i, x1_i, z_i, u) to (i,i)th block
                    blasfeo_dgemm_nn(nx2, nx1, nz, - 1.0, &f_LO_jac_traj[ss], ii * nx2,
                                 2 * nx1 + nu, &dZ_dx1, ii * nz, 0, 0.0,
                                 &aux_G2_x1, ii * nx2, 0, &aux_G2_x1, ii * nx2, 0);
                            // set ith block of aux_G2_x1 to - df_dz(k1_i, x1_i, z_i, u) * dzi_dx1
                    blasfeo_dgemm_nn(nx2, nu, nz, - 1.0, &f_LO_jac_traj[ss], ii * nx2, 2 * nx1 + nu,
                                 &dZ_du, ii * nz, 0, 0.0, &aux_G2_u, ii * nx2, 0,
                                 &aux_G2_u, ii * nx2,  0);
                            // set ith block of aux_G2_u to - df_dz(k1_i, x1_i, z_i, u) * dzi_du
                }

                // BUILD dK2_dwn // dK2_dx1
                blasfeo_dgemm_nn(nK2, nx1, nK1, 1.0, &J_G2_K1, 0, 0, &dK1_dx1, 0, 0, 1.0,
                                 &aux_G2_x1, 0, 0, &aux_G2_x1, 0, 0);

                blasfeo_dgead(nK2, nx1, -1.0, &f_LO_jac_traj[ss], 0, 0, &aux_G2_x1, 0, 0);
                        // add -df_dx1 (at stages concatenated) to aux_G2_x1
                blasfeo_dgemm_nn(nK2, nx1, nK2, -1.0, &M2inv, 0, 0, &aux_G2_x1, 0, 0, 0.0, &dK2_dx1,
                                 0, 0, &dK2_dx1, 0, 0);
                // dK2_du
                blasfeo_dgemm_nn(nK2, nu, nK1, 1.0, &J_G2_K1, 0, 0, &dK1_du, 0, 0, 1.0, &aux_G2_u,
                                 0, 0, &aux_G2_u, 0, 0);
                blasfeo_dgead(nK2, nu, -1.0, &f_LO_jac_traj[ss], 0, 2 * nx1, &aux_G2_u, 0, 0);
                        // add - df_du (at stages concatenated) to aux_G2_u
                blasfeo_dgemm_nn(nK2, nu, nK2, -1.0, &M2inv, 0, 0, &aux_G2_u, 0, 0, 0.0, &dK2_du, 0,
                                 0, &dK2_du, 0, 0);
            }
            // BUILD dxf_dwn
            blasfeo_dgese(nx, nx + nu, 0.0, &dxf_dwn, 0, 0);  // Initialize as unit matrix
            for (int ii = 0; ii < nx; ii++)
            {
                blasfeo_dgein1(1.0, &dxf_dwn, ii, ii);
            }
            for (int ii = 0; ii < num_stages; ii++)
            {
                blasfeo_dgead(nx1, nx1, b_dt[ii], &dK1_dx1, ii * nx1, 0, &dxf_dwn, 0,
                              0);  // derivatives w.r.t. x1
                blasfeo_dgead(nx2, nx1, b_dt[ii], &dK2_dx1, ii * nx2, 0, &dxf_dwn, nx1, 0);

                blasfeo_dgead(nx2, nx2, b_dt[ii], &dK2_dx2, ii * nx2, 0, &dxf_dwn, nx1,
                              nx1);  // derivatives w.r.t. x2

                blasfeo_dgead(nx1, nu, b_dt[ii], &dK1_du, ii * nx1, 0, &dxf_dwn, 0,
                              nx);  // derivatives w.r.t. u
                blasfeo_dgead(nx2, nu, b_dt[ii], &dK2_du, ii * nx2, 0, &dxf_dwn, nx1, nx);
            }
            blasfeo_dgemm_nn(nx, nx, nx, 1.0, &dxf_dwn, 0, 0, &S_forw, 0, 0, 0.0, &S_forw_new, 0, 0,
                             &S_forw_new, 0, 0);
            blasfeo_dgemm_nn(nx, nu, nx, 1.0, &dxf_dwn, 0, 0, &S_forw, 0, nx, 1.0, &dxf_dwn, 0, nx,
                             &S_forw_new, 0, nx);
            blasfeo_dgecp(nx, nx + nu, &S_forw_new, 0, 0, &S_forw, 0, 0);
        }

/* output z and propagate corresponding sensitivities */
        if (ss == 0)
        {
        /* get output value for algebraic states z */
            if (opts->output_z || opts->sens_algebraic)
            {
                for (int ii = 0; ii < nz; ii++)  // ith component of z
                {
                    for (int jj = 0; jj < num_stages; jj++)
                    {
                        Z_work[jj] = blasfeo_dvecex1(&Z_val, nz * jj + ii);
                                // copy values of z_ii in first step, into Z_work
                    }
                    neville_algorithm(0.0, num_stages - 1, mem->c_butcher, Z_work, &out->zn[ii]);
                                // eval polynomial through (c_i, Z_i) at 0.
                }
            }

        /* propagate sensitivities of z */
            if (opts->sens_algebraic)
            {
                for (int ii = 0; ii < nx1; ii++)  // ith component of x0dot_1
                {
                    for (int jj = 0; jj < num_stages; jj++)
                    {
                        Z_work[jj] = blasfeo_dvecex1(&K1_val, nx1 * jj + ii);
                                // copy values of k1_ii in first step, into Z_work
                    }
                    neville_algorithm(0.0, num_stages - 1, mem->c_butcher, Z_work, &out->xn[ii]);
                                // eval polynomial through (c_i, k1_ii) at 0, write in out->xn
                }

                // pack x0dot_1
                blasfeo_pack_dvec(nx1, out->xn, &x0dot_1, 0);
                // pack z0
                blasfeo_pack_dvec(nz, out->zn, &z0, 0);

                // evaluate phi at x0_1, x0_1dot, u, z0;
                // build y_0; // use y_one_stage
                // y_one_stage = Lxdot * x0dot_1
                blasfeo_dgemv_n(ny, nx1, 1.0, &Lxdot, 0, 0, &x0dot_1, 0, 0.0,
                                &y_one_stage, 0, &y_one_stage, 0);
                // y_one_stage += Lx * x0
                blasfeo_dgemv_n(ny, nx1, 1.0, &Lx, 0, 0, &x0_traj, 0, 1.0,
                                &y_one_stage, 0, &y_one_stage, 0);
                // y_one_stage += Lz * z0
                blasfeo_dgemv_n(ny, nz, 1.0, &Lz, 0, 0, &z0, 0, 1.0,
                                &y_one_stage, 0, &y_one_stage, 0);
                // NOTE: alternatively y_one_stage could be computed as
                // yy0 = Y0x * x0_1 + Y0u * u + Y0f * f0,
                //       whereby f0 can be obtained by an interpolation formula.

                // set input for phi
                phi_type_in[0] = BLASFEO_DVEC;
                phi_in[0] = &y_one_stage;

                // set output for phi; store result in first n_out rows of dPHI_dyuhat
                phi_jac_uhat_arg.ai = 0;
                phi_jac_y_arg.ai    = 0;

                acados_tic(&casadi_timer);
                model->phi_jac_y_uhat->evaluate(model->phi_jac_y_uhat, phi_type_in, phi_in,
                                                phi_jac_yuhat_type_out, phi_jac_yuhat_out);
                out->info->ADtime += acados_toc(&casadi_timer);

                /* set up dr0_dxn1u in first rows of J_r_x1u */
                // J_r_x1 = -dphi0_dy * Y0x
                blasfeo_dgemm_nn(n_out, nx1, ny, -1.0, &dPHI_dyuhat, 0, 0, &Y0x, 0, 0, 0.0,
                                    &J_r_x1u, 0, 0, &J_r_x1u, 0, 0);  // w.r.t. x1
                // J_r_u = -dphi0_dy * Y0u
                blasfeo_dgemm_nn(n_out, nu, ny, -1.0, &dPHI_dyuhat, 0, 0, &Y0u, 0, 0, 0.0,
                                    &J_r_x1u, 0, nx1, &J_r_x1u, 0, nx1);

                // J_r_u = J_r_u - dphi0_duhat * Lu
                blasfeo_dgemm_nn(n_out, nu, nuhat, -1.0, &dPHI_dyuhat, 0, ny, &Lu, 0, 0, 1.0,
                                    &J_r_x1u, 0, nx1, &J_r_x1u, 0, nx1);  // w.r.t. u

                /* set up dr0_dff0 = eye(n_out) - dphi_dy * Y0f; */
                blasfeo_dgese(n_out, n_out, 0.0, &dr0_dff0, 0, 0);  // set all to zero
                blasfeo_ddiare(n_out, 1.0, &dr0_dff0, 0, 0);  // dr0_dff0 is unit now
                // dr0_dff0 = dr0_dff0 - dphi_dy * Y0f;
                blasfeo_dgemm_nn(n_out, n_out, ny, -1.0, &dPHI_dyuhat, 0, 0, &Y0f, 0, 0, 1.0,
                        &dr0_dff0, 0, 0, &dr0_dff0, 0, 0);

                /* solve dff0_dxn1u = - dr0_dff0 \ dr0_dxn1u */
                acados_tic(&la_timer);
                blasfeo_dgetrf_rowpivot(n_out, n_out, &dr0_dff0, 0, 0, &dr0_dff0, 0, 0,
                                        ipiv_ff0);        // factorize dr0_dff0

                blasfeo_drowpe(n_out, ipiv_ff0, &J_r_x1u);  // permute also rhs
                blasfeo_dtrsm_llnu(n_out, nx1 + nu, 1.0, &dr0_dff0, 0, 0, &J_r_x1u,
                                     0, 0, &J_r_x1u, 0, 0);
                blasfeo_dtrsm_lunn(n_out, nx1 + nu, 1.0, &dr0_dff0, 0, 0, &J_r_x1u,
                                     0, 0, &J_r_x1u, 0, 0);
                out->info->LAtime += acados_toc(&la_timer);
                // J_r_x1u now contains dff0_dxn1u in first rows

                // copy Z0x, Z0u into dz0_dx1u
                blasfeo_dgecp(nz, nx1, &Z0x, 0, 0, &dz0_dx1u, 0, 0);
                blasfeo_dgecp(nz, nu , &Z0u, 0, 0, &dz0_dx1u, 0, nx1);

                // add Z0f * dff0_dxn1u
                blasfeo_dgemm_nn(nz, nx1 + nu, n_out, -1.0, &Z0f, 0, 0, &J_r_x1u, 0, 0, 1.0,
                                &dz0_dx1u, 0, 0, &dz0_dx1u, 0, 0);

                // extract into out->S_algebraic
                blasfeo_unpack_dmat(nz, nx1, &dz0_dx1u, 0, 0, out->S_algebraic, nz);
                blasfeo_unpack_dmat(nz, nu , &dz0_dx1u, 0, nx1, &out->S_algebraic[nx*nz], nz);

                for (int ii = 0; ii < nx2 * nz; ii++) {  //  dz_dx2_0 = 0
                    out->S_algebraic[nx1 * nz + ii] = 0;
                }

                // reset input for phi
                phi_type_in[0] = BLASFEO_DVEC_ARGS;
                phi_in[0] = &y_in;
            }
        }
    }  // end step loop: ss



/************************************************
 * ADJOINT SENSITIVITY PROPAGATION
 ************************************************/

    if (opts->sens_adj)
    {
        for (int ss = num_steps - 1; ss >= 0; ss--)
        {
            y_in.x = &yy_traj[ss];
            for (int ii = 0; ii < num_stages; ii++)
            {
                blasfeo_dgemm_nn(nx2, nff, nz, -1.0, &f_LO_jac_traj[ss], nx2 * ii, 2 * nx1 + nu,
                        &ZZf, ii * nz, 0, 0.0, &aux_G2_ff, ii * nx2, 0, &aux_G2_ff, ii * nx2, 0);
                        // set aux_G2_ff(i_th block) = - df_dzi * ZZf(ith block);
                blasfeo_dgemm_nn(nx2, nx1, nz, -1.0, &f_LO_jac_traj[ss], nx2 * ii, 2 * nx1 + nu,
                        &ZZx, ii * nz, 0, 0.0, &aux_G2_x1, ii * nx2, 0, &aux_G2_x1, ii * nx2, 0);
                        // set aux_G2_x1(i_th block) = - df_dzi * ZZx(ith block);
                blasfeo_dgemm_nn(nx2, nu, nz, -1.0, &f_LO_jac_traj[ss], nx2 * ii, 2 * nx1 + nu,
                        &ZZu, ii * nz, 0, 0.0, &aux_G2_u, ii * nx2, 0, &aux_G2_u, ii * nx2, 0);
                        // set aux_G2_u(i_th block) = - df_dzi * ZZu(ith block);
                // set up J_G2_K1 (as in forward loop)
                for (int jj = 0; jj < num_stages; jj++)
                {
                    blasfeo_dgecpsc(nx2, nx1, -A_dt[ii + jj * num_stages], &f_LO_jac_traj[ss],
                                 ii * nx2, 0, &J_G2_K1, ii * nx2, jj * nx1);
                        //  copy - T * a_ij * df_dx1(k1_i, x1_i, z_i, u) into (i,j)th block
                }
                blasfeo_dgead(nx2, nx1, -1.0, &f_LO_jac_traj[ss], nx2 * ii, nx1, &J_G2_K1, nx2 * ii,
                              nx1 * ii);
                        //  add -df_dx1dot(k1_i, x1_i, z_i, u) to (i,i)th block
            }
            blasfeo_dgemm_nn(nK2, nff, nK1, 1.0, &J_G2_K1, 0, 0, &KKf, 0, 0, 1.0, &aux_G2_ff, 0, 0,
                             &aux_G2_ff, 0, 0);
            blasfeo_dgemm_nn(nK2, nx1, nK1, 1.0, &J_G2_K1, 0, 0, &KKx, 0, 0, 1.0, &aux_G2_x1, 0, 0,
                             &aux_G2_x1, 0, 0);
            blasfeo_dgemm_nn(nK2, nu, nK1, 1.0, &J_G2_K1, 0, 0, &KKu, 0, 0, 1.0, &aux_G2_u, 0, 0,
                             &aux_G2_u, 0, 0);

            blasfeo_dgead(nK2, nx1, -1.0, &f_LO_jac_traj[ss], 0, 0, &aux_G2_x1, 0, 0);
            blasfeo_dgead(nK2, nu, -1.0, &f_LO_jac_traj[ss], 0, 2 * nx1, &aux_G2_u, 0, 0);

            blasfeo_dgemm_nn(nK2, nff, nK2, -1.0, &M2inv, 0, 0, &aux_G2_ff, 0, 0, 0.0, &dK2_dff, 0,
                             0, &dK2_dff, 0, 0);
            blasfeo_dgemm_nn(nK2, nx1, nK2, -1.0, &M2inv, 0, 0, &aux_G2_x1, 0, 0, 0.0, &dK2_dx1, 0,
                             0, &dK2_dx1, 0, 0);
            blasfeo_dgemm_nn(nK2, nu, nK2, -1.0, &M2inv, 0, 0, &aux_G2_u, 0, 0, 0.0, &dK2_du, 0, 0,
                             &dK2_du, 0, 0);

            blasfeo_dgese(nx, nff, 0.0, &dPsi_dff, 0, 0);  // initialize dPsi_d..
            blasfeo_dgese(nx, nx, 0.0, &dPsi_dx, 0, 0);
            blasfeo_dgese(nx, nu, 0.0, &dPsi_du, 0, 0);
            blasfeo_ddiare(nx, 1.0, &dPsi_dx, 0, 0);  // dPsi_dx is unit now

            // compute dPsi_d..
            for (int ii = 0; ii < num_stages; ii++)
            {
                blasfeo_dgead(nx1, nff, b_dt[ii], &KKf, ii * nx1, 0, &dPsi_dff, 0, 0);
                blasfeo_dgead(nx1, nx1, b_dt[ii], &KKx, ii * nx1, 0, &dPsi_dx, 0, 0);
                blasfeo_dgead(nx1, nu, b_dt[ii], &KKu, ii * nx1, 0, &dPsi_du, 0, 0);

                blasfeo_dgead(nx2, nff, b_dt[ii], &dK2_dff, ii * nx2, 0, &dPsi_dff, nx1, 0);
                blasfeo_dgead(nx2, nx1, b_dt[ii], &dK2_dx1, ii * nx2, 0, &dPsi_dx, nx1, 0);
                blasfeo_dgead(nx2, nx2, b_dt[ii], &dK2_dx2, ii * nx2, 0, &dPsi_dx, nx1, nx1);
                blasfeo_dgead(nx2, nu, b_dt[ii], &dK2_du, ii * nx2, 0, &dPsi_du, nx1, 0);
            }
            // evaluate jacobian of residual function
            // update yy
            blasfeo_dgemv_n(nyy, nff, 1.0, &YYf, 0, 0, &ff_traj[ss], 0, 1.0, &yyss, nyy * ss,
                            &yy_traj[ss], 0);
            // set J_r_ff to unit matrix
            blasfeo_dgese(nff, nff, 0.0, &J_r_ff, 0, 0);
            for (int ii = 0; ii < nff; ii++)
            {
                blasfeo_dgein1(1.0, &J_r_ff, ii, ii);
            }
            for (int ii = 0; ii < num_stages; ii++)
            {
                y_in.xi = ii * ny;  // set input of phi
                phi_jac_uhat_arg.ai = ii * n_out;
                phi_jac_y_arg.ai = ii * n_out;

                acados_tic(&casadi_timer);
                model->phi_jac_y_uhat->evaluate(model->phi_jac_y_uhat, phi_type_in, phi_in,
                                                phi_jac_yuhat_type_out, phi_jac_yuhat_out);
                out->info->ADtime += acados_toc(&casadi_timer);

                // build J_r_ff
                blasfeo_dgemm_nn(n_out, nff, ny, -1.0, &dPHI_dyuhat, ii * n_out, 0, &YYf, ii * ny,
                                 0, 1.0, &J_r_ff, ii * n_out, 0, &J_r_ff, ii * n_out, 0);

                // build J_r_x1u
                blasfeo_dgemm_nn(n_out, nx1, ny, -1.0, &dPHI_dyuhat, ii * n_out, 0, &YYx, ii * ny,
                                 0, 0.0, &J_r_x1u, ii * n_out, 0, &J_r_x1u, ii * n_out,
                                 0);  // w.r.t. x1
                blasfeo_dgemm_nn(n_out, nu, ny, -1.0, &dPHI_dyuhat, ii * n_out, 0, &YYu, ii * ny, 0,
                                 0.0, &J_r_x1u, ii * n_out, nx1, &J_r_x1u, ii * n_out,
                                 nx1);  // w.r.t. u
                blasfeo_dgemm_nn(n_out, nu, nuhat, -1.0, &dPHI_dyuhat, ii * n_out, ny, &Lu, 0, 0,
                                 1.0, &J_r_x1u, ii * n_out, nx1, &J_r_x1u, ii * n_out,
                                 nx1);  // + dPhi_duhat * L_u;
            }
            acados_tic(&la_timer);
            blasfeo_dgetrf_rowpivot(nff, nff, &J_r_ff, 0, 0, &J_r_ff, 0, 0,
                                    ipiv);  // factorize J_r_ff
            out->info->LAtime += acados_toc(&la_timer);

            blasfeo_dgemv_t(nx, nff, 1.0, &dPsi_dff, 0, 0, &lambda, 0, 0.0, &res_val, 0, &res_val,
                            0);  // use res_val to store lambda_ff
            acados_tic(&la_timer);
            blasfeo_dtrsv_utn(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            blasfeo_dtrsv_ltu(nff, &J_r_ff, 0, 0, &res_val, 0, &res_val, 0);
            blasfeo_dvecpei(nff, ipiv, &res_val, 0);  // permute linear syst solution
            out->info->LAtime += acados_toc(&la_timer);

            blasfeo_dveccp(nx + nu, &lambda, 0, &lambda_old, 0);
            blasfeo_dgemv_t(nx, nu, 1.0, &dPsi_du, 0, 0, &lambda_old, 0, 1.0, &lambda_old, nx,
                            &lambda, nx);  // update lambda_u

            blasfeo_dgemv_t(nx, nx, 1.0, &dPsi_dx, 0, 0, &lambda_old, 0, 0.0, &res_val, 0,
                         &lambda, 0);
            blasfeo_dveccp(nx + nu, &lambda, 0, &lambda_old, 0);
            blasfeo_dgemv_t(nff, nx1, -1.0, &J_r_x1u, 0, 0, &res_val, 0, 1.0, &lambda_old, 0,
                            &lambda, 0);
            blasfeo_dgemv_t(nff, nu, -1.0, &J_r_x1u, 0, nx1, &res_val, 0, 1.0, &lambda_old, nx,
                            &lambda, nx);
        }
    }
/* unpack */
    blasfeo_unpack_dvec(nx, &x0_traj, nx * num_steps, out->xn);
    blasfeo_unpack_dmat(nx, nx + nu, &S_forw, 0, 0, out->S_forw, nx);
    blasfeo_unpack_dvec(nx + nu, &lambda, 0, out->S_adj);

    out->info->CPUtime = acados_toc(&tot_timer);
    return 0;
}

void sim_gnsf_config_initialize_default(void *config_)
{
    sim_solver_config *config = config_;
    config->evaluate = &sim_gnsf;
    // opts
    config->opts_calculate_size = &sim_gnsf_opts_calculate_size;
    config->opts_assign = &sim_gnsf_opts_assign;
    config->opts_initialize_default = &sim_gnsf_opts_initialize_default;
    config->opts_update = &sim_gnsf_opts_update;
    // memory & workspace
    config->memory_calculate_size = &sim_gnsf_memory_calculate_size;
    config->memory_assign = &sim_gnsf_memory_assign;
    config->workspace_calculate_size = &sim_gnsf_workspace_calculate_size;
    // model
    config->model_calculate_size = &sim_gnsf_model_calculate_size;
    config->model_assign = &sim_gnsf_model_assign;
    config->model_set_function = &sim_gnsf_model_set_function;
    // dims
    config->dims_calculate_size = &sim_gnsf_dims_calculate_size;
    config->dims_assign = &sim_gnsf_dims_assign;
    config->set_nx = &sim_gnsf_set_nx;
    config->set_nu = &sim_gnsf_set_nu;
    config->set_nz = &sim_gnsf_set_nz;
    config->get_nx = &sim_gnsf_get_nx;
    config->get_nu = &sim_gnsf_get_nu;
    config->get_nz = &sim_gnsf_get_nz;
    return;
}
