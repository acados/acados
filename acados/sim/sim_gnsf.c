/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
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

acados_size_t sim_gnsf_dims_calculate_size()
{
    acados_size_t size = sizeof(sim_gnsf_dims);
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
    dims->nz1 = 0;
    dims->nx1 = 0;
    dims->n_out = 0;
    dims->ny = 0;
    dims->nuhat = 0;

    assert((char *) raw_memory + sim_gnsf_dims_calculate_size() == c_ptr);
    return dims;
}



/************************************************
 * get & set functions
 ************************************************/

void sim_gnsf_dims_set(void *config_, void *dims_, const char *field, const int *value)
{
    sim_gnsf_dims *dims = dims_;

    if (!strcmp(field, "nx"))
    {
        dims->nx = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu = *value;
    }
    else if (!strcmp(field, "nz"))
    {
        dims->nz = *value;
    }
    else if (!strcmp(field, "nx1") || !strcmp(field, "gnsf_nx1"))
    {
        dims->nx1 = *value;
    }
    else if (!strcmp(field, "nz1") || !strcmp(field, "gnsf_nz1"))
    {
        dims->nz1 = *value;
    }
    else if (!strcmp(field, "nout") || !strcmp(field, "gnsf_nout"))
    {
        dims->n_out = *value;
    }
    else if (!strcmp(field, "ny") || !strcmp(field, "gnsf_ny"))
    {
        dims->ny = *value;
    }
    else if (!strcmp(field, "nuhat") || !strcmp(field, "gnsf_nuhat"))
    {
        dims->nuhat = *value;
    }
    else
    {
        printf("\nerror: sim_gnsf_dims_set: field not available: %s\n", field);
        exit(1);
    }
}



void sim_gnsf_dims_get(void *config_, void *dims_, const char *field, int *value)
{
    sim_gnsf_dims *dims = dims_;

    if (!strcmp(field, "nx"))
    {
        *value = dims->nx;
    }
    else if (!strcmp(field, "nu"))
    {
        *value = dims->nu;
    }
    else if (!strcmp(field, "nz"))
    {
        *value = dims->nz;
    }
    else if (!strcmp(field, "nout") || !strcmp(field, "gnsf_nout"))
    {
        *value = dims->n_out;
    }
    else
    {
        printf("\nerror: sim_gnsf_dims_get: field not available: %s\n", field);
        exit(1);
    }
}



/************************************************
 * import function
 ************************************************/

static void sim_gnsf_import_matrices(void *dims_, gnsf_model *model)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    external_function_generic *get_matrices_fun = model->get_gnsf_matrices;

    int nu = dims->nu;
    int nz = dims->nz;
    int nx = dims->nx;
    int nx2 = nx - dims->nx1;
    int nz2 = nz - dims->nz1;

    // ensure compatibility with earlier version without B_LO, ipiv_x, ipiv_z, c_LO
    double tmp_nontriv_f_LO = 0.0;
    double tmp_fully_linear = 0.0;
    for (int ii = 0; ii < (nx2+nz2)*nu; ii++)
        model->B_LO[ii] = 0.0;
    for (int ii = 0; ii < nx2+nz2; ii++)
        model->c_LO[ii] = 0.0;

    // ipiv_x, z
    for (int ii = 0; ii < nx; ii++)
        model->ipiv_x_double[ii] = (double) ii;
    for (int ii = 0; ii < nz; ii++)
        model->ipiv_z_double[ii] = (double) ii;

    // calling the external function
    ext_fun_arg_t ext_fun_type_in[1];
    void *ext_fun_in[1];
    ext_fun_arg_t ext_fun_type_out[17];
    void *ext_fun_out[17];

    ext_fun_type_in[0] = COLMAJ;
    ext_fun_in[0] = model->A;  // just to have some input;

    for (int ii = 0; ii < 17; ii++)
        ext_fun_type_out[ii] = COLMAJ;

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
    ext_fun_out[10] = model->E_LO;
    ext_fun_out[11] = model->B_LO;
    ext_fun_out[12] = &tmp_nontriv_f_LO;
    ext_fun_out[13] = &tmp_fully_linear;
    ext_fun_out[14] = model->ipiv_x_double;
    ext_fun_out[15] = model->ipiv_z_double;
    ext_fun_out[16] = model->c_LO;

    get_matrices_fun->evaluate(get_matrices_fun, ext_fun_type_in, ext_fun_in, ext_fun_type_out, ext_fun_out);

    model->nontrivial_f_LO = (tmp_nontriv_f_LO > 0);
    model->fully_linear = (tmp_fully_linear > 0);
    for (int ii = 0; ii < nx; ii++)
        model->ipiv_x[ii] = (int) model->ipiv_x_double[ii];
    for (int ii = 0; ii < nz; ii++)
        model->ipiv_z[ii] = (int) model->ipiv_z_double[ii];

    // printf("\nimported model matrices\n");
    // printf("\nc_LO =\n");
    // for (int ii = 0; ii < nx2 + nz2; ii++) {
    //     printf("%e \t", model->c_LO[ii]);
    // }

}



/************************************************
 * opts
 ************************************************/

acados_size_t sim_gnsf_opts_calculate_size(void *config_, void *dims)
{
    int ns_max = NS_MAX;

    acados_size_t size = 0;

    size += sizeof(sim_opts);

    size += ns_max * ns_max * sizeof(double);  // A_mat
    size += ns_max * sizeof(double);           // b_vec
    size += ns_max * sizeof(double);           // c_vec

    size += butcher_tableau_work_calculate_size(ns_max);

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}


void *sim_gnsf_opts_assign(void *config_, void *dims, void *raw_memory)
{
    int ns_max = NS_MAX;

    char *c_ptr = (char *) raw_memory;

    sim_opts *opts = (sim_opts *) c_ptr;
    c_ptr += sizeof(sim_opts);

    align_char_to(8, &c_ptr);

    assign_and_advance_double(ns_max * ns_max, &opts->A_mat, &c_ptr);
    assign_and_advance_double(ns_max, &opts->b_vec, &c_ptr);
    assign_and_advance_double(ns_max, &opts->c_vec, &c_ptr);

    // work
    opts->work = c_ptr;
    c_ptr += butcher_tableau_work_calculate_size(ns_max);;

    assert((char *) raw_memory + sim_gnsf_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}



void sim_gnsf_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_opts *opts = opts_;

    // default options
    opts->newton_iter = 3;
    // opts->scheme = NULL;
    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
    opts->jac_reuse = true;
    opts->exact_z_output = false;
    opts->ns = 3;
    opts->collocation_type = GAUSS_LEGENDRE;

    assert(opts->ns <= NS_MAX && "ns > NS_MAX!");

    // butcher tableau
    calculate_butcher_tableau(opts->ns, opts->collocation_type, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);
    // for consistency check
    opts->tableau_size = opts->ns;

    // TODO(oj): check if constr h or cost depend on z, turn on in this case only.
    if (dims->nz > 0)
    {
        opts->output_z = true;
        opts->sens_algebraic = true;
    }
    else
    {
        opts->output_z = false;
        opts->sens_algebraic = false;
    }

    return;
}



void sim_gnsf_opts_update(void *config_, void *dims, void *opts_)
{
    sim_opts *opts = opts_;

    assert(opts->ns <= NS_MAX && "ns > NS_MAX!");

    calculate_butcher_tableau(opts->ns, opts->collocation_type, opts->c_vec, opts->b_vec, opts->A_mat, opts->work);

    opts->tableau_size = opts->ns;

    return;
}



void sim_gnsf_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_set_(opts, field, value);
}



void sim_gnsf_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    sim_opts *opts = (sim_opts *) opts_;
    sim_opts_get_(config_, opts, field, value);
}



/************************************************
 * model
 ************************************************/

acados_size_t sim_gnsf_model_calculate_size(void *config, void *dims_)
{
    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;

    // necessary integers
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int nz      = dims->nz;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    acados_size_t size = 0;
    size += sizeof(gnsf_model);

    // model defining matrices
    size += (nx1 + nz1) * (nx1 + nu + n_out + (nx1 + nz1)) * sizeof(double);  // A, B, C, E

    size += ny * (2 * nx1 + nz1) * sizeof(double);  // L_x, L_xdot, L_z
    size += nuhat * nu * sizeof(double);           // L_u
    size += (nx2 + nz2) * nx2 * sizeof(double);            // A_LO
    size += (nx1 + nz1) * sizeof(double);           // c
    size += (nx2 + nz2) * sizeof(double);           // c_LO
    size += (nx2 + nz2) * (nx2 + nz2) * sizeof(double);  // E_LO
    size += (nx2 + nz2) * nu * sizeof(double); // B_LO

    // ipiv
    size += (nx + nz) * sizeof(int); // ipiv_x, ipiv_z
    size += (nx + nz) * sizeof(double); // ipiv_x_double, ipiv_z_double
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
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int nz      = dims->nz;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    gnsf_model *model = (gnsf_model *) c_ptr;
    c_ptr += sizeof(gnsf_model);

    // set default
    model->auto_import_gnsf = true;

    // assign model matrices
    assign_and_advance_double((nx1 + nz1) * nx1, &model->A, &c_ptr);
    assign_and_advance_double((nx1 + nz1) * nu, &model->B, &c_ptr);
    assign_and_advance_double((nx1 + nz1) * n_out, &model->C, &c_ptr);
    assign_and_advance_double((nx1 + nz1) * (nx1 + nz1), &model->E, &c_ptr);
    assign_and_advance_double(nx1 + nz1 , &model->c, &c_ptr);
    assign_and_advance_double(nx2 + nz2 , &model->c_LO, &c_ptr);

    assign_and_advance_double(ny * nx1, &model->L_x, &c_ptr);
    assign_and_advance_double(ny * nx1, &model->L_xdot, &c_ptr);
    assign_and_advance_double(ny * nz1,  &model->L_z, &c_ptr);

    assign_and_advance_double(nuhat * nu, &model->L_u, &c_ptr);

    assign_and_advance_double((nx2 + nz2) * nx2, &model->A_LO, &c_ptr);
    assign_and_advance_double((nx2 + nz2) * nu, &model->B_LO, &c_ptr);
    assign_and_advance_double((nx2 + nz2) * (nx2 + nz2), &model->E_LO, &c_ptr);

    assign_and_advance_double(nx, &model->ipiv_x_double, &c_ptr);
    assign_and_advance_double(nz, &model->ipiv_z_double, &c_ptr);

    assign_and_advance_int(nx, &model->ipiv_x, &c_ptr);
    assign_and_advance_int(nz, &model->ipiv_z, &c_ptr);

    // initialize with identity permutation
    for (int ii = 0; ii < nx; ii++)
        model->ipiv_x[ii] = ii;
    for (int ii = 0; ii < nz; ii++)
        model->ipiv_z[ii] = ii;

    assert((char *) raw_memory + sim_gnsf_model_calculate_size(config, dims_) >= c_ptr);
    return model;
}



int sim_gnsf_model_set(void *model_, const char *field, void *value)
{
    gnsf_model *model = model_;

    if (!strcmp(field, "phi_fun") || !strcmp(field, "gnsf_phi_fun"))
    {
        model->phi_fun = value;
    }
    else if (!strcmp(field, "phi_fun_jac_y") || !strcmp(field, "gnsf_phi_fun_jac_y"))
    {
        model->phi_fun_jac_y = value;
    }
    else if (!strcmp(field, "phi_jac_y_uhat") || !strcmp(field, "gnsf_phi_jac_y_uhat"))
    {
        model->phi_jac_y_uhat = value;
    }
    else if (!strcmp(field, "f_lo_jac_x1_x1dot_u_z") || !strcmp(field, "gnsf_f_lo_fun_jac_x1k1uz"))
    {
        model->f_lo_fun_jac_x1_x1dot_u_z = value;
    }
    else if (!strcmp(field, "get_gnsf_matrices") || !strcmp(field, "gnsf_get_matrices_fun"))
    {
        model->get_gnsf_matrices = value;
    }
    else
    {
        printf("\nerror: sim_gnsf_model_set: wrong field: %s\n", field);
        exit(1);
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
    sim_opts *opts = (sim_opts *) opts_;

    // int nx      = dims->nx;
    int nu      = dims->nu;
    // int nz      = dims->nz;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    // int nuhat   = dims->nuhat;
    // int nx2     = nx - nx1;
    // int nz2     = nz - nz1;

    int num_stages = opts->ns;

    int nvv = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    // int nK2 = num_stages * (nx2 + nz2);
    int nZ1 = num_stages * nz1;

    char *c_ptr = (char *) raw_memory;
    align_char_to(8, &c_ptr);

    // struct
    gnsf_pre_workspace *work = (gnsf_pre_workspace *) c_ptr;
    c_ptr += sizeof(gnsf_pre_workspace);

    assign_and_advance_int(nK1, &work->ipivEE1, &c_ptr);
    assign_and_advance_int(nZ1, &work->ipivEE2, &c_ptr);
    assign_and_advance_int(nZ1, &work->ipivQQ1, &c_ptr);

    align_char_to(64, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx1, nx1, &work->E11, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx1, nz1, &work->E12, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz1, nx1, &work->E21, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz1, nz1, &work->E22, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx1, nx1, &work->A1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz1, nx1, &work->A2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx1, nu, &work->B1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz1, nu, &work->B2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx1, n_out, &work->C1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz1, n_out, &work->C2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK1, nx1, &work->AA1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nx1, &work->AA2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nu, &work->BB1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nu, &work->BB2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK1, nvv, &work->CC1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nvv, &work->CC2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nZ1, &work->DD1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nK1, &work->DD2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK1, nK1, &work->EE1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nZ1, &work->EE2, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nZ1, nZ1, &work->QQ1, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nyy, nZ1, &work->LLZ, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nx1, &work->LLx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nK1, &work->LLK, &c_ptr);

    // if (opts->sens_algebraic){
    //     // for algebraic sensitivities
    //     assign_and_advance_blasfeo_dmat_mem(nz1,  nz1,  &work->Q1, &c_ptr);
    // }

    assign_and_advance_blasfeo_dvec_mem(nK1, &work->cc1, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nZ1, &work->cc2, &c_ptr);

    assert((char *) raw_memory + sim_gnsf_workspace_calculate_size(config_, dims, opts) >= c_ptr);
    return (void *) work;
}  // cast pre_workspace



int sim_gnsf_precompute(void *config_, sim_in *in, sim_out *out, void *opts_, void *mem_,
                       void *work_)
{
    acados_timer atimer;
    acados_tic(&atimer);

    int status = ACADOS_SUCCESS;

    sim_gnsf_dims *dims = (sim_gnsf_dims *) in->dims;
    sim_opts *opts = opts_;
    gnsf_model *model = in->model;


    if (model->get_gnsf_matrices == NULL && model->auto_import_gnsf)
    {
        printf("sim_gnsf error: get_gnsf_matrices function seems to be unset!\n");
        exit(1);
    }

    if (model->auto_import_gnsf)
        sim_gnsf_import_matrices(dims, model);

    // dimension ints
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nz      = dims->nz;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    int num_stages = opts->ns;
    int num_steps = opts->num_steps;

    int nxz2 = nx2 + nz2;

    int nvv = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nxz2;
    int nZ1 = num_stages * nz1;

    // printf("imported A\n");
    // d_print_exp_mat(nx1, nx1, model->A, nx1);

    // set up precomputation workspace
    gnsf_pre_workspace *work =
        (gnsf_pre_workspace *) gnsf_cast_pre_workspace(config_, dims, opts, work_);
    // set up memory
    sim_gnsf_memory *mem = (sim_gnsf_memory *) mem_;

    double dt = in->T / num_steps;
    if (dt == 0.0)
    {
        printf("sim_gnsf error: simulation time = 0; seems to be unset!\n");
        exit(1);
    }
    mem->dt = dt;

    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
    double *c_vec = opts->c_vec;

    double *c_butcher = mem->c_butcher;
    double *b_dt = mem->b_dt;
    double *A_dt = mem->A_dt;

    double *A_LO = model->A_LO;
    double *B_LO = model->B_LO;

    // workspace stuff
    int *ipivEE1 = work->ipivEE1;
    int *ipivEE2 = work->ipivEE2;
    int *ipivQQ1 = work->ipivQQ1;

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

    // memory - precomputed matrices
    int *ipivM2 = mem->ipivM2;

    struct blasfeo_dmat *KKv = &mem->KKv;
    struct blasfeo_dmat *KKx = &mem->KKx;
    struct blasfeo_dmat *KKu = &mem->KKu;

    struct blasfeo_dmat *YYv = &mem->YYv;
    struct blasfeo_dmat *YYx = &mem->YYx;
    struct blasfeo_dmat *YYu = &mem->YYu;

    struct blasfeo_dmat *ZZv = &mem->ZZv;
    struct blasfeo_dmat *ZZx = &mem->ZZx;
    struct blasfeo_dmat *ZZu = &mem->ZZu;

    struct blasfeo_dmat *ALO = &mem->ALO;
    struct blasfeo_dmat *BLO = &mem->BLO;
    struct blasfeo_dmat *M2_LU = &mem->M2_LU;
    struct blasfeo_dmat *dK2_dx2 = &mem->dK2_dx2;
    struct blasfeo_dmat *dK2_du = &mem->dK2_du;
    struct blasfeo_dmat *dx2f_dx2u = &mem->dx2f_dx2u;

    struct blasfeo_dmat *Lu = &mem->Lu;

    // ONLY for algebraic sensitivity propagation
    // struct blasfeo_dmat *Lx = mem->Lx;
    // struct blasfeo_dmat *Lxdot = mem->Lxdot;
    // struct blasfeo_dmat *Lz = mem->Lz;

    // struct blasfeo_dmat Q1  = work->Q1;

    // struct blasfeo_dmat *K0x = mem->K0x;
    // struct blasfeo_dmat *K0u = mem->K0u;
    // struct blasfeo_dmat *K0v = mem->K0v;

    // struct blasfeo_dmat *Z0x = mem->Z0x;
    // struct blasfeo_dmat *Z0u = mem->Z0u;
    // struct blasfeo_dmat *Z0v = mem->Z0v;

    // struct blasfeo_dmat *Y0x = mem->Y0x;
    // struct blasfeo_dmat *Y0u = mem->Y0u;
    // struct blasfeo_dmat *Y0v = mem->Y0v;

    // struct blasfeo_dmat *ELO_LU = mem->ELO_LU;
    // struct blasfeo_dmat *ELO_inv_ALO = mem->ELO_inv_ALO;
    // int *ipiv_ELO = mem->ipiv_ELO;

    // precomputed vectors
    struct blasfeo_dvec *KK0 = &mem->KK0;
    struct blasfeo_dvec *YY0 = &mem->YY0;
    struct blasfeo_dvec *ZZ0 = &mem->ZZ0;

    // set memory to zeros
    blasfeo_dgese(nK1, nvv, 0.0, &CC1, 0, 0);
    blasfeo_dgese(nZ1, nvv, 0.0, &CC2, 0, 0);
    blasfeo_dgese(nK1, nZ1, 0.0, &DD1, 0, 0);
    blasfeo_dgese(nZ1, nK1, 0.0, &DD2, 0, 0);
    blasfeo_dgese(nK1, nK1, 0.0, &EE1, 0, 0);
    blasfeo_dgese(nZ1, nZ1, 0.0, &EE2, 0, 0);
    blasfeo_dgese(nyy, nK1, 0.0, &LLK, 0, 0);
    blasfeo_dgese(nK2, nK2, 0.0, M2_LU, 0, 0);

    blasfeo_pack_dmat(nx1, nx1, model->E, nx1 + nz1, &E11, 0, 0);

    blasfeo_pack_dmat(nx1, nz1, &model->E[(nx1 + nz1) * nx1], nx1 + nz1, &E12, 0, 0);
    blasfeo_pack_dmat(nz1, nx1, &model->E[nx1], nx1 + nz1, &E21, 0, 0);
    blasfeo_pack_dmat(nz1, nz1, &model->E[nx1 + (nx1 + nz1) * nx1], nx1 + nz1, &E22, 0, 0);

    blasfeo_pack_dmat(nx1, nx1, model->A, nx1 + nz1, &A1, 0, 0);
    blasfeo_pack_dmat(nz1, nx1, &model->A[nx1], nx1 + nz1, &A2, 0, 0);

    blasfeo_pack_dmat(nx1, nu, model->B, nx1 + nz1, &B1, 0, 0);
    blasfeo_pack_dmat(nz1, nu, &model->B[nx1], nx1 + nz1, &B2, 0, 0);

    blasfeo_pack_dmat(nx1, n_out, model->C, nx1 + nz1, &C1, 0, 0);
    blasfeo_pack_dmat(nz1, n_out, &model->C[nx1], nx1 + nz1, &C2, 0, 0);

    blasfeo_pack_dmat(nxz2, nx2, A_LO, nxz2, ALO, 0, 0);
    blasfeo_pack_dmat(nxz2, nu, B_LO, nxz2, BLO, 0, 0);
    blasfeo_pack_dmat(nuhat, nu, model->L_u, nuhat, Lu, 0, 0);

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
        blasfeo_dgecp(nz1, nx1, &A2, 0, 0, &AA2, ii * nz1, 0);
        blasfeo_dgecp(nx1, nu, &B1, 0, 0, &BB1, ii * nx1, 0);
        blasfeo_dgecp(nz1, nu, &B2, 0, 0, &BB2, ii * nz1, 0);

        blasfeo_dgecp(nx1, n_out, &C1, 0, 0, &CC1, ii * nx1, ii * n_out);
        blasfeo_dgecp(nz1, n_out, &C2, 0, 0, &CC2, ii * nz1, ii * n_out);

        blasfeo_dgecpsc(nx1, nz1, -1.0, &E12, 0, 0, &DD1, ii * nx1, ii * nz1);
        blasfeo_dgecpsc(nz1, nx1, -1.0, &E21, 0, 0, &DD2, ii * nz1, ii * nx1);

        blasfeo_dgecp(nx1, nx1, &E11, 0, 0, &EE1, ii * nx1, ii * nx1);
        blasfeo_dgecp(nz1, nz1, &E22, 0, 0, &EE2, ii * nz1, ii * nz1);

        blasfeo_pack_dmat(ny, nz1, model->L_z, ny, &LLZ, ii * ny, ii * nz1);
        blasfeo_pack_dmat(ny, nx1, model->L_x, ny, &LLx, ii * ny, 0);
        blasfeo_pack_dmat(ny, nx1, model->L_xdot, ny, &LLK, ii * ny, ii * nx1);

        blasfeo_pack_dmat(nxz2, nx2, A_LO, nxz2, dK2_dx2, ii * nxz2, 0);
          // dK2_dx2 = repmat(s.ALO,opts.n_stages,1); contains now r.h.s. for dK2_dx2

        blasfeo_pack_dvec(nx1, model->c, 1, &cc1, ii * nx1);
        blasfeo_pack_dvec(nz1 , &model->c[nx1], 1, &cc2, ii * nz1);
    }

    for (int ii = 0; ii < num_stages; ii++)
    {  // perform kronecker product
        for (int jj = 0; jj < num_stages; jj++)
        {
            blasfeo_dgead(nz1, nx1, A_dt[ii * num_stages + jj], &A2, 0, 0,
                             &DD2, jj * nz1, ii * nx1);
            blasfeo_dgead(nx1, nx1, -A_dt[ii * num_stages + jj], &A1, 0, 0, &EE1, jj * nx1,
                          ii * nx1);
            blasfeo_dgead(ny, nx1, A_dt[ii * num_stages + jj], &LLx, 0, 0, &LLK, jj * ny, ii * nx1);
        }
    }

    /************************************************
     *   Compute QQ1 KK*, ZZ* via QQ1
     ************************************************/
    // SOLVE EE1 \ DD1, ... EE1 \ AA1;
    blasfeo_dgetrf_rp(nK1, nK1, &EE1, 0, 0, &EE1, 0, 0, ipivEE1);  // factorize EE1

    blasfeo_drowpe(nK1, ipivEE1, &AA1);  // permute also rhs
    blasfeo_dtrsm_llnu(nK1, nx1, 1.0, &EE1, 0, 0, &AA1, 0, 0, &AA1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nx1, 1.0, &EE1, 0, 0, &AA1, 0, 0, &AA1, 0,
                       0);              // AA1 now contains EE1\AA1
    blasfeo_drowpe(nK1, ipivEE1, &BB1);  // permute also rhs
    blasfeo_dtrsm_llnu(nK1, nu, 1.0, &EE1, 0, 0, &BB1, 0, 0, &BB1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nu, 1.0, &EE1, 0, 0, &BB1, 0, 0, &BB1, 0,
                       0);               // BB1 now contains EE1\BB1
    blasfeo_drowpe(nK1, ipivEE1, &CC1);  // permute also rhs

    blasfeo_dtrsm_llnu(nK1, nvv, 1.0, &EE1, 0, 0, &CC1, 0, 0, &CC1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nvv, 1.0, &EE1, 0, 0, &CC1, 0, 0, &CC1, 0,
                       0);               // CC1 now contains EE1\CC1

    blasfeo_drowpe(nK1, ipivEE1, &DD1);  // permute also rhs
    blasfeo_dtrsm_llnu(nK1, nZ1, 1.0, &EE1, 0, 0, &DD1, 0, 0, &DD1, 0, 0);
    blasfeo_dtrsm_lunn(nK1, nZ1, 1.0, &EE1, 0, 0, &DD1, 0, 0, &DD1, 0,
                       0);  // DD1 now contains EE1\DD1

    // printf("cc1 = (in precompute) \n");
    // blasfeo_print_exp_dvec(nK1, &cc1, 0);

    blasfeo_dvecpe(nK1, ipivEE1, &cc1, 0);  // permute rhs
    blasfeo_dtrsv_lnu(nK1, &EE1, 0, 0, &cc1, 0, &cc1, 0);
    blasfeo_dtrsv_unn(nK1, &EE1, 0, 0, &cc1, 0, &cc1, 0);
    // cc1 now contains EE1\cc1

    // SOLVE EE2 \ DD2, ... EE2 \ AA2;
    blasfeo_dgetrf_rp(nZ1, nZ1, &EE2, 0, 0, &EE2, 0, 0, ipivEE2);  // factorize EE2

    blasfeo_drowpe(nZ1, ipivEE2, &AA2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ1, nx1, 1.0, &EE2, 0, 0, &AA2, 0, 0, &AA2, 0, 0);
    blasfeo_dtrsm_lunn(nZ1, nx1, 1.0, &EE2, 0, 0, &AA2, 0, 0, &AA2, 0, 0);
                                    // AA2 now contains EE2\AA2

    blasfeo_drowpe(nZ1, ipivEE2, &BB2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ1, nu, 1.0, &EE2, 0, 0, &BB2, 0, 0, &BB2, 0, 0);
    blasfeo_dtrsm_lunn(nZ1, nu, 1.0, &EE2, 0, 0, &BB2, 0, 0, &BB2, 0, 0);
                                    // BB2 now contains EE2\BB2

    blasfeo_drowpe(nZ1, ipivEE2, &CC2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ1, nvv, 1.0, &EE2, 0, 0, &CC2, 0, 0, &CC2, 0, 0);
    blasfeo_dtrsm_lunn(nZ1, nvv, 1.0, &EE2, 0, 0, &CC2, 0, 0, &CC2, 0, 0);
                                    // CC2 now contains EE2\CC2

    blasfeo_drowpe(nZ1, ipivEE2, &DD2);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ1, nK1, 1.0, &EE2, 0, 0, &DD2, 0, 0, &DD2, 0, 0);
    blasfeo_dtrsm_lunn(nZ1, nK1, 1.0, &EE2, 0, 0, &DD2, 0, 0, &DD2, 0, 0);
                                    // DD2 now contains EE2\DD2

    blasfeo_dvecpe(nZ1, ipivEE2, &cc2, 0);  // permute rhs
    blasfeo_dtrsv_lnu(nZ1, &EE2, 0, 0, &cc2, 0, &cc2, 0);
    blasfeo_dtrsv_unn(nZ1, &EE2, 0, 0, &cc2, 0, &cc2, 0);
    // cc2 now contains EE2\cc2

    /* Build and factorize QQ1 */
    blasfeo_dgemm_nn(nZ1, nZ1, nK1, -1.0, &DD2, 0, 0, &DD1, 0, 0, 0.0, &QQ1, 0, 0, &QQ1, 0,
                     0);                                               // QQ1 = -DD2*DD1
    blasfeo_ddiare(nZ1, 1.0, &QQ1, 0, 0);                               // add eye(nZ1) to QQ1
    blasfeo_dgetrf_rp(nZ1, nZ1, &QQ1, 0, 0, &QQ1, 0, 0, ipivQQ1);  // factorize QQ1

    /* build ZZv */
    blasfeo_dgemm_nn(nZ1, nvv, nK1, 1.0, &DD2, 0, 0, &CC1, 0, 0, 0.0, ZZv, 0, 0, ZZv, 0,
                     0);                                  // ZZv = DD2 * CC1
    blasfeo_dgead(nZ1, nvv, 1.0, &CC2, 0, 0, ZZv, 0, 0);  // ZZv = ZZv + CC2;
    // solve QQ1\ZZv and store result in ZZv;
    blasfeo_drowpe(nZ1, ipivQQ1, ZZv);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ1, nvv, 1.0, &QQ1, 0, 0, ZZv, 0, 0, ZZv, 0, 0);
    blasfeo_dtrsm_lunn(nZ1, nvv, 1.0, &QQ1, 0, 0, ZZv, 0, 0, ZZv, 0,
                       0);  // ZZv now contains QQ1\ZZv

    /* build ZZu */
    blasfeo_dgemm_nn(nZ1, nu, nK1, 1.0, &DD2, 0, 0, &BB1, 0, 0, 0.0, ZZu, 0, 0, ZZu, 0,
                     0);                                 // ZZu = DD2 * BB1
    blasfeo_dgead(nZ1, nu, 1.0, &BB2, 0, 0, ZZu, 0, 0);  // ZZu = ZZu + BB2;
    // solve QQ1\ZZu and store result in ZZu;
    blasfeo_drowpe(nZ1, ipivQQ1, ZZu);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ1, nu, 1.0, &QQ1, 0, 0, ZZu, 0, 0, ZZu, 0, 0);
    blasfeo_dtrsm_lunn(nZ1, nu, 1.0, &QQ1, 0, 0, ZZu, 0, 0, ZZu, 0,
                       0);  // ZZu now contains QQ1\ZZu

    /*  build ZZx  */
    blasfeo_dgemm_nn(nZ1, nx1, nK1, 1.0, &DD2, 0, 0, &AA1, 0, 0, 0.0, ZZx, 0, 0, ZZx, 0,
                     0);                                  // ZZx = DD2 * AA1;
    blasfeo_dgead(nZ1, nx1, 1.0, &AA2, 0, 0, ZZx, 0, 0);  // ZZx = ZZx + AA2;

    // solve QQ1\ZZx and store result in ZZx;
    blasfeo_drowpe(nZ1, ipivQQ1, ZZu);  // permute also rhs
    blasfeo_dtrsm_llnu(nZ1, nx1, 1.0, &QQ1, 0, 0, ZZx, 0, 0, ZZx, 0, 0);
    blasfeo_dtrsm_lunn(nZ1, nx1, 1.0, &QQ1, 0, 0, ZZx, 0, 0, ZZx, 0,
                       0);  // ZZx now contains QQ1\ZZx

    /* build ZZ0  */
    blasfeo_dgemv_n(nZ1, nK1, 1.0, &DD2, 0, 0, &cc1, 0, 1.0, &cc2, 0, ZZ0, 0);
                                                // ZZ0 = DD2 * cc1 + cc2;
    // solve QQ1\ZZ0 and store result in ZZ0;
    blasfeo_dvecpe(nZ1, ipivQQ1, ZZ0, 0);  // permute rhs
    blasfeo_dtrsv_lnu(nZ1, &QQ1, 0, 0, ZZ0, 0, ZZ0, 0);
    blasfeo_dtrsv_unn(nZ1, &QQ1, 0, 0, ZZ0, 0, ZZ0, 0);

    /* build KKv, KKu, KKx */
    blasfeo_dgemm_nn(nK1, nvv, nZ1, 1.0, &DD1, 0, 0, ZZv, 0, 0, 1.0, &CC1, 0, 0, KKv, 0,
                     0);  // KKv = DD1 * ZZv + CC1
    blasfeo_dgemm_nn(nK1, nu, nZ1, 1.0, &DD1, 0, 0, ZZu, 0, 0, 1.0, &BB1, 0, 0, KKu, 0,
                     0);  // KKu = DD1 * ZZu + BB1
    blasfeo_dgemm_nn(nK1, nx1, nZ1, 1.0, &DD1, 0, 0, ZZx, 0, 0, 1.0, &AA1, 0, 0, KKx, 0,
                     0);  // KKx = DD1 * ZZx + AA1
    blasfeo_dgemv_n(nK1, nZ1, 1.0, &DD1, 0, 0, ZZ0, 0, 1.0, &cc1, 0, KK0, 0);
                          // KK0 = DD1 * ZZ0 + cc1;

    /* build YYx */
    blasfeo_dgemm_nn(nyy, nx1, nK1, 1.0, &LLK, 0, 0, KKx, 0, 0, 1.0, &LLx, 0, 0, YYx, 0, 0);
    blasfeo_dgemm_nn(nyy, nx1, nZ1, 1.0, &LLZ, 0, 0, ZZx, 0, 0, 1.0, YYx, 0, 0, YYx, 0, 0);

    /* build YYu */
    blasfeo_dgemm_nn(nyy, nu, nK1, 1.0, &LLK, 0, 0, KKu, 0, 0, 0.0, YYu, 0, 0, YYu, 0, 0);
    blasfeo_dgemm_nn(nyy, nu, nZ1, 1.0, &LLZ, 0, 0, ZZu, 0, 0, 1.0, YYu, 0, 0, YYu, 0, 0);

    /* build YYv */
    blasfeo_dgemm_nn(nyy, nvv, nK1, 1.0, &LLK, 0, 0, KKv, 0, 0, 0.0, YYv, 0, 0, YYv, 0, 0);
    blasfeo_dgemm_nn(nyy, nvv, nZ1, 1.0, &LLZ, 0, 0, ZZv, 0, 0, 1.0, YYv, 0, 0, YYv, 0, 0);

    /* build YY0 */
    blasfeo_dgemv_n(nyy, nK1, 1.0, &LLK, 0, 0, KK0, 0, 0.0, YY0, 0, YY0, 0);
                                                // YY0 = LLK * KK0;
    blasfeo_dgemv_n(nyy, nZ1, 1.0, &LLZ, 0, 0, ZZ0, 0, 1.0, YY0, 0, YY0, 0);
                                                // YY0 = YY0 + LLZ * ZZ0
    // build M2 in M2_LU
    for (int ii = 0; ii < num_stages; ii++)
    {
        blasfeo_pack_dmat(nxz2, nxz2, model->E_LO, nx2+nz2, M2_LU, ii * nxz2, ii * nxz2);
        for (int jj = 0; jj < num_stages; jj++)
        {
            blasfeo_dgead(nxz2, nx2, - A_dt[ii * num_stages + jj], ALO, 0, 0,
                         M2_LU, jj * nxz2, ii * nxz2);
        }
    }

    // factorize M2 in M2_LU
    blasfeo_dgetrf_rp(nK2, nK2, M2_LU, 0, 0, M2_LU, 0, 0, ipivM2);


    // solve dK2_dx2 = M2 \ dK2_dx2 to obtain dK2_dx2
    blasfeo_drowpe(nK2, ipivM2, dK2_dx2);  // permute rhs dK2_dx2
    blasfeo_dtrsm_llnu(nK2, nx2, 1.0, M2_LU, 0, 0, dK2_dx2, 0, 0, dK2_dx2, 0, 0);
    blasfeo_dtrsm_lunn(nK2, nx2, 1.0, M2_LU, 0, 0, dK2_dx2, 0, 0, dK2_dx2, 0, 0);

    if (!model->nontrivial_f_LO) // precompute dK2_du if f_LO is trivial
    {
        for (int ii = 0; ii < num_stages; ii++)
        {
            blasfeo_dgecp(nxz2, nu, BLO, 0, 0, dK2_du, ii * nxz2,  0);
                    // set ith block of dK2_du to BLO
        }
        blasfeo_drowpe(nK2, ipivM2, dK2_du);  // permute also rhs
        blasfeo_dtrsm_llnu(nK2, nu, 1.0, M2_LU, 0, 0, dK2_du, 0, 0, dK2_du, 0, 0);
        blasfeo_dtrsm_lunn(nK2, nu, 1.0, M2_LU, 0, 0, dK2_du, 0, 0, dK2_du, 0, 0);

        // precompute dx2f_dx2u
        for (int ii = 0; ii < nx2; ii++)
            blasfeo_dgein1(1.0, dx2f_dx2u, ii, ii);

        for (int ii = 0; ii < num_stages; ii++)
        {
            blasfeo_dgead(nx2, nx2, b_dt[ii], dK2_dx2, ii * nxz2, 0, dx2f_dx2u, 0, 0);
            blasfeo_dgead(nx2, nu,  b_dt[ii], dK2_du,  ii * nxz2, 0, dx2f_dx2u, 0, nx2);
        }
    }

    // if (opts->sens_algebraic)
    // {

    //     /* ONLY FOR ALGEBRAIC SENSITIVITIES */
    //     /* Compute K0*, Z0*, Y0*, * \in {f,x,u} */
    //     // pack matrices to memory for algebraic sensitivities
    //     blasfeo_pack_dmat(ny, nx1, model->L_x, ny, Lx, 0, 0);
    //     blasfeo_pack_dmat(ny, nx1, model->L_xdot, ny, Lxdot, 0, 0);
    //     blasfeo_pack_dmat(ny, nz1, model->L_z, ny, Lz, 0, 0);

    //     // SOLVE E11 \ A1,   E11 \ B1,      E11 \ C1,   E11\E12;
    //     blasfeo_dgetrf_rp(nx1, nx1, &E11, 0, 0, &E11, 0, 0, ipivEE1);  // factorize E11

    //     blasfeo_drowpe(nx1, ipivEE1, &A1);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nx1, nx1, 1.0, &E11, 0, 0, &A1, 0, 0, &A1, 0, 0);
    //     blasfeo_dtrsm_lunn(nx1, nx1, 1.0, &E11, 0, 0, &A1, 0, 0, &A1, 0,
    //                     0);              // A1 now contains E11\A1
    //     blasfeo_drowpe(nx1, ipivEE1, &B1);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nx1, nu, 1.0, &E11, 0, 0, &B1, 0, 0, &B1, 0, 0);
    //     blasfeo_dtrsm_lunn(nx1, nu, 1.0, &E11, 0, 0, &B1, 0, 0, &B1, 0,
    //                     0);               // B1 now contains E11\B1

    //     blasfeo_drowpe(nx1, ipivEE1, &C1);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nx1, n_out, 1.0, &E11, 0, 0, &C1, 0, 0, &C1, 0, 0);
    //     blasfeo_dtrsm_lunn(nx1, n_out, 1.0, &E11, 0, 0, &C1, 0, 0, &C1, 0,
    //                     0);               // C1 now contains E11\C1

    //     blasfeo_drowpe(nx1, ipivEE1, &E12);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nx1, nz1, 1.0, &E11, 0, 0, &E12, 0, 0, &E12, 0, 0);
    //     blasfeo_dtrsm_lunn(nx1, nz1, 1.0, &E11, 0, 0, &E12, 0, 0, &E12, 0,
    //                     0);  // E12 now contains E11\E12


    //     // SOLVE E22 \ A2,   E22 \ B2,      E22 \ C2,   E22\E21;
    //     blasfeo_dgetrf_rp(nz1, nz1, &E22, 0, 0, &E22, 0, 0, ipivEE2);  // factorize E22

    //     blasfeo_drowpe(nz1, ipivEE2, &A2);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nz1, nx1, 1.0, &E22, 0, 0, &A2, 0, 0, &A2, 0, 0);
    //     blasfeo_dtrsm_lunn(nz1, nx1, 1.0, &E22, 0, 0, &A2, 0, 0, &A2, 0, 0);
    //                                     // A2 now contains E22\A2

    //     blasfeo_drowpe(nz1, ipivEE2, &B2);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nz1, nu, 1.0, &E22, 0, 0, &B2, 0, 0, &B2, 0, 0);
    //     blasfeo_dtrsm_lunn(nz1, nu, 1.0, &E22, 0, 0, &B2, 0, 0, &B2, 0, 0);
    //                                     // B2 now contains E22\B2

    //     blasfeo_drowpe(nz1, ipivEE2, &C2);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nz1, n_out, 1.0, &E22, 0, 0, &C2, 0, 0, &C2, 0, 0);
    //     blasfeo_dtrsm_lunn(nz1, n_out, 1.0, &E22, 0, 0, &C2, 0, 0, &C2, 0, 0);
    //                                     // C2 now contains E22\C2

    //     blasfeo_drowpe(nz1, ipivEE2, &E21);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nz1, nx1, 1.0, &E22, 0, 0, &E21, 0, 0, &E21, 0, 0);
    //     blasfeo_dtrsm_lunn(nz1, nx1, 1.0, &E22, 0, 0, &E21, 0, 0, &E21, 0, 0);
    //                                     // E21 now contains E22\E21

    //     /* Build and factorize Q1 */
    //     blasfeo_dgemm_nn(nz1, nz1, nx1, -1.0, &E21, 0, 0, &E12, 0, 0, 0.0, &Q1, 0, 0, &Q1, 0,
    //                     0);                                               // Q1 = -E21*E12
    //     blasfeo_ddiare(nz1, 1.0, &Q1, 0, 0);                               // add eye(nz1) to Q1
    //     blasfeo_dgetrf_rp(nz1, nz1, &Q1, 0, 0, &Q1, 0, 0, ipivQQ1);  // factorize Q1

    //     /* build Z0v */
    //     blasfeo_dgemm_nn(nz1, n_out, nx1, 1.0, &E21, 0, 0, &C1, 0, 0, 1.0, &C2, 0, 0, Z0v, 0,
    //                     0);                                  // Z0v = E21 * C1 + C2
    //     // blasfeo_dgead(nz1, n_out, 1.0, &C2, 0, 0, Z0v, 0, 0);  // Z0v = Z0v + C2;

    //     // solve Q1\Z0v and store result in Z0v;
    //     blasfeo_drowpe(nz1, ipivQQ1, Z0v);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nz1, n_out, 1.0, &Q1, 0, 0, Z0v, 0, 0, Z0v, 0, 0);
    //     blasfeo_dtrsm_lunn(nz1, n_out, 1.0, &Q1, 0, 0, Z0v, 0, 0, Z0v, 0,
    //                     0);  // Z0v now contains Q1\Z0v

    //     /* build Z0u */
    //     blasfeo_dgemm_nn(nz1, nu, nx1, 1.0, &E21, 0, 0, &B1, 0, 0, 0.0, Z0u, 0, 0, Z0u, 0,
    //                     0);                                 // Z0u = E21 * B1
    //     blasfeo_dgead(nz1, nu, 1.0, &B2, 0, 0, Z0u, 0, 0);  // Z0u = Z0u + B2;
    //     // solve Q1\Z0u and store result in Z0u;
    //     blasfeo_drowpe(nz1, ipivQQ1, Z0u);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nz1, nu, 1.0, &Q1, 0, 0, Z0u, 0, 0, Z0u, 0, 0);
    //     blasfeo_dtrsm_lunn(nz1, nu, 1.0, &Q1, 0, 0, Z0u, 0, 0, Z0u, 0,
    //                     0);  // Z0u now contains Q1\Z0u

    //     /*  build Z0x  */
    //     blasfeo_dgemm_nn(nz1, nx1, nx1, 1.0, &E21, 0, 0, &A1, 0, 0, 0.0, Z0x, 0, 0, Z0x, 0,
    //                     0);                                  // Z0x = E21 * A1;
    //     blasfeo_dgead(nz1, nx1, 1.0, &A2, 0, 0, Z0x, 0, 0);  // Z0x = Z0x + A2;

    //     // solve Q1\Z0x and store result in Z0x;
    //     blasfeo_drowpe(nz1, ipivQQ1, Z0u);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nz1, nx1, 1.0, &Q1, 0, 0, Z0x, 0, 0, Z0x, 0, 0);
    //     blasfeo_dtrsm_lunn(nz1, nx1, 1.0, &Q1, 0, 0, Z0x, 0, 0, Z0x, 0,
    //                     0);  // Z0x now contains Q1\Z0x

    //     /* build K0v, K0u, K0x */
    //     blasfeo_dgemm_nn(nx1, n_out, nz1, -1.0, &E12, 0, 0, Z0v, 0, 0, 1.0, &C1, 0, 0, K0v, 0,
    //                     0);  // K0v = -E12 * Z0v + C1
    //     blasfeo_dgemm_nn(nx1, nu, nz1, -1.0, &E12, 0, 0, Z0u, 0, 0, 1.0, &B1, 0, 0, K0u, 0,
    //                     0);  // K0u = -E12 * Z0u + B1
    //     blasfeo_dgemm_nn(nx1, nx1, nz1, -1.0, &E12, 0, 0, Z0x, 0, 0, 1.0, &A1, 0, 0, K0x, 0,
    //                     0);  // K0x = -E12 * Z0x + A1

    //     /* build Y0x */
    //     blasfeo_dgemm_nn(ny, nx1, nx1, 1.0, Lxdot, 0, 0, K0x, 0, 0, 1.0, Lx, 0, 0, Y0x, 0, 0);
    //     blasfeo_dgemm_nn(ny, nx1, nz1, 1.0, Lz, 0, 0, Z0x, 0, 0, 1.0, Y0x, 0, 0, Y0x, 0, 0);

    //     /* build Y0u */
    //     blasfeo_dgemm_nn(ny, nu, nx1, 1.0, Lxdot, 0, 0, K0u, 0, 0, 0.0, Y0u, 0, 0, Y0u, 0, 0);
    //     blasfeo_dgemm_nn(ny, nu, nz1, 1.0, Lz, 0, 0, Z0u, 0, 0, 1.0, Y0u, 0, 0, Y0u, 0, 0);

    //     /* build Y0v */
    //     blasfeo_dgemm_nn(ny, n_out, nx1, 1.0, Lxdot, 0, 0, K0v, 0, 0, 0.0, Y0v,
    //                      0, 0, Y0v, 0, 0);
    //     blasfeo_dgemm_nn(ny, n_out, nz1, 1.0, Lz, 0, 0, Z0v, 0, 0, 1.0, Y0v, 0, 0, Y0v, 0, 0);

    //     /** for z2 **/           // factorize ELO
    //     blasfeo_pack_dmat(nxz2, nxz2, model->E_LO, nxz2, ELO_LU, 0, 0);
    //     blasfeo_dgetrf_rp(nxz2, nxz2, ELO_LU, 0, 0, ELO_LU, 0, 0, ipiv_ELO);
    //     //
    //     // blasfeo_print_dmat(nxz2, nxz2, ELO_LU, )
    //     blasfeo_pack_dmat(nxz2, nx2, model->A_LO, nxz2, ELO_inv_ALO, 0, 0);
    //     blasfeo_drowpe(nxz2, ipiv_ELO, ELO_inv_ALO);  // permute also rhs
    //     blasfeo_dtrsm_llnu(nxz2, nx2, 1.0, ELO_LU, 0, 0, ELO_inv_ALO, 0, 0, ELO_inv_ALO, 0, 0);
    //     blasfeo_dtrsm_lunn(nxz2, nx2, 1.0, ELO_LU, 0, 0, ELO_inv_ALO, 0, 0, ELO_inv_ALO, 0, 0);
    // }

    // generate sensitivities
    mem->first_call = true;
    if (model->fully_linear)
    {
        // set identity seed
        for (int jj = 0; jj < nx * (nx + nu); jj++)
            in->S_forw[jj] = 0.0;
        for (int jj = 0; jj < nx; jj++)
            in->S_forw[jj * (nx + 1)] = 1.0;
        in->identity_seed = true;
        sim_gnsf(config_, in, out, opts_, mem_, work_);
    }
    mem->first_call = false;

    return status;
    // double precomputation_time = acados_toc(&atimer) * 1000;
    // printf("time 2 precompute = %f [ms]\n", precomputation_time);
} // sim_gnsf_precompute



/************************************************
 * memory
 ************************************************/

acados_size_t sim_gnsf_memory_calculate_size(void *config, void *dims_, void *opts_)
{
    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_opts *opts = opts_;

    // necessary integers
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nz      = dims->nz;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    int num_stages = opts->ns;

    int nxz2 = nx2 + nz2;
    int nvv = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nxz2;
    int nZ1 = num_stages * nz1;

    acados_size_t size = sizeof(sim_gnsf_memory);

    // scaled butcher table
    size += num_stages * num_stages * sizeof(double);  // A_dt
    size += 2 * num_stages * sizeof(double);           // b_dt, c_butcher;

    size += n_out * sizeof(double); // phi_guess

    size += nK2 * sizeof(int);      // ipivM2

    if (opts->sens_algebraic)
    {
        size += nxz2 * sizeof(int); // ipiv_ELO

        size += 9 * sizeof(struct blasfeo_dmat);  // Z0*, K0*, Y0*, *in {x,u,v}
        size += 2 * sizeof(struct blasfeo_dmat); // ELO_LU, ELO_inv_ALO
        size += 3 * sizeof(struct blasfeo_dmat); // Lx, Lxdot, Lz
    }

    // precomputed matrices
    size += blasfeo_memsize_dmat(nK1, nvv);  // KKv
    size += blasfeo_memsize_dmat(nK1, nx1);  // KKx
    size += blasfeo_memsize_dmat(nK1, nu);   // KKu

    size += blasfeo_memsize_dmat(nyy, nvv);  // YYv
    size += blasfeo_memsize_dmat(nyy, nx1);  // YYx
    size += blasfeo_memsize_dmat(nyy, nu);   // YYu

    size += blasfeo_memsize_dmat(nZ1, nvv);  // ZZv
    size += blasfeo_memsize_dmat(nZ1, nx1);  // ZZx
    size += blasfeo_memsize_dmat(nZ1, nu);   // ZZu

    size += blasfeo_memsize_dmat(nx2 + nz2, nx2);  // ALO
    size += blasfeo_memsize_dmat(nx2 + nz2, nu);  // BLO
    size += blasfeo_memsize_dmat(nK2, nK2);  // M2_LU
    size += blasfeo_memsize_dmat(nK2, nx2);  // dK2_dx2
    size += blasfeo_memsize_dmat(nK2, nu);          // dK2_du
    size += blasfeo_memsize_dmat(nx2, nx2 + nu); // dx2f_dx2u

    size += blasfeo_memsize_dmat(nuhat, nu);  // Lu

    size += blasfeo_memsize_dmat(nx, nx + nu);  // S_forw
    size += blasfeo_memsize_dmat(nz, nx + nu);  // S_algebraic
    // if (opts->sens_algebraic)
    // {
    //     // matrices only needed for algebraic sensitivities
    //     size += blasfeo_memsize_dmat(nz1, nx1);     // Z0x
    //     size += blasfeo_memsize_dmat(nz1, nu);      // Z0u
    //     size += blasfeo_memsize_dmat(nz1, n_out);   // Z0v

    //     size += blasfeo_memsize_dmat(ny, nx1);      // Y0x
    //     size += blasfeo_memsize_dmat(ny, nu);       // Y0u
    //     size += blasfeo_memsize_dmat(ny, n_out);    // Y0v

    //     size += blasfeo_memsize_dmat(nx1, nx1);     // K0x
    //     size += blasfeo_memsize_dmat(nx1, nu);      // K0u
    //     size += blasfeo_memsize_dmat(nx1, n_out);   // K0v

    //     size += blasfeo_memsize_dmat(nxz2, nxz2);   // ELO_LU
    //     size += blasfeo_memsize_dmat(nxz2, nx2);    // ELO_inv_ALO

    //     size += 2 * blasfeo_memsize_dmat(ny, nx1);   // Lx, Lxdot
    //     size += blasfeo_memsize_dmat(ny, nz1);        // Lz
    // }

    size += blasfeo_memsize_dvec(nZ1);  // ZZ0
    size += blasfeo_memsize_dvec(nK1);  // KK0
    size += blasfeo_memsize_dvec(nyy);  // YY0

    size += 1 * 64;  // corresponds to memory alignment
    size += 2 * 8;  // initial memory alignment, alignment for doubles
    make_int_multiple_of(64, &size);

    return size;
}



void *sim_gnsf_memory_assign(void *config, void *dims_, void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_opts *opts = opts_;

    // necessary integers
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nz      = dims->nz;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    int num_stages = opts->ns;

    int nxz2 = nx2 + nz2;
    int nvv = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * nxz2;
    int nZ1 = num_stages * nz1;

    // initial align
    align_char_to(8, &c_ptr);

    // struct
    sim_gnsf_memory *mem = (sim_gnsf_memory *) c_ptr;
    c_ptr += sizeof(sim_gnsf_memory);

    // if (opts->sens_algebraic){
    //     assign_and_advance_int(nxz2, &mem->ipiv_ELO, &c_ptr);
    // }
    assign_and_advance_int(nK2, &mem->ipivM2, &c_ptr);
    align_char_to(8, &c_ptr);

    // assign doubles
    assign_and_advance_double(num_stages * num_stages, &mem->A_dt, &c_ptr);
    assign_and_advance_double(num_stages, &mem->b_dt, &c_ptr);
    assign_and_advance_double(num_stages, &mem->c_butcher, &c_ptr);

    assign_and_advance_double(n_out, &mem->phi_guess, &c_ptr);

    // initialize with zeros for default initialization;
    for (int ii = 0; ii < n_out; ii++) {
        mem->phi_guess[ii] = 0.0;
    }

    // blasfeo_dmat structs
    // if (opts->sens_algebraic){

    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Z0x, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Z0u, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Z0v, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Y0x, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Y0u, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Y0v, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->K0x, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->K0u, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->K0v, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->ELO_LU, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->ELO_inv_ALO, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Lx, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Lxdot, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_structs(1, &mem->Lz, &c_ptr);
    // }
    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dmat_mem
    assign_and_advance_blasfeo_dmat_mem(nK1, nvv, &mem->KKv, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nx1, &mem->KKx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nu, &mem->KKu, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nyy, nvv, &mem->YYv, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nx1, &mem->YYx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nyy, nu, &mem->YYu, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nZ1, nvv, &mem->ZZv, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nx1, &mem->ZZx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nu, &mem->ZZu, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx2 + nz2, nx2, &mem->ALO, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx2 + nz2, nu, &mem->BLO, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nK2, &mem->M2_LU, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nx2, &mem->dK2_dx2, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK2, nu, &mem->dK2_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx2, nx2 + nu, &mem->dx2f_dx2u, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nuhat, nu, &mem->Lu, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &mem->S_forw, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, nx + nu, &mem->S_algebraic, &c_ptr);

    // if (opts->sens_algebraic){
    //     // for algebraic sensitivity propagation
    //     assign_and_advance_blasfeo_dmat_mem(ny, nx1, mem->Lx, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(ny, nx1, mem->Lxdot, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(ny, nz1, mem->Lz, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_mem(nz1,   nx1, mem->Z0x, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(nz1,   nu , mem->Z0u, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(nz1, n_out, mem->Z0v, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_mem(ny,   nx1, mem->Y0x, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(ny,   nu , mem->Y0u, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(ny, n_out, mem->Y0v, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_mem(nx1, nx1  , mem->K0x, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(nx1, nu   , mem->K0u, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(nx1, n_out, mem->K0v, &c_ptr);

    //     assign_and_advance_blasfeo_dmat_mem(nxz2, nxz2, mem->ELO_LU, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(nxz2, nx2 , mem->ELO_inv_ALO, &c_ptr);
    // }

    assign_and_advance_blasfeo_dvec_mem(nZ1, &mem->ZZ0, &c_ptr);  // ZZ0
    assign_and_advance_blasfeo_dvec_mem(nyy, &mem->YY0, &c_ptr);  // YY0
    assign_and_advance_blasfeo_dvec_mem(nK1, &mem->KK0, &c_ptr);  // KK0


    assert((char *) raw_memory + sim_gnsf_memory_calculate_size(config, dims_, opts_) >= c_ptr);
    return mem;
}



int sim_gnsf_memory_set(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    sim_gnsf_memory *mem = (sim_gnsf_memory *) mem_;
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;

    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "phi_guess"))
    {
        double *phi_guess = value;
        for (int ii=0; ii < dims->n_out; ii++)
            mem->phi_guess[ii] = phi_guess[ii];
    }
    else if (!strcmp(field, "guesses_blasfeo"))
    {
        blasfeo_unpack_dvec(dims->n_out, value, 0, mem->phi_guess, 1);
    }
    else
    {
        printf("sim_gnsf_memory_set field %s is not supported! \n", field);
        exit(1);
    }

    return status;
}



int sim_gnsf_memory_set_to_zero(void *config_, void * dims_, void *opts_, void *mem_, const char *field)
{
    sim_gnsf_memory *mem = (sim_gnsf_memory *) mem_;
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;

    int status = ACADOS_SUCCESS;

    if (!strcmp(field, "guesses"))
    {
        for (int ii=0; ii < dims->n_out; ii++)
            mem->phi_guess[ii] = 0.0;
    }
    else
    {
        printf("sim_gnsf_memory_set_to_zero field %s is not supported! \n", field);
        exit(1);
    }

    return status;
}



void sim_gnsf_memory_get(void *config_, void *dims_, void *mem_, const char *field, void *value)
{
    sim_gnsf_memory *mem = mem_;

    if (!strcmp(field, "time_sim"))
    {
		double *ptr = value;
		*ptr = mem->time_sim;
	}
    else if (!strcmp(field, "time_sim_ad"))
    {
		double *ptr = value;
		*ptr = mem->time_ad;
	}
    else if (!strcmp(field, "time_sim_la"))
    {
		double *ptr = value;
		*ptr = mem->time_la;
	}
	else
	{
		printf("sim_gnsf_memory_get field %s is not supported! \n", field);
		exit(1);
	}
}



/************************************************
 * workspace
 ************************************************/

acados_size_t sim_gnsf_workspace_calculate_size(void *config, void *dims_, void *opts_)
{
    // typecast
    sim_gnsf_dims *dims = (sim_gnsf_dims *) dims_;
    sim_opts *opts = opts_;

    // dimension ints
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nz      = dims->nz;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    int num_stages = opts->ns;
    int num_steps = opts->num_steps;

    int nvv = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * (nx2 + nz2);
    int nZ1 = num_stages * nz1;
    int nxz2 = nx2 + nz2;

    /* Calculate workspace size for precompute function */
    acados_size_t pre_size = sizeof(gnsf_pre_workspace);

    make_int_multiple_of(8, &pre_size);
    pre_size += 1 * 8;

    pre_size += (2 * nZ1 + nK1) * sizeof(int);  // ipivEE1, ipivEE2, ipivQQ1

    make_int_multiple_of(64, &pre_size);
    pre_size += 1 * 64;

    pre_size += blasfeo_memsize_dmat(nx1, nx1);     // E11
    pre_size += blasfeo_memsize_dmat(nx1, nz1);     // E12
    pre_size += blasfeo_memsize_dmat(nz1, nx1);     // E21
    pre_size += blasfeo_memsize_dmat(nz1, nz1);     // E22

    pre_size += blasfeo_memsize_dmat(nx1, nx1);     // A1
    pre_size += blasfeo_memsize_dmat(nz1, nx1);     // A2
    pre_size += blasfeo_memsize_dmat(nx1, nu);      // B1
    pre_size += blasfeo_memsize_dmat(nz1, nu);      // B2
    pre_size += blasfeo_memsize_dmat(nx1, n_out);   // C1
    pre_size += blasfeo_memsize_dmat(nz1, n_out);   // C2

    pre_size += blasfeo_memsize_dmat(nK1, nx1);     // AA1
    pre_size += blasfeo_memsize_dmat(nZ1, nx1);     // AA2
    pre_size += blasfeo_memsize_dmat(nK1, nu);      // BB1
    pre_size += blasfeo_memsize_dmat(nZ1, nu);      // BB2

    pre_size += blasfeo_memsize_dmat(nK1, nvv);     // CC1
    pre_size += blasfeo_memsize_dmat(nZ1, nvv);     // CC2
    pre_size += blasfeo_memsize_dmat(nK1, nZ1);     // DD1
    pre_size += blasfeo_memsize_dmat(nZ1, nK1);     // DD2

    pre_size += blasfeo_memsize_dmat(nK1, nK1);     // EE1
    pre_size += blasfeo_memsize_dmat(nZ1, nZ1);     // EE2

    pre_size += blasfeo_memsize_dmat(nZ1, nZ1);     // QQ1

    pre_size += blasfeo_memsize_dmat(nyy, nZ1);     // LLZ
    pre_size += blasfeo_memsize_dmat(nyy, nx1);     // LLx
    pre_size += blasfeo_memsize_dmat(nyy, nK1);     // LLK

    pre_size += blasfeo_memsize_dvec(nK1);  // cc1
    pre_size += blasfeo_memsize_dvec(nZ1);  // cc2

    // if (opts->sens_algebraic){
    //     pre_size += blasfeo_memsize_dmat(nz1, nz1);     // Q1
    // }

    make_int_multiple_of(8, &pre_size);
    pre_size += 1 * 8;

    /* Calculate workspace size for simulation function */
    acados_size_t size = sizeof(gnsf_workspace);
    make_int_multiple_of(8, &size);
    size += 1 * 8;

    size += num_stages * sizeof(double);  // Z_work

    size += 2 * num_steps * sizeof(struct blasfeo_dvec);  // vv_traj, yy_traj
    size += num_steps * sizeof(struct blasfeo_dmat);  // f_LO_jac_traj

    size += nvv * sizeof(int);  // ipiv

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    size += num_steps * blasfeo_memsize_dvec(nvv);      // vv_traj
    size += num_steps * blasfeo_memsize_dvec(nyy);      // yy_traj

    size += blasfeo_memsize_dvec(nZ1);                   // Z1_val
    size += 2 * blasfeo_memsize_dvec(nK1);               // K1_val, x1_stage_val
    size += blasfeo_memsize_dvec(nK2);                   // f_LO_val
    size += blasfeo_memsize_dvec(nK2);                   // K2_val
    size += blasfeo_memsize_dvec(nx * (num_steps + 1));  // x0_traj
    size += blasfeo_memsize_dvec(nvv);                   // res_val
    size += blasfeo_memsize_dvec(nu);                    // u0
    size += blasfeo_memsize_dvec(nx + nu);               // lambda
    size += blasfeo_memsize_dvec(nx + nu);               // lambda_old

    size += blasfeo_memsize_dvec(nyy);              // yyu
    size += blasfeo_memsize_dvec(nyy * num_steps);  // yyss

    size += blasfeo_memsize_dvec(nK1);  // K1u
    size += blasfeo_memsize_dvec(nZ1);   // Zu
    size += 2 * blasfeo_memsize_dvec(nxz2);  // ALOtimesx02, BLOtimesu0

    size += blasfeo_memsize_dvec(nuhat);  // uhat
    size += blasfeo_memsize_dvec(nz);  // z0;

    // if (opts->sens_algebraic){
    //     size += blasfeo_memsize_dvec(nx1);  // x0dot_1;
    //     size += blasfeo_memsize_dvec(ny);  // y_one_stage
    //     size += n_out * sizeof(int);  // ipiv_vv0
    // }

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    size += num_steps * blasfeo_memsize_dmat(nK2, 2 * nx1 + nu + nz1);  // f_LO_jac_traj

    size += blasfeo_memsize_dmat(nvv, nvv);       // J_r_vv
    size += blasfeo_memsize_dmat(nvv, nx1 + nu);  // J_r_x1u

    // if (opts->sens_algebraic)
    // {
    //     size += blasfeo_memsize_dmat(n_out, n_out);  // dr0_dvv0
    //     size += blasfeo_memsize_dmat(nz1, nx1 + nu);  // dz10_dx1u
    //     size += blasfeo_memsize_dmat(nxz2, 2*nx1 + nz1 + nu);  // f_LO_jac0
    //     size += blasfeo_memsize_dmat(nxz2, nx1 + nu);  // sens_z2_rhs
    // }

    size += blasfeo_memsize_dmat(nK1, nx1);  // dK1_dx1
    size += blasfeo_memsize_dmat(nK1, nu);   // dK1_du
    size += blasfeo_memsize_dmat(nZ1, nx1);  // dZ_dx1
    size += blasfeo_memsize_dmat(nZ1, nu);   // dZ_du
    size += blasfeo_memsize_dmat(nK2, nK1);  // J_G2_K1

    size += blasfeo_memsize_dmat(nK2, nx1);     // dK2_dx1
    size += blasfeo_memsize_dmat(nK2, nvv);     // dK2_dvv
    size += blasfeo_memsize_dmat(nx, nx + nu);  // dxf_dwn
    size += blasfeo_memsize_dmat(nx, nx + nu);  // S_forw_new

    size += blasfeo_memsize_dmat(nx, nvv);   // dPsi_dvv
    size += blasfeo_memsize_dmat(nx, nx);    // dPsi_dx
    size += blasfeo_memsize_dmat(nx, nu);    // dPsi_du

    size += blasfeo_memsize_dmat(nvv, ny + nuhat);  // dPHI_dyuhat
    size += blasfeo_memsize_dmat(nz, nx + nu);  // S_algebraic_aux

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
    sim_opts *opts = opts_;

    // dimension ints
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nz      = dims->nz;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    int num_stages = opts->ns;
    int num_steps = opts->num_steps;

    int nvv = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * (nx2 + nz2);
    int nZ1 = num_stages * nz1;
    int nxz2 = nx2 + nz2;

    char *c_ptr = (char *) raw_memory;
    gnsf_workspace *workspace = (gnsf_workspace *) c_ptr;
    c_ptr += sizeof(gnsf_workspace);
    align_char_to(8, &c_ptr);

    assign_and_advance_double(num_stages, &workspace->Z_work, &c_ptr);

    assign_and_advance_int(nvv, &workspace->ipiv, &c_ptr);

    align_char_to(8, &c_ptr);

    assign_and_advance_blasfeo_dmat_structs(num_steps, &workspace->f_LO_jac_traj, &c_ptr);

    assign_and_advance_blasfeo_dvec_structs(num_steps, &workspace->vv_traj, &c_ptr);
    assign_and_advance_blasfeo_dvec_structs(num_steps, &workspace->yy_traj, &c_ptr);

    for (int ii = 0; ii < num_steps; ii++)
    {
        assign_and_advance_blasfeo_dvec_mem(nvv, workspace->vv_traj + ii, &c_ptr);    // vv_traj
        assign_and_advance_blasfeo_dvec_mem(nyy, workspace->yy_traj + ii, &c_ptr);    // yy_traj
    }

    assign_and_advance_blasfeo_dvec_mem(nK2, &workspace->f_LO_val, &c_ptr);  // f_LO_val
    assign_and_advance_blasfeo_dvec_mem(nK1, &workspace->K1_val, &c_ptr);  // K1_val
    assign_and_advance_blasfeo_dvec_mem(nK1, &workspace->x1_stage_val, &c_ptr);  // x1_stage_val
    assign_and_advance_blasfeo_dvec_mem(nZ1, &workspace->Z1_val, &c_ptr);      // Z1_val

    assign_and_advance_blasfeo_dvec_mem(nK2, &workspace->K2_val, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem((num_steps + 1) * nx, &workspace->x0_traj, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nvv, &workspace->res_val, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nu, &workspace->u0, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, &workspace->lambda, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nx + nu, &workspace->lambda_old, &c_ptr);

    assign_and_advance_blasfeo_dvec_mem(nyy, &workspace->yyu, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nyy * num_steps, &workspace->yyss, &c_ptr);

    assign_and_advance_blasfeo_dvec_mem(nK1, &workspace->K1u, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nZ1, &workspace->Zu, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nxz2, &workspace->ALOtimesx02, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nxz2, &workspace->BLOtimesu0, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nuhat, &workspace->uhat, &c_ptr);
    assign_and_advance_blasfeo_dvec_mem(nz, &workspace->z0, &c_ptr);

    // if (opts->sens_algebraic){
        // assign_and_advance_blasfeo_dvec_mem(ny, &workspace->y_one_stage, &c_ptr);
    //     assign_and_advance_blasfeo_dvec_mem(nx1, &workspace->x0dot_1, &c_ptr);

    //     assign_and_advance_int(n_out, &workspace->ipiv_vv0, &c_ptr);
    // }

    // blasfeo_dmat_mem align
    align_char_to(64, &c_ptr);
    for (int ii = 0; ii < num_steps; ii++)
    {
        assign_and_advance_blasfeo_dmat_mem(nK2, 2 * nx1 + nu + nz1, workspace->f_LO_jac_traj + ii,
                                            &c_ptr);  // f_LO_jac_traj
    }

    assign_and_advance_blasfeo_dmat_mem(nvv, nx1 + nu, &workspace->J_r_x1u, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nvv, nvv, &workspace->J_r_vv, &c_ptr);

    // if (opts->sens_algebraic)
    // {
    //     assign_and_advance_blasfeo_dmat_mem(nz1, nx1 + nu, &workspace->dz10_dx1u, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(n_out, n_out, &workspace->dr0_dvv0, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(nxz2, 2*nx1 + nz1 + nu, &workspace->f_LO_jac0, &c_ptr);
    //     assign_and_advance_blasfeo_dmat_mem(nxz2, nx1 + nu, &workspace->sens_z2_rhs, &c_ptr);
    // }

    assign_and_advance_blasfeo_dmat_mem(nK1, nx1, &workspace->dK1_dx1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nK1, nu, &workspace->dK1_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nx1, &workspace->dZ_dx1, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nZ1, nu, &workspace->dZ_du, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK2, nK1, &workspace->J_G2_K1, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nvv, ny + nuhat, &workspace->dPHI_dyuhat, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nK2, nx1, &workspace->dK2_dx1, &c_ptr);
    blasfeo_dgese(nK2, nx1, 0.0, &workspace->dK2_dx1, 0, 0);

    assign_and_advance_blasfeo_dmat_mem(nK2, nvv, &workspace->dK2_dvv, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->dxf_dwn, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx + nu, &workspace->S_forw_new, &c_ptr);

    assign_and_advance_blasfeo_dmat_mem(nx, nvv, &workspace->dPsi_dvv, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nx, &workspace->dPsi_dx, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nx, nu, &workspace->dPsi_du, &c_ptr);
    assign_and_advance_blasfeo_dmat_mem(nz, nx + nu, &workspace->S_algebraic_aux, &c_ptr);

    assert((char *) raw_memory + sim_gnsf_workspace_calculate_size(config, dims_, opts) >= c_ptr);

    return (void *) workspace;
}



int sim_gnsf(void *config, sim_in *in, sim_out *out, void *args, void *mem_, void *work_)
{
    acados_timer tot_timer, casadi_timer, la_timer;
    acados_tic(&tot_timer);

    // typecast
    sim_gnsf_memory *mem = (sim_gnsf_memory *) mem_;
    sim_opts *opts = (sim_opts *) args;
    sim_gnsf_dims *dims = (sim_gnsf_dims *) in->dims;
    gnsf_model *model = in->model;
    gnsf_workspace *workspace =
        (gnsf_workspace *) sim_gnsf_cast_workspace(config, dims, opts, work_);

    if ( opts->ns != opts->tableau_size )
    {
        printf("Error in sim_gnsf: the Butcher tableau size does not match ns.");
        exit(1);
    }
    if (opts->exact_z_output)
    {
        printf("Error in sim_gnsf: option exact_z_output = true not supported.");
        exit(1);
    }

    // necessary integers
    int nx      = dims->nx;
    int nu      = dims->nu;
    int nz      = dims->nz;
    int nx1     = dims->nx1;
    int nz1     = dims->nz1;
    int n_out   = dims->n_out;
    int ny      = dims->ny;
    int nuhat   = dims->nuhat;
    int nx2     = nx - nx1;
    int nz2     = nz - nz1;

    int num_stages = opts->ns;
    int num_steps  = opts->num_steps;
    int newton_iter = opts->newton_iter;

    int nvv = num_stages * n_out;
    int nyy = num_stages * ny;
    int nK1 = num_stages * nx1;
    int nK2 = num_stages * (nx2 + nz2);
    int nZ1 = num_stages * nz1;

    int nxz2 = nx2 + nz2;

    // assert - only use supported features
    if (mem->dt != in->T / opts->num_steps)
    {
        printf("ERROR sim_gnsf: mem->dt n!= in->T/opts->num_steps, check initialization\n");
        exit(1);
    }

    // assign variables from workspace
    double *Z_work = workspace->Z_work;

    struct blasfeo_dmat *J_r_vv =
        &workspace->J_r_vv;  // store the the jacobian of the residual w.r.t. vv
    int *ipiv = workspace->ipiv;
    struct blasfeo_dmat *J_r_x1u = &workspace->J_r_x1u;  // needed for sensitivity propagation

    struct blasfeo_dmat *dK1_dx1 = &workspace->dK1_dx1;
    struct blasfeo_dmat *dK1_du = &workspace->dK1_du;
    struct blasfeo_dmat *dZ_dx1 = &workspace->dZ_dx1;
    struct blasfeo_dmat *dZ_du = &workspace->dZ_du;
    struct blasfeo_dmat *J_G2_K1 = &workspace->J_G2_K1;
    struct blasfeo_dmat *dK2_dx1 = &workspace->dK2_dx1;
    struct blasfeo_dmat *dK2_dvv = &workspace->dK2_dvv;
    struct blasfeo_dmat *dxf_dwn = &workspace->dxf_dwn;
    struct blasfeo_dmat *S_forw_new = &workspace->S_forw_new;  // used to avoid side effects
    struct blasfeo_dmat *S_forw = &mem->S_forw;
    struct blasfeo_dmat *S_algebraic = &mem->S_algebraic; // to store
    struct blasfeo_dmat *S_algebraic_aux = &workspace->S_algebraic_aux; // to perform permutation

    struct blasfeo_dmat *f_LO_jac_traj = workspace->f_LO_jac_traj;

    struct blasfeo_dvec *vv_traj = workspace->vv_traj;
    struct blasfeo_dvec *yy_traj = workspace->yy_traj;

    struct blasfeo_dvec *f_LO_val = &workspace->f_LO_val;
    struct blasfeo_dvec *K1_val = &workspace->K1_val;
    struct blasfeo_dvec *x1_stage_val = &workspace->x1_stage_val;
    struct blasfeo_dvec *Z1_val = &workspace->Z1_val;

    struct blasfeo_dvec *yyu = &workspace->yyu;
    struct blasfeo_dvec *yyss = &workspace->yyss;
    struct blasfeo_dmat *dPHI_dyuhat = &workspace->dPHI_dyuhat;

    struct blasfeo_dvec *K2_val = &workspace->K2_val;
    struct blasfeo_dvec *x0_traj = &workspace->x0_traj;
    struct blasfeo_dvec *res_val = &workspace->res_val;
    struct blasfeo_dvec *u0 = &workspace->u0;
    struct blasfeo_dvec *lambda = &workspace->lambda;
    struct blasfeo_dvec *lambda_old = &workspace->lambda_old;

    struct blasfeo_dmat *dPsi_dvv = &workspace->dPsi_dvv;
    struct blasfeo_dmat *dPsi_dx = &workspace->dPsi_dx;
    struct blasfeo_dmat *dPsi_du = &workspace->dPsi_du;

    struct blasfeo_dvec *K1u = &workspace->K1u;
    struct blasfeo_dvec *Zu = &workspace->Zu;
    struct blasfeo_dvec *ALOtimesx02 = &workspace->ALOtimesx02;
    struct blasfeo_dvec *BLOtimesu0 = &workspace->BLOtimesu0;
    struct blasfeo_dvec *uhat = &workspace->uhat;
    struct blasfeo_dvec *z0 = &workspace->z0;

    int *ipiv_x = model->ipiv_x;
    int *ipiv_z = model->ipiv_z;

    // memory only available if (opts->sens_algebraic)
    // struct blasfeo_dvec *y_one_stage = &workspace->y_one_stage;
    // struct blasfeo_dmat *dz10_dx1u = &workspace->dz10_dx1u;
    // struct blasfeo_dmat *f_LO_jac0 = &workspace->f_LO_jac0;
    // struct blasfeo_dmat *sens_z2_rhs = &workspace->sens_z2_rhs;
    // struct blasfeo_dvec *x0dot_1 = &workspace->x0dot_1;

    // memory - precomputed matrices
    double *A_dt = mem->A_dt;
    double *b_dt = mem->b_dt;

    int *ipivM2 = mem->ipivM2;

    struct blasfeo_dmat *KKv = &mem->KKv;
    struct blasfeo_dmat *KKx = &mem->KKx;
    struct blasfeo_dmat *KKu = &mem->KKu;

    struct blasfeo_dmat *YYv = &mem->YYv;
    struct blasfeo_dmat *YYx = &mem->YYx;
    struct blasfeo_dmat *YYu = &mem->YYu;

    struct blasfeo_dmat *ZZv = &mem->ZZv;
    struct blasfeo_dmat *ZZx = &mem->ZZx;
    struct blasfeo_dmat *ZZu = &mem->ZZu;

    struct blasfeo_dmat *ALO = &mem->ALO;
    struct blasfeo_dmat *BLO = &mem->BLO;
    struct blasfeo_dmat *M2_LU = &mem->M2_LU;
    struct blasfeo_dmat *dK2_dx2 = &mem->dK2_dx2;
    struct blasfeo_dmat *dK2_du = &mem->dK2_du;
    struct blasfeo_dmat *dx2f_dx2u = &mem->dx2f_dx2u;
    struct blasfeo_dmat *Lu = &mem->Lu;

    // precomputed vectors
    struct blasfeo_dvec *KK0 = &mem->KK0;
    struct blasfeo_dvec *YY0 = &mem->YY0;
    struct blasfeo_dvec *ZZ0 = &mem->ZZ0;

    // ONLY available for algebraic sensitivity propagation
    // struct blasfeo_dmat *Z0x = mem->Z0x;
    // struct blasfeo_dmat *Z0u = mem->Z0u;
    // struct blasfeo_dmat *Z0v = mem->Z0v;

    // struct blasfeo_dmat *Y0x = mem->Y0x;
    // struct blasfeo_dmat *Y0u = mem->Y0u;
    // struct blasfeo_dmat *Y0v = mem->Y0v;

    // struct blasfeo_dmat *K0x = mem->K0x;
    // struct blasfeo_dmat *K0u = mem->K0u;
    // struct blasfeo_dmat *K0v = mem->K0v;

    // struct blasfeo_dmat *Lx = mem->Lx;
    // struct blasfeo_dmat *Lxdot = mem->Lxdot;
    // struct blasfeo_dmat *Lz = mem->Lz;

    // struct blasfeo_dmat *ELO_LU = mem->ELO_LU;

    // int *ipiv_ELO = mem->ipiv_ELO;
    // int *ipiv_vv0 = workspace->ipiv_vv0;

    // struct blasfeo_dmat *dr0_dvv0 = &workspace->dr0_dvv0;

    // transform inputs to blasfeo and apply permutation ipiv_x
    blasfeo_pack_dvec(nu, in->u, 1, u0, 0);
    blasfeo_pack_dvec(nx, &in->x[0], 1, x0_traj, 0);
    blasfeo_dvecpe(nx, ipiv_x, x0_traj, 0);
    blasfeo_pack_dvec(nx + nu, &in->S_adj[0], 1, lambda, 0);
    blasfeo_pack_dvec(nx + nu, &in->S_adj[0], 1, lambda_old, 0);
    blasfeo_dvecpe(nx, ipiv_x, lambda, 0);
    blasfeo_dvecpe(nx, ipiv_x, lambda_old, 0);


    if (model->fully_linear && !mem->first_call)
    {
        // xf = x_0 + S_forw_x * x0 + S_forw_u * u0; 
        blasfeo_dgemv_n(nx, nx, 1.0, S_forw, 0, 0, x0_traj, 0, 0.0,
                        x0_traj, 0, x0_traj, nx * num_steps);
        blasfeo_dgemv_n(nx, nu, 1.0, S_forw, 0, nx, u0, 0, 1.0,
                        x0_traj, nx * num_steps, x0_traj, nx * num_steps);

        // sensitivities are constant and already computed - only multiply with seed
        // dont change S_forw!
        if (!in->identity_seed)
        {
            // pack seed into dxf_dwn, permute
            blasfeo_pack_dmat(nx, nx + nu, &in->S_forw[0], nx, dxf_dwn, 0, 0);
            blasfeo_drowpe(nx, ipiv_x, dxf_dwn);
            blasfeo_dcolpe(nx, ipiv_x, dxf_dwn);


            // S_forw_new_x = S_forw_x * dxf_dx
            blasfeo_dgemm_nn(nx, nx, nx, 1.0, S_forw, 0, 0, dxf_dwn, 0, 0, 0.0, S_forw_new, 0, 0,
                                    S_forw_new, 0, 0);
            // S_forw_new_u = S_forw_x * S_forw_new_u + dxf_du
            blasfeo_dgemm_nn(nx, nu, nx, 1.0, S_forw, 0, 0, dxf_dwn, 0, nx, 1.0, S_forw_new, 0, nx,
                                S_forw_new, 0, nx);
        }
        else
        {
            // copy into S_forw_new and unpack from there later
            blasfeo_dgecp(nx, nx + nu, S_forw, 0, 0, S_forw_new, 0, 0);
        }
        

        // adjoint
        blasfeo_dgemv_t(nx, nx+nu, 1.0, S_forw, 0, 0, lambda_old, 0, 0.0, lambda_old, 0, lambda, 0);

        // z0 = S_algebraic_x * x0
        blasfeo_dgemv_n(nz, nx, 1.0, S_algebraic, 0, 0, x0_traj, 0, 0.0, z0, 0, z0, 0);
        // z0 += S_algebraic_u * u;
        blasfeo_dgemv_n(nz, nu, 1.0, S_algebraic, 0, nx, u0, 0, 1.0, z0, 0 , z0, 0);

        // S_algebraic available already..
    }
    else
    {
        // pack seed into S_forw, permute
        blasfeo_pack_dmat(nx, nx + nu, &in->S_forw[0], nx, S_forw, 0, 0);
        blasfeo_drowpe(nx, ipiv_x, S_forw);
        blasfeo_dcolpe(nx, ipiv_x, S_forw);

        // initialize vv for first step, for further steps initialize with last vv value in step loop
        for (int i = 0; i < num_stages; i++)
            blasfeo_pack_dvec(n_out, mem->phi_guess, 1, &vv_traj[0], i * n_out);

        // printf("GNSF_IRK initial guess:\n");
        // blasfeo_print_dvec(num_stages * n_out, &vv_traj[0], 0);
        // exit(1);

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
        phi_in[1] = uhat;

        // set output for phi_fun
        phi_fun_type_out[0] = BLASFEO_DVEC_ARGS;
        phi_fun_out[0] = &phi_fun_val_arg;
        phi_fun_val_arg.x = res_val;

        // set output for phi_fun_jac_y
        phi_fun_jac_y_type_out[0] = BLASFEO_DVEC_ARGS;
        phi_fun_jac_y_out[0] = &phi_fun_val_arg;

        phi_fun_jac_y_type_out[1] = BLASFEO_DMAT_ARGS;  // jac_y
        phi_jac_y_arg.A = dPHI_dyuhat;
        phi_jac_y_arg.aj = 0;  // never changes
        phi_fun_jac_y_out[1] = &phi_jac_y_arg;

        // set output for phi_jac_yuhat // jac_y
        phi_jac_yuhat_type_out[0] = BLASFEO_DMAT_ARGS;
        phi_jac_yuhat_out[0] = &phi_jac_y_arg;
        // jac_uhat
        phi_jac_yuhat_type_out[1] = BLASFEO_DMAT_ARGS;
        phi_jac_uhat_arg.A = dPHI_dyuhat;
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
        struct blasfeo_dvec_args f_lo_in_z1;
        f_lo_fun_in[2] = &f_lo_in_z1;
        // f_lo_in[3]: u;
        f_lo_fun_in[3] = u0;

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
        if (nx1 > 0 || nz1 > 0)
        {
            blasfeo_dgemv_n(nyy, nu, 1.0, YYu, 0, 0, u0, 0, 1.0, YY0, 0, yyu, 0);
            blasfeo_dgemv_n(nK1, nu, 1.0, KKu, 0, 0, u0, 0, 1.0, KK0, 0, K1u, 0);
            blasfeo_dgemv_n(nZ1, nu, 1.0, ZZu, 0, 0, u0, 0, 1.0, ZZ0, 0, Zu , 0);
        }

        // uhat = Lu * u0;
        blasfeo_dgemv_n(nuhat, nu, 1.0, Lu, 0, 0, u0, 0, 0.0, uhat, 0, uhat, 0);
        // BLOtimesu0 = BLO * u0 + c_LO
        blasfeo_pack_dvec(nxz2, model->c_LO, 1, BLOtimesu0, 0);
        blasfeo_dgemv_n(nxz2, nu, 1.0, BLO, 0, 0, u0, 0, 1.0, BLOtimesu0, 0, BLOtimesu0, 0);

        /************************************************
         * FORWARD LOOP
         ************************************************/
        // printf("GNSF: nx %d, nz %d, nu %d, nx1 %d, nz1 %d\n", nx, nz, nu, nx1, nz1);

        for (int ss = 0; ss < num_steps; ss++)
        {
            // STEP LOOP
            // initialize lifted variables vv with solution of previous step
            if (ss > 0)
                blasfeo_dveccp(nvv, &vv_traj[ss-1], 0, &vv_traj[ss], 0);

            if (nx1 > 0 || nz1 > 0)
            {
                // yyss = YY0 + YYu * u + YYx * x0(at_stage)
                blasfeo_dgemv_n(nyy, nx1, 1.0, YYx, 0, 0, x0_traj, ss*nx, 1.0, yyu, 0, yyss, nyy*ss);

                y_in.x = &yy_traj[ss];
                for (int iter = 0; iter < newton_iter; iter++)
                {  // NEWTON-ITERATION
                    /* EVALUATE RESIDUAL FUNCTION & JACOBIAN */

                    blasfeo_dgemv_n(nyy, nvv, 1.0, YYv, 0, 0, &vv_traj[ss], 0, 1.0, yyss, nyy * ss,
                                    &yy_traj[ss], 0);
                    // printf("yy =  \n");
                    // blasfeo_print_exp_dvec(nyy, &yy_traj[ss], 0);
                    if ((opts->jac_reuse && (ss == 0) && (iter == 0)) || (!opts->jac_reuse))
                    {
                        // set J_r_vv to unit matrix
                        blasfeo_dgese(nvv, nvv, 0.0, J_r_vv, 0, 0);
                        for (int ii = 0; ii < nvv; ii++)
                        {
                            blasfeo_dgein1(1.0, J_r_vv, ii, ii);
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

                            // build jacobian J_r_vv
                            blasfeo_dgemm_nn(n_out, nvv, ny, -1.0, dPHI_dyuhat, ii * n_out, 0, YYv,
                                            ii*ny, 0, 1.0, J_r_vv, ii*n_out, 0, J_r_vv, ii*n_out, 0);
                        }
                        else
                        {
                            acados_tic(&casadi_timer);
                            model->phi_fun->evaluate(model->phi_fun, phi_type_in, phi_in, phi_fun_type_out,
                                                    phi_fun_out);
                            out->info->ADtime += acados_toc(&casadi_timer);
                        }
                        // printf("\ngnsf: phi residual for newton %d, stage %d\n", iter, ii);
                        // blasfeo_print_dvec(n_out, res_val, ii*n_out);
                    }

                    // set up res_val
                    blasfeo_dveccpsc(nvv, -1.0, res_val, 0, res_val,
                                    0);  // set res_val = - res_val; NOW: res_val = -PHI_val;
                    blasfeo_dvecad(nvv, 1.0, &vv_traj[ss], 0, res_val, 0);
                            // set res_val = res_val + vv_traj;
                            // this is the actual value of the residual function!
                    acados_tic(&la_timer);
                    // factorize J_r_vv
                    if ((opts->jac_reuse && (ss == 0) & (iter == 0)) || (!opts->jac_reuse))
                    {
                        blasfeo_dgetrf_rp(nvv, nvv, J_r_vv, 0, 0, J_r_vv, 0, 0, ipiv);
                    }

                    /* Solve linear system and update vv */
                    blasfeo_dvecpe(nvv, ipiv, res_val, 0);  // permute r.h.s.
                    blasfeo_dtrsv_lnu(nvv, J_r_vv, 0, 0, res_val, 0, res_val, 0);
                    blasfeo_dtrsv_unn(nvv, J_r_vv, 0, 0, res_val, 0, res_val, 0);
                    out->info->LAtime += acados_toc(&la_timer);

                    blasfeo_daxpy(nvv, -1.0, res_val, 0, &vv_traj[ss], 0, &vv_traj[ss], 0);

                }  // END NEWTON-ITERATION

                // compute K1 and Z values
                blasfeo_dgemv_n(nK1, nvv, 1.0, KKv, 0, 0, &vv_traj[ss], 0, 1.0, K1u, 0, K1_val,
                                0);  // K1u contains KKu * u0 + KK0;
                blasfeo_dgemv_n(nK1, nx1, 1.0, KKx, 0, 0, x0_traj, ss * nx, 1.0, K1_val, 0,
                                K1_val, 0);
                if (nz1)
                {
                    blasfeo_dgemv_n(nZ1, nvv, 1.0, ZZv, 0, 0, &vv_traj[ss], 0, 1.0, Zu, 0, Z1_val,
                                    0);  // Zu contains ZZu * u0 + ZZ0;
                    blasfeo_dgemv_n(nZ1, nx1, 1.0, ZZx, 0, 0, x0_traj, ss * nx, 1.0, Z1_val, 0,
                                    Z1_val, 0);
                }
                // build x1 stage values
                for (int ii = 0; ii < num_stages; ii++)
                {
                    blasfeo_dveccp(nx1, x0_traj, ss * nx, x1_stage_val, nx1 * ii);
                    for (int jj = 0; jj < num_stages; jj++)
                    {
                        blasfeo_daxpy(nx1, A_dt[ii + num_stages * jj], K1_val, nx1 * jj, x1_stage_val,
                                    nx1 * ii, x1_stage_val, nx1 * ii);
                    }
                }
            } // if (nx1 > 0 || nz1 > 0)

            if (nxz2)
            {

                /* SIMULATE LINEAR OUTPUT SYSTEM */
                blasfeo_dgemv_n(nxz2, nx2, 1.0, ALO, 0, 0, x0_traj, ss * nx + nx1, 1.0, BLOtimesu0,
                                0, ALOtimesx02, 0);
                        // ALOtimesx02 = ALO * x2part(x0_traj(step)) + BLO * u0 + c_LO;

                if (model->nontrivial_f_LO) // Evaluate f_LO + jacobians
                {
                    f_lo_in_x1.x = x1_stage_val;
                    f_lo_in_k1.x = K1_val;
                    f_lo_in_z1.x = Z1_val;

                    f_lo_val_out.x = f_LO_val;
                    f_lo_jac_out.A = &f_LO_jac_traj[ss];
                    for (int ii = 0; ii < num_stages; ii++)
                    {
                        f_lo_in_x1.xi = ii * nx1;
                        f_lo_in_k1.xi = ii * nx1;
                        f_lo_in_z1.xi = ii * nz1;

                        f_lo_val_out.xi = ii * nxz2;
                        f_lo_jac_out.ai = ii * nxz2;

                        acados_tic(&casadi_timer);
                        model->f_lo_fun_jac_x1_x1dot_u_z->evaluate(model->f_lo_fun_jac_x1_x1dot_u_z,
                                                                f_lo_fun_type_in, f_lo_fun_in,
                                                                f_lo_fun_type_out, f_lo_fun_out);
                        out->info->ADtime += acados_toc(&casadi_timer);
                        blasfeo_dvecsc(nxz2, -1.0, f_LO_val, nxz2 * ii);  // f_LO_val = - f_LO_val
                        blasfeo_dvecad(nxz2, -1.0, ALOtimesx02, 0, f_LO_val, nxz2 * ii);
                        // f_LO_val = - ALOtimesx02 (includes BLO * u + c_LO) (actual rhs)
                    }
                }
                else
                    for (int ii = 0; ii < num_stages; ii++)
                    {
                        // f_LO_val = - ALOtimesx02 (includes BLO * u + c_LO) (actual rhs)
                        blasfeo_dveccpsc(nxz2, -1.0, ALOtimesx02, 0, f_LO_val, nxz2 * ii);
                    }

                acados_tic(&la_timer);
                /* Solve linear system and update K2 */
                blasfeo_dvecpe(nK2, ipivM2, f_LO_val, 0);  // permute r.h.s.
                blasfeo_dtrsv_lnu(nK2, M2_LU, 0, 0, f_LO_val, 0, f_LO_val, 0);
                blasfeo_dtrsv_unn(nK2, M2_LU, 0, 0, f_LO_val, 0, f_LO_val, 0);
                blasfeo_dveccpsc(nK2, -1.0, f_LO_val, 0, K2_val, 0);
                out->info->LAtime += acados_toc(&la_timer);
            }

            /* Get simulation result */
            blasfeo_dveccp(nx, x0_traj, nx * ss, x0_traj, nx * (ss + 1));
            for (int ii = 0; ii < num_stages; ii++)
            {
                blasfeo_daxpy(nx1, b_dt[ii], K1_val, ii * nx1, x0_traj, nx * (ss + 1), x0_traj,
                            nx * (ss + 1));
                blasfeo_daxpy(nx2, b_dt[ii], K2_val, ii * nxz2, x0_traj, nx1 + nx * (ss + 1),
                            x0_traj, nx1 + nx * (ss + 1));
            }

            // Forward Sensitivities (via IND)
            if (opts->sens_forw)
            {
                if (nx1 > 0 || nz1 > 0)
                {
                    // evaluate jacobian of residual function
                    // update yy
                    blasfeo_dgemv_n(nyy, nvv, 1.0, YYv, 0, 0, &vv_traj[ss], 0, 1.0, yyss, nyy * ss,
                                    &yy_traj[ss], 0);
                    // set J_r_vv to unit matrix
                    blasfeo_dgese(nvv, nvv, 0.0, J_r_vv, 0, 0);
                    for (int ii = 0; ii < nvv; ii++)
                    {
                        blasfeo_dgein1(1.0, J_r_vv, ii, ii);
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

                        // build J_r_vv
                        blasfeo_dgemm_nn(n_out, nvv, ny, -1.0, dPHI_dyuhat, ii * n_out, 0, YYv, ii * ny,
                                        0, 1.0, J_r_vv, ii * n_out, 0, J_r_vv, ii * n_out, 0);
                        // build J_r_x1u
                        blasfeo_dgemm_nn(n_out, nx1, ny, -1.0, dPHI_dyuhat, ii * n_out, 0, YYx, ii * ny,
                                        0, 0.0, J_r_x1u, ii * n_out, 0, J_r_x1u, ii * n_out,
                                        0);  // w.r.t. x1
                        blasfeo_dgemm_nn(n_out, nu, ny, -1.0, dPHI_dyuhat, ii * n_out, 0, YYu, ii * ny, 0,
                                        0.0, J_r_x1u, ii * n_out, nx1, J_r_x1u, ii * n_out,
                                        nx1);  // w.r.t. u
                        blasfeo_dgemm_nn(n_out, nu, nuhat, -1.0, dPHI_dyuhat, ii * n_out, ny, Lu, 0, 0,
                                        1.0, J_r_x1u, ii * n_out, nx1, J_r_x1u, ii * n_out,
                                        nx1);  // + dPhi_duhat * L_u;
                    }
                    acados_tic(&la_timer);
                    blasfeo_dgetrf_rp(nvv, nvv, J_r_vv, 0, 0, J_r_vv, 0, 0,
                                            ipiv);        // factorize J_r_vv
                    // printf("dPHI_dyuhat = (forward, ss = %d) \n", ss);
                    // blasfeo_print_exp_dmat(nvv, ny+nuhat, dPHI_dyuhat, 0, 0);

                    blasfeo_drowpe(nvv, ipiv, J_r_x1u);  // permute also rhs
                    blasfeo_dtrsm_llnu(nvv, nx1 + nu, 1.0, J_r_vv, 0, 0, J_r_x1u, 0, 0, J_r_x1u, 0, 0);
                    blasfeo_dtrsm_lunn(nvv, nx1 + nu, 1.0, J_r_vv, 0, 0, J_r_x1u, 0, 0, J_r_x1u, 0, 0);
                    out->info->LAtime += acados_toc(&la_timer);

                    blasfeo_dgemm_nn(nK1, nx1, nvv, -1.0, KKv, 0, 0, J_r_x1u, 0, 0, 1.0, KKx, 0, 0,
                                    dK1_dx1, 0, 0);
                    blasfeo_dgemm_nn(nK1, nu, nvv, -1.0, KKv, 0, 0, J_r_x1u, 0, nx1, 1.0, KKu, 0, 0,
                                    dK1_du, 0, 0);  // Blasfeo HP & Reference divver here
                    blasfeo_dgemm_nn(nZ1, nx1, nvv, -1.0, ZZv, 0, 0, J_r_x1u, 0, 0, 1.0, ZZx, 0, 0,
                                    dZ_dx1, 0, 0);
                    blasfeo_dgemm_nn(nZ1, nu, nvv, -1.0, ZZv, 0, 0, J_r_x1u, 0, nx1, 1.0, ZZu, 0, 0,
                                    dZ_du, 0, 0);
                }

                if (nxz2)
                {
                    // BUILD J_G2_wn, J_G2_K1 and right hand sides for linear systems dK2_dx1, dK2_du
                    if (model->nontrivial_f_LO)
                    {
                        for (int ii = 0; ii < num_stages; ii++)
                        {
                            for (int jj = 0; jj < num_stages; jj++)
                            {
                                blasfeo_dgecpsc(nxz2, nx1, -A_dt[ii + jj * num_stages], &f_LO_jac_traj[ss],
                                                ii * nxz2, 0, J_G2_K1, ii * nxz2, jj * nx1);
                                    //  copy - T * a_ij * df_dx1(k1_i, x1_i, z_i, u) into (i,j)th block
                            }
                            blasfeo_dgead(nxz2, nx1, -1.0, &f_LO_jac_traj[ss], ii * nxz2, nx1,
                                        J_G2_K1, ii * nxz2, ii * nx1);
                                    //  add - df_dx1dot(k1_i, x1_i, z_i, u) to (i,i)th block
                            blasfeo_dgemm_nn(nxz2, nx1, nz1, 1.0, &f_LO_jac_traj[ss], ii * nxz2,
                                        2 * nx1 + nu, dZ_dx1, ii * nz1, 0, 0.0,
                                        dK2_dx1, ii * nxz2, 0, dK2_dx1, ii * nxz2, 0);
                                    // set ith block of dK2_dx1 to df_dz(k1_i, x1_i, z_i, u) * dzi_dx1
                            blasfeo_dgemm_nn(nxz2, nu, nz1, 1.0, &f_LO_jac_traj[ss], ii * nxz2,
                                        2 * nx1 + nu, dZ_du, ii * nz1, 0, 1.0, BLO,
                                        0, 0, dK2_du, ii * nxz2,  0);
                                    // set ith block of dK2_du to df_dz(k1_i, x1_i, z_i, u) * dzi_du +  BLO
                        }
                        // dK2_dx1
                        blasfeo_dgemm_nn(nK2, nx1, nK1, -1.0, J_G2_K1, 0, 0, dK1_dx1, 0, 0, 1.0,
                                        dK2_dx1, 0, 0, dK2_dx1, 0, 0);
                            // dK2_dx1 = dK2_dx1 - JG2_K1 * dK1_dx1
                        blasfeo_dgead(nK2, nx1, 1.0, &f_LO_jac_traj[ss], 0, 0, dK2_dx1, 0, 0);
                            // add df_dx1 (at stages concatenated) to dK2_dx1

                        // solve for dK2_dx1
                        acados_tic(&la_timer);
                        blasfeo_drowpe(nK2, ipivM2, dK2_dx1);  // permute also rhs
                        blasfeo_dtrsm_llnu(nK2, nx1, 1.0, M2_LU, 0, 0,
                                        dK2_dx1, 0, 0, dK2_dx1, 0, 0);
                        blasfeo_dtrsm_lunn(nK2, nx1, 1.0, M2_LU, 0, 0, dK2_dx1, 0, 0, dK2_dx1, 0, 0);
                        out->info->LAtime += acados_toc(&la_timer);

                        // dK2_du
                        blasfeo_dgemm_nn(nK2, nu, nK1, -1.0, J_G2_K1, 0, 0, dK1_du, 0, 0, 1.0, dK2_du,
                                        0, 0, dK2_du, 0, 0);  // dK2_du += J_G2_K1 * dK1_du;

                        // BUILD dK2_du
                        blasfeo_dgead(nK2, nu, 1.0, &f_LO_jac_traj[ss], 0, 2 * nx1, dK2_du, 0, 0);
                            // add - df_du (at stages concatenated) to dK2_du
                        // solve for dK2_du
                        acados_tic(&la_timer);
                        blasfeo_drowpe(nK2, ipivM2, dK2_du);  // permute also rhs
                        blasfeo_dtrsm_llnu(nK2, nu, 1.0, M2_LU, 0, 0, dK2_du, 0, 0, dK2_du, 0, 0);
                        blasfeo_dtrsm_lunn(nK2, nu, 1.0, M2_LU, 0, 0, dK2_du, 0, 0, dK2_du, 0, 0);
                        out->info->LAtime += acados_toc(&la_timer);
                    }
                }
                // BUILD dxf_dwn
                blasfeo_dgese(nx, nx + nu, 0.0, dxf_dwn, 0, 0);  // Initialize as unit matrix
                for (int ii = 0; ii < nx; ii++)
                {
                    blasfeo_dgein1(1.0, dxf_dwn, ii, ii);
                }
                for (int ii = 0; ii < num_stages; ii++)
                {
                    // dK1_dw
                    blasfeo_dgead(nx1, nx1, b_dt[ii], dK1_dx1, ii * nx1, 0, dxf_dwn, 0, 0);
                    blasfeo_dgead(nx1, nu , b_dt[ii], dK1_du , ii * nx1, 0, dxf_dwn, 0, nx);
                }
                // dK2_dw
                if (model->nontrivial_f_LO)
                {
                    for (int ii = 0; ii < num_stages; ii++)
                    {
                        blasfeo_dgead(nx2, nx1, b_dt[ii], dK2_dx1, ii * nxz2, 0, dxf_dwn, nx1, 0);
                        blasfeo_dgead(nx2, nx2, b_dt[ii], dK2_dx2, ii * nxz2, 0, dxf_dwn, nx1, nx1);
                        blasfeo_dgead(nx2, nu,  b_dt[ii], dK2_du,  ii * nxz2, 0, dxf_dwn, nx1, nx);
                    }
                }
                else
                {
                    // copy from precomputed dx2f_dx2u
                    blasfeo_dgecp(nx2, nx2 + nu, dx2f_dx2u, 0, 0, dxf_dwn, nx1, nx1);
                }
                // omit matrix multiplication for identity seed
                if (in->identity_seed && ss == 0)
                {
                    blasfeo_dgecp(nx, nx + nu, dxf_dwn, 0, 0, S_forw, 0, 0);
                    blasfeo_dgecp(nx, nx + nu, dxf_dwn, 0, 0, S_forw_new, 0, 0);
                }
                else
                {
                    blasfeo_dgemm_nn(nx, nx, nx, 1.0, dxf_dwn, 0, 0, S_forw, 0, 0, 0.0, S_forw_new, 0, 0,
                                    S_forw_new, 0, 0);
                    blasfeo_dgemm_nn(nx, nu, nx, 1.0, dxf_dwn, 0, 0, S_forw, 0, nx, 1.0, dxf_dwn, 0, nx,
                                    S_forw_new, 0, nx);
                    blasfeo_dgecp(nx, nx + nu, S_forw_new, 0, 0, S_forw, 0, 0);
                }
            }

    /* output z and propagate corresponding sensitivities */
            if (ss == 0)
            {
            /* get output value for algebraic states z */
                if (opts->output_z || opts->sens_algebraic)
                {
                    for (int ii = 0; ii < nz1; ii++)  // ith component of z1
                    {
                        for (int jj = 0; jj < num_stages; jj++)
                        {
                            Z_work[jj] = blasfeo_dvecex1(Z1_val, nz1 * jj + ii);
                                    // copy values of z_ii in first step, into Z_work
                        }
                        neville_algorithm(0.0, num_stages - 1, mem->c_butcher, Z_work, &out->zn[ii]);
                                    // eval polynomial through (c_i, Z_i) at 0.
                    }
                    for (int ii = 0; ii < nz2; ii++)  // ith component of z2
                    {
                        for (int jj = 0; jj < num_stages; jj++)
                        {
                            Z_work[jj] = blasfeo_dvecex1(K2_val, nxz2 * jj + nx2 + ii);
                                    // copy values of z_ii in first step, into Z_work
                        }
                        neville_algorithm(0.0, num_stages-1, mem->c_butcher, Z_work, &out->zn[ii+nz1]);
                                    // eval polynomial through (c_i, Z_i) at 0.
                    }
                }
                // pack z0
                blasfeo_pack_dvec(nz, out->zn, 1, z0, 0);

                if (opts->sens_algebraic)
                {
                    if (!opts->sens_forw)
                    {
                        printf("\nsim_gnsf: algebraic sensitivities only supported with forward sensitivities\
                            \nplease set sens_forw to true");
                        exit(1);
                    }

                    // dz1_du
                    double interpolated_value;
                    for (int jj = 0; jj < nu; jj++)
                    {
                        for (int ii = 0; ii < nz1; ii++)
                        {
                            for (int kk = 0; kk < num_stages; kk++)
                            {
                                Z_work[kk] = blasfeo_dgeex1(dZ_du, kk*nz1+ii, jj);
                            }
                            neville_algorithm(0.0, num_stages - 1, opts->c_vec, Z_work, &interpolated_value);
                                        // eval polynomial through (c_kk, k_kk) at 0.
                            out->S_algebraic[ii + (jj+nx)*nz] = interpolated_value;
                            // printf("\ndz[ii=%d]_dxu[jj=%d] = %e\n", ii, jj, interpolated_value);
                            // blasfeo_pack_dvec(1, &interpolated_value, 1, xtdot, ii);
                        }
                    }
                    // dz1_dx1
                    for (int jj = 0; jj < nx1; jj++)
                    {
                        for (int ii = 0; ii < nz1; ii++)
                        {
                            for (int kk = 0; kk < num_stages; kk++)
                            {
                                Z_work[kk] = blasfeo_dgeex1(dZ_dx1, kk*nz1+ii, jj);
                            }
                            neville_algorithm(0.0, num_stages - 1, opts->c_vec, Z_work, &interpolated_value);
                                        // eval polynomial at 0.
                            out->S_algebraic[ii + jj*nz] = interpolated_value;
                        }
                    }
                    // dz1_dx2 = 0
                    for (int jj = 0; jj < nx2; jj++)
                    {
                        for (int ii = 0; ii < nz1; ii++)
                        {
                            out->S_algebraic[ii + (jj+nx1)*nz] = 0.0;
                        }
                    }
                    // dz2_dx1
                    for (int jj = 0; jj < nx1; jj++)
                    {
                        for (int ii = 0; ii < nz2; ii++)
                        {
                            for (int kk = 0; kk < num_stages; kk++)
                            {
                                Z_work[kk] = blasfeo_dgeex1(dK2_dx1, kk * nxz2 + nx2 + ii, jj);
                            }
                            neville_algorithm(0.0, num_stages - 1, opts->c_vec, Z_work, &interpolated_value);
                                        // eval polynomial at 0.
                            out->S_algebraic[ii+nz1+jj*nz] = interpolated_value;
                        }
                    }
                    // dz2_dx2
                    for (int jj = 0; jj < nx2; jj++)
                    {
                        for (int ii = 0; ii < nz2; ii++)
                        {
                            for (int kk = 0; kk < num_stages; kk++)
                            {
                                Z_work[kk] = blasfeo_dgeex1(dK2_dx2, kk*nxz2+nx2+ii, jj);
                            }
                            neville_algorithm(0.0, num_stages - 1, opts->c_vec, Z_work, &interpolated_value);
                                        // eval polynomial at 0.
                            out->S_algebraic[ii + nz1 + (jj+nx1)*nz] = interpolated_value;
                        }
                    }
                    // dz2_du
                    for (int jj = 0; jj < nu; jj++)
                    {
                        for (int ii = 0; ii < nz2; ii++)
                        {
                            for (int kk = 0; kk < num_stages; kk++)
                            {
                                Z_work[kk] = blasfeo_dgeex1(dK2_du, kk*nxz2+nx2+ii, jj);
                            }
                            neville_algorithm(0.0, num_stages - 1, opts->c_vec, Z_work, &interpolated_value);
                                        // eval polynomial at 0.
                            out->S_algebraic[ii + nz1 + (jj+nx)*nz] = interpolated_value;
                        }
                    }
                }
                blasfeo_pack_dmat(nz, nx + nu, out->S_algebraic, nz, S_algebraic, 0, 0);
                blasfeo_dgecp(nz, nx+nu, S_algebraic, 0, 0, S_algebraic_aux, 0, 0);

            // /* propagate sensitivities of z1 */
            //     if (opts->sens_algebraic)
            //     {
            //         for (int ii = 0; ii < nx1; ii++)  // ith component of x0dot_1
            //         {
            //             for (int jj = 0; jj < num_stages; jj++)
            //             {
            //                 Z_work[jj] = blasfeo_dvecex1(K1_val, nx1 * jj + ii);
            //                         // copy values of k1_ii in first step, into Z_work
            //             }
            //             neville_algorithm(0.0, num_stages - 1, mem->c_butcher, Z_work, &out->xn[ii]);
            //                         // eval polynomial through (c_i, k1_ii) at 0, write in out->xn
            //         }

            //         // pack x0dot_1
            //         blasfeo_pack_dvec(nx1, out->xn, 1, x0dot_1, 0);

            //         // evaluate phi at x0_1, x0_1dot, u, z0;
            //         // build y_0; // use y_one_stage
            //         // y_one_stage = Lxdot * x0dot_1
            //         blasfeo_dgemv_n(ny, nx1, 1.0, Lxdot, 0, 0, x0dot_1, 0, 0.0,
            //                         y_one_stage, 0, y_one_stage, 0);
            //         // y_one_stage += Lx * x0
            //         blasfeo_dgemv_n(ny, nx1, 1.0, Lx, 0, 0, x0_traj, 0, 1.0,
            //                         y_one_stage, 0, y_one_stage, 0);
            //         // y_one_stage += Lz * z0
            //         blasfeo_dgemv_n(ny, nz1, 1.0, Lz, 0, 0, z0, 0, 1.0,
            //                         y_one_stage, 0, y_one_stage, 0);
            //         // NOTE: alternatively y_one_stage could be computed as
            //         // yy0 = Y0x * x0_1 + Y0u * u + Y0v * f0,
            //         //       whereby f0 can be obtained by an interpolation formula.

            //         // set input for phi
            //         phi_type_in[0] = BLASFEO_DVEC;
            //         phi_in[0] = y_one_stage;

            //         // set output for phi; store result in first n_out rows of dPHI_dyuhat
            //         phi_jac_uhat_arg.ai = 0;
            //         phi_jac_y_arg.ai    = 0;

            //         acados_tic(&casadi_timer);
            //         model->phi_jac_y_uhat->evaluate(model->phi_jac_y_uhat, phi_type_in, phi_in,
            //                                         phi_jac_yuhat_type_out, phi_jac_yuhat_out);
            //         out->info->ADtime += acados_toc(&casadi_timer);

            //         /* set up dr0_dxn1u in first rows of J_r_x1u */
            //         // J_r_x1 = -dphi0_dy * Y0x
            //         blasfeo_dgemm_nn(n_out, nx1, ny, -1.0, dPHI_dyuhat, 0, 0, Y0x, 0, 0, 0.0,
            //                             J_r_x1u, 0, 0, J_r_x1u, 0, 0);  // w.r.t. x1
            //         // J_r_u = -dphi0_dy * Y0u
            //         blasfeo_dgemm_nn(n_out, nu, ny, -1.0, dPHI_dyuhat, 0, 0, Y0u, 0, 0, 0.0,
            //                             J_r_x1u, 0, nx1, J_r_x1u, 0, nx1);

            //         // J_r_u = J_r_u - dphi0_duhat * Lu
            //         blasfeo_dgemm_nn(n_out, nu, nuhat, -1.0, dPHI_dyuhat, 0, ny, Lu, 0, 0, 1.0,
            //                             J_r_x1u, 0, nx1, J_r_x1u, 0, nx1);  // w.r.t. u

            //         /* set up dr0_dvv0 = eye(n_out) - dphi_dy * Y0v; */
            //         blasfeo_dgese(n_out, n_out, 0.0, dr0_dvv0, 0, 0);  // set all to zero
            //         blasfeo_ddiare(n_out, 1.0, dr0_dvv0, 0, 0);  // dr0_dvv0 is unit now
            //         // dr0_dvv0 = dr0_dvv0 - dphi_dy * Y0v;
            //         blasfeo_dgemm_nn(n_out, n_out, ny, -1.0, dPHI_dyuhat, 0, 0, Y0v, 0, 0, 1.0,
            //                 dr0_dvv0, 0, 0, dr0_dvv0, 0, 0);

            //         /* solve dvv0_dxn1u = - dr0_dvv0 \ dr0_dxn1u */
            //         acados_tic(&la_timer);
            //         blasfeo_dgetrf_rp(n_out, n_out, dr0_dvv0, 0, 0, dr0_dvv0, 0, 0,
            //                                 ipiv_vv0);        // factorize dr0_dvv0

            //         blasfeo_drowpe(n_out, ipiv_vv0, J_r_x1u);  // permute also rhs
            //         blasfeo_dtrsm_llnu(n_out, nx1 + nu, 1.0, dr0_dvv0, 0, 0, J_r_x1u,
            //                             0, 0, J_r_x1u, 0, 0);
            //         blasfeo_dtrsm_lunn(n_out, nx1 + nu, 1.0, dr0_dvv0, 0, 0, J_r_x1u,
            //                             0, 0, J_r_x1u, 0, 0);
            //         out->info->LAtime += acados_toc(&la_timer);
            //         // J_r_x1u now contains dvv0_dxn1u in first rows

            //         // copy Z0x, Z0u into dz10_dx1u
            //         blasfeo_dgecp(nz1, nx1, Z0x, 0, 0, dz10_dx1u, 0, 0);
            //         blasfeo_dgecp(nz1, nu , Z0u, 0, 0, dz10_dx1u, 0, nx1);

            //         // add Z0v * dvv0_dxn1u
            //         blasfeo_dgemm_nn(nz1, nx1 + nu, n_out, -1.0, Z0v, 0, 0, J_r_x1u, 0, 0, 1.0,
            //                         dz10_dx1u, 0, 0, dz10_dx1u, 0, 0);

            //         // extract into out->S_algebraic
            //         blasfeo_unpack_dmat(nz1, nx1, dz10_dx1u, 0, 0, out->S_algebraic, nz);
            //         blasfeo_unpack_dmat(nz1, nu , dz10_dx1u, 0, nx1, &out->S_algebraic[nx*nz], nz);

            //         for (int ii = 0; ii < nx2 * nz1; ii++) {  //  dz_dx2_0 = 0
            //             out->S_algebraic[nx1 * nz1 + ii] = 0;
            //         }

            //         // reset input for phi
            //         phi_type_in[0] = BLASFEO_DVEC_ARGS;
            //         phi_in[0] = &y_in;
            //     }
            // /* propagate sensitivities of z2 */
            //     if (opts->sens_algebraic && nz2)
            //     {
            //         // eval f_lo
            //         f_lo_in_x1.x  = x0_traj;
            //         f_lo_in_x1.xi = 0;
            //         f_lo_in_k1.x  = x0dot_1;
            //         f_lo_in_k1.xi = 0;
            //         f_lo_in_z1.x  = z0;
            //         f_lo_in_z1.xi = 0;

            //         f_lo_val_out.xi = 0;
            //         f_lo_jac_out.ai = 0;
            //         f_lo_jac_out.A = f_LO_jac0;

            //         acados_tic(&casadi_timer);
            //         model->f_lo_fun_jac_x1_x1dot_u_z->evaluate(model->f_lo_fun_jac_x1_x1dot_u_z,
            //                                                 f_lo_fun_type_in, f_lo_fun_in,
            //                                                 f_lo_fun_type_out, f_lo_fun_out);
            //         out->info->ADtime += acados_toc(&casadi_timer);

            //         // reset f_lo in & output
            //         f_lo_in_x1.x = x1_stage_val;
            //         f_lo_in_k1.x = K1_val;
            //         f_lo_in_z1.x = Z1_val;
            //         f_lo_jac_out.A = &f_LO_jac_traj[ss];

            //         // get dk01_dx01u, note: J_r_x1u contains dvv0_dxn1u in first rows
            //         // use dK1_dx1, dK1_du (first rows)
            //         blasfeo_dgemm_nn(nx1, nx1, n_out, 1.0, K0v, 0, 0, J_r_x1u, 0, 0, 1.0,
            //                 K0x, 0, 0, dK1_dx1, 0, 0);
            //         blasfeo_dgemm_nn(nx1, nu, n_out, 1.0, K0v, 0, 0, J_r_x1u, 0, nx1, 1.0,
            //                 K0u, 0, 0, dK1_du, 0, 0);

            //         /* set up right hand side: sens_z2_rhs */
            //         // sens_z2_rhs = dg2_dk1 * dk1_dx1u + dg2_dx1u
            //         blasfeo_dgemm_nn(nxz2, nx1, nx1, 1.0, f_LO_jac0, nx1, 0, dK1_dx1,
            //                 0, 0, 1.0, f_LO_jac0, 0, 0, sens_z2_rhs, 0, 0);
            //         blasfeo_dgemm_nn(nxz2, nu, nx1, 1.0, f_LO_jac0, nx1, 0, dK1_du,
            //                 0, 0, 1.0, f_LO_jac0, 0, 2*nx1, sens_z2_rhs, 0, nx1);

            //         // sens_z2_rhs += dg2_dz1 * dz1_dx1u
            //         blasfeo_dgemm_nn(nxz2, nx1 + nu, nz1, 1.0, f_LO_jac0, 2*nx1+nu, 0, dz10_dx1u,
            //                 0, 0, 1.0, sens_z2_rhs, 0, 0, sens_z2_rhs, 0, 0);

            //         /* solve E_LO \ sens_z2_rhs */
            //         acados_tic(&la_timer);
            //         blasfeo_drowpe(nxz2, ipiv_ELO, sens_z2_rhs);  // permute also rhs
            //         blasfeo_dtrsm_llnu(nxz2, nx1 + nu, 1.0, ELO_LU, 0, 0, sens_z2_rhs,
            //                             0, 0, sens_z2_rhs, 0, 0);
            //         blasfeo_dtrsm_lunn(nxz2, nx1 + nu, 1.0, ELO_LU, 0, 0, sens_z2_rhs,
            //                             0, 0, sens_z2_rhs, 0, 0);
            //         out->info->LAtime += acados_toc(&la_timer);

            //         // extract into out->S_algebraic
            //         blasfeo_unpack_dmat(nz2, nx1, sens_z2_rhs, nx2, 0, &out->S_algebraic[nz1], nz);
            //         blasfeo_unpack_dmat(nz2, nu , sens_z2_rhs, nx2, nx1,
            //                             &out->S_algebraic[nx*nz + nz1], nz);
            //         blasfeo_unpack_dmat(nz2, nx2, mem->ELO_inv_ALO, nx2, 0,
            //                             &out->S_algebraic[nx1*nz + nz1], nz);
            //        }

            } // if ss == 0;
            if (ss == num_steps-1)
            {
                // store last vv values for next initialization
                blasfeo_unpack_dvec(n_out, &vv_traj[ss], (num_stages-1) * n_out, mem->phi_guess, 1);
            }
        }  // end step loop: ss



    /************************************************
     * ADJOINT SENSITIVITY PROPAGATION
     ************************************************/

        if (opts->sens_adj)
        {
            for (int ss = num_steps - 1; ss >= 0; ss--)
            {
                /*  SET UP Right Hand Sides for LINEAR SYSTEMS and J_G2_K1 */
                y_in.x = &yy_traj[ss];
                if (model->nontrivial_f_LO)
                {
                    for (int ii = 0; ii < num_stages; ii++)
                    {
                        blasfeo_dgemm_nn(nxz2, nvv, nz1, 1.0, &f_LO_jac_traj[ss], nxz2 * ii, 2 * nx1 + nu,
                                ZZv, ii * nz1, 0, 0.0, dK2_dvv, ii * nxz2, 0, dK2_dvv, ii * nxz2, 0);
                                // set dK2_dvv(i_th block) = df_dzi * ZZv(ith block);
                        blasfeo_dgemm_nn(nxz2, nx1, nz1, 1.0, &f_LO_jac_traj[ss], nxz2 * ii, 2 * nx1 + nu,
                                ZZx, ii * nz1, 0, 0.0, dK2_dx1, ii * nxz2, 0, dK2_dx1, ii * nxz2, 0);
                                // set dK2_dx1(i_th block) = df_dzi * ZZx(ith block);
                        blasfeo_dgemm_nn(nxz2, nu, nz1, 1.0, &f_LO_jac_traj[ss], nxz2 * ii, 2 * nx1 + nu,
                                ZZu, ii * nz1, 0, 1.0, BLO, 0, 0, dK2_du, ii * nxz2, 0);
                                // set dK2_du(i_th block) = df_dzi * ZZu(ith block) + BLO;

                        // set up J_G2_K1 (as in forward loop)
                        for (int jj = 0; jj < num_stages; jj++)
                        {
                            blasfeo_dgecpsc(nxz2, nx1, -A_dt[ii + jj * num_stages], &f_LO_jac_traj[ss],
                                        ii * nxz2, 0, J_G2_K1, ii * nxz2, jj * nx1);
                                //  copy - T * a_ij * df_dx1(k1_i, x1_i, z_i, u) into (i,j)th block
                        }
                        blasfeo_dgead(nxz2, nx1, -1.0, &f_LO_jac_traj[ss], nxz2 * ii, nx1,
                                J_G2_K1, nxz2 * ii, nx1 * ii);
                                //  add -df_dx1dot(k1_i, x1_i, z_i, u) to (i,i)th block
                    }
                    blasfeo_dgemm_nn(nK2, nvv, nK1, -1.0, J_G2_K1, 0, 0, KKv, 0, 0, 1.0, dK2_dvv, 0, 0,
                                    dK2_dvv, 0, 0);  // dK2_dvv += -J_G2_K1 * KKv
                    blasfeo_dgemm_nn(nK2, nx1, nK1, -1.0, J_G2_K1, 0, 0, KKx, 0, 0, 1.0, dK2_dx1, 0, 0,
                                    dK2_dx1, 0, 0);  // dK2_dx1 += -JG2_K1 * KKx
                    blasfeo_dgemm_nn(nK2, nu, nK1, -1.0, J_G2_K1, 0, 0, KKu, 0, 0, 1.0, dK2_du, 0, 0,
                                    dK2_du, 0, 0);  // dK2_du += -JG2_K1 * KKu

                    blasfeo_dgead(nK2, nx1, 1.0, &f_LO_jac_traj[ss], 0, 0, dK2_dx1, 0, 0);
                    blasfeo_dgead(nK2, nu,  1.0, &f_LO_jac_traj[ss], 0, 2 * nx1, dK2_du, 0, 0);

                    // solve dK2_du = M2 \ dK2_du
                    blasfeo_drowpe(nK2, ipivM2, dK2_du);  // permute also rhs
                    blasfeo_dtrsm_llnu(nK2, nu, 1.0, M2_LU, 0, 0, dK2_du, 0, 0, dK2_du, 0, 0);
                    blasfeo_dtrsm_lunn(nK2, nu, 1.0, M2_LU, 0, 0, dK2_du, 0, 0, dK2_du, 0, 0);
                    out->info->LAtime += acados_toc(&la_timer);
                }
                

                /*  SOLVE LINEAR SYSTEMS  */
                acados_tic(&la_timer);
                // solve dK2_dvv = M2 \ dK2_dvv
                blasfeo_drowpe(nK2, ipivM2, dK2_dvv);  // permute also rhs
                blasfeo_dtrsm_llnu(nK2, nvv, 1.0, M2_LU, 0, 0, dK2_dvv, 0, 0, dK2_dvv, 0, 0);
                blasfeo_dtrsm_lunn(nK2, nvv, 1.0, M2_LU, 0, 0, dK2_dvv, 0, 0, dK2_dvv, 0, 0);

                // solve dK2_dx1 = M2 \ dK2_dx1
                blasfeo_drowpe(nK2, ipivM2, dK2_dx1);  // permute also rhs
                blasfeo_dtrsm_llnu(nK2, nx1, 1.0, M2_LU, 0, 0, dK2_dx1, 0, 0, dK2_dx1, 0, 0);
                blasfeo_dtrsm_lunn(nK2, nx1, 1.0, M2_LU, 0, 0, dK2_dx1, 0, 0, dK2_dx1, 0, 0);


                blasfeo_dgese(nx, nvv, 0.0, dPsi_dvv, 0, 0);  // initialize dPsi_d..
                blasfeo_dgese(nx, nx, 0.0, dPsi_dx, 0, 0);
                blasfeo_dgese(nx, nu, 0.0, dPsi_du, 0, 0);
                blasfeo_ddiare(nx, 1.0, dPsi_dx, 0, 0);  // dPsi_dx is unit now

                // compute dPsi_d..
                for (int ii = 0; ii < num_stages; ii++)
                {
                    blasfeo_dgead(nx1, nvv, b_dt[ii], KKv, ii * nx1, 0, dPsi_dvv, 0, 0);
                    blasfeo_dgead(nx1, nx1, b_dt[ii], KKx, ii * nx1, 0, dPsi_dx, 0, 0);
                    blasfeo_dgead(nx1, nu, b_dt[ii], KKu, ii * nx1, 0, dPsi_du, 0, 0);

                    blasfeo_dgead(nx2, nvv, b_dt[ii], dK2_dvv, ii * nxz2, 0, dPsi_dvv, nx1, 0);
                    blasfeo_dgead(nx2, nx1, b_dt[ii], dK2_dx1, ii * nxz2, 0, dPsi_dx, nx1, 0);
                    blasfeo_dgead(nx2, nx2, b_dt[ii], dK2_dx2, ii * nxz2, 0, dPsi_dx, nx1, nx1);
                    blasfeo_dgead(nx2, nu, b_dt[ii], dK2_du, ii * nxz2, 0, dPsi_du, nx1, 0);
                }

                if (nx1 > 0 || nz1 > 0)
                {
                    // evaluate jacobian of residual function
                    // update yy
                    blasfeo_dgemv_n(nyy, nvv, 1.0, YYv, 0, 0, &vv_traj[ss], 0, 1.0, yyss, nyy * ss,
                                    &yy_traj[ss], 0);
                    // set J_r_vv to unit matrix
                    blasfeo_dgese(nvv, nvv, 0.0, J_r_vv, 0, 0);
                    for (int ii = 0; ii < nvv; ii++)
                    {
                        blasfeo_dgein1(1.0, J_r_vv, ii, ii);
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

                        // build J_r_vv
                        blasfeo_dgemm_nn(n_out, nvv, ny, -1.0, dPHI_dyuhat, ii * n_out, 0, YYv, ii * ny,
                                        0, 1.0, J_r_vv, ii * n_out, 0, J_r_vv, ii * n_out, 0);

                        // build J_r_x1u
                        blasfeo_dgemm_nn(n_out, nx1, ny, -1.0, dPHI_dyuhat, ii * n_out, 0, YYx, ii * ny,
                                        0, 0.0, J_r_x1u, ii * n_out, 0, J_r_x1u, ii * n_out,
                                        0);  // w.r.t. x1
                        blasfeo_dgemm_nn(n_out, nu, ny, -1.0, dPHI_dyuhat, ii * n_out, 0, YYu, ii * ny, 0,
                                        0.0, J_r_x1u, ii * n_out, nx1, J_r_x1u, ii * n_out,
                                        nx1);  // w.r.t. u
                        blasfeo_dgemm_nn(n_out, nu, nuhat, -1.0, dPHI_dyuhat, ii * n_out, ny, Lu, 0, 0,
                                        1.0, J_r_x1u, ii * n_out, nx1, J_r_x1u, ii * n_out,
                                        nx1);  // + dPhi_duhat * L_u;
                    }
                    acados_tic(&la_timer);
                    blasfeo_dgetrf_rp(nvv, nvv, J_r_vv, 0, 0, J_r_vv, 0, 0,
                                            ipiv);  // factorize J_r_vv
                    out->info->LAtime += acados_toc(&la_timer);

                    blasfeo_dgemv_t(nx, nvv, 1.0, dPsi_dvv, 0, 0, lambda, 0, 0.0, res_val, 0, res_val,
                                    0);  // use res_val to store lambda_vv
                    acados_tic(&la_timer);
                    blasfeo_dtrsv_utn(nvv, J_r_vv, 0, 0, res_val, 0, res_val, 0);
                    blasfeo_dtrsv_ltu(nvv, J_r_vv, 0, 0, res_val, 0, res_val, 0);
                    blasfeo_dvecpei(nvv, ipiv, res_val, 0);  // permute linear syst solution
                    out->info->LAtime += acados_toc(&la_timer);
                }

                blasfeo_dveccp(nx + nu, lambda, 0, lambda_old, 0);
                blasfeo_dgemv_t(nx, nu, 1.0, dPsi_du, 0, 0, lambda_old, 0, 1.0, lambda_old, nx,
                                lambda, nx);  // update lambda_u

                blasfeo_dgemv_t(nx, nx, 1.0, dPsi_dx, 0, 0, lambda_old, 0, 0.0, res_val, 0,
                            lambda, 0);
                blasfeo_dveccp(nx + nu, lambda, 0, lambda_old, 0);
                blasfeo_dgemv_t(nvv, nx1, -1.0, J_r_x1u, 0, 0, res_val, 0, 1.0, lambda_old, 0,
                                lambda, 0);
                blasfeo_dgemv_t(nvv, nu, -1.0, J_r_x1u, 0, nx1, res_val, 0, 1.0, lambda_old, nx,
                                lambda, nx);
            }
        }
    }
/* unpack */
    // printf("GNSF: x before permutation\n");
    // blasfeo_print_exp_dvec(nx, x0_traj, nx * num_steps);

    blasfeo_dvecpei(nx, ipiv_x, x0_traj, nx * num_steps);
    blasfeo_unpack_dvec(nx, x0_traj, nx * num_steps, out->xn, 1);

    if (opts->sens_forw)
    {
        blasfeo_drowpei(nx, ipiv_x, S_forw_new);
        blasfeo_dcolpei(nx, ipiv_x, S_forw_new);
        blasfeo_unpack_dmat(nx, nx + nu, S_forw_new, 0, 0, out->S_forw, nx);
    }
    if (opts->sens_adj)
    {
        blasfeo_dvecpei(nx, ipiv_x, lambda, 0);
        blasfeo_unpack_dvec(nx + nu, lambda, 0, out->S_adj, 1);
    }
    if (opts->sens_algebraic)
    {
        // permute rows and cols
        blasfeo_dgecp(nz, nx+nu, S_algebraic, 0, 0, S_algebraic_aux, 0, 0);
        blasfeo_drowpei(nz, ipiv_z, S_algebraic_aux);
        blasfeo_dcolpei(nx, ipiv_x, S_algebraic_aux);
        blasfeo_unpack_dmat(nz, nx+nu, S_algebraic_aux, 0, 0, out->S_algebraic, nz);
    }
    if (opts->output_z)
    {
        blasfeo_dvecpei(nz, ipiv_z, z0, 0);
        blasfeo_unpack_dvec(nz, z0, 0, out->zn, 1);
    }

    out->info->CPUtime = acados_toc(&tot_timer);

	mem->time_sim = out->info->CPUtime;
	mem->time_ad = out->info->ADtime;
	mem->time_la = out->info->LAtime;

    return ACADOS_SUCCESS;
}



void sim_gnsf_config_initialize_default(void *config_)
{
    sim_config *config = config_;
    config->evaluate = &sim_gnsf;
    config->precompute = &sim_gnsf_precompute;
    // opts
    config->opts_calculate_size = &sim_gnsf_opts_calculate_size;
    config->opts_assign = &sim_gnsf_opts_assign;
    config->opts_initialize_default = &sim_gnsf_opts_initialize_default;
    config->opts_update = &sim_gnsf_opts_update;
    config->opts_set = &sim_gnsf_opts_set;
    config->opts_get = &sim_gnsf_opts_get;
    // memory & workspace
    config->memory_calculate_size = &sim_gnsf_memory_calculate_size;
    config->memory_assign = &sim_gnsf_memory_assign;
    config->memory_set = &sim_gnsf_memory_set;
    config->memory_set_to_zero = &sim_gnsf_memory_set_to_zero;
    config->memory_get = &sim_gnsf_memory_get;
    config->workspace_calculate_size = &sim_gnsf_workspace_calculate_size;
    // model
    config->model_calculate_size = &sim_gnsf_model_calculate_size;
    config->model_assign = &sim_gnsf_model_assign;
    config->model_set = &sim_gnsf_model_set;
    // dims
    config->dims_calculate_size = &sim_gnsf_dims_calculate_size;
    config->dims_assign = &sim_gnsf_dims_assign;
    config->dims_set = &sim_gnsf_dims_set;
    config->dims_get = &sim_gnsf_dims_get;
    return;
}
