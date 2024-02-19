/*
 * Copyright (c) The acados authors.
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


#include "acados/ocp_nlp/ocp_nlp_cost_nls.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_nlp_cost_nls_dims_calculate_size(void *config_)
{
    acados_size_t size = sizeof(ocp_nlp_cost_nls_dims);

    return size;
}



void *ocp_nlp_cost_nls_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_nls_dims *dims = (ocp_nlp_cost_nls_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_nls_dims);

    assert((char *) raw_memory + ocp_nlp_cost_nls_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void ocp_nlp_cost_nls_dims_set(void *config_, void *dims_, const char *field, int* value)
{
    ocp_nlp_cost_nls_dims *dims = (ocp_nlp_cost_nls_dims *) dims_;
    if (!strcmp(field, "nx"))
    {
        dims->nx = *value;
    }
    else if (!strcmp(field, "nz"))
    {
        dims->nz = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu = *value;
    }
    else if (!strcmp(field, "ny"))
    {
        dims->ny = *value;
    }
    else if (!strcmp(field, "ns"))
    {
        dims->ns = *value;
    }
    else
    {
        printf("\nerror: dimension type: %s not available in module\n", field);
        exit(1);
    }
}



/* dimension getters */
static void ocp_nlp_cost_nls_get_ny(void *config_, void *dims_, int* value)
{
    ocp_nlp_cost_nls_dims *dims = (ocp_nlp_cost_nls_dims *) dims_;
    *value = dims->ny;
}



void ocp_nlp_cost_nls_dims_get(void *config_, void *dims_, const char *field, int* value)
{
    if (!strcmp(field, "ny"))
    {
        ocp_nlp_cost_nls_get_ny(config_, dims_, value);
    }
    else
    {
        printf("error: ocp_nlp_cost_nls_dims_get: attempt to get dimensions of non-existing field %s\n", field);
        exit(1);
    }
}



/************************************************
 * model
 ************************************************/

acados_size_t ocp_nlp_cost_nls_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_nls_model);

    size += 64;  // blasfeo_mem align

    size += 1 * blasfeo_memsize_dmat(ny, ny);  // W
    size += 1 * blasfeo_memsize_dvec(ny);      // y_ref
    size += 2 * blasfeo_memsize_dvec(2 * ns);  // Z, z

    return size;
}



void *ocp_nlp_cost_nls_model_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_cost_nls_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // int nx = dims->nx;
    // int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    // struct
    ocp_nlp_cost_nls_model *model = (ocp_nlp_cost_nls_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_nls_model);

    model->nls_y_fun = NULL;
    model->nls_y_fun_jac = NULL;
    model->nls_y_hess = NULL;

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dmat
    // W
    assign_and_advance_blasfeo_dmat_mem(ny, ny, &model->W, &c_ptr);

    // blasfeo_dvec
    // y_ref
    assign_and_advance_blasfeo_dvec_mem(ny, &model->y_ref, &c_ptr);
    blasfeo_dvecse(ny, 0.0, &model->y_ref, 0);

    // Z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->Z, &c_ptr);
    // z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->z, &c_ptr);

    // default initialization
    model->scaling = 1.0;
    model->t = 0.0;

    // initialize to 1 to update factorization of W in precompute
    model->W_changed = 1;

    // assert
    assert((char *) raw_memory + ocp_nlp_cost_nls_model_calculate_size(config_, dims) >= c_ptr);

    return model;
}



int ocp_nlp_cost_nls_model_set(void *config_, void *dims_, void *model_,
                                         const char *field, void *value_)
{
    int status = ACADOS_SUCCESS;

    if ( !config_ || !dims_ || !model_ || !value_ )
    {
        printf("ocp_nlp_cost_nls_model_set: got NULL pointer \n");
        exit(1);
    }

    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;

    // int nx = dims->nx;
    // int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    if (!strcmp(field, "W"))
    {
        double *W_col_maj = (double *) value_;
        blasfeo_pack_dmat(ny, ny, W_col_maj, ny, &model->W, 0, 0);
        model->W_changed = 1;
    }
    else if (!strcmp(field, "y_ref") || !strcmp(field, "yref"))
    {
        double *y_ref = (double *) value_;
        blasfeo_pack_dvec(ny, y_ref, 1, &model->y_ref, 0);
    }
    else if (!strcmp(field, "Z"))
    {
        double *Z = (double *) value_;
        blasfeo_pack_dvec(ns, Z, 1, &model->Z, 0);
        blasfeo_pack_dvec(ns, Z, 1, &model->Z, ns);
    }
    else if (!strcmp(field, "Zl"))
    {
        double *Zl = (double *) value_;
        blasfeo_pack_dvec(ns, Zl, 1, &model->Z, 0);
    }
    else if (!strcmp(field, "Zu"))
    {
        double *Zu = (double *) value_;
        blasfeo_pack_dvec(ns, Zu, 1, &model->Z, ns);
    }
    else if (!strcmp(field, "z"))
    {
        double *z = (double *) value_;
        blasfeo_pack_dvec(ns, z, 1, &model->z, 0);
        blasfeo_pack_dvec(ns, z, 1, &model->z, ns);
    }
    else if (!strcmp(field, "zl"))
    {
        double *zl = (double *) value_;
        blasfeo_pack_dvec(ns, zl, 1, &model->z, 0);
    }
    else if (!strcmp(field, "zu"))
    {
        double *zu = (double *) value_;
        blasfeo_pack_dvec(ns, zu, 1, &model->z, ns);
    }
    else if (!strcmp(field, "nls_y_fun") || !strcmp(field, "nls_res"))
    {
        model->nls_y_fun = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "nls_y_fun_jac") || !strcmp(field, "nls_res_jac"))
    {
        model->nls_y_fun_jac = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "nls_y_hess") || !strcmp(field, "nls_hess"))
    {
        model->nls_y_hess = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "scaling"))
    {
        double *scaling_ptr = (double *) value_;
        model->scaling = *scaling_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_cost_nls_model_set\n", field);
        exit(1);
    }
    return status;
}



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_cost_nls_opts_calculate_size(void *config_, void *dims_)
{
    // ocp_nlp_cost_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_nls_opts);
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_cost_nls_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    // ocp_nlp_cost_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_nls_opts *opts = (ocp_nlp_cost_nls_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_nls_opts);

    assert((char *) raw_memory + ocp_nlp_cost_nls_opts_calculate_size(config_, dims_) >= c_ptr);

    return opts;
}



void ocp_nlp_cost_nls_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_opts *opts = opts_;

    opts->gauss_newton_hess = 1;

    return;
}



void ocp_nlp_cost_nls_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    // ocp_nlp_cost_nls_opts *opts = opts_;

    return;
}



void ocp_nlp_cost_nls_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_opts *opts = opts_;

    if(!strcmp(field, "gauss_newton_hess"))
    {
        int *int_ptr = value;
        opts->gauss_newton_hess = *int_ptr;
    }
    else if(!strcmp(field, "exact_hess"))
    {
        int *int_ptr = value;
        if(*int_ptr==0)
        {
            opts->gauss_newton_hess = 1;
        }
        else
        {
            opts->gauss_newton_hess = 0;
        }
    }
    else if(!strcmp(field, "integrator_cost"))
    {
        int *opt_val = (int *) value;
        opts->integrator_cost = *opt_val;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_cost_nls_opts_set\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_cost_nls_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_dims *dims = dims_;
    // ocp_nlp_cost_nls_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_nls_memory);

    size += 1 * blasfeo_memsize_dmat(ny, ny);            // W_chol
    size += 1 * blasfeo_memsize_dmat(nu + nx, ny);       // Jt
    size += 1 * blasfeo_memsize_dvec(ny);                // res
    size += 1 * blasfeo_memsize_dvec(nu + nx + 2 * ns);  // grad

    size += 64;  // blasfeo_mem align

    return size;
}



void *ocp_nlp_cost_nls_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_dims *dims = dims_;
    // ocp_nlp_cost_nls_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    // struct
    ocp_nlp_cost_nls_memory *memory = (ocp_nlp_cost_nls_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_nls_memory);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // W_chol
    assign_and_advance_blasfeo_dmat_mem(ny, ny, &memory->W_chol, &c_ptr);
    // Jt
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ny, &memory->Jt, &c_ptr);
    // res
    assign_and_advance_blasfeo_dvec_mem(ny, &memory->res, &c_ptr);
    // grad
    assign_and_advance_blasfeo_dvec_mem(nu + nx + 2 * ns, &memory->grad, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_cost_nls_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}



double *ocp_nlp_cost_nls_memory_get_fun_ptr(void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    return &memory->fun;
}


struct blasfeo_dmat *ocp_nlp_cost_nls_memory_get_W_chol_ptr(void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    return &memory->W_chol;
}


struct blasfeo_dvec *ocp_nlp_cost_nls_model_get_y_ref_ptr(void *in_)
{
    ocp_nlp_cost_nls_model *model = in_;

    return &model->y_ref;
}

struct blasfeo_dvec *ocp_nlp_cost_nls_memory_get_grad_ptr(void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    return &memory->grad;
}



void ocp_nlp_cost_nls_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    memory->RSQrq = RSQrq;

    return;
}



void ocp_nlp_cost_nls_memory_set_Z_ptr(struct blasfeo_dvec *Z, void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    memory->Z = Z;

    return;
}



void ocp_nlp_cost_nls_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    memory->ux = ux;

    return;
}



void ocp_nlp_cost_nls_memory_set_z_alg_ptr(struct blasfeo_dvec *z_alg, void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    memory->z_alg = z_alg;
}



void ocp_nlp_cost_nls_memory_set_dzdux_tran_ptr(struct blasfeo_dmat *dzdux_tran, void *memory_)
{
    ocp_nlp_cost_nls_memory *memory = memory_;

    memory->dzdux_tran = dzdux_tran;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_cost_nls_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_dims *dims = dims_;
    // ocp_nlp_cost_nls_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_nls_workspace);

    size += 1 * blasfeo_memsize_dmat(nu + nx, ny);       // tmp_nv_ny
    size += 1 * blasfeo_memsize_dmat(nu + nx, nu + nx);  // tmp_nv_nv
    size += 1 * blasfeo_memsize_dmat(nu + nx, ny);  // Cyt_tilde
    size += 1 * blasfeo_memsize_dmat(nz, ny);           // Vz
    size += 1 * blasfeo_memsize_dvec(ny);                // tmp_ny
    size += 1 * blasfeo_memsize_dvec(2*ns);              // tmp_2ns
    size += 1 * blasfeo_memsize_dvec(nz);           // tmp_nz

    size += 64;  // blasfeo_mem align
//    size += 8;

    return size;
}



static void ocp_nlp_cost_nls_cast_workspace(void *config_, void *dims_, void *opts_, void *work_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_workspace *work = work_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;
    int nz = dims->nz;
    int ns = dims->ns;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_cost_nls_workspace);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // tmp_nv_ny
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ny, &work->tmp_nv_ny, &c_ptr);

    // tmp_nv_nv
    assign_and_advance_blasfeo_dmat_mem(nu + nx, nu + nx, &work->tmp_nv_nv, &c_ptr);

    // Vz
    assign_and_advance_blasfeo_dmat_mem(nz, ny, &work->Vz, &c_ptr);

    // Cyt_tilde
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ny, &work->Cyt_tilde, &c_ptr);

    // tmp_ny
    assign_and_advance_blasfeo_dvec_mem(ny, &work->tmp_ny, &c_ptr);

    // tmp_2ns
    assign_and_advance_blasfeo_dvec_mem(2*ns, &work->tmp_2ns, &c_ptr);

    // tmp_nz
    assign_and_advance_blasfeo_dvec_mem(nz, &work->tmp_nz, &c_ptr);

    assert((char *) work + ocp_nlp_cost_nls_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

static void ocp_nlp_cost_nls_update_W_factorization(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;
    ocp_nlp_cost_nls_memory *memory = memory_;

    ocp_nlp_cost_nls_cast_workspace(config_, dims_, opts_, work_);

    int ny = dims->ny;

    if (model->W_changed)
    {
        blasfeo_dpotrf_l(ny, &model->W, 0, 0, &memory->W_chol, 0, 0);
        model->W_changed = 0;
    }
    return;
}



void ocp_nlp_cost_nls_precompute(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{
    ocp_nlp_cost_nls_model *model = model_;
    model->W_changed = 1;
    ocp_nlp_cost_nls_update_W_factorization(config_, dims_, model_, opts_, memory_, work_);
    return;
}



void ocp_nlp_cost_nls_initialize(void *config_, void *dims_, void *model_, void *opts_,
                                 void *memory_, void *work_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;
    ocp_nlp_cost_nls_memory *memory = memory_;

    ocp_nlp_cost_nls_cast_workspace(config_, dims_, opts_, work_);
    ocp_nlp_cost_nls_update_W_factorization(config_, dims_, model_, opts_, memory_, work_);

    int ns = dims->ns;
    blasfeo_dveccpsc(2*ns, model->scaling, &model->Z, 0, memory->Z, 0);

    return;
}



void ocp_nlp_cost_nls_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_,
                                         void *memory_, void *work_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;
    ocp_nlp_cost_nls_opts *opts = opts_;
    ocp_nlp_cost_nls_memory *memory = memory_;
    ocp_nlp_cost_nls_workspace *work = work_;

    ocp_nlp_cost_nls_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    ext_fun_arg_t ext_fun_type_in[5];
    void *ext_fun_in[5];
    ext_fun_arg_t ext_fun_type_out[3];
    void *ext_fun_out[3];

    struct blasfeo_dvec_args x_in;  // input x of external fun;
    struct blasfeo_dvec_args u_in;  // input u of external fun;

    if (opts->integrator_cost == 0)
    {
        x_in.x = memory->ux;
        u_in.x = memory->ux;

        x_in.xi = nu;
        u_in.xi = 0;

        ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
        ext_fun_in[0] = &x_in;
        ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
        ext_fun_in[1] = &u_in;
        ext_fun_type_in[2] = BLASFEO_DVEC;
        ext_fun_in[2] = memory->z_alg;
        ext_fun_type_in[3] = COLMAJ;
        ext_fun_in[3] = &model->t;

        ext_fun_type_out[0] = BLASFEO_DVEC;
        ext_fun_out[0] = &memory->res;  // fun: ny
        ext_fun_type_out[1] = BLASFEO_DMAT;
        ext_fun_out[1] = &memory->Jt;  // jac': (nu+nx) * ny
        ext_fun_type_out[2] = BLASFEO_DMAT;
        ext_fun_out[2] = &work->Vz;  // jac_yexpr_z:  ny * nz

        // evaluate external function
        model->nls_y_fun_jac->evaluate(model->nls_y_fun_jac, ext_fun_type_in, ext_fun_in,
                                    ext_fun_type_out, ext_fun_out);

        /* gradient */
        // res = res - y_ref
        blasfeo_daxpy(ny, -1.0, &model->y_ref, 0, &memory->res, 0, &memory->res, 0);

        // printf("W\n");
        // blasfeo_print_dmat(ny, ny, &model->W, 0, 0);

        // printf("res\n");
        // blasfeo_print_dvec(ny, &memory->res, 0);

        // tmp_ny = W * res
        blasfeo_dsymv_l(ny, 1.0, &model->W, 0, 0, &memory->res, 0,
                        0.0, &model->y_ref, 0, &work->tmp_ny, 0);

        if (nz > 0)
        {
            // Cy_tilde = Jt + dzdux_tran*Vz^T
            blasfeo_dgemm_nt(nu + nx, ny, nz, 1.0, memory->dzdux_tran, 0, 0,
                    &work->Vz, 0, 0, 1.0, &memory->Jt, 0, 0, &work->Cyt_tilde, 0, 0);

            // grad = Cyt_tilde * tmp_ny
            blasfeo_dgemv_n(nu+nx, ny, 1.0, &work->Cyt_tilde, 0, 0, &work->tmp_ny, 0,
                            0.0, &memory->grad, 0, &memory->grad, 0);


            // gauss-newton component update
            // tmp_nv_ny = W_chol * Cyt_tilde
            blasfeo_dtrmm_rlnn(nu + nx, ny, 1.0, &memory->W_chol, 0, 0,
                            &work->Cyt_tilde, 0, 0, &work->tmp_nv_ny, 0, 0);
        }
        else
        {
            // grad = Jt * tmp_ny
            blasfeo_dgemv_n(nu+nx, ny, 1.0, &memory->Jt, 0, 0, &work->tmp_ny, 0,
                            0.0, &memory->grad, 0, &memory->grad, 0);
            // gauss-newton component update
            // tmp_nv_ny = Jt * W_chol, where W_chol is lower triangular
            blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, &memory->W_chol, 0, 0, &memory->Jt, 0, 0,
                                &work->tmp_nv_ny, 0, 0);
        }

        // function
        memory->fun = 0.5 * blasfeo_ddot(ny, &work->tmp_ny, 0, &memory->res, 0);
        // printf("tmp_ny\n");
        // blasfeo_print_dvec(ny, &work->tmp_ny, 0);

        // printf("W_chol\n");
        // blasfeo_print_dmat(ny, ny, &memory->W_chol, 0, 0);

        // printf("Jt\n");
        // blasfeo_print_dmat(nu+nx, ny, &memory->Jt, 0, 0);


        /* hessian */
        if (opts->gauss_newton_hess)
        {
            // RSQrq += scaling * tmp_nv_ny * tmp_nv_ny^T
            blasfeo_dsyrk_ln(nu+nx, ny, model->scaling, &work->tmp_nv_ny, 0, 0, &work->tmp_nv_ny, 0, 0,
                            1.0, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);
        }
        else
        {
            if (nz > 0)
            {
                printf("\nocp_nlp_cost_nls_update_qp_matrices: nz > 0 only implemented for gauss_newton_hess.\n");
                exit(1);
            }
            // NOTE(oj): this should add the non-Gauss-Newton term to RSQrq,
            // the product < r, d2_d[x,u] r >, where the cost is 0.5 * norm2(r(x,u))^2
            // exact hessian of ls cost

            // ext_fun_[type_]in 0, 1, 2 are the same as before.
            ext_fun_type_in[3] = BLASFEO_DVEC;
            ext_fun_in[3] = &work->tmp_ny;  // fun: ny
            ext_fun_type_in[4] = COLMAJ;
            ext_fun_in[4] = &model->t;

            ext_fun_type_out[0] = BLASFEO_DMAT;
            ext_fun_out[0] = &work->tmp_nv_nv;   // hess*fun: (nu+nx) * (nu+nx)

            // evaluate external function
            model->nls_y_hess->evaluate(model->nls_y_hess, ext_fun_type_in, ext_fun_in,
                                    ext_fun_type_out, ext_fun_out);

            // RSQrq += scaling * (tmp_nv_nv + tmp_nv_ny * tmp_nv_ny^T)
            blasfeo_dsyrk_ln(nu+nx, ny, model->scaling, &work->tmp_nv_ny, 0, 0, &work->tmp_nv_ny, 0, 0,
                            1.0, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);
            blasfeo_dgead(nu+nx, nu+nx, model->scaling, &work->tmp_nv_nv, 0, 0, memory->RSQrq, 0, 0);
        }
    } // end if (opts->integrator_cost == 0)

    // slack update gradient
    blasfeo_dveccp(2*ns, &model->z, 0, &memory->grad, nu+nx);
    blasfeo_dvecmulacc(2*ns, &model->Z, 0, memory->ux, nu+nx, &memory->grad, nu+nx);

    // slack update function value
    // tmp_2ns = 2 * z + Z .* slack
    blasfeo_dveccpsc(2*ns, 2.0, &model->z, 0, &work->tmp_2ns, 0);
    blasfeo_dvecmulacc(2*ns, &model->Z, 0, memory->ux, nu+nx, &work->tmp_2ns, 0);
    // fun += .5 * (tmp_2ns .* slack)
    memory->fun += 0.5 * blasfeo_ddot(2*ns, &work->tmp_2ns, 0, memory->ux, nu+nx);

    // scale
    if (model->scaling!=1.0)
    {
        blasfeo_dvecsc(nu+nx+2*ns, model->scaling, &memory->grad, 0);
        memory->fun *= model->scaling;
    }
    // printf("cost_fun: %e\n", memory->fun);

    // printf("cost grad\n");
    // blasfeo_print_exp_tran_dvec(nu+nx+2*ns, &memory->grad, 0);

    // printf("cost hess\n");
    // blasfeo_print_dmat(nu+nx, nu+nx, memory->RSQrq, 0, 0);
    // blasfeo_print_tran_dvec(2*ns, memory->Z, 0);
    // blasfeo_print_tran_dvec(nu+nx+2*ns, &memory->grad, 0);
    // exit(1);

    return;
}



void ocp_nlp_cost_nls_compute_fun(void *config_, void *dims_, void *model_,
                                  void *opts_, void *memory_, void *work_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;
    ocp_nlp_cost_nls_opts *opts = opts_;
    ocp_nlp_cost_nls_memory *memory = memory_;
    ocp_nlp_cost_nls_workspace *work = work_;

    struct blasfeo_dvec *ux = memory->ux;

    ocp_nlp_cost_nls_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    if (opts->integrator_cost == 0)
    {
        ext_fun_arg_t nls_y_fun_type_in[4];
        void *nls_y_fun_in[4];
        ext_fun_arg_t nls_y_fun_type_out[1];
        void *nls_y_fun_out[1];

        struct blasfeo_dvec_args x_in;  // input x of external fun;
        struct blasfeo_dvec_args u_in;  // input u of external fun;

        x_in.x = ux;
        x_in.xi = nu;

        u_in.x = ux;
        u_in.xi = 0;

        nls_y_fun_type_in[0] = BLASFEO_DVEC_ARGS;
        nls_y_fun_in[0] = &x_in;

        nls_y_fun_type_in[1] = BLASFEO_DVEC_ARGS;
        nls_y_fun_in[1] = &u_in;

        nls_y_fun_type_in[2] = BLASFEO_DVEC;
        nls_y_fun_in[2] = memory->z_alg;

        nls_y_fun_type_in[3] = COLMAJ;
        nls_y_fun_in[3] = &model->t;

        nls_y_fun_type_out[0] = BLASFEO_DVEC;
        nls_y_fun_out[0] = &memory->res;  // fun: ny

        if (model->nls_y_fun == 0)
        {
            printf("ocp_nlp_cost_nls_compute_fun: nls_y_fun is not provided. Exiting.\n");
            exit(1);
        }
        // evaluate external function
        model->nls_y_fun->evaluate(model->nls_y_fun, nls_y_fun_type_in, nls_y_fun_in,
                                nls_y_fun_type_out, nls_y_fun_out);

        // res = res - y_ref
        blasfeo_daxpy(ny, -1.0, &model->y_ref, 0, &memory->res, 0, &memory->res, 0);
        // tmp_ny = W_chol * nls_res
        blasfeo_dtrmv_ltn(ny, &memory->W_chol, 0, 0, &memory->res, 0, &work->tmp_ny, 0);

        memory->fun = 0.5 * blasfeo_ddot(ny, &work->tmp_ny, 0, &work->tmp_ny, 0);
    }

    // slack update function value
    blasfeo_dveccpsc(2*ns, 2.0, &model->z, 0, &work->tmp_2ns, 0);
    blasfeo_dvecmulacc(2*ns, &model->Z, 0, ux, nu+nx, &work->tmp_2ns, 0);
    memory->fun += 0.5 * blasfeo_ddot(2*ns, &work->tmp_2ns, 0, ux, nu+nx);

    // scale
    if (model->scaling!=1.0)
    {
        memory->fun *= model->scaling;
    }

    return;

}



void ocp_nlp_cost_nls_config_initialize_default(void *config_)
{
    ocp_nlp_cost_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_cost_nls_dims_calculate_size;
    config->dims_assign = &ocp_nlp_cost_nls_dims_assign;
    config->dims_set = &ocp_nlp_cost_nls_dims_set;
    config->dims_get = &ocp_nlp_cost_nls_dims_get;
    config->model_calculate_size = &ocp_nlp_cost_nls_model_calculate_size;
    config->model_assign = &ocp_nlp_cost_nls_model_assign;
    config->model_set = &ocp_nlp_cost_nls_model_set;
    config->opts_calculate_size = &ocp_nlp_cost_nls_opts_calculate_size;
    config->opts_assign = &ocp_nlp_cost_nls_opts_assign;
    config->opts_initialize_default = &ocp_nlp_cost_nls_opts_initialize_default;
    config->opts_update = &ocp_nlp_cost_nls_opts_update;
    config->opts_set = &ocp_nlp_cost_nls_opts_set;
    config->memory_calculate_size = &ocp_nlp_cost_nls_memory_calculate_size;
    config->memory_assign = &ocp_nlp_cost_nls_memory_assign;
    config->memory_get_fun_ptr = &ocp_nlp_cost_nls_memory_get_fun_ptr;
    config->memory_get_grad_ptr = &ocp_nlp_cost_nls_memory_get_grad_ptr;
    config->memory_get_W_chol_ptr = &ocp_nlp_cost_nls_memory_get_W_chol_ptr;
    config->model_get_y_ref_ptr = &ocp_nlp_cost_nls_model_get_y_ref_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_cost_nls_memory_set_ux_ptr;
    config->memory_set_z_alg_ptr = &ocp_nlp_cost_nls_memory_set_z_alg_ptr;
    config->memory_set_dzdux_tran_ptr = &ocp_nlp_cost_nls_memory_set_dzdux_tran_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_cost_nls_memory_set_RSQrq_ptr;
    config->memory_set_Z_ptr = &ocp_nlp_cost_nls_memory_set_Z_ptr;
    config->workspace_calculate_size = &ocp_nlp_cost_nls_workspace_calculate_size;
    config->initialize = &ocp_nlp_cost_nls_initialize;
    config->update_qp_matrices = &ocp_nlp_cost_nls_update_qp_matrices;
    config->compute_fun = &ocp_nlp_cost_nls_compute_fun;
    config->config_initialize_default = &ocp_nlp_cost_nls_config_initialize_default;
    config->precompute = &ocp_nlp_cost_nls_precompute;

    return;
}
