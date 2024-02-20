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


#include "acados/ocp_nlp/ocp_nlp_cost_conl.h"
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

acados_size_t ocp_nlp_cost_conl_dims_calculate_size(void *config_)
{
    acados_size_t size = sizeof(ocp_nlp_cost_conl_dims);

    return size;
}



void *ocp_nlp_cost_conl_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_conl_dims *dims = (ocp_nlp_cost_conl_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_conl_dims);

    assert((char *) raw_memory + ocp_nlp_cost_conl_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void ocp_nlp_cost_conl_dims_set(void *config_, void *dims_, const char *field, int* value)
{
    ocp_nlp_cost_conl_dims *dims = (ocp_nlp_cost_conl_dims *) dims_;
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
static void ocp_nlp_cost_conl_get_ny(void *config_, void *dims_, int* value)
{
    ocp_nlp_cost_conl_dims *dims = (ocp_nlp_cost_conl_dims *) dims_;
    *value = dims->ny;
}



void ocp_nlp_cost_conl_dims_get(void *config_, void *dims_, const char *field, int* value)
{
    if (!strcmp(field, "ny"))
    {
        ocp_nlp_cost_conl_get_ny(config_, dims_, value);
    }
    else
    {
        printf("error: ocp_nlp_cost_conl_dims_get: attempt to get dimensions of non-existing field %s\n", field);
        exit(1);
    }
}



/************************************************
 * model
 ************************************************/

acados_size_t ocp_nlp_cost_conl_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_cost_conl_dims *dims = dims_;

    int ny = dims->ny;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_conl_model);

    size += 64;  // blasfeo_mem align

    size += 1 * blasfeo_memsize_dvec(ny);      // y_ref
    size += 2 * blasfeo_memsize_dvec(2 * ns);  // Z, z

    return size;
}



void *ocp_nlp_cost_conl_model_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_cost_conl_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    int ny = dims->ny;
    int ns = dims->ns;

    // struct
    ocp_nlp_cost_conl_model *model = (ocp_nlp_cost_conl_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_conl_model);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

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

    // assert
    assert((char *) raw_memory + ocp_nlp_cost_conl_model_calculate_size(config_, dims) >= c_ptr);

    return model;
}



int ocp_nlp_cost_conl_model_set(void *config_, void *dims_, void *model_,
                                         const char *field, void *value_)
{
    int status = ACADOS_SUCCESS;

    if ( !config_ || !dims_ || !model_ || !value_ )
    {
        printf("ocp_nlp_cost_conl_model_set: got NULL pointer \n");
        exit(1);
    }

    ocp_nlp_cost_conl_dims *dims = dims_;
    ocp_nlp_cost_conl_model *model = model_;

    int ny = dims->ny;
    int ns = dims->ns;

    if (!strcmp(field, "y_ref") || !strcmp(field, "yref"))
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
    else if (!strcmp(field, "conl_cost_fun_jac_hess"))
    {
        model->conl_cost_fun_jac_hess = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "conl_cost_fun"))
    {
        model->conl_cost_fun = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "scaling"))
    {
        double *scaling_ptr = (double *) value_;
        model->scaling = *scaling_ptr;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_cost_conl_model_set\n", field);
        exit(1);
    }
    return status;
}



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_cost_conl_opts_calculate_size(void *config_, void *dims_)
{

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_conl_opts);
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_cost_conl_opts_assign(void *config_, void *dims_, void *raw_memory)
{

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_conl_opts *opts = (ocp_nlp_cost_conl_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_conl_opts);

    assert((char *) raw_memory + ocp_nlp_cost_conl_opts_calculate_size(config_, dims_) >= c_ptr);

    return opts;
}



void ocp_nlp_cost_conl_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_cost_conl_opts *opts = opts_;

    opts->gauss_newton_hess = 1;

    return;
}



void ocp_nlp_cost_conl_opts_update(void *config_, void *dims_, void *opts_)
{
    return;
}



void ocp_nlp_cost_conl_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_cost_conl_opts *opts = opts_;

    if(!strcmp(field, "gauss_newton_hess"))
    {
        int *int_ptr = value;
        opts->gauss_newton_hess = *int_ptr;  // NOTE: we always use a Gauss-Newton Hessian
    }
    else if(!strcmp(field, "integrator_cost"))
    {
        int *opt_val = (int *) value;
        opts->integrator_cost = *opt_val;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_cost_conl_opts_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_cost_conl_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_cost_conl_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;
    int ny = dims->ny;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_conl_memory);

    size += 1 * blasfeo_memsize_dvec(nu + nx + 2 * ns);  // grad
    size += 64;  // blasfeo_mem align

    size += 1 * blasfeo_memsize_dmat(ny, ny);            // W_chol

    return size;
}



void *ocp_nlp_cost_conl_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_cost_conl_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;
    int ny = dims->ny;

    // struct
    ocp_nlp_cost_conl_memory *memory = (ocp_nlp_cost_conl_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_conl_memory);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // W_chol
    assign_and_advance_blasfeo_dmat_mem(ny, ny, &memory->W_chol, &c_ptr);

    // grad
    assign_and_advance_blasfeo_dvec_mem(nu + nx + 2 * ns, &memory->grad, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_cost_conl_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}



double *ocp_nlp_cost_conl_memory_get_fun_ptr(void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;
    return &memory->fun;
}

struct blasfeo_dmat *ocp_nlp_cost_conl_memory_get_W_chol_ptr(void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;
    return &memory->W_chol;
}



struct blasfeo_dvec *ocp_nlp_cost_conl_memory_get_grad_ptr(void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;

    return &memory->grad;
}



struct blasfeo_dvec *ocp_nlp_cost_conl_model_get_y_ref_ptr(void *in_)
{
    ocp_nlp_cost_conl_model *model = in_;

    return &model->y_ref;
}


void ocp_nlp_cost_conl_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;

    memory->RSQrq = RSQrq;

    return;
}



void ocp_nlp_cost_conl_memory_set_Z_ptr(struct blasfeo_dvec *Z, void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;

    memory->Z = Z;

    return;
}



void ocp_nlp_cost_conl_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;

    memory->ux = ux;

    return;
}



void ocp_nlp_cost_conl_memory_set_z_alg_ptr(struct blasfeo_dvec *z_alg, void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;

    memory->z_alg = z_alg;
}



void ocp_nlp_cost_conl_memory_set_dzdux_tran_ptr(struct blasfeo_dmat *dzdux_tran, void *memory_)
{
    ocp_nlp_cost_conl_memory *memory = memory_;

    memory->dzdux_tran = dzdux_tran;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_cost_conl_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_cost_conl_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_conl_workspace);

    size += 1 * blasfeo_memsize_dmat(ny, ny);            // W
    size += 1 * blasfeo_memsize_dmat(nu + nx, ny);       // Jt_ux
    size += 1 * blasfeo_memsize_dmat(nu + nx, ny);       // Jt_ux_tilde
    size += 1 * blasfeo_memsize_dmat(nz, ny);            // Jt_z
    size += 1 * blasfeo_memsize_dmat(nu + nx, ny);       // tmp_nv_ny
    size += 1 * blasfeo_memsize_dvec(ny);                // tmp_ny
    size += 1 * blasfeo_memsize_dvec(2*ns);              // tmp_2ns

    size += 64;  // blasfeo_mem align

    return size;
}



static void ocp_nlp_cost_conl_cast_workspace(void *config_, void *dims_, void *opts_, void *work_)
{
    ocp_nlp_cost_conl_dims *dims = dims_;
    ocp_nlp_cost_conl_workspace *work = work_;

    // extract dims
    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_cost_conl_workspace);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // W
    assign_and_advance_blasfeo_dmat_mem(ny, ny, &work->W, &c_ptr);

    // Jt_ux
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ny, &work->Jt_ux, &c_ptr);

    // Jt_ux_tilde
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ny, &work->Jt_ux_tilde, &c_ptr);

    // Jt_z
    assign_and_advance_blasfeo_dmat_mem(nz, ny, &work->Jt_z, &c_ptr);

    // tmp_nv_ny
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ny, &work->tmp_nv_ny, &c_ptr);

    // tmp_ny
    assign_and_advance_blasfeo_dvec_mem(ny, &work->tmp_ny, &c_ptr);

    // tmp_2ns
    assign_and_advance_blasfeo_dvec_mem(2*ns, &work->tmp_2ns, &c_ptr);

    assert((char *) work + ocp_nlp_cost_conl_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_cost_conl_precompute(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{
    return;
}

void ocp_nlp_cost_conl_initialize(void *config_, void *dims_, void *model_, void *opts_,
                                  void *memory_, void *work_)
{
    ocp_nlp_cost_conl_dims *dims = dims_;
    ocp_nlp_cost_conl_model *model = model_;
    ocp_nlp_cost_conl_memory *memory = memory_;

    ocp_nlp_cost_conl_cast_workspace(config_, dims, opts_, work_);

    int ns = dims->ns;

    blasfeo_dveccpsc(2*ns, model->scaling, &model->Z, 0, memory->Z, 0);

    return;
}

// NOTE: it implementing without W_chol should have the same computational complexity.


void ocp_nlp_cost_conl_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_,
                                         void *memory_, void *work_)
{
    // NOTE: We assume that opts->gauss_newton_hess is True (this is checked in the interface)

    ocp_nlp_cost_conl_dims *dims = dims_;
    ocp_nlp_cost_conl_model *model = model_;
    ocp_nlp_cost_conl_memory *memory = memory_;
    ocp_nlp_cost_conl_workspace *work = work_;
    ocp_nlp_cost_conl_opts *opts = (ocp_nlp_cost_conl_opts *) opts_;

    ocp_nlp_cost_conl_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    if (opts->integrator_cost == 0)
    {
        ext_fun_arg_t conl_fun_jac_hess_type_in[5];
        void *conl_fun_jac_hess_in[5];
        ext_fun_arg_t conl_fun_jac_hess_type_out[5];
        void *conl_fun_jac_hess_out[5];

        // INPUT
        struct blasfeo_dvec_args x_in;  // input x
        x_in.x = memory->ux;
        x_in.xi = nu;

        struct blasfeo_dvec_args u_in;  // input u
        u_in.x = memory->ux;
        u_in.xi = 0;

        conl_fun_jac_hess_type_in[0] = BLASFEO_DVEC_ARGS;
        conl_fun_jac_hess_in[0] = &x_in;
        conl_fun_jac_hess_type_in[1] = BLASFEO_DVEC_ARGS;
        conl_fun_jac_hess_in[1] = &u_in;

        conl_fun_jac_hess_type_in[2] = BLASFEO_DVEC;
        conl_fun_jac_hess_in[2] = memory->z_alg;
        conl_fun_jac_hess_type_in[3] = BLASFEO_DVEC;
        conl_fun_jac_hess_in[3] = &model->y_ref;

        conl_fun_jac_hess_type_in[4] = COLMAJ;
        conl_fun_jac_hess_in[4] = &model->t;

        // OUTPUT
        conl_fun_jac_hess_type_out[0] = COLMAJ;
        conl_fun_jac_hess_out[0] = &memory->fun;         // fun: scalar
        conl_fun_jac_hess_type_out[1] = BLASFEO_DVEC;
        conl_fun_jac_hess_out[1] = &work->tmp_ny;        // grad of outer loss wrt residual, ny
        conl_fun_jac_hess_type_out[2] = BLASFEO_DMAT;
        conl_fun_jac_hess_out[2] = &work->Jt_ux;         // inner Jacobian wrt ux, transposed, (nu+nx) x ny
        conl_fun_jac_hess_type_out[3] = BLASFEO_DMAT;
        conl_fun_jac_hess_out[3] = &work->Jt_z;          // inner Jacobian wrt z, transposed, nz x ny
        conl_fun_jac_hess_type_out[4] = BLASFEO_DMAT;
        conl_fun_jac_hess_out[4] = &work->W;             // outer hessian: ny x ny

        // evaluate external function
        model->conl_cost_fun_jac_hess->evaluate(model->conl_cost_fun_jac_hess, conl_fun_jac_hess_type_in,
                                                conl_fun_jac_hess_in, conl_fun_jac_hess_type_out, conl_fun_jac_hess_out);

        // hessian of outer loss function
        blasfeo_dpotrf_l(ny, &work->W, 0, 0, &memory->W_chol, 0, 0);

        if (nz > 0)
        {
            // Jt_ux_tilde = Jt_ux + dzdux_tran*Jt_z
            blasfeo_dgemm_nn(nu + nx, ny, nz, 1.0, memory->dzdux_tran, 0, 0,
                    &work->Jt_z, 0, 0, 1.0, &work->Jt_ux, 0, 0, &work->Jt_ux_tilde, 0, 0);

            // grad = Jt_ux_tilde * tmp_ny
            blasfeo_dgemv_n(nu+nx, ny, 1.0, &work->Jt_ux_tilde, 0, 0, &work->tmp_ny, 0,
                            0.0, &memory->grad, 0, &memory->grad, 0);

            // tmp_nv_ny = Jt_ux_tilde * W_chol
            blasfeo_dtrmm_rlnn(nu + nx, ny, 1.0, &memory->W_chol, 0, 0,
                            &work->Jt_ux_tilde, 0, 0, &work->tmp_nv_ny, 0, 0);
        }
        else
        {
            // grad = Jt_ux * tmp_ny
            blasfeo_dgemv_n(nu+nx, ny, 1.0, &work->Jt_ux, 0, 0, &work->tmp_ny, 0,
                            0.0, &memory->grad, 0, &memory->grad, 0);

            // tmp_nv_ny = Jt_ux * W_chol, where W_chol is lower triangular
            blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, &memory->W_chol, 0, 0, &work->Jt_ux, 0, 0,
                                &work->tmp_nv_ny, 0, 0);

        }
        // RSQrq += scaling * tmp_nv_ny * tmp_nv_ny^T
        blasfeo_dsyrk_ln(nu+nx, ny, model->scaling, &work->tmp_nv_ny, 0, 0, &work->tmp_nv_ny, 0, 0,
                        1.0, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);
    }

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
    if(model->scaling!=1.0)
    {
        blasfeo_dvecsc(nu+nx+2*ns, model->scaling, &memory->grad, 0);
        memory->fun *= model->scaling;
    }

    return;
}



void ocp_nlp_cost_conl_compute_fun(void *config_, void *dims_, void *model_,
                                  void *opts_, void *memory_, void *work_)
{
    ocp_nlp_cost_conl_dims *dims = dims_;
    ocp_nlp_cost_conl_model *model = model_;
    ocp_nlp_cost_conl_memory *memory = memory_;
    ocp_nlp_cost_conl_workspace *work = work_;
    ocp_nlp_cost_conl_opts *opts = opts_;

    ocp_nlp_cost_conl_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    struct blasfeo_dvec *ux = memory->ux;

    if (opts->integrator_cost == 0)
    {
        /* specify input types and pointers for external cost function */
        ext_fun_arg_t ext_fun_type_in[5];
        void *ext_fun_in[5];
        ext_fun_arg_t ext_fun_type_out[1];
        void *ext_fun_out[1];

        // INPUT
        struct blasfeo_dvec_args x_in;  // input x
        x_in.x = ux;
        x_in.xi = nu;

        struct blasfeo_dvec_args u_in;  // input u
        u_in.x = ux;
        u_in.xi = 0;

        ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
        ext_fun_in[0] = &x_in;
        ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
        ext_fun_in[1] = &u_in;
        ext_fun_type_in[2] = BLASFEO_DVEC;
        ext_fun_in[2] = memory->z_alg;
        ext_fun_type_in[3] = BLASFEO_DVEC;
        ext_fun_in[3] = &model->y_ref;
        ext_fun_type_in[4] = COLMAJ;
        ext_fun_in[4] = &model->t;

        // OUTPUT
        ext_fun_type_out[0] = COLMAJ;
        ext_fun_out[0] = &memory->fun;  // function: scalar

        model->conl_cost_fun->evaluate(model->conl_cost_fun, ext_fun_type_in, ext_fun_in,
                                    ext_fun_type_out, ext_fun_out);
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



void ocp_nlp_cost_conl_config_initialize_default(void *config_)
{
    ocp_nlp_cost_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_cost_conl_dims_calculate_size;
    config->dims_assign = &ocp_nlp_cost_conl_dims_assign;
    config->dims_set = &ocp_nlp_cost_conl_dims_set;
    config->dims_get = &ocp_nlp_cost_conl_dims_get;
    config->model_calculate_size = &ocp_nlp_cost_conl_model_calculate_size;
    config->model_assign = &ocp_nlp_cost_conl_model_assign;
    config->model_set = &ocp_nlp_cost_conl_model_set;
    config->opts_calculate_size = &ocp_nlp_cost_conl_opts_calculate_size;
    config->opts_assign = &ocp_nlp_cost_conl_opts_assign;
    config->opts_initialize_default = &ocp_nlp_cost_conl_opts_initialize_default;
    config->opts_update = &ocp_nlp_cost_conl_opts_update;
    config->opts_set = &ocp_nlp_cost_conl_opts_set;
    config->memory_calculate_size = &ocp_nlp_cost_conl_memory_calculate_size;
    config->memory_assign = &ocp_nlp_cost_conl_memory_assign;
    config->memory_get_fun_ptr = &ocp_nlp_cost_conl_memory_get_fun_ptr;
    config->memory_get_grad_ptr = &ocp_nlp_cost_conl_memory_get_grad_ptr;
    config->memory_get_W_chol_ptr = &ocp_nlp_cost_conl_memory_get_W_chol_ptr;
    config->model_get_y_ref_ptr = &ocp_nlp_cost_conl_model_get_y_ref_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_cost_conl_memory_set_ux_ptr;
    config->memory_set_z_alg_ptr = &ocp_nlp_cost_conl_memory_set_z_alg_ptr;
    config->memory_set_dzdux_tran_ptr = &ocp_nlp_cost_conl_memory_set_dzdux_tran_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_cost_conl_memory_set_RSQrq_ptr;
    config->memory_set_Z_ptr = &ocp_nlp_cost_conl_memory_set_Z_ptr;
    config->workspace_calculate_size = &ocp_nlp_cost_conl_workspace_calculate_size;
    config->initialize = &ocp_nlp_cost_conl_initialize;
    config->update_qp_matrices = &ocp_nlp_cost_conl_update_qp_matrices;
    config->compute_fun = &ocp_nlp_cost_conl_compute_fun;
    config->config_initialize_default = &ocp_nlp_cost_conl_config_initialize_default;
    config->precompute = &ocp_nlp_cost_conl_precompute;

    return;
}
