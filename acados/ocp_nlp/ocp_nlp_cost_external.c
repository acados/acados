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


#include "acados/ocp_nlp/ocp_nlp_cost_external.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"



/************************************************
 * model
 ************************************************/

acados_size_t ocp_nlp_cost_external_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_cost_dims *dims = dims_;

    int nx = dims->nx;
    int nu = dims->nu;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_external_model);

    size += 1 * 64;  // blasfeo_mem align
    size += blasfeo_memsize_dmat(nx+nu, nx+nu);

    size += ocp_nlp_cost_common_model_calculate_size(dims);  // common

    return size;
}



void *ocp_nlp_cost_external_model_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_cost_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    int nx = dims->nx;
    int nu = dims->nu;

    // struct
    ocp_nlp_cost_external_model *model = (ocp_nlp_cost_external_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_model);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);
    // numerical_hessian
    assign_and_advance_blasfeo_dmat_mem(nx+nu, nx+nu, &model->numerical_hessian, &c_ptr);

    // common
    model->common = ocp_nlp_cost_common_model_assign(dims, &c_ptr);

    // assert
    assert((char *) raw_memory + ocp_nlp_cost_external_model_calculate_size(config_, dims_) >=
           c_ptr);

    return model;
}



int ocp_nlp_cost_external_model_set(void *config_, void *dims_, void *model_,
                                         const char *field, void *value_)
{
    int status = ACADOS_SUCCESS;

    if ( !config_ || !dims_ || !model_ || !value_ )
    {
        printf("ocp_nlp_cost_external_model_set: got Null pointer \n");
        exit(1);
    }

    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;

    int nx = dims->nx;
    int nu = dims->nu;

    if (!strcmp(field, "ext_cost_fun"))
    {
        model->ext_cost_fun = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_fun_jac_hes") || !strcmp(field, "ext_cost_fun_jac_hess"))
    {
        model->ext_cost_fun_jac_hess = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_fun_jac"))
    {
        model->ext_cost_fun_jac = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_hess_xu_p"))
    {
        model->ext_cost_hess_xu_p = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_adj_ux_pdiff"))
    {
        model->ext_cost_adj_ux_pdiff = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_grad_p"))
    {
        model->ext_cost_grad_p = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_num_hess"))
    {
        double *numerical_hessian = (double *) value_;
        blasfeo_pack_dmat(nx+nu, nx+nu, numerical_hessian, nx+nu, &model->numerical_hessian, 0, 0);
    }
    else if (ocp_nlp_cost_common_model_set(dims, model->common, field, value_))
    {
        // Zl, Zu, z, zl, zu, scaling handled by common setter
    }
    else
    {
        printf("\nerror: %s not available in module ocp_nlp_cost_external_model_set\n", field);
        exit(1);
    }
    return status;
}



int ocp_nlp_cost_external_model_get(void *config_, void *dims_, void *model_,
                                         const char *field, void *value_)
{
    int status = ACADOS_SUCCESS;

    if ( !config_ || !dims_ || !model_ || !value_ )
    {
        printf("ocp_nlp_cost_external_model_set: got Null pointer \n");
        exit(1);
    }

    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;

    int nx = dims->nx;
    int nu = dims->nu;

    double * value = (double *) value_;

    if (!strcmp(field, "ext_cost_num_hess"))
    {
        blasfeo_unpack_dmat(nx+nu, nx+nu, &model->numerical_hessian, 0, 0, value, nx+nu);
    }
    else if (ocp_nlp_cost_common_model_get(dims, model->common, field, value_))
    {
        // Zl, Zu, zl, zu, scaling handled by common getter
    }
    else
    {
        printf("\nerror: %s not available in module ocp_nlp_cost_external_model_get\n", field);
        exit(1);
    }
    return status;
}


/************************************************
 * options
 ************************************************/

void ocp_nlp_cost_external_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_cost_external_opts *opts = opts_;

    // NOTE: the exact hessian is always computed if no custom hessian is provided,
    // ignore "exact_hess" option

    return;
}

/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_cost_external_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_dims *dims = dims_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_external_memory);

    size += ocp_nlp_cost_common_memory_calculate_size(dims);

    return size;
}



void *ocp_nlp_cost_external_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // struct
    ocp_nlp_cost_external_memory *memory = (ocp_nlp_cost_external_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_memory);

    // common
    memory->common = ocp_nlp_cost_common_memory_assign(dims, &c_ptr);

    assert((char *) raw_memory +
               ocp_nlp_cost_external_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}



void *ocp_nlp_cost_external_memory_get(void *memory_, const char *field)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    void *out = ocp_nlp_cost_common_memory_get(memory->common, field);
    if (out)
    {
        return out;
    }

    printf("\nerror: field %s not available in ocp_nlp_cost_external_memory_get\n", field);
    exit(1);
}

void ocp_nlp_cost_external_memory_set(void *config_, void *dims_, void *memory_, const char *field, void *value)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    if (ocp_nlp_cost_common_memory_set(memory->common, field, value))
    {
        // most handled by common setter
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_cost_external_memory_set\n", field);
        exit(1);
    }
}


/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_cost_external_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_common_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ns = dims->ns;
    int np_global = dims->np_global;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_external_workspace);

    if (opts->with_solution_sens_wrt_params_forw)
    {
        size += 1 * blasfeo_memsize_dmat(nu + nx, np_global);  // cost_grad_params_jac
    }
    if (opts->with_solution_sens_wrt_params_adj)
    {
        size += 1 * blasfeo_memsize_dvec(np_global);  // adj_cost_ux_pdiff
    }
    size += 1 * blasfeo_memsize_dmat(nu+nx, nu+nx);  // tmp_nunx_nunx
    size += 1 * blasfeo_memsize_dmat(nz, nz);  // tmp_nz_nz
    size += 1 * blasfeo_memsize_dmat(nz, nu+nx);  // tmp_nz_nunx
    size += 1 * blasfeo_memsize_dvec(nu+nx+nz);  // tmp_nunxnz

    size += 1 * blasfeo_memsize_dvec(2*ns);  // tmp_2ns

    size += 64;  // blasfeo_mem align

    return size;
}



static void ocp_nlp_cost_external_cast_workspace(void *config_, void *dims_, void *opts_,
                                                 void *work_)
{
    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_workspace *work = work_;
    ocp_nlp_cost_common_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;
    int ns = dims->ns;
    int np_global = dims->np_global;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_cost_external_workspace);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);


    if (opts->with_solution_sens_wrt_params_forw)
    {
        assign_and_advance_blasfeo_dmat_mem(nu + nx, np_global, &work->cost_grad_params_jac, &c_ptr);
    }

    // tmp_nunx_nunx
    assign_and_advance_blasfeo_dmat_mem(nu + nx, nu + nx, &work->tmp_nunx_nunx, &c_ptr);

    // tmp_nz_nz
    assign_and_advance_blasfeo_dmat_mem(nz, nz, &work->tmp_nz_nz, &c_ptr);

    // tmp_nz_nunx
    assign_and_advance_blasfeo_dmat_mem(nz, nu+nx, &work->tmp_nz_nunx, &c_ptr);

    // tmp_nunxnz
    assign_and_advance_blasfeo_dvec_mem(nu + nx + nz, &work->tmp_nunxnz, &c_ptr);

    // tmp_2ns
    assign_and_advance_blasfeo_dvec_mem(2*ns, &work->tmp_2ns, &c_ptr);

    if (opts->with_solution_sens_wrt_params_adj)
    {
        assign_and_advance_blasfeo_dvec_mem(np_global, &work->adj_cost_ux_pdiff, &c_ptr);
    }

    assert((char *) work_ + ocp_nlp_cost_external_workspace_calculate_size(config_, dims_, opts_) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_cost_external_precompute(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{
    return;
}



void ocp_nlp_cost_external_initialize(void *config_, void *dims_, void *model_, void *opts_,
                                      void *memory_, void *work_)
{
    ocp_nlp_cost_external_model *model = model_;
    ocp_nlp_cost_external_memory *memory = memory_;

    ocp_nlp_cost_common_initialize(dims_, model->common, memory->common);

    return;
}



void ocp_nlp_cost_external_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_,
                                              void *memory_, void *work_)
{
    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    ocp_nlp_cost_common_opts *opts = opts_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_common_memory *mem_common = memory->common;
    ocp_nlp_cost_external_workspace *work = work_;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;

    /* specify input types and pointers for external cost function */
    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[5];
    void *ext_fun_out[5];

    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = mem_common->ux;
    u_in.xi = 0;
    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = mem_common->ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;
    ext_fun_type_in[2] = BLASFEO_DVEC;
    ext_fun_in[2] = mem_common->z_alg;

    // OUTPUT
    ext_fun_type_out[0] = COLMAJ;
    ext_fun_out[0] = &mem_common->fun;  // fun: scalar

    ext_fun_type_out[1] = BLASFEO_DVEC;
    ext_fun_out[1] = &work->tmp_nunxnz;  // tmp_nunxnz: nu+nx+nz

    if (opts->use_numerical_hessian > 0)
    {
        // evaluate external function
        model->ext_cost_fun_jac->evaluate(model->ext_cost_fun_jac, ext_fun_type_in,
                                            ext_fun_in, ext_fun_type_out, ext_fun_out);
        // custom hessian
        if (opts->add_hess_contribution)
        {
            blasfeo_dgead(nx+nu, nx+nu, model->common->scaling, &model->numerical_hessian, 0, 0, mem_common->RSQrq, 0, 0);
        }
        else
        {
            blasfeo_dgecpsc(nx+nu, nx+nu, model->common->scaling, &model->numerical_hessian, 0, 0, mem_common->RSQrq, 0, 0);
        }
    }
    else
    {
        // additional output
        ext_fun_type_out[2] = BLASFEO_DMAT;
        ext_fun_out[2] = &work->tmp_nunx_nunx;   // hess: (nu+nx) * (nu+nx)
        ext_fun_type_out[3] = BLASFEO_DMAT;
        ext_fun_out[3] = &work->tmp_nz_nz;       // hess_z: nz x nz
        ext_fun_type_out[4] = BLASFEO_DMAT;
        ext_fun_out[4] = &work->tmp_nz_nunx;    // hess_z_nunx: nz x nu+nx

        // evaluate external function
        model->ext_cost_fun_jac_hess->evaluate(model->ext_cost_fun_jac_hess, ext_fun_type_in,
                                            ext_fun_in, ext_fun_type_out, ext_fun_out);

        // hessian contribution from xu with scaling
        if (opts->add_hess_contribution)
        {
            // add to RSQrq
            blasfeo_dgead(nx+nu, nx+nu, model->common->scaling, &work->tmp_nunx_nunx, 0, 0, mem_common->RSQrq, 0, 0);
        }
        else
        {
            // copy to RSQrq
            blasfeo_dgecpsc(nx+nu, nx+nu, model->common->scaling, &work->tmp_nunx_nunx, 0, 0, mem_common->RSQrq, 0, 0);
        }

        if (nz > 0)
        {
            // NOTE: we compute the Hessian as follows:
            // H = d2l_dxu2 + dz_dux.T * d2l_dz2 * dz_dux + d2l_dux_dz * dz_dux + (d2l_dux_dz * dz_dux).T
            // the term d2z_dux2 is dropped!

            // compute and add cross terms (NOTE: only lower triangular is computed)
            blasfeo_dsyr2k_ln(nu+nx, nz, model->common->scaling, mem_common->dzdux_tran, 0, 0, &work->tmp_nz_nunx, 0, 0, 1., mem_common->RSQrq, 0, 0, mem_common->RSQrq, 0, 0);

            // hessian contribution from z
            blasfeo_dgemm_nt(nz, nu+nx, nz, 1., &work->tmp_nz_nz, 0, 0, mem_common->dzdux_tran, 0, 0, 0.0, &work->tmp_nz_nunx, 0, 0, &work->tmp_nz_nunx, 0, 0);
            blasfeo_dgemm_nn(nu+nx, nu+nx, nz, model->common->scaling, mem_common->dzdux_tran, 0, 0, &work->tmp_nz_nunx, 0, 0, 1., mem_common->RSQrq, 0, 0, mem_common->RSQrq, 0, 0);
        }
    }

    // gradient
    blasfeo_dveccp(nu+nx, &work->tmp_nunxnz, 0, &mem_common->grad, 0);
    if (nz > 0)
    {
        blasfeo_dgemv_n(nu+nx, nz, 1.0, mem_common->dzdux_tran, 0, 0, &work->tmp_nunxnz, nu+nx, 1., &mem_common->grad, 0, &mem_common->grad, 0);
    }

    // slack update gradient and function value, and scale
    cost_common_update_gradient_with_slacks_and_scale(dims, model->common, memory->common);
    cost_common_add_slack_contributions_to_fun_and_scale(dims, model->common, memory->common, &work->tmp_2ns);

    return;
}



void ocp_nlp_cost_external_compute_gradient(void *config_, void *dims_, void *model_, void *opts_,
                                 void *memory_, void *work_)
{
    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    // ocp_nlp_cost_external_opts *opts = opts_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_external_workspace *work = work_;
    ocp_nlp_cost_common_memory *mem_common = memory->common;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nz = dims->nz;
    int nu = dims->nu;

    /* specify input types and pointers for external cost function */
    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[2];
    void *ext_fun_out[2];

    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = mem_common->ux;
    u_in.xi = 0;
    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = mem_common->ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;
    ext_fun_type_in[2] = BLASFEO_DVEC;
    ext_fun_in[2] = mem_common->z_alg;

    // OUTPUT
    ext_fun_type_out[0] = COLMAJ;
    ext_fun_out[0] = &mem_common->fun;  // fun: scalar

    ext_fun_type_out[1] = BLASFEO_DVEC;
    ext_fun_out[1] = &work->tmp_nunxnz;  // tmp_nunxnz: nu+nx+nz

    // evaluate external function
    model->ext_cost_fun_jac->evaluate(model->ext_cost_fun_jac, ext_fun_type_in,
                                        ext_fun_in, ext_fun_type_out, ext_fun_out);

    // gradient
    blasfeo_dveccp(nu+nx, &work->tmp_nunxnz, 0, &mem_common->grad, 0);
    if (nz > 0)
    {
        blasfeo_dgemv_n(nu+nx, nz, 1.0, mem_common->dzdux_tran, 0, 0, &work->tmp_nunxnz, nu+nx, 1., &mem_common->grad, 0, &mem_common->grad, 0);
    }

    // slack update gradient
    cost_common_update_gradient_with_slacks_and_scale(dims, model->common, memory->common);
}



void ocp_nlp_cost_external_compute_fun(void *config_, void *dims_, void *model_,
                                       void *opts_, void *memory_, void *work_)
{

    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    // ocp_nlp_cost_external_opts *opts = opts_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_external_workspace *work = work_;
    ocp_nlp_cost_common_memory *mem_common = memory->common;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    struct blasfeo_dvec *ux = mem_common->ux;

    int nu = dims->nu;

    /* specify input types and pointers for external cost function */
    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[1];
    void *ext_fun_out[1];

    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = ux;
    u_in.xi = 0;

    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;
    ext_fun_type_in[2] = BLASFEO_DVEC;
    ext_fun_in[2] = mem_common->z_alg;
    // OUTPUT
    ext_fun_type_out[0] = COLMAJ;
    ext_fun_out[0] = &mem_common->fun;  // function: scalar

    // evaluate external function
    if (model->ext_cost_fun == 0)
    {
        printf("ocp_nlp_cost_external_compute_fun: ext_cost_fun is not provided. Exiting.\n");
        exit(1);
    }
    model->ext_cost_fun->evaluate(model->ext_cost_fun, ext_fun_type_in, ext_fun_in,
                                  ext_fun_type_out, ext_fun_out);

    // slack update function value and scale
    cost_common_add_slack_contributions_to_fun_and_scale(dims, model->common, memory->common, &work->tmp_2ns);

    return;
}

void ocp_nlp_cost_external_compute_jac_p(void *config_, void *dims_, void *model_,
                                       void *opts_, void *memory_, void *work_)
{
    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    ocp_nlp_cost_external_workspace *work = work_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_common_memory *mem_common = memory->common;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    struct blasfeo_dvec *ux = mem_common->ux;

    int nu = dims->nu;
    int nx = dims->nx;
    int np_global = dims->np_global;

    /* specify input types and pointers for external cost function */
    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[1];
    void *ext_fun_out[1];

    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = ux;
    u_in.xi = 0;

    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;
    ext_fun_type_in[2] = BLASFEO_DVEC;
    ext_fun_in[2] = mem_common->z_alg;

    // OUTPUT
    ext_fun_type_out[0] = BLASFEO_DMAT;
    ext_fun_out[0] = &work->cost_grad_params_jac;

    // evaluate external function
    if (model->ext_cost_hess_xu_p == 0)
    {
        printf("ocp_nlp_cost_external_compute_jac_p: ext_cost_hess_xu_p is not provided. Exiting.\n");
        exit(1);
    }
    model->ext_cost_hess_xu_p->evaluate(model->ext_cost_hess_xu_p, ext_fun_type_in, ext_fun_in,
                                  ext_fun_type_out, ext_fun_out);

    // add contribution to stationarity jacobian:
    // jac_lag_stat_p_global += scaling * cost_grad_params_jac
    blasfeo_dgead(nu+nx, np_global, model->common->scaling, &work->cost_grad_params_jac, 0, 0, mem_common->jac_lag_stat_p_global, 0, 0);

    return;
}



void ocp_nlp_cost_external_compute_adj_sol_sens_pdiff(void *config_, void *dims_, void *model_,
                                       void *opts_, void *memory_, void *work_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_external_workspace *work = work_;
    ocp_nlp_cost_common_memory *mem_common = memory->common;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    struct blasfeo_dvec *ux = mem_common->ux;

    int nu = dims->nu;
    int np_global = dims->np_global;

    /* specify input types and pointers for external cost function */
    ext_fun_arg_t ext_fun_type_in[4];
    void *ext_fun_in[4];
    ext_fun_arg_t ext_fun_type_out[1];
    void *ext_fun_out[1];

    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = ux;
    u_in.xi = 0;

    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;
    ext_fun_type_in[2] = BLASFEO_DVEC;
    ext_fun_in[2] = mem_common->z_alg;

    // seed_ux
    ext_fun_type_in[3] = BLASFEO_DVEC;
    ext_fun_in[3] = mem_common->seed_ux;

    // OUTPUT
    ext_fun_type_out[0] = BLASFEO_DVEC;
    ext_fun_out[0] = &work->adj_cost_ux_pdiff;

    // evaluate external function
    if (model->ext_cost_adj_ux_pdiff == 0)
    {
        printf("ocp_nlp_cost_external_compute_jac_p: ext_cost_adj_ux_pdiff is not provided. Exiting.\n");
        exit(1);
    }
    model->ext_cost_adj_ux_pdiff->evaluate(model->ext_cost_adj_ux_pdiff, ext_fun_type_in, ext_fun_in,
                                  ext_fun_type_out, ext_fun_out);
    blasfeo_dvecad(np_global, model->common->scaling, &work->adj_cost_ux_pdiff, 0, mem_common->adj_lag_p_global, 0);
    return;
}

void ocp_nlp_cost_external_eval_grad_p(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_, struct blasfeo_dvec *out)
{
    ocp_nlp_cost_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_common_memory *mem_common = memory->common;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    struct blasfeo_dvec *ux = mem_common->ux;

    int nu = dims->nu;
    int np_global = dims->np_global;

    /* specify input types and pointers for external cost function */
    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[1];
    void *ext_fun_out[1];

    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = ux;
    u_in.xi = 0;

    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;
    ext_fun_type_in[2] = BLASFEO_DVEC;
    ext_fun_in[2] = mem_common->z_alg;

    // OUTPUT
    ext_fun_type_out[0] = BLASFEO_DVEC;
    ext_fun_out[0] = out;

    // evaluate external function
    model->ext_cost_grad_p->evaluate(model->ext_cost_grad_p, ext_fun_type_in, ext_fun_in,
                                  ext_fun_type_out, ext_fun_out);

    // scale
    if(model->common->scaling != 1.0)
    {
        blasfeo_dvecsc(np_global, model->common->scaling, out, 0);
    }

    return;
}


size_t ocp_nlp_cost_external_get_external_fun_workspace_requirement(void *config_, void *dims_, void *opts_, void *model_)
{
    ocp_nlp_cost_external_model *model = model_;

    size_t size = 0;
    size_t tmp_size;

    tmp_size = external_function_get_workspace_requirement_if_defined(model->ext_cost_fun);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->ext_cost_fun_jac);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->ext_cost_fun_jac_hess);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->ext_cost_grad_p);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->ext_cost_hess_xu_p);
    size = size > tmp_size ? size : tmp_size;
    tmp_size = external_function_get_workspace_requirement_if_defined(model->ext_cost_adj_ux_pdiff);
    size = size > tmp_size ? size : tmp_size;

    return size;
}


void ocp_nlp_cost_external_set_external_fun_workspaces(void *config_, void *dims_, void *opts_, void *model_, void *workspace_)
{
    ocp_nlp_cost_external_model *model = model_;
    external_function_set_fun_workspace_if_defined(model->ext_cost_fun, workspace_);
    external_function_set_fun_workspace_if_defined(model->ext_cost_fun_jac, workspace_);
    external_function_set_fun_workspace_if_defined(model->ext_cost_fun_jac_hess, workspace_);
    external_function_set_fun_workspace_if_defined(model->ext_cost_grad_p, workspace_);
    external_function_set_fun_workspace_if_defined(model->ext_cost_hess_xu_p, workspace_);
    external_function_set_fun_workspace_if_defined(model->ext_cost_adj_ux_pdiff, workspace_);
}



/* config */

void ocp_nlp_cost_external_config_initialize_default(void *config_, int stage)
{
    ocp_nlp_cost_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_cost_dims_calculate_size;
    config->dims_assign = &ocp_nlp_cost_dims_assign;
    config->dims_set = &ocp_nlp_cost_dims_set;
    config->dims_get = &ocp_nlp_cost_dims_get;
    config->model_calculate_size = &ocp_nlp_cost_external_model_calculate_size;
    config->model_assign = &ocp_nlp_cost_external_model_assign;
    config->model_set = &ocp_nlp_cost_external_model_set;
    config->model_get = &ocp_nlp_cost_external_model_get;
    config->opts_calculate_size = &ocp_nlp_cost_common_opts_calculate_size;
    config->opts_assign = &ocp_nlp_cost_common_opts_assign;
    config->opts_initialize_default = &ocp_nlp_cost_common_opts_initialize_default;
    config->opts_update = &ocp_nlp_cost_external_opts_update;
    config->opts_set = &ocp_nlp_cost_common_opts_set;
    config->opts_get = &ocp_nlp_cost_common_opts_get;
    config->memory_calculate_size = &ocp_nlp_cost_external_memory_calculate_size;
    config->memory_assign = &ocp_nlp_cost_external_memory_assign;
    config->memory_get = &ocp_nlp_cost_external_memory_get;
    config->memory_set = &ocp_nlp_cost_external_memory_set;
    config->workspace_calculate_size = &ocp_nlp_cost_external_workspace_calculate_size;
    config->get_external_fun_workspace_requirement = &ocp_nlp_cost_external_get_external_fun_workspace_requirement;
    config->set_external_fun_workspaces = &ocp_nlp_cost_external_set_external_fun_workspaces;
    config->initialize = &ocp_nlp_cost_external_initialize;
    config->update_qp_matrices = &ocp_nlp_cost_external_update_qp_matrices;
    config->compute_fun = &ocp_nlp_cost_external_compute_fun;
    config->compute_jac_p = &ocp_nlp_cost_external_compute_jac_p;
    config->compute_adj_sol_sens_pdiff = &ocp_nlp_cost_external_compute_adj_sol_sens_pdiff;
    config->compute_gradient = &ocp_nlp_cost_external_compute_gradient;
    config->eval_grad_p = &ocp_nlp_cost_external_eval_grad_p;
    config->config_initialize_default = &ocp_nlp_cost_external_config_initialize_default;
    config->precompute = &ocp_nlp_cost_external_precompute;
    config->stage = stage;

    return;
}
