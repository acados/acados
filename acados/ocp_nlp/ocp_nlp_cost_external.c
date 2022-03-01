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


#include "acados/ocp_nlp/ocp_nlp_cost_external.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_nlp_cost_external_dims_calculate_size(void *config_)
{
    acados_size_t size = sizeof(ocp_nlp_cost_external_dims);

    return size;
}



void *ocp_nlp_cost_external_dims_assign(void *config_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_external_dims *dims = (ocp_nlp_cost_external_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_dims);

    assert((char *) raw_memory + ocp_nlp_cost_external_dims_calculate_size(config_) >= c_ptr);

    return dims;
}



void ocp_nlp_cost_external_dims_initialize(void *config_, void *dims_, int nx, int nu, int ny, int ns, int nz)
{
    ocp_nlp_cost_external_dims *dims = dims_;

    dims->nx = nx;
    dims->nu = nu;
    dims->ns = ns;

    // TODO nz

    return;
}



void ocp_nlp_cost_external_dims_set(void *config_, void *dims_, const char *field, int* value)
{
    ocp_nlp_cost_external_dims *dims = (ocp_nlp_cost_external_dims *) dims_;

    if (!strcmp(field, "nx"))
    {
        dims->nx = *value;
    }
    else if (!strcmp(field, "nz"))
    {
        // do nothing
        // TODO(oj): implement constraints with daes
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu = *value;
    }
    else if (!strcmp(field, "ns"))
    {
        dims->ns = *value;
    }
    else
    {
        printf("\nerror: ocp_nlp_cost_external_dims_set: dimension type %s not available.\n", field);
        exit(1);
    }

    return;
}



void ocp_nlp_cost_external_dims_get(void *config_, void *dims_, const char *field, int* value)
{
        printf("error: ocp_nlp_cost_external_dims_get: attempt to get dimensions of non-existing field %s\n", field);
        exit(1);
}



/************************************************
 * model
 ************************************************/

acados_size_t ocp_nlp_cost_external_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_cost_external_dims *dims = dims_;

    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_external_model);

    size += 1 * 64;  // blasfeo_mem align
    size += blasfeo_memsize_dmat(nx+nu, nx+nu);

    size += 2 * blasfeo_memsize_dvec(2 * ns);  // Z, z

    return size;
}



void *ocp_nlp_cost_external_model_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_cost_external_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    // struct
    ocp_nlp_cost_external_model *model = (ocp_nlp_cost_external_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_model);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);
    // numerical_hessian
    assign_and_advance_blasfeo_dmat_mem(nx+nu, nx+nu, &model->numerical_hessian, &c_ptr);

    // blasfeo_dvec
    // Z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->Z, &c_ptr);
    // z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->z, &c_ptr);

    // default initialization
    model->scaling = 1.0;

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

    ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;

    int ns = dims->ns;
    int nx = dims->nx;
    int nu = dims->nu;

    if (!strcmp(field, "ext_cost_fun"))
    {
        model->ext_cost_fun = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_fun_jac_hes") | !strcmp(field, "ext_cost_fun_jac_hess"))
    {
        model->ext_cost_fun_jac_hess = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_fun_jac"))
    {
        model->ext_cost_fun_jac = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "ext_cost_num_hess"))
    {
        double *numerical_hessian = (double *) value_;
        blasfeo_pack_dmat(nx+nu, nx+nu, numerical_hessian, nx+nu, &model->numerical_hessian, 0, 0);
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
    else if (!strcmp(field, "scaling"))
    {
        double *scaling_ptr = (double *) value_;
        model->scaling = *scaling_ptr;
    }
    else
    {
        printf("\nerror: model entry: %s not available in module ocp_nlp_cost_external\n", field);
        exit(1);
    }
    return status;
}



/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_cost_external_opts_calculate_size(void *config_, void *dims_)
{
    // ocp_nlp_cost_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_external_opts);
    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_cost_external_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    // ocp_nlp_cost_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_cost_external_opts *opts = (ocp_nlp_cost_external_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_opts);

    assert((char *) raw_memory + ocp_nlp_cost_external_opts_calculate_size(config_, dims_) >=
           c_ptr);

    return opts;
}



void ocp_nlp_cost_external_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_external_opts *opts = opts_;

    opts->use_numerical_hessian = 0;

    return;
}



void ocp_nlp_cost_external_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    // ocp_nlp_cost_external_opts *opts = opts_;

    // opts->gauss_newton_hess = 1;

    return;
}



void ocp_nlp_cost_external_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_external_opts *opts = opts_;

    if(!strcmp(field, "exact_hess"))
    {
        // do nothing: the exact hessian is always computed if no custom hessian is provided
    }
    else if(!strcmp(field, "numerical_hessian"))
    {
        int *opt_val = (int *) value;
        opts->use_numerical_hessian = *opt_val;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_cost_external_opts_set\n", field);
        exit(1);
    }

    return;

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_cost_external_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_external_dims *dims = dims_;
    // ocp_nlp_cost_external_opts *opts = opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_external_memory);

    size += 1 * blasfeo_memsize_dvec(nu + nx + 2 * ns);  // grad

    size += 64;  // blasfeo_mem align

    return size;
}



void *ocp_nlp_cost_external_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_external_dims *dims = dims_;
    // ocp_nlp_cost_external_opts *opts = opts_;

    char *c_ptr = (char *) raw_memory;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    // struct
    ocp_nlp_cost_external_memory *memory = (ocp_nlp_cost_external_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_memory);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // grad
    assign_and_advance_blasfeo_dvec_mem(nu + nx + 2 * ns, &memory->grad, &c_ptr);

    assert((char *) raw_memory +
               ocp_nlp_cost_external_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return memory;
}



double *ocp_nlp_cost_external_memory_get_fun_ptr(void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    return &memory->fun;
}



struct blasfeo_dvec *ocp_nlp_cost_external_memory_get_grad_ptr(void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    return &memory->grad;
}



void ocp_nlp_cost_external_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    memory->RSQrq = RSQrq;

    return;
}



void ocp_nlp_cost_external_memory_set_Z_ptr(struct blasfeo_dvec *Z, void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    memory->Z = Z;
}



void ocp_nlp_cost_external_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    memory->ux = ux;

    return;
}



void ocp_nlp_cost_external_memory_set_tmp_ux_ptr(struct blasfeo_dvec *tmp_ux, void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    memory->tmp_ux = tmp_ux;

    return;
}



void ocp_nlp_cost_external_memory_set_z_alg_ptr(struct blasfeo_dvec *z_alg, void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    memory->z_alg = z_alg;
}



void ocp_nlp_cost_external_memory_set_dzdux_tran_ptr(struct blasfeo_dmat *dzdux_tran, void *memory_)
{
    ocp_nlp_cost_external_memory *memory = memory_;

    memory->dzdux_tran = dzdux_tran;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_cost_external_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_cost_external_dims *dims = dims_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_cost_external_workspace);

    size += 1 * blasfeo_memsize_dmat(nu+nx, nu+nx);  // tmp_nv_nv

    size += 1 * blasfeo_memsize_dvec(2*ns);  // tmp_2ns

    size += 64;  // blasfeo_mem align
    
    return size;
}



static void ocp_nlp_cost_external_cast_workspace(void *config_, void *dims_, void *opts_,
                                                 void *work_)
{
    ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_workspace *work = work_;
    
    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;
    
    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_cost_external_workspace);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);
    
    // tmp_nv_nv
    assign_and_advance_blasfeo_dmat_mem(nu + nx, nu + nx, &work->tmp_nv_nv, &c_ptr);

    // tmp_2ns
    assign_and_advance_blasfeo_dvec_mem(2*ns, &work->tmp_2ns, &c_ptr);

    assert((char *) work_ + ocp_nlp_cost_external_workspace_calculate_size(config_, dims_, opts_) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_cost_external_initialize(void *config_, void *dims_, void *model_, void *opts_,
                                      void *memory_, void *work_)
{
    ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    ocp_nlp_cost_external_memory *memory = memory_;
    // ocp_nlp_cost_external_workspace *work= work_;

    // ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    int ns = dims->ns;

    blasfeo_dveccpsc(2*ns, model->scaling, &model->Z, 0, memory->Z, 0);

    return;
}



void ocp_nlp_cost_external_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_,
                                              void *memory_, void *work_)
{
    ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    ocp_nlp_cost_external_opts *opts = opts_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_external_workspace *work = work_;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    /* specify input types and pointers for external cost function */
    // TODO(oj): add z
    ext_fun_arg_t ext_fun_type_in[2];
    void *ext_fun_in[2];
    ext_fun_arg_t ext_fun_type_out[3];
    void *ext_fun_out[3];


    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = memory->ux;
    u_in.xi = 0;
    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = memory->ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;

    // OUTPUT
    ext_fun_type_out[0] = COLMAJ;
    ext_fun_out[0] = &memory->fun;  // fun: scalar
    ext_fun_type_out[1] = BLASFEO_DVEC;
    ext_fun_out[1] = &memory->grad;  // grad: nu+nx

    if (opts->use_numerical_hessian > 0)
    {
        // evaluate external function
        model->ext_cost_fun_jac->evaluate(model->ext_cost_fun_jac, ext_fun_type_in,
                                            ext_fun_in, ext_fun_type_out, ext_fun_out);
        // custom hessian
        blasfeo_dgead(nx+nu, nx+nu, model->scaling, &model->numerical_hessian, 0, 0, memory->RSQrq, 0, 0);
    }
    else
    {
        // additional output
        ext_fun_type_out[2] = BLASFEO_DMAT;
        ext_fun_out[2] = &work->tmp_nv_nv;   // hess: (nu+nx) * (nu+nx)
        // evaluate external function
        model->ext_cost_fun_jac_hess->evaluate(model->ext_cost_fun_jac_hess, ext_fun_type_in,
                                            ext_fun_in, ext_fun_type_out, ext_fun_out);
        // hessian contribution
        blasfeo_dgead(nx+nu, nx+nu, model->scaling, &work->tmp_nv_nv, 0, 0, memory->RSQrq, 0, 0);
    }

    // slack update gradient
    blasfeo_dveccp(2*ns, &model->z, 0, &memory->grad, nu+nx);
    blasfeo_dvecmulacc(2*ns, &model->Z, 0, memory->ux, nu+nx, &memory->grad, nu+nx);

    // slack update function value
    blasfeo_dveccpsc(2*ns, 2.0, &model->z, 0, &work->tmp_2ns, 0);
    blasfeo_dvecmulacc(2*ns, &model->Z, 0, memory->ux, nu+nx, &work->tmp_2ns, 0);
    memory->fun += 0.5 * blasfeo_ddot(2*ns, &work->tmp_2ns, 0, memory->ux, nu+nx);

    // scale
    if (model->scaling!=1.0)
    {
        blasfeo_dvecsc(nu+nx+2*ns, model->scaling, &memory->grad, 0);
        memory->fun *= model->scaling;
    }

    // blasfeo_print_dmat(nu+nx, nu+nx, memory->RSQrq, 0, 0);
    // blasfeo_print_tran_dvec(2*ns, memory->Z, 0);
    // blasfeo_print_tran_dvec(nu+nx+2*ns, &memory->grad, 0);

    return;
}



void ocp_nlp_cost_external_compute_fun(void *config_, void *dims_, void *model_,
                                       void *opts_, void *memory_, void *work_)
{
//    printf("\nerror: ocp_external_cost_nls_compute_fun: not implemented yet\n");
//    exit(1);

    ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    // ocp_nlp_cost_external_opts *opts = opts_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_external_workspace *work = work_;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    /* specify input types and pointers for external cost function */
    // TODO(oj): add z
    ext_fun_arg_t ext_fun_type_in[2];
    void *ext_fun_in[2];
    ext_fun_arg_t ext_fun_type_out[1];
    void *ext_fun_out[1];

    // INPUT
    struct blasfeo_dvec_args u_in;  // input u
    u_in.x = memory->tmp_ux;
    u_in.xi = 0;

    struct blasfeo_dvec_args x_in;  // input x
    x_in.x = memory->tmp_ux;
    x_in.xi = nu;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;
    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;

    // OUTPUT
    ext_fun_type_out[0] = COLMAJ;
    ext_fun_out[0] = &memory->fun;  // function: scalar

    // evaluate external function
    if (model->ext_cost_fun == 0)
    {
        printf("ocp_nlp_cost_external_compute_fun: ext_cost_fun is not provided. Exiting.\n");
        exit(1);
    }
    model->ext_cost_fun->evaluate(model->ext_cost_fun, ext_fun_type_in, ext_fun_in,
                                  ext_fun_type_out, ext_fun_out);

    // slack update function value
    blasfeo_dveccpsc(2*ns, 2.0, &model->z, 0, &work->tmp_2ns, 0);
    blasfeo_dvecmulacc(2*ns, &model->Z, 0, memory->tmp_ux, nu+nx, &work->tmp_2ns, 0);
    memory->fun += 0.5 * blasfeo_ddot(2*ns, &work->tmp_2ns, 0, memory->tmp_ux, nu+nx);

    // scale
    if(model->scaling!=1.0)
    {
        memory->fun *= model->scaling;
    }

    return;
}



/* config */

void ocp_nlp_cost_external_config_initialize_default(void *config_)
{
    ocp_nlp_cost_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_cost_external_dims_calculate_size;
    config->dims_assign = &ocp_nlp_cost_external_dims_assign;
    config->dims_initialize = &ocp_nlp_cost_external_dims_initialize;
    config->dims_set = &ocp_nlp_cost_external_dims_set;
    config->dims_get = &ocp_nlp_cost_external_dims_get;
    config->model_calculate_size = &ocp_nlp_cost_external_model_calculate_size;
    config->model_assign = &ocp_nlp_cost_external_model_assign;
    config->model_set = &ocp_nlp_cost_external_model_set;
    config->opts_calculate_size = &ocp_nlp_cost_external_opts_calculate_size;
    config->opts_assign = &ocp_nlp_cost_external_opts_assign;
    config->opts_initialize_default = &ocp_nlp_cost_external_opts_initialize_default;
    config->opts_update = &ocp_nlp_cost_external_opts_update;
    config->opts_set = &ocp_nlp_cost_external_opts_set;
    config->memory_calculate_size = &ocp_nlp_cost_external_memory_calculate_size;
    config->memory_assign = &ocp_nlp_cost_external_memory_assign;
    config->memory_get_fun_ptr = &ocp_nlp_cost_external_memory_get_fun_ptr;
    config->memory_get_grad_ptr = &ocp_nlp_cost_external_memory_get_grad_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_cost_external_memory_set_ux_ptr;
    config->memory_set_tmp_ux_ptr = &ocp_nlp_cost_external_memory_set_tmp_ux_ptr;
    config->memory_set_z_alg_ptr = &ocp_nlp_cost_external_memory_set_z_alg_ptr;
    config->memory_set_dzdux_tran_ptr = &ocp_nlp_cost_external_memory_set_dzdux_tran_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_cost_external_memory_set_RSQrq_ptr;
    config->memory_set_Z_ptr = &ocp_nlp_cost_external_memory_set_Z_ptr;
    config->workspace_calculate_size = &ocp_nlp_cost_external_workspace_calculate_size;
    config->initialize = &ocp_nlp_cost_external_initialize;
    config->update_qp_matrices = &ocp_nlp_cost_external_update_qp_matrices;
    config->compute_fun = &ocp_nlp_cost_external_compute_fun;
    config->config_initialize_default = &ocp_nlp_cost_external_config_initialize_default;

    return;
}
