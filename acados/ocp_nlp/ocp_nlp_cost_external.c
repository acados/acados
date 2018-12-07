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

/************************************************
 * dims
 ************************************************/

int ocp_nlp_cost_external_dims_calculate_size(void *config_)
{
    int size = sizeof(ocp_nlp_cost_external_dims);

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

void ocp_nlp_cost_external_dims_initialize(void *config_, void *dims_, int nx, int nu, int ny,
                                           int ns)
{
    ocp_nlp_cost_external_dims *dims = dims_;

    dims->nx = nx;
    dims->nu = nu;
    dims->ns = ns;

    return;
}

static void ocp_nlp_cost_external_set_nx(void *config_, void *dims_, int *nx)
{
    ocp_nlp_cost_external_dims *dims = (ocp_nlp_cost_external_dims *) dims_;
    dims->nx = *nx;
}


static void ocp_nlp_cost_external_set_nu(void *config_, void *dims_, int *nu)
{
    ocp_nlp_cost_external_dims *dims = (ocp_nlp_cost_external_dims *) dims_;
    dims->nu = *nu;
}


static void ocp_nlp_cost_external_set_ns(void *config_, void *dims_, int *ns)
{
    ocp_nlp_cost_external_dims *dims = (ocp_nlp_cost_external_dims *) dims_;
    dims->ns = *ns;
}



void ocp_nlp_cost_external_dims_set(void *config_, void *dims_, const char *field, int* value)
{
    if (!strcmp(field, "nx"))
    {
        ocp_nlp_cost_external_set_nx(config_, dims_, value);
    }
    else if (!strcmp(field, "nz"))
    {
        // do nothing
        // TODO(oj): implement constraints with daes
    }
    else if (!strcmp(field, "nu"))
    {
        ocp_nlp_cost_external_set_nu(config_, dims_, value);
    }
    else if (!strcmp(field, "ns"))
    {
        ocp_nlp_cost_external_set_ns(config_, dims_, value);
    }
    else
    {
        printf("\nerror: dimension type not available in module\n");
        exit(1);
    }
}


/************************************************
 * model
 ************************************************/

int ocp_nlp_cost_external_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_cost_external_dims *dims = dims_;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;
    // int ny = dims->ny;
    int ns = dims->ns;

    int size = 0;

    size += sizeof(ocp_nlp_cost_external_model);

    size += 1 * 64;  // blasfeo_mem align

    size += 2 * blasfeo_memsize_dvec(2 * ns);  // Z, z

    return size;
}

void *ocp_nlp_cost_external_model_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_cost_external_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

    // int nx = dims->nx;
    // int nu = dims->nu;
    // int ny = dims->ny;
    int ns = dims->ns;

    // struct
    ocp_nlp_cost_external_model *model = (ocp_nlp_cost_external_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_model);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dvec
    // Z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->Z, &c_ptr);
    // z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->z, &c_ptr);

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
        status = ACADOS_FAILURE;

    ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;

    if (!strcmp(field, "ext_cost"))
    {
        model->ext_cost = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "Z "))
    {
        double *Z = (double *) value_;
        blasfeo_pack_dvec(2 * dims->ns, Z, &model->Z, 0);
    }
    else if (!strcmp(field, "z"))
    {
        double *z = (double *) value_;
        blasfeo_pack_dvec(2 * dims->ns, z, &model->z, 0);
    }
    else
    {
        printf("\nerror: model entry: %s not available in module ocp_nlp_cost_nls\n", field);
        status = ACADOS_FAILURE;
    }
    return status;
}

/************************************************
 * options
 ************************************************/

int ocp_nlp_cost_external_opts_calculate_size(void *config_, void *dims_)
{
    // ocp_nlp_cost_config *config = config_;

    int size = 0;

    size += sizeof(ocp_nlp_cost_external_opts);

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
    // ocp_nlp_cost_external_opts *opts = opts_;

    // opts->gauss_newton_hess = 1;

    return;
}

void ocp_nlp_cost_external_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    // ocp_nlp_cost_external_opts *opts = opts_;

    // opts->gauss_newton_hess = 1;

    return;
}

/************************************************
 * memory
 ************************************************/

int ocp_nlp_cost_external_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_external_dims *dims = dims_;
    // ocp_nlp_cost_external_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    // int ny = dims->ny;
    int ns = dims->ns;

    int size = 0;

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
    // int ny = dims->ny;
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

/************************************************
 * workspace
 ************************************************/

int ocp_nlp_cost_external_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_external_dims *dims = dims_;
    // ocp_nlp_cost_external_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;

    int size = 0;

    size += sizeof(ocp_nlp_cost_external_workspace);

    size += 1 * (nu + nx) * sizeof(double);                            // ext_cost_in
    size += 1 * ((nu + nx) + (nu + nx) * (nu + nx)) * sizeof(double);  // ext_cost_out

    size += 8;

    return size;
}

static void ocp_nlp_cost_external_cast_workspace(void *config_, void *dims_, void *opts_,
                                                 void *work_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_external_dims *dims = dims_;
    // ocp_nlp_cost_external_opts *opts = opts_;
    ocp_nlp_cost_external_workspace *work = work_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_cost_external_workspace);

    align_char_to(8, &c_ptr);

    // external_jac_in
    assign_and_advance_double(nu + nx, &work->ext_cost_in, &c_ptr);
    // external_jac_out
    assign_and_advance_double((nu + nx) + (nu + nx) * (nu + nx), &work->ext_cost_out, &c_ptr);

    assert((char *) work + ocp_nlp_cost_external_workspace_calculate_size(config_, dims, opts_) >=
           c_ptr);

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

    blasfeo_dveccp(2 * ns, &model->Z, 0, memory->Z, 0);

    return;
}

void ocp_nlp_cost_external_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_,
                                              void *memory_, void *work_)
{
    ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    // ocp_nlp_cost_external_opts *opts = opts_;
    ocp_nlp_cost_external_memory *memory = memory_;
    ocp_nlp_cost_external_workspace *work = work_;

    ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

    int nx = dims->nx;
    int nu = dims->nu;
    int ns = dims->ns;

    // XXX large enough ?
    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[3];
    void *ext_fun_out[3];

    // unpack ls cost input
    blasfeo_unpack_dvec(nu + nx, memory->ux, 0, work->ext_cost_in);

    ext_fun_type_in[0] = COLMAJ;
    ext_fun_in[0] = work->ext_cost_in + 0;  // ux: nu+nx

    ext_fun_type_out[0] = COLMAJ;
    ext_fun_out[0] = work->ext_cost_out + 0;  // grad: nu+nx
    ext_fun_type_out[1] = COLMAJ;
    ext_fun_out[1] = work->ext_cost_out + nx + nu;  // hess: (nu+nx) * (nu+nx)

    // evaluate external function
    model->ext_cost->evaluate(model->ext_cost, ext_fun_type_in, ext_fun_in, ext_fun_type_out,
                              ext_fun_out);

    // pack gradient
    blasfeo_pack_dvec(nu + nx, ext_fun_out[0], &memory->grad, 0);  // rq

    // pack Hessian
    blasfeo_pack_dmat(nu + nx, nu + nx, ext_fun_out[1], nx + nu, memory->RSQrq, 0, 0);  // RSQrq

    // slacks
    blasfeo_dveccp(2 * ns, &model->z, 0, &memory->grad, nu + nx);
    blasfeo_dvecmulacc(2 * ns, &model->Z, 0, memory->ux, nu + nx, &memory->grad, nu + nx);

    // blasfeo_print_dmat(nu+nx, nu+nx, memory->RSQrq, 0, 0);
    // blasfeo_print_tran_dvec(2*ns, memory->Z, 0);
    // blasfeo_print_tran_dvec(nu+nx+2*ns, &memory->grad, 0);
    // exit(1);

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
    config->model_calculate_size = &ocp_nlp_cost_external_model_calculate_size;
    config->model_assign = &ocp_nlp_cost_external_model_assign;
    config->model_set = &ocp_nlp_cost_external_model_set;
    config->opts_calculate_size = &ocp_nlp_cost_external_opts_calculate_size;
    config->opts_assign = &ocp_nlp_cost_external_opts_assign;
    config->opts_initialize_default = &ocp_nlp_cost_external_opts_initialize_default;
    config->opts_update = &ocp_nlp_cost_external_opts_update;
    config->memory_calculate_size = &ocp_nlp_cost_external_memory_calculate_size;
    config->memory_assign = &ocp_nlp_cost_external_memory_assign;
    config->memory_get_grad_ptr = &ocp_nlp_cost_external_memory_get_grad_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_cost_external_memory_set_ux_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_cost_external_memory_set_RSQrq_ptr;
    config->memory_set_Z_ptr = &ocp_nlp_cost_external_memory_set_Z_ptr;
    config->workspace_calculate_size = &ocp_nlp_cost_external_workspace_calculate_size;
    config->initialize = &ocp_nlp_cost_external_initialize;
    config->update_qp_matrices = &ocp_nlp_cost_external_update_qp_matrices;
    config->config_initialize_default = &ocp_nlp_cost_external_config_initialize_default;

    return;
}
