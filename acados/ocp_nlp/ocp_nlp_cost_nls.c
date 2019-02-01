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

int ocp_nlp_cost_nls_dims_calculate_size(void *config_)
{
    int size = sizeof(ocp_nlp_cost_nls_dims);

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



void ocp_nlp_cost_nls_dims_initialize(void *config_, void *dims_, int nx, int nu, int ny, int ns)
{
    ocp_nlp_cost_nls_dims *dims = dims_;

    dims->nx = nx;
    dims->nu = nu;
    dims->ny = ny;
    dims->ns = ns;

    return;
}



static void ocp_nlp_cost_nls_set_nx(void *config_, void *dims_, int *nx)
{
    ocp_nlp_cost_nls_dims *dims = (ocp_nlp_cost_nls_dims *) dims_;
    dims->nx = *nx;
}



static void ocp_nlp_cost_nls_set_nu(void *config_, void *dims_, int *nu)
{
    ocp_nlp_cost_nls_dims *dims = (ocp_nlp_cost_nls_dims *) dims_;
    dims->nu = *nu;
}



static void ocp_nlp_cost_nls_set_ny(void *config_, void *dims_, int *ny)
{
    ocp_nlp_cost_nls_dims *dims = (ocp_nlp_cost_nls_dims *) dims_;
    dims->ny = *ny;
}



static void ocp_nlp_cost_nls_set_ns(void *config_, void *dims_, int *ns)
{
    ocp_nlp_cost_nls_dims *dims = (ocp_nlp_cost_nls_dims *) dims_;
    dims->ns = *ns;
}



void ocp_nlp_cost_nls_dims_set(void *config_, void *dims_, const char *field, int* value)
{
    if (!strcmp(field, "nx"))
    {
        ocp_nlp_cost_nls_set_nx(config_, dims_, value);
    }
    else if (!strcmp(field, "nz"))
    {
        // do nothing
        // TODO(oj): implement cost with daes
    }
    else if (!strcmp(field, "nu"))
    {
        ocp_nlp_cost_nls_set_nu(config_, dims_, value);
    }
    else if (!strcmp(field, "ny"))
    {
        ocp_nlp_cost_nls_set_ny(config_, dims_, value);
    }
    else if (!strcmp(field, "ns"))
    {
        ocp_nlp_cost_nls_set_ns(config_, dims_, value);
    }
    else
    {
        printf("\nerror: dimension type: %s not available in module\n", field);
        exit(1);
    }
}



/************************************************
 * model
 ************************************************/

int ocp_nlp_cost_nls_model_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;

    // extract dims
    // int nx = dims->nx;
    // int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    int size = 0;

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

    model->nls_hess = NULL;
    model->nls_res_jac  = NULL;

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // blasfeo_dmat
    // W
    assign_and_advance_blasfeo_dmat_mem(ny, ny, &model->W, &c_ptr);

    // blasfeo_dvec
    // y_ref
    assign_and_advance_blasfeo_dvec_mem(ny, &model->y_ref, &c_ptr);
    // Z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->Z, &c_ptr);
    // z
    assign_and_advance_blasfeo_dvec_mem(2 * ns, &model->z, &c_ptr);

	// default initialization
	model->scaling = 1.0;

    // assert
    assert((char *) raw_memory + ocp_nlp_cost_nls_model_calculate_size(config_, dims) >= c_ptr);

    return model;
}



int ocp_nlp_cost_nls_model_set(void *config_, void *dims_, void *model_,
                                         const char *field, void *value_)
{
    int status = ACADOS_SUCCESS;

    if ( !config_ || !dims_ || !model_ || !value_ )
        status = ACADOS_FAILURE;

    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;

    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    if (!strcmp(field, "W"))
    {
        double *W_col_maj = (double *) value_;
        blasfeo_pack_dmat(ny, ny, W_col_maj, ny, &model->W, 0, 0);
    }
    else if (!strcmp(field, "y_ref") || !strcmp(field, "yref"))
    {
        double *y_ref = (double *) value_;
        blasfeo_pack_dvec(ny, y_ref, &model->y_ref, 0);
    }
    else if (!strcmp(field, "Z"))
    {
        double *Z = (double *) value_;
        blasfeo_pack_dvec(ns, Z, &model->Z, 0);
        blasfeo_pack_dvec(ns, Z, &model->Z, ns);
    }
    else if (!strcmp(field, "Zl"))
    {
        double *Zl = (double *) value_;
        blasfeo_pack_dvec(ns, Zl, &model->Z, 0);
    }
    else if (!strcmp(field, "Zu"))
    {
        double *Zu = (double *) value_;
        blasfeo_pack_dvec(ns, Zu, &model->Z, ns);
    }
    else if (!strcmp(field, "z"))
    {
        double *z = (double *) value_;
        blasfeo_pack_dvec(ns, z, &model->z, 0);
        blasfeo_pack_dvec(ns, z, &model->z, ns);
    }
    else if (!strcmp(field, "zl"))
    {
        double *zl = (double *) value_;
        blasfeo_pack_dvec(ns, zl, &model->z, 0);
    }
    else if (!strcmp(field, "zu"))
    {
        double *zu = (double *) value_;
        blasfeo_pack_dvec(ns, zu, &model->z, ns);
    }
    else if (!strcmp(field, "nls_res_jac"))
    {
        model->nls_res_jac = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "nls_hess"))
    {
        model->nls_hess = (external_function_generic *) value_;
    }
    else if (!strcmp(field, "scaling"))
    {
        double *scaling_ptr = (double *) value_;
        model->scaling = *scaling_ptr;
    }
    else
    {
        printf("\nerror: model entry: %s not available in module ocp_nlp_cost_nls\n", field);
		exit(1);
//        status = ACADOS_FAILURE;
    }
    return status;
}



/************************************************
 * options
 ************************************************/

int ocp_nlp_cost_nls_opts_calculate_size(void *config_, void *dims_)
{
    // ocp_nlp_cost_config *config = config_;

    int size = 0;

    size += sizeof(ocp_nlp_cost_nls_opts);

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



/************************************************
 * memory
 ************************************************/

int ocp_nlp_cost_nls_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_dims *dims = dims_;
    // ocp_nlp_cost_nls_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    int size = 0;

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



/************************************************
 * workspace
 ************************************************/

int ocp_nlp_cost_nls_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_dims *dims = dims_;
    // ocp_nlp_cost_nls_opts *opts = opts_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;

    int size = 0;

    size += sizeof(ocp_nlp_cost_nls_workspace);

    size += 1 * blasfeo_memsize_dmat(nu + nx, ny);  // tmp_nv_ny
    size += 1 * blasfeo_memsize_dvec(ny);           // tmp_ny

    size += 64;  // blasfeo_mem align
    size += 8;
    return size;
}



static void ocp_nlp_cost_nls_cast_workspace(void *config_, void *dims_, void *opts_, void *work_)
{
    // ocp_nlp_cost_config *config = config_;
    ocp_nlp_cost_nls_dims *dims = dims_;
    // ocp_nlp_cost_nls_opts *opts = opts_;
    ocp_nlp_cost_nls_workspace *work = work_;

    // extract dims
    int nx = dims->nx;
    int nu = dims->nu;
    int ny = dims->ny;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_cost_nls_workspace);

    align_char_to(8, &c_ptr);

    // blasfeo_mem align
    align_char_to(64, &c_ptr);

    // tmp_nv_ny
    assign_and_advance_blasfeo_dmat_mem(nu + nx, ny, &work->tmp_nv_ny, &c_ptr);

    // tmp_ny
    assign_and_advance_blasfeo_dvec_mem(ny, &work->tmp_ny, &c_ptr);

    assert((char *) work + ocp_nlp_cost_nls_workspace_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

// TODO move factorization of W into pre-compute???
void ocp_nlp_cost_nls_initialize(void *config_, void *dims_, void *model_, void *opts_,
                                 void *memory_, void *work_)
{
    ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;
    ocp_nlp_cost_nls_memory *memory = memory_;
    // ocp_nlp_cost_nls_workspace *work= work_;

    ocp_nlp_cost_nls_cast_workspace(config_, dims, opts_, work_);

    // int nx = dims->nx;
    // int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    // TODO(all): recompute factorization only if W are re-tuned ???
    blasfeo_dpotrf_l(ny, &model->W, 0, 0, &memory->W_chol, 0, 0);

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
    int nu = dims->nu;
    int ny = dims->ny;
    int ns = dims->ns;

    ext_fun_arg_t ext_fun_type_in[3];
    void *ext_fun_in[3];
    ext_fun_arg_t ext_fun_type_out[3];
    void *ext_fun_out[3];

    struct blasfeo_dvec_args x_in;  // input x of external fun;
    struct blasfeo_dvec_args u_in;  // input u of external fun;

    x_in.x = memory->ux;
    u_in.x = memory->ux;

    x_in.xi = nu;
    u_in.xi = 0;

    ext_fun_type_in[0] = BLASFEO_DVEC_ARGS;
    ext_fun_in[0] = &x_in;

    ext_fun_type_in[1] = BLASFEO_DVEC_ARGS;
    ext_fun_in[1] = &u_in;

    ext_fun_type_out[0] = BLASFEO_DVEC;
    ext_fun_out[0] = &memory->res;  // fun: ny
    ext_fun_type_out[1] = BLASFEO_DMAT;
    ext_fun_out[1] = &memory->Jt;  // jac': (nu+nx) * ny

    // evaluate external function
    model->nls_res_jac->evaluate(model->nls_res_jac, ext_fun_type_in, ext_fun_in, ext_fun_type_out,
                             ext_fun_out);

    /* gradient */
    // res = res - y_ref
    blasfeo_daxpy(ny, -1.0, &model->y_ref, 0, &memory->res, 0, &memory->res, 0);

    // printf("W\n");
    // blasfeo_print_dmat(ny, ny, &model->W, 0, 0);

    // printf("res\n");
    // blasfeo_print_dvec(ny, &memory->res, 0);

    // TODO(all): use lower triangular chol of W to save n_y^2 flops
    // tmp_ny = W * res
    blasfeo_dsymv_l(ny, ny, 1.0, &model->W, 0, 0, &memory->res, 0, 0.0, &work->tmp_ny, 0,
                    &work->tmp_ny, 0);
    // grad = Jt * tmp_ny
    blasfeo_dgemv_n(nu + nx, ny, 1.0, &memory->Jt, 0, 0, &work->tmp_ny, 0, 0.0, &memory->grad, 0,
                    &memory->grad, 0);

    // printf("tmp_ny\n");
    // blasfeo_print_dvec(ny, &work->tmp_ny, 0);

    // printf("W_chol\n");
    // blasfeo_print_dmat(ny, ny, &memory->W_chol, 0, 0);

    // printf("Jt\n");
    // blasfeo_print_dmat(nu+nx, ny, &memory->Jt, 0, 0);


    /* hessian */

    if (opts->gauss_newton_hess)
    {
        // gauss-newton approximation of hessian of ls cost
        // tmp_nv_ny = Jt * W_chol,     where W_chol is lower triangular
        blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, &memory->W_chol, 0, 0, &memory->Jt, 0, 0,
                           &work->tmp_nv_ny, 0, 0);

        // blasfeo_print_dmat(nu + nx, ny, &work->tmp_nv_ny, 0, 0);

        // RSQrq = tmp_nv_ny * tmp_nv_ny^T
        blasfeo_dsyrk_ln(nu+nx, ny, model->scaling, &work->tmp_nv_ny, 0, 0, &work->tmp_nv_ny, 0, 0, 0.0,
                         memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);

        // blasfeo_print_dmat(nu+nx, nu+nx, memory->RSQrq, 0, 0);
    }
    else
    {
        // TODO(oj): test this
        // NOTE(oj): this should add the non-Gauss-Newton term to RSQrq, the product < r, d2_d[x,u] r >,
        //  where the cost is 0.5 * norm2(r(x,u))^2
        // exact hessian of ls cost

        // ext_fun_[type_]in 0,1 are the same as before.
        ext_fun_type_in[2] = BLASFEO_DVEC;
        ext_fun_in[2] = &work->tmp_ny;  // fun: ny

        ext_fun_type_out[0] = BLASFEO_DMAT;
        ext_fun_out[0] = memory->RSQrq;  // hess: (nu+nx) * (nu+nx)

        // evaluate external function
        model->nls_hess->evaluate(model->nls_hess, ext_fun_type_in, ext_fun_in, ext_fun_type_out,
                                  ext_fun_out);

        // gauss-newton component update
        // tmp_nv_ny = Jt * W_chol, where W_chol is lower triangular
        blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, &memory->W_chol, 0, 0, &memory->Jt, 0, 0,
                           &work->tmp_nv_ny, 0, 0);
        // RSQrq = RSQrq + tmp_nv_ny * tmp_nv_ny^T
        blasfeo_dsyrk_ln(nu+nx, ny, model->scaling, &work->tmp_nv_ny, 0, 0, &work->tmp_nv_ny, 0, 0,
                         model->scaling, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);
    }

    // slacks
    blasfeo_dveccp(2*ns, &model->z, 0, &memory->grad, nu+nx);
    blasfeo_dvecmulacc(2*ns, &model->Z, 0, memory->ux, nu+nx, &memory->grad, nu+nx);

	// scale
	if(model->scaling!=1.0)
	{
		blasfeo_dvecsc(nu+nx+2*ns, model->scaling, &memory->grad, 0);
	}

    // blasfeo_print_dmat(nu+nx, nu+nx, memory->RSQrq, 0, 0);
    // blasfeo_print_tran_dvec(2*ns, memory->Z, 0);
    // blasfeo_print_tran_dvec(nu+nx+2*ns, &memory->grad, 0);
    // exit(1);

    return;
}

void ocp_nlp_cost_nls_config_initialize_default(void *config_)
{
    ocp_nlp_cost_config *config = config_;

    config->dims_calculate_size = &ocp_nlp_cost_nls_dims_calculate_size;
    config->dims_assign = &ocp_nlp_cost_nls_dims_assign;
    config->dims_initialize = &ocp_nlp_cost_nls_dims_initialize;
    config->dims_set = &ocp_nlp_cost_nls_dims_set;
    config->model_calculate_size = &ocp_nlp_cost_nls_model_calculate_size;
    config->model_assign = &ocp_nlp_cost_nls_model_assign;
    config->model_set = &ocp_nlp_cost_nls_model_set;
    config->opts_calculate_size = &ocp_nlp_cost_nls_opts_calculate_size;
    config->opts_assign = &ocp_nlp_cost_nls_opts_assign;
    config->opts_initialize_default = &ocp_nlp_cost_nls_opts_initialize_default;
    config->opts_update = &ocp_nlp_cost_nls_opts_update;
    config->memory_calculate_size = &ocp_nlp_cost_nls_memory_calculate_size;
    config->memory_assign = &ocp_nlp_cost_nls_memory_assign;
    config->memory_get_grad_ptr = &ocp_nlp_cost_nls_memory_get_grad_ptr;
    config->memory_set_ux_ptr = &ocp_nlp_cost_nls_memory_set_ux_ptr;
    config->memory_set_RSQrq_ptr = &ocp_nlp_cost_nls_memory_set_RSQrq_ptr;
    config->memory_set_Z_ptr = &ocp_nlp_cost_nls_memory_set_Z_ptr;
    config->workspace_calculate_size = &ocp_nlp_cost_nls_workspace_calculate_size;
    config->initialize = &ocp_nlp_cost_nls_initialize;
    config->update_qp_matrices = &ocp_nlp_cost_nls_update_qp_matrices;
    config->config_initialize_default = &ocp_nlp_cost_nls_config_initialize_default;

    return;
}
