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

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
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



void ocp_nlp_cost_nls_dims_initialize(void *config_, void *dims_, int nx, int nu, int ny)
{
	ocp_nlp_cost_nls_dims *dims = dims_;

	dims->nx = nx;
	dims->nu = nu;
	dims->ny = ny;

	return;
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

	int size = 0;

	size += sizeof(ocp_nlp_cost_nls_model);

	size += 64; // blasfeo_mem align

	size += blasfeo_memsize_dmat(ny, ny); // W
	size += blasfeo_memsize_dvec(ny); // y_ref

	return size;

}



void *ocp_nlp_cost_nls_model_assign(void *config_, void *dims_, void *raw_memory)
{
	ocp_nlp_cost_nls_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;

	// int nx = dims->nx;
	// int nu = dims->nu;
	int ny = dims->ny;

	// struct
    ocp_nlp_cost_nls_model *model = (ocp_nlp_cost_nls_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_nls_model);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// blasfeo_dmat
	// W
	assign_and_advance_blasfeo_dmat_mem(ny, ny, &model->W, &c_ptr);

	// blasfeo_dvec
	// y_ref
	assign_and_advance_blasfeo_dvec_mem(ny, &model->y_ref, &c_ptr);

	// assert
    assert((char *) raw_memory + ocp_nlp_cost_nls_model_calculate_size(config_, dims) >= c_ptr);

	return model;

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

    assert((char*)raw_memory + ocp_nlp_cost_nls_opts_calculate_size(config_, dims_) >= c_ptr);

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
//	ocp_nlp_cost_nls_opts *opts = opts_;

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

	int size = 0;

    size += sizeof(ocp_nlp_cost_nls_memory);

	size += 1*blasfeo_memsize_dmat(ny, ny); // W_chol
	size += 1*blasfeo_memsize_dmat(nu+nx, ny); // Jt
	size += 1*blasfeo_memsize_dvec(ny); // res
	size += 1*blasfeo_memsize_dvec(nu+nx); // grad

	size += 64; // blasfeo_mem align

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

	// struct
    ocp_nlp_cost_nls_memory *memory = (ocp_nlp_cost_nls_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_nls_memory);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// W_chol
	assign_and_advance_blasfeo_dmat_mem(ny, ny, &memory->W_chol, &c_ptr);
	// Jt
	assign_and_advance_blasfeo_dmat_mem(nu+nx, ny, &memory->Jt, &c_ptr);
	// res
	assign_and_advance_blasfeo_dvec_mem(ny, &memory->res, &c_ptr);
	// grad
	assign_and_advance_blasfeo_dvec_mem(nu+nx, &memory->grad, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_cost_nls_memory_calculate_size(config_, dims, opts_) >= c_ptr);

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
	ocp_nlp_cost_nls_opts *opts = opts_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int ny = dims->ny;

	int size = 0;

    size += sizeof(ocp_nlp_cost_nls_workspace);

	size += 1*blasfeo_memsize_dmat(nu+nx, ny); // tmp_nv_ny
	size += 1*blasfeo_memsize_dvec(ny); // tmp_ny

	size += 1*(nu+nx)*sizeof(double); // nls_jac_in
	size += 1*(ny+ny*(nu+nx))*sizeof(double); // nls_jac_out
	if (!opts->gauss_newton_hess)
	{
		size += 1*(nu+nx+ny)*sizeof(double); // nls_hess_in
		size += 1*((nu+nx)*(nu+nx))*sizeof(double); // nls_hess_out
	}

	size += 64; // blasfeo_mem align
	size += 8;
	return size;

}



static void ocp_nlp_cost_nls_cast_workspace(void *config_, void *dims_, void *opts_, void *work_)
{

	// ocp_nlp_cost_config *config = config_;
	ocp_nlp_cost_nls_dims *dims = dims_;
	ocp_nlp_cost_nls_opts *opts = opts_;
	ocp_nlp_cost_nls_workspace *work = work_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int ny = dims->ny;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_cost_nls_workspace);

	align_char_to(8, &c_ptr);

	// nls_jac_in
	assign_and_advance_double(nu+nx, &work->nls_jac_in, &c_ptr);
	// nls_jac_out
	assign_and_advance_double(ny+ny*(nu+nx), &work->nls_jac_out, &c_ptr);
	if (!opts->gauss_newton_hess)
	{
		// nls_hess_in
		assign_and_advance_double(nu+nx+ny, &work->nls_hess_in, &c_ptr);
		// nls_hess_out
		assign_and_advance_double((nu+nx)*(nu+nx), &work->nls_hess_out, &c_ptr);
	}

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// tmp_nv_ny
	assign_and_advance_blasfeo_dmat_mem(nu+nx, ny, &work->tmp_nv_ny, &c_ptr);

	// tmp_ny
	assign_and_advance_blasfeo_dvec_mem(ny, &work->tmp_ny, &c_ptr);

    assert((char *)work + ocp_nlp_cost_nls_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

	return;
}



/************************************************
* functions
************************************************/

void ocp_nlp_cost_nls_initialize(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{

	ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;
    ocp_nlp_cost_nls_memory *memory= memory_;
    // ocp_nlp_cost_nls_workspace *work= work_;

	ocp_nlp_cost_nls_cast_workspace(config_, dims, opts_, work_);

	// int nx = dims->nx;
	// int nu = dims->nu;
	int ny = dims->ny;

	// TODO recompute factorization only if W are re-tuned ???
	blasfeo_dpotrf_l(ny, &model->W, 0, 0, &memory->W_chol, 0, 0);

	return;

}



void ocp_nlp_cost_nls_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{

	ocp_nlp_cost_nls_dims *dims = dims_;
    ocp_nlp_cost_nls_model *model = model_;
    ocp_nlp_cost_nls_opts *opts = opts_;
    ocp_nlp_cost_nls_memory *memory= memory_;
    ocp_nlp_cost_nls_workspace *work= work_;

	ocp_nlp_cost_nls_cast_workspace(config_, dims, opts_, work_);

	int nx = dims->nx;
	int nu = dims->nu;
	int ny = dims->ny;

	double *ext_fun_in[3]; // XXX large enough ?
	double *ext_fun_out[3]; // XXX large enough ?

	// unpack input
	blasfeo_unpack_dvec(nu, memory->ux, 0, work->nls_jac_in+nx);
	blasfeo_unpack_dvec(nx, memory->ux, nu, work->nls_jac_in);

	ext_fun_in[0] = work->nls_jac_in+0; // x: nx
	ext_fun_in[1] = work->nls_jac_in+nx; // u: nu

	ext_fun_out[0] = work->nls_jac_out+0; // fun: ny
	ext_fun_out[1] = work->nls_jac_out+ny; // jac: ny*(nx+nu)

	// evaluate external function (that assumes variables stacked as [x; u] )
	model->nls_jac->evaluate(model->nls_jac, ext_fun_in, ext_fun_out);

	// pack residuals into res
	blasfeo_pack_dvec(ny, work->nls_jac_out, &memory->res, 0);
	// pack jacobian into Jt
	blasfeo_pack_tran_dmat(ny, nx, work->nls_jac_out+ny, ny, &memory->Jt, nu, 0);
	blasfeo_pack_tran_dmat(ny, nu, work->nls_jac_out+ny+ny*nx, ny, &memory->Jt, 0, 0);

	/* gradient */

	blasfeo_daxpy(ny, -1.0, &model->y_ref, 0, &memory->res, 0, &memory->res, 0);

	// TODO use lower triangular chol of W to save n_y^2 flops
	blasfeo_dsymv_l(ny, ny, 1.0, &model->W, 0, 0, &memory->res, 0, 0.0, &work->tmp_ny, 0, &work->tmp_ny, 0);
	blasfeo_dgemv_n(nu+nx, ny, 1.0, &memory->Jt, 0, 0, &work->tmp_ny, 0, 0.0, &memory->grad, 0, &memory->grad, 0);

	/* hessian */

	if (opts->gauss_newton_hess)
	{

		// gauss-newton approximation of hessian of ls cost

		blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, &memory->W_chol, 0, 0, &memory->Jt, 0, 0, &work->tmp_nv_ny, 0, 0);
		blasfeo_dsyrk_ln(nu+nx, ny, 1.0, &work->tmp_nv_ny, 0, 0, &work->tmp_nv_ny, 0, 0, 0.0, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);

	}
	else
	{
		// exact hessian of ls cost

		// unpack input
		blasfeo_unpack_dvec(nu, memory->ux, 0, work->nls_hess_in+nx);
		blasfeo_unpack_dvec(nx, memory->ux, nu, work->nls_hess_in);
		blasfeo_unpack_dvec(ny, &work->tmp_ny, 0, work->nls_hess_in+nx+nu);

		ext_fun_in[0] = work->nls_hess_in+0; // x: nx
		ext_fun_in[1] = work->nls_hess_in+nx; // u: nu
		ext_fun_in[2] = work->nls_hess_in+nx+nu; // fun: ny

		ext_fun_out[0] = work->nls_hess_out+0; // hess: (nx+nu)*(nx+nu)

		// evaluate external function (that assumes variables stacked as [x; u] )
		model->nls_hess->evaluate(model->nls_hess, ext_fun_in, ext_fun_out);

		// pack hessian
		blasfeo_pack_dmat(nx, nx, work->nls_hess_out, nx+nu, memory->RSQrq, nu, nu); // Q
		blasfeo_pack_tran_dmat(nu, nx, work->nls_hess_out+nx, nx+nu, memory->RSQrq, nu, 0); // S
		blasfeo_pack_dmat(nu, nu, work->nls_hess_out+nx+nx*(nx+nu), nx+nu, memory->RSQrq, 0, 0); // R

		// gauss-newton component update
		blasfeo_dtrmm_rlnn(nu+nx, ny, 1.0, &memory->W_chol, 0, 0, &memory->Jt, 0, 0, &work->tmp_nv_ny, 0, 0);
		blasfeo_dsyrk_ln(nu+nx, ny, 1.0, &work->tmp_nv_ny, 0, 0, &work->tmp_nv_ny, 0, 0, 1.0, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);

	}


	return;

}



void ocp_nlp_cost_nls_config_initialize_default(void *config_)
{
	ocp_nlp_cost_config *config = config_;

	config->dims_calculate_size = &ocp_nlp_cost_nls_dims_calculate_size;
	config->dims_assign = &ocp_nlp_cost_nls_dims_assign;
	config->dims_initialize = &ocp_nlp_cost_nls_dims_initialize;
	config->model_calculate_size = &ocp_nlp_cost_nls_model_calculate_size;
	config->model_assign = &ocp_nlp_cost_nls_model_assign;
	config->opts_calculate_size = &ocp_nlp_cost_nls_opts_calculate_size;
	config->opts_assign = &ocp_nlp_cost_nls_opts_assign;
	config->opts_initialize_default = &ocp_nlp_cost_nls_opts_initialize_default;
	config->opts_update = &ocp_nlp_cost_nls_opts_update;
	config->memory_calculate_size = &ocp_nlp_cost_nls_memory_calculate_size;
	config->memory_assign = &ocp_nlp_cost_nls_memory_assign;
	config->memory_get_grad_ptr = &ocp_nlp_cost_nls_memory_get_grad_ptr;
	config->memory_set_ux_ptr = &ocp_nlp_cost_nls_memory_set_ux_ptr;
	config->memory_set_RSQrq_ptr = &ocp_nlp_cost_nls_memory_set_RSQrq_ptr;
	config->workspace_calculate_size = &ocp_nlp_cost_nls_workspace_calculate_size;
	config->initialize = &ocp_nlp_cost_nls_initialize;
	config->update_qp_matrices = &ocp_nlp_cost_nls_update_qp_matrices;
	config->config_initialize_default = &ocp_nlp_cost_nls_config_initialize_default;

	return;

}

