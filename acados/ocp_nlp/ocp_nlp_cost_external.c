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



void ocp_nlp_cost_external_dims_initialize(void *config_, void *dims_, int nx, int nu, int ny)
{
	ocp_nlp_cost_external_dims *dims = dims_;

	dims->nx = nx;
	dims->nu = nu;

	return;
}



/************************************************
* model
************************************************/

int ocp_nlp_cost_external_model_calculate_size(void *config_, void *dims_)
{
	// extract dims
	// int nx = dims->nx;
	// int nu = dims->nu;
	// int ny = dims->ny;

	int size = 0;

	size += sizeof(ocp_nlp_cost_external_model);

	return size;

}



void *ocp_nlp_cost_external_model_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

	// int nx = dims->nx;
	// int nu = dims->nu;
	// int ny = dims->ny;

	// struct
    ocp_nlp_cost_external_model *model = (ocp_nlp_cost_external_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_model);

	// assert
    assert((char *) raw_memory + ocp_nlp_cost_external_model_calculate_size(config_, dims_) >= c_ptr);

	return model;

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

    assert((char*)raw_memory + ocp_nlp_cost_external_opts_calculate_size(config_, dims_) >= c_ptr);

    return opts;
}



void ocp_nlp_cost_external_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
	// ocp_nlp_cost_config *config = config_;
//	ocp_nlp_cost_external_opts *opts = opts_;

//	opts->gauss_newton_hess = 1;

	return;

}



void ocp_nlp_cost_external_opts_update(void *config_, void *dims_, void *opts_)
{
	// ocp_nlp_cost_config *config = config_;
//	ocp_nlp_cost_external_opts *opts = opts_;

//	opts->gauss_newton_hess = 1;

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
//	int ny = dims->ny;

	int size = 0;

    size += sizeof(ocp_nlp_cost_external_memory);

	size += 1*blasfeo_memsize_dvec(nu+nx); // grad

	size += 64; // blasfeo_mem align

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
//	int ny = dims->ny;

	// struct
    ocp_nlp_cost_external_memory *memory = (ocp_nlp_cost_external_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_cost_external_memory);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// grad
	assign_and_advance_blasfeo_dvec_mem(nu+nx, &memory->grad, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_cost_external_memory_calculate_size(config_, dims, opts_) >= c_ptr);

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

	size += 1*(nu+nx)*sizeof(double); // ext_cost_in
	size += 1*((nu+nx)+(nu+nx)*(nu+nx))*sizeof(double); // ext_cost_out

	return size;

}



static void ocp_nlp_cost_external_cast_workspace(void *config_, void *dims_, void *opts_, void *work_)
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

	// external_jac_in
	assign_and_advance_double(nu+nx, &work->ext_cost_in, &c_ptr);
	// external_jac_out
	assign_and_advance_double((nu+nx)+(nu+nx)*(nu+nx), &work->ext_cost_out, &c_ptr);

    assert((char *)work + ocp_nlp_cost_external_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

	return;
}



/************************************************
* functions
************************************************/

void ocp_nlp_cost_external_initialize(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{

//    ocp_nlp_cost_external_model *model = model_;
//    ocp_nlp_cost_external_memory *memory= memory_;
    // ocp_nlp_cost_external_workspace *work= work_;

//	ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

	return;

}



void ocp_nlp_cost_external_update_qp_matrices(void *config_, void *dims_, void *model_, void *opts_, void *memory_, void *work_)
{

	ocp_nlp_cost_external_dims *dims = dims_;
    ocp_nlp_cost_external_model *model = model_;
    // ocp_nlp_cost_external_opts *opts = opts_;
    ocp_nlp_cost_external_memory *memory= memory_;
    ocp_nlp_cost_external_workspace *work= work_;

	ocp_nlp_cost_external_cast_workspace(config_, dims, opts_, work_);

	int nx = dims->nx;
	int nu = dims->nu;

	// unpack ls cost input
	blasfeo_unpack_dvec(nu, memory->ux, 0, work->ext_cost_in+nx);
	blasfeo_unpack_dvec(nx, memory->ux, nu, work->ext_cost_in);

	// evaluate external function (that assumes variables stacked as [x; u] )
	model->ext_cost->evaluate(model->ext_cost, work->ext_cost_in, work->ext_cost_out);

	// pack gradient
	blasfeo_pack_dvec(nx, work->ext_cost_out, &memory->grad, nu); // q
	blasfeo_pack_dvec(nu, work->ext_cost_out+nx, &memory->grad, 0); // r

	// pack Hessian
	double *d_ptr = work->ext_cost_out+nu+nx;
	blasfeo_pack_dmat(nx, nx, d_ptr, nx+nu, memory->RSQrq, nu, nu); // Q
	blasfeo_pack_tran_dmat(nu, nx, d_ptr+nx, nx+nu, memory->RSQrq, nu, 0); // S
	blasfeo_pack_dmat(nu, nu, d_ptr+nx*(nu+nx+1), nx+nu, memory->RSQrq, 0, 0); // R

	return;

}



/* config */

void ocp_nlp_cost_external_config_initialize_default(void *config_)
{
	ocp_nlp_cost_config *config = config_;

	config->dims_calculate_size = &ocp_nlp_cost_external_dims_calculate_size;
	config->dims_assign = &ocp_nlp_cost_external_dims_assign;
	config->dims_initialize = &ocp_nlp_cost_external_dims_initialize;
	config->model_calculate_size = &ocp_nlp_cost_external_model_calculate_size;
	config->model_assign = &ocp_nlp_cost_external_model_assign;
	config->opts_calculate_size = &ocp_nlp_cost_external_opts_calculate_size;
	config->opts_assign = &ocp_nlp_cost_external_opts_assign;
	config->opts_initialize_default = &ocp_nlp_cost_external_opts_initialize_default;
	config->opts_update = &ocp_nlp_cost_external_opts_update;
	config->memory_calculate_size = &ocp_nlp_cost_external_memory_calculate_size;
	config->memory_assign = &ocp_nlp_cost_external_memory_assign;
	config->memory_get_grad_ptr = &ocp_nlp_cost_external_memory_get_grad_ptr;
	config->memory_set_ux_ptr = &ocp_nlp_cost_external_memory_set_ux_ptr;
	config->memory_set_RSQrq_ptr = &ocp_nlp_cost_external_memory_set_RSQrq_ptr;
	config->workspace_calculate_size = &ocp_nlp_cost_external_workspace_calculate_size;
	config->initialize = &ocp_nlp_cost_external_initialize;
	config->update_qp_matrices = &ocp_nlp_cost_external_update_qp_matrices;
	config->config_initialize_default = &ocp_nlp_cost_external_config_initialize_default;

	return;

}


