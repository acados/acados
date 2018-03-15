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

#include "acados/ocp_nlp/ocp_nlp_constraints.h"

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
#include "acados/ocp_qp/ocp_qp_common.h"



/************************************************
* config
************************************************/

int ocp_nlp_constraints_config_calculate_size()
{

	int size = 0;

	size += sizeof(ocp_nlp_constraints_config);

	return size;

}



ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory)
{

	char *c_ptr = raw_memory;

	ocp_nlp_constraints_config *config = (ocp_nlp_constraints_config *) c_ptr;
	c_ptr += sizeof(ocp_nlp_constraints_config);

	return config;

}



/************************************************
* dims
************************************************/

int ocp_nlp_constraints_dims_calculate_size()
{
    int size = sizeof(ocp_nlp_constraints_dims);

    return size;
}



ocp_nlp_constraints_dims *ocp_nlp_constraints_dims_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_dims *dims = (ocp_nlp_constraints_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_dims);

    assert((char *) raw_memory + ocp_nlp_constraints_dims_calculate_size() >= c_ptr);

    return dims;
}




/************************************************
* linear constraints
************************************************/

/* model */

int ocp_nlp_constraints_calculate_size(void *config, ocp_nlp_constraints_dims *dims)
{
	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;

	int size = 0;

    size += sizeof(ocp_nlp_constraints_model);

	size += sizeof(int)*nb;  // idxb
	size += blasfeo_memsize_dvec(2*nb+2*ng); // d
	size += blasfeo_memsize_dmat(nu+nx, ng); // DCt

	size += 64; // blasfeo_mem align

	return size;
}



void *ocp_nlp_constraints_assign(void *config, ocp_nlp_constraints_dims *dims, void *raw_memory)
{
	char *c_ptr = (char *) raw_memory;

	// extract sizes
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;

	// struct
    ocp_nlp_constraints_model *model = (ocp_nlp_constraints_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_model);

	// dims
	model->dims = dims;

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// blasfeo_dmat
	// DCt
	assign_and_advance_blasfeo_dmat_mem(nu+nx, ng, &model->DCt, &c_ptr);

	// blasfeo_dvec
	// d
	assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng, &model->d, &c_ptr);

    // idxb
    assign_and_advance_int(dims->nbx+dims->nbu, &model->idxb, &c_ptr);

	// h
	model->h = NULL;

	// assert
    assert((char *) raw_memory + ocp_nlp_constraints_calculate_size(config, dims) >= c_ptr);

	return model;
}



/* options */

int ocp_nlp_constraints_opts_calculate_size(void *config_, ocp_nlp_constraints_dims *dims)
{
    int size = 0;

    size += sizeof(ocp_nlp_constraints_opts);

    return size;
}



void *ocp_nlp_constraints_opts_assign(void *config_, ocp_nlp_constraints_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_constraints_opts *opts = (ocp_nlp_constraints_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_opts);

    assert((char*)raw_memory + ocp_nlp_constraints_opts_calculate_size(config_, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_constraints_opts_initialize_default(void *config_, ocp_nlp_constraints_dims *dims, void *opts_)
{
//	ocp_nlp_constraints_opts *opts = opts_;

	return;

}



/* memory */

int ocp_nlp_constraints_memory_calculate_size(void *config_, ocp_nlp_constraints_dims *dims, void *opts_)
{
	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;

	int size = 0;

    size += sizeof(ocp_nlp_constraints_memory);

	size += 1*blasfeo_memsize_dvec(2*nb+2*ng+2*nh);  // fun
	size += 1*blasfeo_memsize_dvec(nu+nx);  // adj

	size += 1*64;  // blasfeo_mem align

	return size;
}



void *ocp_nlp_constraints_memory_assign(void *config_, ocp_nlp_constraints_dims *dims, void *opts_, void *raw_memory)
{
	char *c_ptr = (char *) raw_memory;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;

	// struct
    ocp_nlp_constraints_memory *memory = (ocp_nlp_constraints_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_constraints_memory);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// fun
	assign_and_advance_blasfeo_dvec_mem(2*nb+2*ng+2*nh, &memory->fun, &c_ptr);
	// adj
	assign_and_advance_blasfeo_dvec_mem(nu+nx, &memory->adj, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_constraints_memory_calculate_size(config_, dims, opts_) >= c_ptr);

	return memory;
}



struct blasfeo_dvec *ocp_nlp_constraints_memory_get_fun_ptr(void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	return &memory->fun;
}



struct blasfeo_dvec *ocp_nlp_constraints_memory_get_adj_ptr(void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	return &memory->adj;
}



void ocp_nlp_constraints_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->ux = ux;
}



void ocp_nlp_constraints_memory_set_lam_ptr(struct blasfeo_dvec *lam, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->lam = lam;
}



void ocp_nlp_constraints_memory_set_DCt_ptr(struct blasfeo_dmat *DCt, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->DCt = DCt;
}



void ocp_nlp_constraints_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->RSQrq = RSQrq;
}



void ocp_nlp_constraints_memory_set_idxb_ptr(int *idxb, void *memory_)
{
	ocp_nlp_constraints_memory *memory = memory_;

	memory->idxb = idxb;
}



/* workspace */

int ocp_nlp_constraints_workspace_calculate_size(void *config_, ocp_nlp_constraints_dims *dims, void *opts_)
{
	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;
	int nq = dims->nq;

	int size = 0;

    size += sizeof(ocp_nlp_constraints_workspace);

	size += 1*blasfeo_memsize_dvec(nb+ng+nh);  // tmp_nbg
	size += (nx+nu)*sizeof(double);  // nl_constraint_input
	size += nh*(1+nx+nu)*sizeof(double);  // nl_constraint_output
	if (nq > 0) {
		size += nq*(nx+nu)*sizeof(double);
		size += 1*blasfeo_memsize_dmat(nx+nu, nq);
	}

	size += 2*64;  // blasfeo_mem align

	return size;

}



static void ocp_nlp_constraints_cast_workspace(void *config_, ocp_nlp_constraints_dims *dims, void *opts_, void *work_)
{
	ocp_nlp_constraints_workspace *work = work_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;
	int nq = dims->nq;

    char *c_ptr = (char *) work_;
    c_ptr += sizeof(ocp_nlp_constraints_workspace);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// tmp_nbg
	assign_and_advance_blasfeo_dvec_mem(nb+ng+nh, &work->tmp_nbg, &c_ptr);
	work->nl_constraint_input = (double *) c_ptr;
	c_ptr += (nx+nu)*sizeof(double);
	work->nl_constraint_output = (double *) c_ptr;
	c_ptr += nh*(1+nx+nu)*sizeof(double);
	if (nq > 0) {
		c_ptr += nq*(nx+nu)*sizeof(double);
		align_char_to(64, &c_ptr);
		assign_and_advance_blasfeo_dmat_mem(nx+nu, nq, &work->jacobian_quadratic, &c_ptr);	
	}

    assert((char *)work + ocp_nlp_constraints_workspace_calculate_size(config_, dims, opts_) >= c_ptr);

	return;
}



/* functions */

void ocp_nlp_constraints_initialize_qp(void *config_, ocp_nlp_constraints_dims *dims, void *model_, void *opts, void *memory_, void *work_)
{

	ocp_nlp_constraints_model *model = model_;
	ocp_nlp_constraints_memory *memory = memory_;

	// loop index
	int j;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;

	// initialize idxb
	for (j=0; j<nb; j++)
	{
		memory->idxb[j] = model->idxb[j];
	}

	// initialize general constraints matrix
	blasfeo_dgecp(nu+nx, ng, &model->DCt, 0, 0, memory->DCt, 0, 0);

	return;

}


void ocp_nlp_constraints_update_qp_matrices(void *config_, ocp_nlp_constraints_dims *dims, void *model_, void *opts_, void *memory_, void *work_)
{

	ocp_nlp_constraints_model *model = model_;
	ocp_nlp_constraints_memory *memory = memory_;
	ocp_nlp_constraints_workspace *work = work_;

	ocp_nlp_constraints_cast_workspace(config_, dims, opts_, work_);

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nb = dims->nb;
	int ng = dims->ng;
	int nh = dims->nh;
	int nq = dims->nq;

	blasfeo_dvecex_sp(nb, 1.0, model->idxb, memory->ux, 0, &work->tmp_nbg, 0);
	blasfeo_dgemv_t(nu+nx, ng, 1.0, memory->DCt, 0, 0, memory->ux, 0, 0.0, &work->tmp_nbg, nb, &work->tmp_nbg, nb);

	if (nh>0)
	{
		blasfeo_unpack_dvec(nx, memory->ux, nu, work->nl_constraint_input);
		blasfeo_unpack_dvec(nu, memory->ux, 0, work->nl_constraint_input+nx);

		model->h->evaluate(model->h, work->nl_constraint_input, work->nl_constraint_output);

		blasfeo_pack_dvec(nh, work->nl_constraint_output, &work->tmp_nbg, nb+ng);
		blasfeo_pack_tran_dmat(nh, nx, work->nl_constraint_output+nh, nh, memory->DCt, ng, nu);
		blasfeo_pack_tran_dmat(nh, nu, work->nl_constraint_output+nh+nx*nh, nh, memory->DCt, ng, 0);
	}

	if (nq>0)
	{
		if (nh != 1) {
			printf("Not implemented");
			exit(1);
		}
		double lam = blasfeo_dvecex1(memory->lam, 2*(nb+ng)+nh) - blasfeo_dvecex1(memory->lam, nb+ng);
		blasfeo_pack_tran_dmat(nq, nx+nu, work->nl_constraint_output+nh*(1+nx+nu), nq, &work->jacobian_quadratic, 0, 0);
		blasfeo_dsyrk_ln(nx+nu, nq, 2*lam, &work->jacobian_quadratic, 0, 0, &work->jacobian_quadratic, 0, 0,
			0.0, memory->RSQrq, 0, 0, memory->RSQrq, 0, 0);
	}

	blasfeo_daxpy(nb+ng+nh, -1.0, &work->tmp_nbg, 0, &model->d, 0, &memory->fun, 0);
	blasfeo_daxpy(nb+ng+nh, -1.0, &model->d, nb+ng+nh, &work->tmp_nbg, 0, &memory->fun, nb+ng+nh);

	// nlp_mem: ineq_adj
	blasfeo_dvecse(nu+nx, 0.0, &memory->adj, 0);
	blasfeo_daxpy(nb+ng+nh, -1.0, memory->lam, nb+ng+nh, memory->lam, 0, &work->tmp_nbg, 0);
	blasfeo_dvecad_sp(nb, 1.0, &work->tmp_nbg, 0, model->idxb, &memory->adj, 0);
	blasfeo_dgemv_n(nu+nx, ng+nh, 1.0, memory->DCt, 0, 0, &work->tmp_nbg, nb, 1.0, &memory->adj, 0, &memory->adj, 0);

	return;

}



void ocp_nlp_constraints_config_initialize_default(void *config_)
{
	ocp_nlp_constraints_config *config = config_;

	config->model_calculate_size = &ocp_nlp_constraints_calculate_size;
	config->model_assign = &ocp_nlp_constraints_assign;
	config->opts_calculate_size = &ocp_nlp_constraints_opts_calculate_size;
	config->opts_assign = &ocp_nlp_constraints_opts_assign;
	config->opts_initialize_default = &ocp_nlp_constraints_opts_initialize_default;
	config->memory_calculate_size = &ocp_nlp_constraints_memory_calculate_size;
	config->memory_assign = &ocp_nlp_constraints_memory_assign;
	config->memory_get_fun_ptr = &ocp_nlp_constraints_memory_get_fun_ptr;
	config->memory_get_adj_ptr = &ocp_nlp_constraints_memory_get_adj_ptr;
	config->memory_set_ux_ptr = &ocp_nlp_constraints_memory_set_ux_ptr;
	config->memory_set_lam_ptr = &ocp_nlp_constraints_memory_set_lam_ptr;
	config->memory_set_DCt_ptr = &ocp_nlp_constraints_memory_set_DCt_ptr;
	config->memory_set_RSQrq_ptr = &ocp_nlp_constraints_memory_set_RSQrq_ptr;
	config->memory_set_idxb_ptr = &ocp_nlp_constraints_memory_set_idxb_ptr;
	config->workspace_calculate_size = &ocp_nlp_constraints_workspace_calculate_size;
	config->initialize_qp = &ocp_nlp_constraints_initialize_qp;
	config->update_qp_matrices = &ocp_nlp_constraints_update_qp_matrices;
	config->config_initialize_default = &ocp_nlp_constraints_config_initialize_default;

	return;

}




