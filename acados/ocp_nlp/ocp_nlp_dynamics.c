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

#include "acados/ocp_nlp/ocp_nlp_dynamics.h"

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
* config
************************************************/

int ocp_nlp_dynamics_config_calculate_size()
{

	int size = 0;

	size += sizeof(ocp_nlp_dynamics_config);

	size += sim_solver_config_calculate_size();

	return size;

}



ocp_nlp_dynamics_config *ocp_nlp_dynamics_config_assign(void *raw_memory)
{

	char *c_ptr = raw_memory;

	ocp_nlp_dynamics_config *config = (ocp_nlp_dynamics_config *) c_ptr;
	c_ptr += sizeof(ocp_nlp_dynamics_config);

	config->sim_solver = sim_solver_config_assign(c_ptr);
	c_ptr += sim_solver_config_calculate_size();

	return config;

}



/************************************************
* dims
************************************************/

int ocp_nlp_dynamics_dims_calculate_size()
{
    int size = 0;

	size += sizeof(ocp_nlp_dynamics_dims);

	size += sim_dims_calculate_size();

    return size;
}



ocp_nlp_dynamics_dims *ocp_nlp_dynamics_dims_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_nlp_dynamics_dims *dims = (ocp_nlp_dynamics_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_dims);

	dims->sim = sim_dims_assign(c_ptr);
	c_ptr += sim_dims_calculate_size();

    assert((char *) raw_memory + ocp_nlp_dynamics_dims_calculate_size() >= c_ptr);

    return dims;
}



/************************************************
* options
************************************************/

int ocp_nlp_dynamics_opts_calculate_size(void *config_, ocp_nlp_dynamics_dims *dims)
{
	ocp_nlp_dynamics_config *config = config_;

    int size = 0;

    size += sizeof(ocp_nlp_dynamics_opts);

    size += config->sim_solver->opts_calculate_size(config->sim_solver, dims->sim);

    return size;
}



void *ocp_nlp_dynamics_opts_assign(void *config_, ocp_nlp_dynamics_dims *dims, void *raw_memory)
{
	ocp_nlp_dynamics_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_dynamics_opts *opts = (ocp_nlp_dynamics_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_opts);

	opts->sim_solver = config->sim_solver->opts_assign(config->sim_solver, dims->sim, c_ptr);
	c_ptr += config->sim_solver->opts_calculate_size(config->sim_solver, dims->sim);

    assert((char*)raw_memory + ocp_nlp_dynamics_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_dynamics_opts_initialize_default(void *config_, ocp_nlp_dynamics_dims *dims, void *opts_)
{
	ocp_nlp_dynamics_config *config = config_;
	ocp_nlp_dynamics_opts *opts = opts_;

	config->sim_solver->opts_initialize_default(config->sim_solver, dims->sim, opts->sim_solver);

	return;

}



/************************************************
* memory
************************************************/

int ocp_nlp_dynamics_memory_calculate_size(void *config_, ocp_nlp_dynamics_dims *dims)
{
	// ocp_nlp_dynamics_config *config = config_;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nx1 = dims->nx1;

	int size = 0;

    size += sizeof(ocp_nlp_dynamics_memory);

	size += 1*blasfeo_memsize_dvec(nu+nx+nx1); // adj
	size += 1*blasfeo_memsize_dvec(nx1); // fun

	size += 64; // blasfeo_mem align

	return size;
}



void *ocp_nlp_dynamics_memory_assign(void *config_, ocp_nlp_dynamics_dims *dims, void *raw_memory)
{
	ocp_nlp_dynamics_config *config = config_;

	char *c_ptr = (char *) raw_memory;

	// extract dims
	int nx = dims->nx;
	int nu = dims->nu;
	int nx1 = dims->nx1;

	// struct
    ocp_nlp_dynamics_memory *memory = (ocp_nlp_dynamics_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_memory);

	// blasfeo_mem align
	align_char_to(64, &c_ptr);

	// adj
	assign_blasfeo_dvec_mem(nu+nx+nx1, &memory->adj, &c_ptr);
	// fun
	assign_blasfeo_dvec_mem(nx1, &memory->fun, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_dynamics_memory_calculate_size(config, dims) >= c_ptr);

	return memory;
}



/************************************************
* workspace
************************************************/

int ocp_nlp_dynamics_workspace_calculate_size(void *config, ocp_nlp_dynamics_dims *dims)
{

	int size = 0;

	return size;

}



/************************************************
* dynamics
************************************************/

int ocp_nlp_dynamics_model_calculate_size(void *config_, ocp_nlp_dynamics_dims *dims)
{
	ocp_nlp_dynamics_config *config = config_;

	// extract dims
	// int nx = dims->nx;
	// int nu = dims->nu;

	int size = 0;

    size += sizeof(ocp_nlp_dynamics_model);

	size += config->sim_solver->model_calculate_size(config->sim_solver, dims->sim);

	return size;
}



void *ocp_nlp_dynamics_model_assign(void *config_, ocp_nlp_dynamics_dims *dims, void *raw_memory)
{
	ocp_nlp_dynamics_config *config = config_;

	char *c_ptr = (char *) raw_memory;

	// extract dims
	// int nx = dims->nx;
	// int nu = dims->nu;

	// struct
    ocp_nlp_dynamics_model *model = (ocp_nlp_dynamics_model *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_model);

	// dims
	model->dims = dims;

	model->sim_model = config->sim_solver->model_assign(config->sim_solver, dims->sim, c_ptr);
	c_ptr += config->sim_solver->model_calculate_size(config->sim_solver, dims->sim);

    assert((char *) raw_memory + ocp_nlp_dynamics_model_calculate_size(config, dims) >= c_ptr);

	return model;
}



void ocp_nlp_dynamics_update_qp_matrices(void *config_, ocp_nlp_dynamics_dims *dims, ocp_nlp_dynamics_memory *mem, sim_in *sim_in, sim_out *sim_out, void *sim_opts, void *sim_mem, void *sim_work)
{

	ocp_nlp_dynamics_config *config = config_;

	int nx = dims->nx;
	int nu = dims->nu;
	int nx1 = dims->nx1;
	int nu1 = dims->nu1;

	// pass state and control to integrator
	blasfeo_unpack_dvec(nu, mem->ux, 0, sim_in->u);
	blasfeo_unpack_dvec(nx, mem->ux, nu, sim_in->x);

	// call integrator
	config->sim_solver->evaluate(config->sim_solver, sim_in, sim_out, sim_opts, sim_mem, sim_work);

	// TODO(rien): transition functions for changing dimensions not yet implemented!

	// B
	blasfeo_pack_tran_dmat(nx1, nu, sim_out->S_forw+nx1*nx, nx1, mem->BAbt, 0, 0);
	// A
	blasfeo_pack_tran_dmat(nx1, nx, sim_out->S_forw+0, nx1, mem->BAbt, nu, 0);

	// fun
	blasfeo_pack_dvec(nx1, sim_out->xn, &mem->fun, 0);
	blasfeo_daxpy(nx1, -1.0, mem->ux1, nu1, &mem->fun, 0, &mem->fun, 0);

	// adj TODO if not computed by the integrator
	blasfeo_dgemv_n(nu+nx, nx1, -1.0, mem->BAbt, 0, 0, mem->pi, 0, 0.0, &mem->adj, 0, &mem->adj, 0);
	blasfeo_dveccp(nx1, mem->pi, 0, &mem->adj, nu+nx);

	return;
}



void ocp_nlp_dynamics_config_initialize_default(void *config_)
{
	ocp_nlp_dynamics_config *config = config_;

	config->model_calculate_size = &ocp_nlp_dynamics_model_calculate_size;
	config->model_assign = &ocp_nlp_dynamics_model_assign;
	config->opts_calculate_size = &ocp_nlp_dynamics_opts_calculate_size;
	config->opts_assign = &ocp_nlp_dynamics_opts_assign;
	config->opts_initialize_default = &ocp_nlp_dynamics_opts_initialize_default;
	config->memory_calculate_size = &ocp_nlp_dynamics_memory_calculate_size;
	config->memory_assign = &ocp_nlp_dynamics_memory_assign;
	config->config_initialize_default = &ocp_nlp_dynamics_config_initialize_default;

	return;

}

