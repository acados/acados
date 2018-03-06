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



void ocp_nlp_dynamics_config_initialize_default(void *config_)
{
	ocp_nlp_dynamics_config *config = config_;

	config->model_calculate_size = &ocp_nlp_dynamics_model_calculate_size;
	config->model_assign = &ocp_nlp_dynamics_model_assign;
	config->config_initialize_default = &ocp_nlp_dynamics_config_initialize_default;

	return;

}

