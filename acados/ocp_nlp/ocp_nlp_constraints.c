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
* constraints
************************************************/

int ocp_nlp_constraints_model_calculate_size(void *config, ocp_nlp_constraints_dims *dims)
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



void *ocp_nlp_constraints_model_assign(void *config, ocp_nlp_constraints_dims *dims, void *raw_memory)
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
	assign_blasfeo_dmat_mem(nu+nx, ng, &model->DCt, &c_ptr);

	// blasfeo_dvec
	// d
	assign_blasfeo_dvec_mem(2*nb+2*ng, &model->d, &c_ptr);

    // idxb
    assign_int(dims->nbx+dims->nbu, &model->idxb, &c_ptr);

	// assert
    assert((char *) raw_memory + ocp_nlp_constraints_model_calculate_size(config, dims) >= c_ptr);

	return model;
}



// TODO mem and work if needed
void ocp_nlp_constraints_initialize_qp(void *config, ocp_nlp_constraints_dims *dims, ocp_nlp_constraints_model *model, int *idxb, struct blasfeo_dmat *DCt, void *mem, void *work)
{

	// loop index
	int j;

	int nx, nu, nb, ng;

	nx = dims->nx;
	nu = dims->nu;
	nb = dims->nb;
	ng = dims->ng;

	// initialize idxb
	for (j=0; j<nb; j++)
	{
		idxb[j] = model->idxb[j];
	}

	// initialize general constraints matrix
	blasfeo_dgecp(nu+nx, ng, &model->DCt, 0, 0, DCt, 0, 0);

	return;

}



void ocp_nlp_constraints_config_initialize_default(void *config_)
{
	ocp_nlp_constraints_config *config = config_;

	config->model_calculate_size = &ocp_nlp_constraints_model_calculate_size;
	config->model_assign = &ocp_nlp_constraints_model_assign;
	config->initialize_qp = &ocp_nlp_constraints_initialize_qp;
	config->config_initialize_default = &ocp_nlp_constraints_config_initialize_default;

	return;

}

