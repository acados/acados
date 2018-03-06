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
 *    Author: Jonathan Frey
 */

// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_gnsf.h"
#include "acados/sim/sim_gnsf_casadi_wrapper.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "blasfeo/include/blasfeo_d_kernel.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux.h"


int gnsf2_dims_calculate_size()
{
    int size = sizeof(gnsf2_dims);
    return size;
}


gnsf2_dims *gnsf2_dims_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;
    gnsf2_dims *dims = (gnsf2_dims *) c_ptr;
    c_ptr += sizeof(gnsf2_dims);
    assert((char *) raw_memory + gnsf_dims_calculate_size() == c_ptr);
    return dims;
}

void gnsf2_get_dims( gnsf2_dims *dims, casadi_function_t get_ints_fun)
{
    double *ints_out;
    ints_out = (double*) calloc(9,sizeof(double));
    export_from_ML_wrapped(ints_out, ints_out, get_ints_fun);

    dims->nx = (int) ints_out[0];
    dims->nu = (int) ints_out[1]; 
    dims->nz = (int) ints_out[2];
    dims->nx1 = (int) ints_out[3];
    dims->nx2 = (int) ints_out[4];
    dims->num_stages = (int) ints_out[5];
    dims->num_steps = (int) ints_out[6];
    dims->n_out = (int) ints_out[7];
    dims->n_in = (int) ints_out[8];
    free(ints_out);
}

void sim_gnsf2_config_initialize_default(void *config_)
{
	sim_solver_config *config = config_;
// TODO:!!
	// config->evaluate = &gnsf2_simulate;
	// config->opts_calculate_size = &gnsf2_opts_calculate_size;
	// config->opts_assign = &gnsf2_opts_assign;
	// config->opts_initialize_default = &sim_irk_opts_initialize_default; TODO
	// config->memory_calculate_size = &sim_irk_memory_calculate_size; TODO
	// config->memory_assign = &sim_irk_memory_assign;  TODO
	config->workspace_calculate_size = &gnsf_calculate_workspace_size;
	// config->model_calculate_size = &gnsf2_model_calculate_size;
	// config->model_assign = &gnsf2_model_assign;

	return;
}