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

#ifndef ACADOS_OCP_NLP_OCP_NLP_DYNAMICSS_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_DYNAMICSS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
// acados
#include "acados/utils/types.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"
//#include "acados/ocp_nlp/ocp_nlp_common.h" // XXX loop of include !!!



/************************************************
* dims
************************************************/

typedef struct
{
	sim_dims *sim;
    int nx; // number of states at the current stage
    int nu; // number of inputs at the current stage
	int nx1; // number of states at the next stage
	int nu1; // number of inputes at the next stage
} ocp_nlp_dynamics_dims;

//
int ocp_nlp_dynamics_dims_calculate_size();
//
ocp_nlp_dynamics_dims *ocp_nlp_dynamics_dims_assign(void *raw_memory);



/************************************************
* config
************************************************/

typedef struct
{
	int (*model_calculate_size) (void *config, ocp_nlp_dynamics_dims *dims);
	void *(*model_assign) (void *config, ocp_nlp_dynamics_dims *dims, void *raw_memory);
	void (*model_set_T) (double T, void *model);
	int (*opts_calculate_size) (void *config, ocp_nlp_dynamics_dims *dims);
	void *(*opts_assign) (void *config, ocp_nlp_dynamics_dims *dims, void *raw_memory);
	void (*opts_initialize_default) (void *config, ocp_nlp_dynamics_dims *dims, void *opts);
	void (*opts_update) (void *config, ocp_nlp_dynamics_dims *dims, void *opts);
	int (*memory_calculate_size) (void *config, ocp_nlp_dynamics_dims *dims, void *opts);
	void *(*memory_assign) (void *config, ocp_nlp_dynamics_dims *dims, void *opts, void *raw_memory);
	struct blasfeo_dvec *(*memory_get_fun_ptr) (void *memory_);
	struct blasfeo_dvec *(*memory_get_adj_ptr) (void *memory_);
	void (*memory_set_ux_ptr) (struct blasfeo_dvec *ux, void *memory_);
	void (*memory_set_ux1_ptr) (struct blasfeo_dvec *ux1, void *memory_);
	void (*memory_set_pi_ptr) (struct blasfeo_dvec *pi, void *memory_);
	void (*memory_set_BAbt_ptr) (struct blasfeo_dmat *BAbt, void *memory_);
	int (*workspace_calculate_size) (void *config, ocp_nlp_dynamics_dims *dims, void *opts);
	void (*initialize) (void *config_, ocp_nlp_dynamics_dims *dims, void *model_, void *opts_, void *mem_, void *work_);
	void (*update_qp_matrices) (void *config_, ocp_nlp_dynamics_dims *dims, void *model_, void *opts_, void *mem_, void *work_);
	void (*config_initialize_default) (void *config);
    sim_solver_config *sim_solver;
} ocp_nlp_dynamics_config;

//
int ocp_nlp_dynamics_config_calculate_size();
//
ocp_nlp_dynamics_config *ocp_nlp_dynamics_config_assign(void *raw_memory);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_COMMON_H_
