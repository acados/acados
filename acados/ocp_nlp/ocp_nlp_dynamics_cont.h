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

#ifndef ACADOS_OCP_NLP_OCP_NLP_DYNAMICSS_CONT_H_
#define ACADOS_OCP_NLP_OCP_NLP_DYNAMICSS_CONT_H_

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
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"



/************************************************
* options
************************************************/

typedef struct
{
    void *sim_solver;
} ocp_nlp_dynamics_cont_opts;

//
int ocp_nlp_dynamics_cont_opts_calculate_size(void *config, ocp_nlp_dynamics_dims *dims);
//
void *ocp_nlp_dynamics_cont_opts_assign(void *config, ocp_nlp_dynamics_dims *dims, void *raw_memory);
//
void ocp_nlp_dynamics_cont_opts_initialize_default(void *config, ocp_nlp_dynamics_dims *dims, void *opts);



/************************************************
* memory
************************************************/

typedef struct
{
	struct blasfeo_dvec fun;
	struct blasfeo_dvec adj;
	struct blasfeo_dvec *ux; // pointer to ux in nlp_out at current stage
	struct blasfeo_dvec *ux1; // pointer to ux in nlp_out at next stage
	struct blasfeo_dvec *pi; // pointer to pi in nlp_out at current stage
	struct blasfeo_dmat *BAbt; // pointer to BAbt in qp_in
	void *sim_solver; // sim solver memory
} ocp_nlp_dynamics_cont_memory;

//
int ocp_nlp_dynamics_cont_memory_calculate_size(void *config, ocp_nlp_dynamics_dims *dims, void *opts);
//
void *ocp_nlp_dynamics_cont_memory_assign(void *config, ocp_nlp_dynamics_dims *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_cont_memory_get_fun_ptr(void *memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_cont_memory_get_adj_ptr(void *memory);
//
void ocp_nlp_dynamics_cont_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory);
//
void ocp_nlp_dynamics_cont_memory_set_ux1_ptr(struct blasfeo_dvec *ux1, void *memory);
//
void ocp_nlp_dynamics_cont_memory_set_pi_ptr(struct blasfeo_dvec *pi, void *memory);
//
void ocp_nlp_dynamics_cont_memory_set_BAbt_ptr(struct blasfeo_dmat *BAbt, void *memory);



/************************************************
* workspace
************************************************/

typedef struct
{
    sim_in *sim_in;
    sim_out *sim_out;
    void *sim_solver; // sim solver workspace
} ocp_nlp_dynamics_cont_workspace;

int ocp_nlp_dynamics_cont_workspace_calculate_size(void *config, ocp_nlp_dynamics_dims *dims, void *opts);


/************************************************
* model
************************************************/

typedef struct
{
	ocp_nlp_dynamics_dims *dims;
	void *sim_model;
//	double *state_transition; // TODO
	double T; // simulation time
} ocp_nlp_dynamics_cont_model;

//
int ocp_nlp_dynamics_cont_model_calculate_size(void *config, ocp_nlp_dynamics_dims *dims);
//
void *ocp_nlp_dynamics_cont_model_assign(void *config, ocp_nlp_dynamics_dims *dims, void *raw_memory);
//
void ocp_nlp_dynamics_cont_model_set_T(double T, void *model);



/************************************************
* functions
************************************************/

//
void ocp_nlp_dynamics_cont_config_initialize_default(void *config);
//
void ocp_nlp_dynamics_cont_update_qp_matrices(void *config_, ocp_nlp_dynamics_dims *dims, void *model_, void *opts, void *mem, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_CONT_H_
