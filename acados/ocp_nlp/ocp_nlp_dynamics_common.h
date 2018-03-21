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
* options
************************************************/

typedef struct
{
    void *sim_solver;
} ocp_nlp_dynamics_opts;

//
int ocp_nlp_dynamics_opts_calculate_size(void *config, ocp_nlp_dynamics_dims *dims);
//
void *ocp_nlp_dynamics_opts_assign(void *config, ocp_nlp_dynamics_dims *dims, void *raw_memory);
//
void ocp_nlp_dynamics_opts_initialize_default(void *config, ocp_nlp_dynamics_dims *dims, void *opts);



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
} ocp_nlp_dynamics_memory;

//
int ocp_nlp_dynamics_memory_calculate_size(void *config, ocp_nlp_dynamics_dims *dims, void *opts);
//
void *ocp_nlp_dynamics_memory_assign(void *config, ocp_nlp_dynamics_dims *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_memory_get_fun_ptr(void *memory);
//
struct blasfeo_dvec *ocp_nlp_dynamics_memory_get_adj_ptr(void *memory);
//
void ocp_nlp_dynamics_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory);
//
void ocp_nlp_dynamics_memory_set_ux1_ptr(struct blasfeo_dvec *ux1, void *memory);
//
void ocp_nlp_dynamics_memory_set_pi_ptr(struct blasfeo_dvec *pi, void *memory);
//
void ocp_nlp_dynamics_memory_set_BAbt_ptr(struct blasfeo_dmat *BAbt, void *memory);



/************************************************
* workspace
************************************************/

typedef struct
{
    sim_in *sim_in;
    sim_out *sim_out;
    void *sim_solver; // sim solver workspace
} ocp_nlp_dynamics_workspace;

int ocp_nlp_dynamics_workspace_calculate_size(void *config, ocp_nlp_dynamics_dims *dims, void *opts);


/************************************************
* dynamics
************************************************/

typedef struct
{
	ocp_nlp_dynamics_dims *dims;
	void *sim_model;
//	double *state_transition; // TODO
	double T; // simulation time
} ocp_nlp_dynamics_model;

//
int ocp_nlp_dynamics_model_calculate_size(void *config, ocp_nlp_dynamics_dims *dims);
//
void *ocp_nlp_dynamics_model_assign(void *config, ocp_nlp_dynamics_dims *dims, void *raw_memory);
//
void ocp_nlp_dynamics_model_set_T(double T, void *model);
//
void ocp_nlp_dynamics_update_qp_matrices(void *config_, ocp_nlp_dynamics_dims *dims, void *model_, void *opts, void *mem, void *work_);



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
	int (*memory_calculate_size) (void *config, ocp_nlp_dynamics_dims *dims, void *opts);
	void *(*memory_assign) (void *config, ocp_nlp_dynamics_dims *dims, void *opts, void *raw_memory);
	struct blasfeo_dvec *(*memory_get_fun_ptr) (void *memory_);
	struct blasfeo_dvec *(*memory_get_adj_ptr) (void *memory_);
	void (*memory_set_ux_ptr) (struct blasfeo_dvec *ux, void *memory_);
	void (*memory_set_ux1_ptr) (struct blasfeo_dvec *ux1, void *memory_);
	void (*memory_set_pi_ptr) (struct blasfeo_dvec *pi, void *memory_);
	void (*memory_set_BAbt_ptr) (struct blasfeo_dmat *BAbt, void *memory_);
	int (*workspace_calculate_size) (void *config, ocp_nlp_dynamics_dims *dims, void *opts);
	void (*update_qp_matrices) (void *config_, ocp_nlp_dynamics_dims *dims, void *model_, void *opts_, void *mem_, void *work_);
	void (*config_initialize_default) (void *config);
    sim_solver_config *sim_solver;
} ocp_nlp_dynamics_config;

//
int ocp_nlp_dynamics_config_calculate_size();
//
ocp_nlp_dynamics_config *ocp_nlp_dynamics_config_assign(void *raw_memory);
//
void ocp_nlp_dynamics_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // ACADOS_OCP_NLP_OCP_NLP_DYNAMICS_COMMON_H_
