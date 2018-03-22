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

#ifndef ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_constraints.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"



/************************************************
* structures
************************************************/

/************************************************
* dims
************************************************/

typedef struct
{
	// NOTE(giaf) add some (arrays of) int to work outside (and independent of) sub-modules? e.g. nv, ne, nc (variables, equality-constraints, inequality-constraints
	ocp_nlp_cost_dims **cost;
	ocp_nlp_dynamics_dims **dynamics;
	ocp_nlp_constraints_dims **constraints;
	ocp_qp_dims *qp_solver; // xcond_solver instrad ???
    int N;
} ocp_nlp_dims;



/************************************************
* in
************************************************/

typedef struct
{
	ocp_nlp_dims *dims;  // pointer to nlp dimensions

	double *Ts; // length of sampling intervals

    void **cost;
	void **dynamics;
	void **constraints;

    // TODO(rien): what about invariants, e.g., algebraic constraints?

    bool freezeSens;  // TODO(dimitris): shouldn't this be in the integrator args?
} ocp_nlp_in;



/************************************************
* out
************************************************/

typedef struct
{
	struct blasfeo_dvec *ux;
	struct blasfeo_dvec *pi;
	struct blasfeo_dvec *lam;
	struct blasfeo_dvec *t;
} ocp_nlp_out;



/************************************************
* memory TODO move to sqp ???
************************************************/

typedef struct
{
	struct blasfeo_dvec *cost_grad;
	struct blasfeo_dvec *ineq_fun;
	struct blasfeo_dvec *ineq_adj;
	struct blasfeo_dvec *dyn_fun;
	struct blasfeo_dvec *dyn_adj;
} ocp_nlp_memory;



/************************************************
* residuals
************************************************/

typedef struct
{
	struct blasfeo_dvec *res_g; // stationarity
	struct blasfeo_dvec *res_b; // dynamics
	struct blasfeo_dvec *res_d; // inequality constraints
	struct blasfeo_dvec *res_m; // complementarity
	double inf_norm_res_g;
	double inf_norm_res_b;
	double inf_norm_res_d;
	double inf_norm_res_m;
	int memsize;
} ocp_nlp_res;



/************************************************
* config
************************************************/

typedef struct
{
	int N; // number of stages

	// all the others
    int (*evaluate)(void *config, ocp_nlp_dims *dims, ocp_nlp_in *qp_in, ocp_nlp_out *qp_out, void *opts_, void *mem, void *work);
    int (*opts_calculate_size)(void *config, ocp_nlp_dims *dims);
    void *(*opts_assign)(void *config, ocp_nlp_dims *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, ocp_nlp_dims *dims, void *opts_);
    int (*memory_calculate_size)(void *config, ocp_nlp_dims *dims, void *opts_);
    void *(*memory_assign)(void *config, ocp_nlp_dims *dims, void *opts_, void *raw_memory);
    int (*workspace_calculate_size)(void *config, ocp_nlp_dims *dims, void *opts_);
    ocp_qp_xcond_solver_config *qp_solver;
//    sim_solver_config **sim_solvers;
    ocp_nlp_dynamics_config **dynamics;
	ocp_nlp_cost_config **cost;
	ocp_nlp_constraints_config **constraints;
} ocp_nlp_solver_config;





/************************************************
* headers
************************************************/

/************************************************
* config
************************************************/

//
int ocp_nlp_solver_config_calculate_size(int N);
//
ocp_nlp_solver_config *ocp_nlp_solver_config_assign(int N, void *raw_memory);

/************************************************
* dims
************************************************/

//
int ocp_nlp_dims_calculate_size(int N);
//
ocp_nlp_dims *ocp_nlp_dims_assign(int N, void *raw_memory);
//
void ocp_nlp_dims_initialize(int *nx, int *nu, int *ny, int *nbx, int *nbu, int *ng, int *nh, int *nq, int *ns, ocp_nlp_dims *dims);

/************************************************
* in
************************************************/

//
int ocp_nlp_in_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory);

/************************************************
* out
************************************************/

//
int ocp_nlp_out_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory);

/************************************************
* memory
************************************************/

//
int ocp_nlp_memory_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory);

/************************************************
* residuals
************************************************/

//
int ocp_nlp_res_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_res *ocp_nlp_res_assign(ocp_nlp_dims *dims, void *raw_memory);
//
void ocp_nlp_res_compute(ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_res *res, ocp_nlp_memory *mem); //ocp_nlp_res_workspace *work);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
