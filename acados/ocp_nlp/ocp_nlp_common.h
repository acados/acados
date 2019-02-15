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

#include "acados/ocp_nlp/ocp_nlp_constraints_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * config
 ************************************************/

typedef struct
{
    int N;  // number of stages

    // solver-specific implementations of memory management functions
    int (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts_);
    void (*opts_update)(void *config, void *dims, void *opts_);
    int (*memory_calculate_size)(void *config, void *dims, void *opts_);
    void *(*memory_assign)(void *config, void *dims, void *opts_, void *raw_memory);
    int (*workspace_calculate_size)(void *config, void *dims, void *opts_);
    void (*opts_set)(void *config_, void *opts_, const char *field, const void* value);
    int (*dynamics_opts_set)(void *config, void *opts_, int stage,
                                     const char *field, void *value);
    // evaluate solver
    int (*evaluate)(void *config, void *dims, void *qp_in, void *qp_out,
                    void *opts_, void *mem, void *work);
    // prepare memory
    int (*precompute)(void *config, void *dims, void *qp_in, void *qp_out,
                void *opts_, void *mem, void *work);
    // initalize this struct with default values
    void (*config_initialize_default)(void *config);
    // general getter
    void (*get)(void *config_, void *mem_, const char *field, void *return_value_);
    // config structs of submodules
    ocp_qp_xcond_solver_config *qp_solver;
    ocp_nlp_dynamics_config **dynamics;
    ocp_nlp_cost_config **cost;
    ocp_nlp_constraints_config **constraints;
    ocp_nlp_reg_config *regularization;

} ocp_nlp_config;

//
int ocp_nlp_config_calculate_size(int N);
//
ocp_nlp_config *ocp_nlp_config_assign(int N, void *raw_memory);



/************************************************
 * dims
 ************************************************/

typedef struct
{
    void **cost;
    void **dynamics;
    void **constraints;
    ocp_qp_dims *qp_solver;  // xcond solver instead ??

    int *nv;  // number of primal variables (states+controls+slacks)
    int *nx;  // number of differential states
    int *nu;  // number of inputs
    int *ni;  // number of two-sided inequality constraints TODO make one-sided ???
    int *nz;  // number of algebraic variables
    int *ns;  // number of slack variables
    int N;    // number of shooting nodes
} ocp_nlp_dims;

//
int ocp_nlp_dims_calculate_size_self(int N);
//
int ocp_nlp_dims_calculate_size(void *config);
//
ocp_nlp_dims *ocp_nlp_dims_assign_self(int N, void *raw_memory);
//
ocp_nlp_dims *ocp_nlp_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_dims_set_opt_vars(void *config_, void *dims_,
                               const char *field, const void* value_array);
//
void ocp_nlp_dims_set_constraints(void *config_, void *dims_, int stage,
                                  const char *field, const void* value_field);
//
void ocp_nlp_dims_set_cost(void *config_, void *dims_, int stage, const char *field,
                           const void* value_field);
//
void ocp_nlp_dims_set_dynamics(void *config_, void *dims_, int stage, const char *field,
                               const void* value);

/************************************************
 * in
 ************************************************/

typedef struct
{
    double *Ts;  // length of sampling intervals

    void **cost;
    void **dynamics;
    void **constraints;
} ocp_nlp_in;

//
int ocp_nlp_in_calculate_size_self(int N);
//
int ocp_nlp_in_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_in *ocp_nlp_in_assign_self(int N, void *raw_memory);
//
ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory);


/************************************************
 * out
 ************************************************/

typedef struct
{
    struct blasfeo_dvec *ux;
    struct blasfeo_dvec *z;
    struct blasfeo_dvec *pi;
    struct blasfeo_dvec *lam;
    struct blasfeo_dvec *t;

    int sqp_iter;
    int qp_iter;
    double inf_norm_res;
    double total_time;
} ocp_nlp_out;

//
int ocp_nlp_out_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                void *raw_memory);



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

//
int ocp_nlp_memory_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                      void *raw_memory);



/************************************************
 * residuals
 ************************************************/

typedef struct
{
    struct blasfeo_dvec *res_g;  // stationarity
    struct blasfeo_dvec *res_b;  // dynamics
    struct blasfeo_dvec *res_d;  // inequality constraints
    struct blasfeo_dvec *res_m;  // complementarity
    double inf_norm_res_g;
    double inf_norm_res_b;
    double inf_norm_res_d;
    double inf_norm_res_m;
    int memsize;
} ocp_nlp_res;

//
int ocp_nlp_res_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_res *ocp_nlp_res_assign(ocp_nlp_dims *dims, void *raw_memory);
//
void ocp_nlp_res_compute(ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_res *res,
                         ocp_nlp_memory *mem);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
