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

#ifndef INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_
#define INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/sim/sim_gnsf.h"
// acados_c
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/sim_interface.h"


/// Solution methods for optimal control problems.
typedef enum
{
    SQP,
    SQP_RTI,
    INVALID_NLP_SOLVER,
} ocp_nlp_solver_t;


/// Types of the cost function.
typedef enum
{
    LINEAR_LS,
    NONLINEAR_LS,
    EXTERNALLY_PROVIDED,
    INVALID_COST,
} ocp_nlp_cost_t;


/// Types of the system dynamics, discrete or continuous time.
typedef enum
{
    CONTINUOUS_MODEL,
    DISCRETE_MODEL,
    INVALID_DYNAMICS,
} ocp_nlp_dynamics_t;


/// Bound types
typedef enum
{
    /// Comprises simple bounds, polytopic constraints,
    /// general non-linear constraints.
    BGH,

    /// Comprises simple bounds, polytopic constraints,
    /// general non-linear constraints, and positive definite constraints.
    BGHP,

    INVALID_CONSTRAINT,
} ocp_nlp_constraints_t;


/// Regularization types
typedef enum
{
    NO_REGULARIZE,
    MIRROR,
    PROJECT,
    CONVEXIFY,
    INVALID_REGULARIZE,
} ocp_nlp_reg_t;


/// Structure to store the configuration of a non-linear program
typedef struct
{
    /// Qp solver configuration.
    ocp_qp_solver_plan ocp_qp_solver_plan;

    /// Simulation solver configuration for each stage.
    sim_solver_plan *sim_solver_plan;

    /// Nlp solver type.
    ocp_nlp_solver_t nlp_solver;

    /// Regularization type, defaults to no regularization.
    ocp_nlp_reg_t regularization;

    /// Cost type for each stage.
    ocp_nlp_cost_t *nlp_cost;

    /// Dynamics type for each stage.
    ocp_nlp_dynamics_t *nlp_dynamics;

    /// Constraints type for each stage.
    ocp_nlp_constraints_t *nlp_constraints;

    /// Horizon length.
    int N;

} ocp_nlp_plan;


/// Structure to store state of the non-linear programming solver
typedef struct
{
    ocp_nlp_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} ocp_nlp_solver;


/// Construct an empty plan object (nlp configuration), all fields are set to a
/// default/invalid state.
ocp_nlp_plan *ocp_nlp_plan_create(int N);

/// Destroys a plan object, frees memory.
void ocp_nlp_plan_destroy(void* plan_);

/* config */
ocp_nlp_config *ocp_nlp_config_create(ocp_nlp_plan plan);
//
void ocp_nlp_config_destroy(void *config_);

/* dims */
ocp_nlp_dims *ocp_nlp_dims_create(void *config_);
//
void ocp_nlp_dims_destroy(void *dims_);
//

/* in */
ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void ocp_nlp_in_destroy(void *in);
//
void ocp_nlp_in_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage,
        const char *field, void *value);
//
int ocp_nlp_dynamics_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                           int stage, const char *fun_type, void *fun_ptr);
//
// TODO remove and use ocp_nlp_dynamics_model_set instead !!!
int nlp_set_discrete_model_in_stage(ocp_nlp_config *config, ocp_nlp_in *in, int stage,
                                    void *fun_ptr);
//
int ocp_nlp_cost_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                           int stage, const char *field, void *value);
//
int ocp_nlp_constraints_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                     ocp_nlp_in *in, int stage, const char *field, void *value);

/* out */
//
ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void ocp_nlp_out_destroy(void *out);
//
void ocp_nlp_out_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
                     int stage, const char *field, void *value);
//
void ocp_nlp_out_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
                     int stage, const char *field, void *value);

/* opts */
//
void *ocp_nlp_opts_create(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void ocp_nlp_opts_destroy(void *opts);
//
void ocp_nlp_opts_set(ocp_nlp_config *config, void *opts_,
                      const char *field, const void* value);
//
void ocp_nlp_dynamics_opts_set(ocp_nlp_config *config, void *opts_, int stage,
                                         const char *field, void *value);
//
void ocp_nlp_opts_update(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_);


/* solver */
//
ocp_nlp_solver *ocp_nlp_solver_create(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_);
//
void ocp_nlp_solver_destroy(void *solver);
//
int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out);
//
int ocp_nlp_precompute(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out);

/* get */
void ocp_nlp_get(ocp_nlp_config *config, ocp_nlp_solver *solver,
                 const char *field, void *return_value_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_
