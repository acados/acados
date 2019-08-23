/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
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


/// Constraint types
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
    PROJECT_REDUC_HESS,
    CONVEXIFY,
    INVALID_REGULARIZE,
} ocp_nlp_reg_t;


/// Structure to store the configuration of a non-linear program
typedef struct
{
    /// QP solver configuration.
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


/// Structure to store the state/configuration for the non-linear programming solver
typedef struct
{
    ocp_nlp_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} ocp_nlp_solver;


/// Constructs an empty plan struct (user nlp configuration), all fields are set to a
/// default/invalid state.
///
/// \param N Horizon length
ocp_nlp_plan *ocp_nlp_plan_create(int N);

/// Destructor for plan struct, frees memory.
///
/// \param plan_ The plan struct to destroy.
void ocp_nlp_plan_destroy(void* plan_);


/// Constructs an nlp configuration struct from a plan.
///
/// \param plan The plan (user nlp configuration).
ocp_nlp_config *ocp_nlp_config_create(ocp_nlp_plan plan);

/// Desctructor of the nlp configuration.
///
/// \param config_ The configuration struct.
void ocp_nlp_config_destroy(void *config_);


/// Constructs an struct that contains the dimensions of the variables.
///
/// \param config_ The configuration struct.
ocp_nlp_dims *ocp_nlp_dims_create(void *config_);

/// Destructor of the dimensions struct.
///
/// \param dims_ The dimensions struct.
void ocp_nlp_dims_destroy(void *dims_);


/// Constructs an input struct for a non-linear programs.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_config *config, ocp_nlp_dims *dims);

/// Destructor of the inputs struct.
///
/// \param in The inputs struct.
void ocp_nlp_in_destroy(void *in);


/// Sets the sampling times for the given stage.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param in The inputs struct.
/// \param stage Stage number.
/// \param field Has to be "Ts" (TBC other options).
/// \param value The sampling times (floating point).
void ocp_nlp_in_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage,
        const char *field, void *value);


/// Sets the function pointers to the dynamics functions for the given stage.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param in The inputs struct.
/// \param stage Stage number.
/// \param fun_type The name of the function type, either impl_ode_fun,
///     impl_ode_fun_jac_x_xdot, impl_ode_jac_x_xdot_u (TBC)
/// \param fun_ptr Function pointer to the dynamics function.
int ocp_nlp_dynamics_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
		int stage, const char *fun_type, void *fun_ptr);


/// Sets the function pointers to the cost functions for the given stage.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param in The inputs struct.
/// \param stage Stage number.
/// \param field The name of the field, either nls_res_jac,
///     y_ref, W (others TBC)
/// \param value Cost values.
int ocp_nlp_cost_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
		int stage, const char *field, void *value);


/// Sets the function pointers to the constraints functions for the given stage.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param in The inputs struct.
/// \param stage Stage number.
/// \param field The name of the field, either lb, ub (others TBC)
/// \param value Constraints function or values.
int ocp_nlp_constraints_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
		ocp_nlp_in *in, int stage, const char *field, void *value);

/* out */

/// Constructs an outputs struct for the non-linear program.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_config *config, ocp_nlp_dims *dims);

/// Destructor of the outputs struct.
///
/// \param out The outputs struct.
void ocp_nlp_out_destroy(void *out);


/// TBD
void ocp_nlp_out_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
		int stage, const char *field, void *value);

/// TBD
void ocp_nlp_out_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
		int stage, const char *field, void *value);
//
int ocp_nlp_dims_get_from_attr(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
		int stage, const char *field);

/* opts */

/// Creates an options struct for the non-linear program.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
void *ocp_nlp_opts_create(ocp_nlp_config *config, ocp_nlp_dims *dims);

/// Destructor of the options.
///
/// \param opts The options struct.
void ocp_nlp_opts_destroy(void *opts);

/// Sets an option.
///
/// \param config The configuration struct.
/// \param opts_ The options struct.
/// \param field Name of the option.
/// \param value Value of the option.
void ocp_nlp_opts_set(ocp_nlp_config *config, void *opts_, const char *field, void* value);

/// TBC
/// Set the option for the dynamics in a given stage.
///
/// \param config The configuration struct.
/// \param opts_ The options struct.
/// \param stage Stage number.
/// \param field Name of the option.
/// \param value Value of the option.
void ocp_nlp_dynamics_opts_set(ocp_nlp_config *config, void *opts_, int stage,
		const char *field, void *value);

/// TBC
/// Set the option for the cost in a given stage.
///
/// \param config The configuration struct.
/// \param opts_ The options struct.
/// \param stage Stage number.
/// \param field Name of the option.
/// \param value Value of the option.
void ocp_nlp_cost_opts_set(ocp_nlp_config *config, void *opts_, int stage,
		const char *field, void *value);

/// TBC
/// Set the option for the constraints in a given stage.
///
/// \param config The configuration struct.
/// \param opts_ The options struct.
/// \param stage Stage number.
/// \param field Name of the option.
/// \param value Value of the option.
void ocp_nlp_constraints_opts_set(ocp_nlp_config *config, void *opts_, int stage,
		const char *field, void *value);

/// TBC
/// Updates the options.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param opts_ The options struct.
void ocp_nlp_opts_update(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_);


/* solver */

/// Creates an ocp solver.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param opts_ The options struct.
/// \return The solver.
ocp_nlp_solver *ocp_nlp_solver_create(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_);

/// Destructor of the solver.
///
/// \param solver The solver struct.
void ocp_nlp_solver_destroy(void *solver);

/// Solves the optimal control problem. Call ocp_nlp_precompute before
/// calling this functions (TBC).
///
/// \param solver The solver struct.
/// \param nlp_in The inputs struct.
/// \param nlp_out The outputs struct.
int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out);

/// Performs precomputations for the solver. Needs to be called before
/// ocl_nlp_solve (TBC).
///
/// \param solver The solver struct.
/// \param nlp_in The inputs struct.
/// \param nlp_out The outputs struct.
int ocp_nlp_precompute(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out);

//
void ocp_nlp_eval_param_sens(ocp_nlp_solver *solver, char *field, int stage, int index, ocp_nlp_out *sens_nlp_out);

/* get */
/// TBD
void ocp_nlp_get(ocp_nlp_config *config, ocp_nlp_solver *solver,
		const char *field, void *return_value_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_
