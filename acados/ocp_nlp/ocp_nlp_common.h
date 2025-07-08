/*
 * Copyright (c) The acados authors.
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



/// \defgroup ocp_nlp ocp_nlp
/// @{
/// @}

/// \defgroup ocp_nlp_solver ocp_nlp_solver
/// @{
/// @}

/// \ingroup ocp_nlp
/// @{

/// \ingroup ocp_nlp_solver
/// @{

/// \defgroup ocp_nlp_common ocp_nlp_common
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_nlp/ocp_nlp_constraints_common.h"
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/ocp_nlp/ocp_nlp_qpscaling.h"
#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_xcond_solver.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"


/************************************************
 * config
 ************************************************/
// NOTE: in here only void* arguments, as ocp_nlp_in etc are defined based on config.
typedef struct ocp_nlp_config
{
    int N;  // number of stages
    int with_feasible_qp;

    // solver-specific implementations of memory management functions
    acados_size_t (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts_);
    void (*opts_update)(void *config, void *dims, void *opts_);
    acados_size_t (*memory_calculate_size)(void *config, void *dims, void *opts_, void *in);
    void *(*memory_assign)(void *config, void *dims, void *opts_, void *in, void *raw_memory);
    acados_size_t (*workspace_calculate_size)(void *config, void *dims, void *opts_, void *in);
    void (*opts_set)(void *config_, void *opts_, const char *field, void* value);
    void (*opts_set_at_stage)(void *config_, void *opts_, size_t stage, const char *field, void* value);
    // evaluate solver // TODO rename into solve
    int (*evaluate)(void *config, void *dims, void *nlp_in, void *nlp_out, void *opts_, void *mem, void *work);
    int (*setup_qp_matrices_and_factorize)(void *config, void *dims, void *nlp_in, void *nlp_out, void *opts_, void *mem, void *work);

    void (*eval_kkt_residual)(void *config, void *dims, void *nlp_in, void *nlp_out, void *opts_, void *mem, void *work);
    void (*eval_param_sens)(void *config, void *dims, void *opts_, void *mem, void *work,
                            char *field, int stage, int index, void *sens_nlp_out);
    void (*eval_lagr_grad_p)(void *config, void *dims, void *nlp_in, void *opts_, void *mem, void *work,
                            const char *field, void *grad_p);
    void (*eval_solution_sens_adj_p)(void *config_, void *dims_,
                        void *opts_, void *mem_, void *work_, void *sens_nlp_out,
                        const char *field, int stage, void *grad_p);
    void (*step_update)(void *config, void *dims, void *in,
            void *out_start, void *opts, void *mem, void *work,
            void *out_destination, void* solver_mem, double alpha, bool full_step_dual);
    // prepare memory
    int (*precompute)(void *config, void *dims, void *nlp_in, void *nlp_out, void *opts_, void *mem, void *work);
    void (*memory_reset_qp_solver)(void *config, void *dims, void *nlp_in, void *nlp_out, void *opts_, void *mem, void *work);
    // initialize this struct with default values
    void (*config_initialize_default)(void *config);
    // general getter
    void (*get)(void *config_, void *dims, void *mem_, const char *field, void *return_value_);
    void (*opts_get)(void *config_, void *dims, void *opts_, const char *field, void *return_value_);
    void (*work_get)(void *config_, void *dims, void *work_, const char *field, void *return_value_);
    //
    void (*terminate)(void *config, void *mem, void *work);

    bool (*is_real_time_algorithm)();


    // config structs of submodules
    ocp_qp_xcond_solver_config *qp_solver; // TODO rename xcond_solver
    ocp_qp_xcond_solver_config *relaxed_qp_solver;
    ocp_nlp_dynamics_config **dynamics;
    ocp_nlp_cost_config **cost;
    ocp_nlp_constraints_config **constraints;
    ocp_nlp_reg_config *regularize;
    ocp_nlp_globalization_config *globalization;

} ocp_nlp_config;

//
acados_size_t ocp_nlp_config_calculate_size(int N);
//
ocp_nlp_config *ocp_nlp_config_assign(int N, void *raw_memory);



/************************************************
 * dims
 ************************************************/

/// Structure to store dimensions/number of variables.
typedef struct ocp_nlp_dims
{
    void **cost;
    void **dynamics;
    void **constraints;
    ocp_qp_xcond_solver_dims *qp_solver;  // xcond solver instead ??
    ocp_qp_xcond_solver_dims *relaxed_qp_solver;  // xcond solver instead ??
    ocp_nlp_reg_dims *regularize;
    ocp_nlp_qpscaling_dims *qpscaling;
    ocp_nlp_qpscaling_dims *relaxed_qpscaling;

    int *nv;  // number of primal variables (states+controls+slacks)
    int *nx;  // number of differential states
    int *nu;  // number of inputs
    int *nz;  // number of algebraic variables
    int *ns;  // number of slack variables
    int *np;  // number of parameters
    // constraints
    int *ni;  // number of two-sided inequality constraints: nb+ng+nh+ns+nphi
    int *nb;  // number of two-sided bounds
    int *ng;  // number of two-sided general linear constraints
    int *ni_nl;  // number of two-sided nonlinear inequalities

    int np_global;  // number of global parameters
    int n_global_data;  // size of global_data; expressions that only depend on p_global; detected automatically during code generation
    int N;    // number of shooting nodes

    // total dimensions
    int nx_total;  // total number of states
    int nu_total;  // total number of controls
    int ns_total;  // total number of slack variables
    int nz_total;  // total number of algebraic variables
    // int npi_total; // = nx_total - nx_0
    int ni_total;  // total number of inequalities
    int np_total;  // total number of parameters
    int nbx_total;  // total number of state bounds
    int nbu_total;  // total number of control bounds
    int ng_total;  // total number of general linear constraints
    int nh_total;  // total number of nonlinear inequalities
    int nphi_total;  // total number of nonlinear inequalities

    void *raw_memory; // Pointer to allocated memory, to be used for freeing
} ocp_nlp_dims;

//
acados_size_t ocp_nlp_dims_calculate_size(void *config);
//
ocp_nlp_dims *ocp_nlp_dims_assign(void *config, void *raw_memory);

/// Sets the dimension of optimization variables
/// (states, controls, algebraic variables, slack variables).
///
/// \param config_ The configuration struct.
/// \param dims_ The dimension struct.
/// \param field The type of optimization variables, either nx, nu, nz, or ns.
/// \param value_array Number of variables for each stage.
void ocp_nlp_dims_set_opt_vars(void *config_, void *dims_,
                               const char *field, const void* value_array);

/// Sets the dimensions of constraints functions for a stage
/// (bounds on states, bounds on controls, equality constraints,
/// inequality constraints).
///
/// \param config_ The configuration struct.
/// \param dims_ The dimension struct.
/// \param stage Stage number.
/// \param field The type of constraint/bound, either nbx, nbu, ng, or nh.
/// \param value_field Number of constraints/bounds for the given stage.
void ocp_nlp_dims_set_constraints(void *config_, void *dims_, int stage,
                                  const char *field, const void* value_field);

/// Sets the dimensions of the cost terms for a stage.
///
/// \param config_ The configuration struct.
/// \param dims_ The dimension struct.
/// \param stage Stage number.
/// \param field Type of cost dimension
/// \param value_field Number of cost terms/residuals for the given stage.
void ocp_nlp_dims_set_cost(void *config_, void *dims_, int stage, const char *field,
                           const void* value_field);

/// Sets the dimensions of the dynamics for a stage.
///
/// \param config_ The configuration struct.
/// \param dims_ The dimension struct.
/// \param stage Stage number.
/// \param field TBD
/// \param value TBD
void ocp_nlp_dims_set_dynamics(void *config_, void *dims_, int stage, const char *field,
                               const void* value);

void ocp_nlp_dims_set_global(void *config_, void *dims_, const char *field, int value_field);


/************************************************
 * Inputs
 ************************************************/

/// Struct for storing the inputs of an OCP NLP solver
typedef struct ocp_nlp_in
{
    /// Timesteps.
    double *Ts;

    /// Parameter values.
    double **parameter_values;

    /// Global data
    double *global_data;

    /// Constraint mask
    struct blasfeo_dvec *dmask;

    /// Pointers to cost functions (TBC).
    void **cost;

    /// Pointers to dynamics functions (TBC).
    void **dynamics;

    /// Pointers to constraints functions (TBC).
    void **constraints;

    /// Pointer to allocated memory, to be used for freeing.
    void *raw_memory;

} ocp_nlp_in;

//
acados_size_t ocp_nlp_in_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_in *ocp_nlp_in_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory);


/************************************************
 * out
 ************************************************/

typedef struct ocp_nlp_out
{
    struct blasfeo_dvec *ux;  // NOTE: this contains [u; x; s_l; s_u]! - rename to uxs?
    struct blasfeo_dvec *z;  // algebraic variables
    struct blasfeo_dvec *pi;  // multipliers for dynamics
    struct blasfeo_dvec *lam;  // inequality multipliers

    // NOTE: the inequalities are internally organized in the following order:
    // [ lbu lbx lg lh lphi ubu ubx ug uh uphi; lsbu lsbx lsg lsh lsphi usbu usbx usg ush usphi]
    double inf_norm_res;

    void *raw_memory; // Pointer to allocated memory, to be used for freeing

} ocp_nlp_out;

//
acados_size_t ocp_nlp_out_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_out *ocp_nlp_out_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                void *raw_memory);



/************************************************
 * options
 ************************************************/
typedef struct ocp_nlp_opts
{
    ocp_qp_xcond_solver_opts *qp_solver_opts; // xcond solver opts instead ???
    void *regularize;
    void *qpscaling;
    void *globalization;  // globalization_opts
    void **dynamics;     // dynamics_opts
    void **cost;         // cost_opts
    void **constraints;  // constraints_opts
    double levenberg_marquardt;  // LM factor to be added to the hessian before regularization
    int reuse_workspace;
    int num_threads;
    int print_level;
    int fixed_hess;
    int log_primal_step_norm; // compute and log the max norm of the primal steps
    int log_dual_step_norm; // compute and log the max norm of the dual steps
    int max_iter; // maximum number of (SQP/DDP) iterations
    int qp_iter_max; // maximum iter of QP solver, stored to remember.
    double tau_min;  // minimum value of the barrier parameter, for IPMs

    // Flag for usage of adaptive levenberg marquardt strategy
    bool with_adaptive_levenberg_marquardt;
    double adaptive_levenberg_marquardt_lam;
    double adaptive_levenberg_marquardt_mu_min;
    double adaptive_levenberg_marquardt_mu0;
    double adaptive_levenberg_marquardt_obj_scalar;

    int with_solution_sens_wrt_params;
    int with_value_sens_wrt_params;
    double solution_sens_qp_t_lam_min;

    int ext_qp_res;
    int qp_warm_start;

    bool warm_start_first_qp; // to set qp_warm_start in first iteration
    bool warm_start_first_qp_from_nlp;  // if True first QP will be initialized using values from NLP iterate, otherwise from previous QP solution.
    bool eval_residual_at_max_iter; // if convergence should be checked after last iterations or only throw max_iter reached
    bool store_iterates; // flag indicating whether intermediate iterates should be stored

    bool with_anderson_acceleration;

    // termination tolerances
    double tol_stat;     // exit tolerance on stationarity condition
    double tol_eq;       // exit tolerance on equality constraints
    double tol_ineq;     // exit tolerance on inequality constraints
    double tol_comp;     // exit tolerance on complementarity condition
    double tol_unbounded; // exit threshold when objective function seems to be unbounded
    double tol_min_step_norm; // exit tolerance for small step

} ocp_nlp_opts;

//
acados_size_t ocp_nlp_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_opts_set(void *config_, void *opts_, const char *field, void* value);
//
void ocp_nlp_opts_set_at_stage(void *config, void *opts, int stage, const char *field, void *value);


/************************************************
 * residuals
 ************************************************/

typedef struct ocp_nlp_res
{
    struct blasfeo_dvec *res_stat;  // stationarity
    struct blasfeo_dvec *res_eq;  // dynamics
    struct blasfeo_dvec *res_ineq;  // inequality constraints
    struct blasfeo_dvec *res_comp;  // complementarity
    struct blasfeo_dvec tmp;  // tmp
    double inf_norm_res_stat;
    double inf_norm_res_eq;
    double inf_norm_res_ineq;
    double inf_norm_res_comp;
    acados_size_t memsize;
} ocp_nlp_res;

//
acados_size_t ocp_nlp_res_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_res *ocp_nlp_res_assign(ocp_nlp_dims *dims, void *raw_memory);
//
void ocp_nlp_res_get_inf_norm(ocp_nlp_res *res, double *out);

/************************************************
 * timings
 ************************************************/

typedef struct ocp_nlp_timings
{
    // these timers are reset at every solver call
    double time_qp_sol;
    double time_qp_solver_call;
    double time_qp_xcond;
    double time_qpscaling;
    double time_lin;
    double time_reg;
    double time_tot;
    double time_glob;
    double time_sim;
    double time_sim_la;
    double time_sim_ad;
    // these are not
    double time_solution_sensitivities;
    double time_feedback;
    double time_preparation;
} ocp_nlp_timings;


void ocp_nlp_timings_get(ocp_nlp_config *config, ocp_nlp_timings *timings, const char *field, void *return_value_);

void ocp_nlp_timings_reset(ocp_nlp_timings *timings);


/************************************************
 * memory
 ************************************************/

typedef struct ocp_nlp_memory
{
//    void *qp_solver_mem; // xcond solver mem instead ???
    ocp_qp_xcond_solver_memory *qp_solver_mem; // xcond solver mem instead ???
    void *regularize_mem;
    void *qpscaling;
    void *globalization; // globalization memory
    void **dynamics;     // dynamics memory
    void **cost;         // cost memory
    void **constraints;  // constraints memory

    // intermediate iterates
    struct ocp_nlp_out ** iterates;

    // residuals
    ocp_nlp_res *nlp_res;

    // timings
    ocp_nlp_timings *nlp_timings;

    // qp in & out
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;

    // scaled qp in & out (just pointers)
    ocp_qp_in *scaled_qp_in;
    ocp_qp_out *scaled_qp_out;

    // for Anderson acceleration
    ocp_qp_out *prev_qp_out;
    ocp_qp_out *anderson_step;

    // QP stuff not entering the qp_in struct
    struct blasfeo_dmat *dzduxt; // dzdux transposed
    struct blasfeo_dvec *z_alg; // z_alg, output algebraic variables

    struct blasfeo_dvec *cost_grad;
    struct blasfeo_dvec *ineq_fun;
    struct blasfeo_dvec *ineq_adj;
    struct blasfeo_dvec *dyn_fun;
    struct blasfeo_dvec *dyn_adj;

    // optimal value gradient wrt params
    struct blasfeo_dmat *jac_lag_stat_p_global;  // jacobian of stationarity condition wrt p_global (nv, np_global)
    struct blasfeo_dmat *jac_ineq_p_global;  // jacobian of nonlinear inequalities wrt p_global (ni_nl, np_global)
    struct blasfeo_dmat *jac_dyn_p_global;  // jacobian of dynamics wrt p_global (nx_next, np_global)
    struct blasfeo_dvec out_np_global;

    double cost_value;
    double qp_cost_value;
    double predicted_infeasibility_reduction; // used for funnel globalization
    double predicted_optimality_reduction; // // used for funnel globalization
    double objective_multiplier; // used for funnel globalization
    int compute_hess;

    int status;
    int iter;

    double adaptive_levenberg_marquardt_mu;
    double adaptive_levenberg_marquardt_mu_bar;

    bool *set_sim_guess; // indicate if there is new explicitly provided guess for integration variables
    double *primal_step_norm;
    double *dual_step_norm;

    struct blasfeo_dvec *sim_guess;
    acados_size_t workspace_size;

} ocp_nlp_memory;

//
acados_size_t ocp_nlp_memory_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts, ocp_nlp_in *in);
//
ocp_nlp_memory *ocp_nlp_memory_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                      ocp_nlp_opts *opts, ocp_nlp_in *in, void *raw_memory);
//
void ocp_nlp_memory_get(ocp_nlp_config *config, ocp_nlp_memory *nlp_mem, const char *field, void *return_value_);
//
void ocp_nlp_memory_get_at_stage(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_memory *nlp_mem, int stage, const char *field, void *return_value_);

/************************************************
 * workspace
 ************************************************/

typedef struct ocp_nlp_workspace
{

    void *qp_work;
    void **dynamics;     // dynamics_workspace
    void **cost;         // cost_workspace
    void **constraints;  // constraints_workspace

    // temp QP out
    ocp_qp_out *tmp_qp_out;

    // qp residuals
    ocp_qp_res *qp_res;
    ocp_qp_res_ws *qp_res_ws;

    // qp seed (for solution sensitivities)
    ocp_qp_seed *qp_seed;

    // for globalization: -> move to module?!
    ocp_nlp_out *tmp_nlp_out;
    ocp_nlp_out *weight_merit_fun;
    struct blasfeo_dvec tmp_nv;
    struct blasfeo_dvec tmp_2ni;
    struct blasfeo_dvec dxnext_dy;

    // optimal value gradient wrt params
    struct blasfeo_dvec tmp_np_global;
    // AS-RTI
    double *tmp_nv_double;

    int *tmp_nins;

} ocp_nlp_workspace;

//
acados_size_t ocp_nlp_workspace_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *opts, ocp_nlp_in *nlp_in);
//
ocp_nlp_workspace *ocp_nlp_workspace_assign(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                ocp_nlp_opts *opts, ocp_nlp_in *nlp_in, ocp_nlp_memory *mem, void *raw_memory);



/************************************************
 * function
 ************************************************/
//
void ocp_nlp_alias_memory_to_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_initialize_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_set_primal_variable_pointers_in_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in,
                                                       ocp_nlp_out *nlp_out, ocp_nlp_memory *nlp_mem);
//
void ocp_nlp_approximate_qp_matrices(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
             ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_approximate_qp_vectors_sqp(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                 ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_zero_order_qp_update(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_level_c_update(ocp_nlp_config *config,
    ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts,
    ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_update_variables_sqp(void *config_, void *dims_,
            void *in_, void *out_, void *opts_, void *mem_, void *work_,
            void *out_destination_, void *solver_mem, double alpha, bool full_step_dual);
//
void ocp_nlp_convert_primaldelta_absdual_step_to_delta_step(ocp_nlp_config *config, ocp_nlp_dims *dims,
        ocp_nlp_out *out, ocp_qp_out *step);
//
double ocp_nlp_compute_anderson_gamma(ocp_nlp_workspace *work, ocp_qp_out *new_qp_step, ocp_qp_out *new_minus_old_qp_step);
//
int ocp_nlp_precompute_common(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);

//
void ocp_nlp_initialize_qp_from_nlp(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_qp_in *qp_in,
            ocp_nlp_out *out, ocp_qp_out *qp_out);

//
void ocp_nlp_res_compute(ocp_nlp_dims *dims, ocp_nlp_opts *opts, ocp_nlp_in *in, ocp_nlp_out *out,
                         ocp_nlp_res *res, ocp_nlp_memory *mem, ocp_nlp_workspace *work);

double ocp_nlp_compute_delta_dual_norm_inf(ocp_nlp_dims *dims, ocp_nlp_workspace *work, ocp_nlp_out *nlp_out, ocp_qp_out *qp_out);
//
void copy_ocp_nlp_out(ocp_nlp_dims *dims, ocp_nlp_out *from, ocp_nlp_out *to);

//
void ocp_nlp_cost_compute(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_eval_constraints_common(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
void ocp_nlp_get_cost_value_from_submodules(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);
//
int ocp_nlp_common_setup_qp_matrices_and_factorize(ocp_nlp_config *config, ocp_nlp_dims *dims_, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out,
                ocp_nlp_opts *nlp_opts, ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work);
//
void ocp_nlp_params_jac_compute(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work);

void ocp_nlp_common_eval_param_sens(ocp_nlp_config *config, ocp_nlp_dims *dims,
                        ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
                        char *field, int stage, int index, ocp_nlp_out *sens_nlp_out);
//
void ocp_nlp_common_eval_lagr_grad_p(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                        ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
                        const char *field, void *grad_p);
//
void ocp_nlp_common_eval_solution_sens_adj_p(ocp_nlp_config *config, ocp_nlp_dims *dims,
                        ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work,
                        ocp_nlp_out *sens_nlp_out, const char *field, int stage, void *grad_p);
//
void ocp_nlp_add_levenberg_marquardt_term(ocp_nlp_config *config, ocp_nlp_dims *dims,
    ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem,
    ocp_nlp_workspace *work, double alpha, int iter, ocp_qp_in *qp_in);
//
double ocp_nlp_compute_dual_pi_norm_inf(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out);
//
double ocp_nlp_compute_dual_lam_norm_inf(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out);
//
double ocp_nlp_get_l1_infeasibility(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_memory *nlp_mem);
//
int ocp_nlp_perform_second_order_correction(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                            ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_opts *nlp_opts,
                                            ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work,
                                            ocp_qp_in *qp_in, ocp_qp_out *qp_out);
//
int ocp_nlp_solve_qp_and_correct_dual(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *nlp_opts,
                     ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work,
                     bool precondensed_lhs, ocp_qp_in *scaled_qp_in_, ocp_qp_in *qp_in_, ocp_qp_out *scaled_qp_out_,ocp_qp_out *qp_out_,
                     ocp_qp_xcond_solver *xcond_solver);
//
int ocp_nlp_solve_qp(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_opts *nlp_opts,
    ocp_nlp_memory *nlp_mem, ocp_nlp_workspace *nlp_work,
    ocp_qp_in *qp_in_, ocp_qp_out *qp_out_,
    ocp_qp_xcond_solver *xcond_solver);
//
double ocp_nlp_compute_qp_objective_value(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out, ocp_nlp_workspace *nlp_work);
//
double ocp_nlp_compute_gradient_directional_derivative(ocp_nlp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out);

// print / debug functionality
void ocp_nlp_dump_qp_out_to_file(ocp_qp_out *qp_out, int sqp_iter, int soc);
void ocp_nlp_dump_qp_in_to_file(ocp_qp_in *qp_in, int sqp_iter, int soc);
void ocp_nlp_common_print_iteration_header();
void ocp_nlp_common_print_iteration(int iter_count, ocp_nlp_res *nlp_res);

void ocp_nlp_update_variables_sqp_delta_primal_dual(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
            ocp_nlp_out *out, ocp_nlp_opts *opts, ocp_nlp_memory *mem, ocp_nlp_workspace *work, double alpha, ocp_qp_out *step);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
/// @}
/// @}
/// @}
