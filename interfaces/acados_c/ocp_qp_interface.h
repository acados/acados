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


#ifndef INTERFACES_ACADOS_C_OCP_QP_INTERFACE_H_
#define INTERFACES_ACADOS_C_OCP_QP_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_xcond_solver.h"


/// QP solver types (Enumeration).
///
/// Full list of fields:
///   PARTIAL_CONDENSING_HPIPM
///   PARTIAL_CONDENSING_HPMPC
///   PARTIAL_CONDENSING_OOQP
///   PARTIAL_CONDENSING_OSQP
///   PARTIAL_CONDENSING_QPDUNES
///   FULL_CONDENSING_HPIPM
///   FULL_CONDENSING_QPOASES
///   FULL_CONDENSING_QORE
///   FULL_CONDENSING_OOQP
///   INVALID_QP_SOLVER
///
/// Note: In this enumeration the partial condensing solvers have to be
///       specified before the full condensing solvers.
typedef enum {
    PARTIAL_CONDENSING_HPIPM,
#ifdef ACADOS_WITH_HPMPC
    PARTIAL_CONDENSING_HPMPC,
#endif
#ifdef ACADOS_WITH_OOQP
    PARTIAL_CONDENSING_OOQP,
#endif
#ifdef ACADOS_WITH_OSQP
    PARTIAL_CONDENSING_OSQP,
#endif
#ifdef ACADOS_WITH_QPDUNES
    PARTIAL_CONDENSING_QPDUNES,
#endif
    FULL_CONDENSING_HPIPM,
#ifdef ACADOS_WITH_QPOASES
    FULL_CONDENSING_QPOASES,
#endif
#ifdef ACADOS_WITH_QORE
    FULL_CONDENSING_QORE,
#endif
#ifdef ACADOS_WITH_OOQP
    FULL_CONDENSING_OOQP,
#endif
    INVALID_QP_SOLVER,
} ocp_qp_solver_t;


/// Struct containing qp solver
typedef struct
{
    ocp_qp_solver_t qp_solver;
} ocp_qp_solver_plan;


/// Linear ocp configuration.
typedef struct
{
    ocp_qp_xcond_solver_config *config;
    ocp_qp_xcond_solver_dims *dims;
    void *opts;
    void *mem;
    void *work;
} ocp_qp_solver;


/// Initializes the qp solver configuration.
/// TBC should this be private/static?
void ocp_qp_xcond_solver_config_initialize_from_plan(ocp_qp_solver_t solver_name, ocp_qp_xcond_solver_config *solver_config);

/// Constructs a qp solver config and Initializes with default values.
///
/// \param plan The qp solver plan struct.
ocp_qp_xcond_solver_config *ocp_qp_xcond_solver_config_create(ocp_qp_solver_plan plan);

/// Destructor for config struct, frees memory.
///
/// \param config_ The config object to destroy.
void ocp_qp_xcond_solver_config_free(ocp_qp_xcond_solver_config *config);


/// Constructs a struct that contains the dimensions for the variables of the qp.
///
/// \param N The number of variables.
ocp_qp_dims *ocp_qp_dims_create(int N);

/// Destructor of the dimensions struct.
///
/// \param dims_ The dimensions struct.
void ocp_qp_dims_free(void *dims);


//
ocp_qp_xcond_solver_dims *ocp_qp_xcond_solver_dims_create(ocp_qp_xcond_solver_config *config, int N);
//
void ocp_qp_xcond_solver_dims_free(ocp_qp_xcond_solver_dims *dims_);

/// Constructs an input object for the qp.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
ocp_qp_in *ocp_qp_in_create(ocp_qp_dims *dims);

/// Destructor of the inputs struct.
///
/// \param in_ The inputs struct.
void ocp_qp_in_free(void *in_);


/// Constructs an outputs object for the qp.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
ocp_qp_out *ocp_qp_out_create(ocp_qp_dims *dims);

/// Destructor of the outputs struct.
///
/// \param out_ The outputs struct.
void ocp_qp_out_free(void *out_);


/// Constructs an options object for the qp.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
void *ocp_qp_xcond_solver_opts_create(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims);

/// Destructor of the options struct.
///
/// \param opts_ The options struct to destroy.
void ocp_qp_xcond_solver_opts_free(ocp_qp_xcond_solver_opts *opts);


/// TBC Should be private/static?
int ocp_qp_calculate_size(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims, void *opts_);


/// TBC Reserves memory? TBC Should this be private?
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param opts_ The options struct.
/// \param raw_memory The TBD.
ocp_qp_solver *ocp_qp_assign(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims, void *opts_,
                             void *raw_memory);

/// Creates a qp solver. Reserves memory.
///
/// \param config The configuration struct.
/// \param dims The dimensions struct.
/// \param opts_ The options struct.
ocp_qp_solver *ocp_qp_create(ocp_qp_xcond_solver_config *config, ocp_qp_xcond_solver_dims *dims, void *opts_);

/// Solves the qp.
///
/// \param solver The solver.
/// \param qp_in The inputs struct.
/// \param qp_out The outputs struct.
int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out);


/// Calculates the infinity norm of the residuals.
///
/// \param solver The solver.
/// \param qp_in The inputs struct.
/// \param qp_out The outputs struct.
/// \param res Output array for the residuals.
void ocp_qp_inf_norm_residuals(ocp_qp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                               double *res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_OCP_QP_INTERFACE_H_
