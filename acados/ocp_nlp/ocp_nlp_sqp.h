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


/// \addtogroup ocp_nlp
/// @{
/// \addtogroup ocp_nlp_solver
/// @{
/// \addtogroup ocp_nlp_sqp ocp_nlp_sqp
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_SQP_H_
#define ACADOS_OCP_NLP_OCP_NLP_SQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"



/************************************************
 * options
 ************************************************/

typedef struct
{
    void *qp_solver_opts;
    void *regularize;
    void **dynamics;     // dynamics_opts
    void **cost;         // cost_opts
    void **constraints;  // constraints_opts
    double tol_stat;     // exit tolerance on stationarity condition
    double tol_eq;       // exit tolerance on equality constraints
    double tol_ineq;     // exit tolerance on inequality constraints
    double tol_comp;     // exit tolerance on complemetarity condition
    int max_iter;
    int reuse_workspace;
    int num_threads;
	int ext_qp_res;      // compute external QP residuals (i.e. at SQP level) at each SQP iteration (for debugging)
	int qp_warm_start;
} ocp_nlp_sqp_opts;

//
int ocp_nlp_sqp_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_sqp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_sqp_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_sqp_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_sqp_opts_set(void *config_, void *opts_, const char *field, void* value);
//
void ocp_nlp_sqp_dyanimcs_opts_set(void *config, void *opts, int stage, const char *field, void *value);



/************************************************
 * memory
 ************************************************/

typedef struct
{
	// qp in & out
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;
	// QP stuff not entering the qp_in struct
    struct blasfeo_dmat *dzduxt; // dzdux transposed
    struct blasfeo_dvec *z_alg; // z_alg

    //    ocp_nlp_dims *dims;
    void *qp_solver_mem;

    void *regularize_mem;

    void **dynamics;     // dynamics memory
    void **cost;         // cost memory
    void **constraints;  // constraints memory

    // residuals
    ocp_nlp_res *nlp_res;

    // nlp memory
    ocp_nlp_memory *nlp_mem;

    int status;

    int sqp_iter;

    double time_qp_sol;
    double time_lin;
    double time_reg;
    double time_tot;

	double *stat;
	int stat_m;
	int stat_n;
} ocp_nlp_sqp_memory;

//
int ocp_nlp_sqp_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_nlp_sqp_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);



/************************************************
 * workspace
 ************************************************/

typedef struct
{
	// temp QP in & out (to be used as workspace)
    ocp_qp_in *tmp_qp_in;
    ocp_qp_out *tmp_qp_out;

    // QP solver
    void *qp_work;
	ocp_qp_res *qp_res;
	ocp_qp_res_ws *qp_res_ws;

    void **dynamics;     // dynamics_workspace
    void **cost;         // cost_workspace
    void **constraints;  // constraints_workspace

} ocp_nlp_sqp_work;

//
int ocp_nlp_sqp_workspace_calculate_size(void *config, void *dims, void *opts_);



/************************************************
 * functions
 ************************************************/

//
int ocp_nlp_sqp(void *config, void *dims, void *nlp_in, void *nlp_out,
                void *args, void *mem, void *work_);
//
void ocp_nlp_sqp_config_initialize_default(void *config_);
//
int ocp_nlp_sqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_SQP_H_
/// @}
/// @}
/// @}
