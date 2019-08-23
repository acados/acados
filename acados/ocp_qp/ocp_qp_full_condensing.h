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


#ifndef ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_H_
#define ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"



typedef struct
{
	ocp_qp_dims *orig_dims;
	dense_qp_dims *fcond_dims;
} ocp_qp_full_condensing_dims;



typedef struct ocp_qp_full_condensing_opts_
{
    struct d_cond_qp_arg *hpipm_opts;
    dense_qp_dims *fcond_dims;  // TODO(all): move to dims
    int cond_hess; // 0 cond only rhs, 1 cond hess + rhs
    int expand_dual_sol; // 0 primal sol only, 1 primal + dual sol
	int ric_alg;
	int mem_qp_in; // allocate qp_in in memory
} ocp_qp_full_condensing_opts;



typedef struct ocp_qp_full_condensing_memory_
{
    struct d_cond_qp_ws *hpipm_workspace;
	// in memory
	dense_qp_in *fcond_qp_in;
	dense_qp_out *fcond_qp_out;
	// only pointer
    ocp_qp_in *ptr_qp_in;
	qp_info *qp_out_info; // info in fcond_qp_in
} ocp_qp_full_condensing_memory;



//
int ocp_qp_full_condensing_opts_calculate_size(void *dims);
//
void *ocp_qp_full_condensing_opts_assign(void *dims, void *raw_memory);
//
void ocp_qp_full_condensing_opts_initialize_default(void *dims, void *opts_);
//
void ocp_qp_full_condensing_opts_update(void *dims, void *opts_);
//
void ocp_qp_full_condensing_opts_set(void *opts_, const char *field, void* value);
//
int ocp_qp_full_condensing_memory_calculate_size(void *dims, void *opts_);
//
void *ocp_qp_full_condensing_memory_assign(void *dims, void *opts_, void *raw_memory);
//
int ocp_qp_full_condensing_workspace_calculate_size(void *dims, void *opts_);
//
int ocp_qp_full_condensing(void *in, void *out, void *opts, void *mem, void *work);
//
int ocp_qp_full_expansion(void *in, void *out, void *opts, void *mem, void *work);
//
void ocp_qp_full_condensing_config_initialize_default(void *config_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_H_
