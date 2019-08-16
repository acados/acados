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


#ifndef ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_SOLVER_H_
#define ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/utils/types.h"



typedef struct ocp_qp_full_condensing_solver_opts_
{
    ocp_qp_full_condensing_opts *cond_opts;
    void *qp_solver_opts;
} ocp_qp_full_condensing_solver_opts;



typedef struct ocp_qp_full_condensing_solver_memory_
{
    ocp_qp_full_condensing_memory *cond_memory;
    void *solver_memory;
    dense_qp_in *qpd_in;
    dense_qp_out *qpd_out;
} ocp_qp_full_condensing_solver_memory;



typedef struct ocp_qp_full_condensing_solver_workspace_
{
    void *cond_work;
    void *solver_workspace;
    // TODO(dimitris): move from memory to workspace
    // dense_qp_in *qpd_in;
    // dense_qp_out *qpd_out;
} ocp_qp_full_condensing_solver_workspace;



//
int ocp_qp_full_condensing_solver_opts_calculate_size(void *config, ocp_qp_dims *dims);
//
void *ocp_qp_full_condensing_solver_opts_assign(void *config, ocp_qp_dims *dims, void *raw_memory);
//
void ocp_qp_full_condensing_solver_opts_initialize_default(void *config, ocp_qp_dims *dims,
                                                           void *opts_);
//
void ocp_qp_full_condensing_solver_opts_update(void *config, ocp_qp_dims *dims, void *opts_);
//
void ocp_qp_full_condensing_solver_opts_set(void *config_, void *opts_, const char *field, void* value);
//
int ocp_qp_full_condensing_solver_memory_calculate_size(void *config, ocp_qp_dims *dims,
                                                        void *opts_);
//
void *ocp_qp_full_condensing_solver_memory_assign(void *config, ocp_qp_dims *dims, void *opts_,
                                                  void *raw_memory);
//
int ocp_qp_full_condensing_solver_workspace_calculate_size(void *config, ocp_qp_dims *dims,
                                                           void *opts_);
//
int ocp_qp_full_condensing_solver(void *config, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_,
                                  void *mem_, void *work_);
//
void ocp_qp_full_condensing_solver_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_FULL_CONDENSING_SOLVER_H_
