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


#ifndef ACADOS_OCP_NLP_OCP_NLP_COST_EXTERNAL_H_
#define ACADOS_OCP_NLP_OCP_NLP_COST_EXTERNAL_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo/include/blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_cost_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"

/************************************************
 * dims
 ************************************************/

typedef struct
{
    int nx;  // number of states
    int nu;  // number of inputs
    int ns;  // number of slacks
} ocp_nlp_cost_external_dims;

//
int ocp_nlp_cost_external_dims_calculate_size(void *config);
//
void *ocp_nlp_cost_external_dims_assign(void *config, void *raw_memory);
//
void ocp_nlp_cost_external_dims_initialize(void *config, void *dims, int nx, int nu, int ny, int ns, int nz);
void ocp_nlp_cost_external_dims_set(void *config_, void *dims_, const char *field, int* value);
//
void ocp_nlp_cost_external_dims_get(void *config_, void *dims_, const char *field, int* value);

/************************************************
 * model
 ************************************************/

typedef struct
{
    external_function_generic *ext_cost;  // gradient and hessian
    struct blasfeo_dvec Z;
    struct blasfeo_dvec z;
    double scaling;
} ocp_nlp_cost_external_model;

//
int ocp_nlp_cost_external_model_calculate_size(void *config, void *dims);
//
void *ocp_nlp_cost_external_model_assign(void *config, void *dims, void *raw_memory);



/************************************************
 * options
 ************************************************/

typedef struct
{
    int dummy; // struct can't be void
} ocp_nlp_cost_external_opts;

//
int ocp_nlp_cost_external_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_cost_external_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_cost_external_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_cost_external_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_cost_external_opts_set(void *config, void *opts, const char *field, void *value);



/************************************************
 * memory
 ************************************************/

typedef struct
{
    struct blasfeo_dvec grad;    // gradient of cost function
    struct blasfeo_dvec *ux;     // pointer to ux in nlp_out
    struct blasfeo_dmat *RSQrq;  // pointer to RSQrq in qp_in
    struct blasfeo_dvec *Z;      // pointer to Z in qp_in
    struct blasfeo_dvec *z_alg;         ///< pointer to z in sim_out
    struct blasfeo_dmat *dzdux_tran;    ///< pointer to sensitivity of a wrt ux in sim_out
} ocp_nlp_cost_external_memory;

//
int ocp_nlp_cost_external_memory_calculate_size(void *config, void *dims, void *opts);
//
void *ocp_nlp_cost_external_memory_assign(void *config, void *dims, void *opts, void *raw_memory);
//
struct blasfeo_dvec *ocp_nlp_cost_external_memory_get_grad_ptr(void *memory_);
//
void ocp_nlp_cost_external_memory_set_RSQrq_ptr(struct blasfeo_dmat *RSQrq, void *memory);
//
void ocp_nlp_cost_ls_memory_set_Z_ptr(struct blasfeo_dvec *Z, void *memory);
//
void ocp_nlp_cost_external_memory_set_ux_ptr(struct blasfeo_dvec *ux, void *memory_);
//
void ocp_nlp_cost_external_memory_set_z_alg_ptr(struct blasfeo_dvec *z_alg, void *memory_);
//
void ocp_nlp_cost_external_memory_set_dzdux_tran_ptr(struct blasfeo_dmat *dzdux_tran, void *memory_);

/************************************************
 * workspace
 ************************************************/

typedef struct
{
    struct blasfeo_dmat tmp_nv_nv;
} ocp_nlp_cost_external_workspace;

//
int ocp_nlp_cost_external_workspace_calculate_size(void *config, void *dims, void *opts);

/************************************************
 * functions
 ************************************************/

//
void ocp_nlp_cost_external_config_initialize_default(void *config);
//
void ocp_nlp_cost_external_initialize(void *config_, void *dims, void *model_, void *opts_,
                                      void *mem_, void *work_);
//
void ocp_nlp_cost_external_update_qp_matrices(void *config_, void *dims, void *model_, void *opts_,
                                              void *memory_, void *work_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COST_EXTERNAL_H_
