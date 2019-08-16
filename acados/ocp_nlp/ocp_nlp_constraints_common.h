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


/// \ingroup ocp_nlp
/// @{

/// \defgroup ocp_nlp_constraints ocp_nlp_constraints
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/external_function_generic.h"
#include "acados/utils/types.h"



/************************************************
 * config
 ************************************************/

typedef struct
{
    int (*dims_calculate_size)(void *config);
    void *(*dims_assign)(void *config, void *raw_memory);
    void (*dims_initialize)(void *config, void *dims, int nx, int nu, int nbx, int nbu, int ng,
                            int nh, int nq, int ns);
    int (*model_calculate_size)(void *config, void *dims);
    void *(*model_assign)(void *config, void *dims, void *raw_memory);
    int (*model_set)(void *config_, void *dims_, void *model_, const char *field, void *value);
    int (*opts_calculate_size)(void *config, void *dims);
    void *(*opts_assign)(void *config, void *dims, void *raw_memory);
    void (*opts_initialize_default)(void *config, void *dims, void *opts);
    void (*opts_update)(void *config, void *dims, void *opts);
    void (*opts_set)(void *config, void *opts, char *field, void *value);
    int (*memory_calculate_size)(void *config, void *dims, void *opts);
    struct blasfeo_dvec *(*memory_get_fun_ptr)(void *memory);
    struct blasfeo_dvec *(*memory_get_adj_ptr)(void *memory);
    void (*memory_set_ux_ptr)(struct blasfeo_dvec *ux, void *memory);
    void (*memory_set_lam_ptr)(struct blasfeo_dvec *lam, void *memory);
    void (*memory_set_DCt_ptr)(struct blasfeo_dmat *DCt, void *memory);
    void (*memory_set_RSQrq_ptr)(struct blasfeo_dmat *RSQrq, void *memory);
    void (*memory_set_idxb_ptr)(int *idxb, void *memory);
    void (*memory_set_idxs_ptr)(int *idxs, void *memory);
    void *(*memory_assign)(void *config, void *dims, void *opts, void *raw_memory);
    int (*workspace_calculate_size)(void *config, void *dims, void *opts);
    void (*initialize)(void *config, void *dims, void *model, void *opts, void *mem, void *work);
    void (*update_qp_matrices)(void *config, void *dims, void *model, void *opts, void *mem,
                               void *work);
    void (*config_initialize_default)(void *config);
    // dimension setters
    void (*dims_set)(void *config_, void *dims_, const char *field, const int *value);
    void (*dims_get)(void *config_, void *dims_, const char *field, int* value);
} ocp_nlp_constraints_config;

//
int ocp_nlp_constraints_config_calculate_size();
//
ocp_nlp_constraints_config *ocp_nlp_constraints_config_assign(void *raw_memory);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_CONSTRAINTS_COMMON_H_
/// @}
/// @}
