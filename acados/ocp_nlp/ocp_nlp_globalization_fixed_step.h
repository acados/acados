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


/// \addtogroup ocp_nlp
/// @{
/// \addtogroup ocp_nlp_globalization
/// @{

#ifndef ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FIXED_STEP_H_
#define ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FIXED_STEP_H_

#ifdef __cplusplus
extern "C" {
#endif

// blasfeo
#include "blasfeo/include/blasfeo_common.h"

// acados
#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/utils/types.h"

/************************************************
 * options
 ************************************************/
typedef struct
{
    ocp_nlp_globalization_opts *globalization_opts;

} ocp_nlp_globalization_fixed_step_opts;

//
acados_size_t ocp_nlp_globalization_fixed_step_opts_calculate_size(void *config, void *dims);
//
void *ocp_nlp_globalization_fixed_step_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_nlp_globalization_fixed_step_opts_initialize_default(void *config, void *dims, void *opts);
//
void ocp_nlp_globalization_fixed_step_opts_update(void *config, void *dims, void *opts);
//
void ocp_nlp_globalization_fixed_step_opts_set(void *config, void *opts, const char *field, void* value);


/************************************************
 * memory
 ************************************************/

// TODO: remove stufF!
typedef struct
{
    double step_norm;
    double alpha;
} ocp_nlp_globalization_fixed_step_memory;

//
acados_size_t ocp_nlp_globalization_fixed_step_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_nlp_globalization_fixed_step_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//

/************************************************
 * functions
 ************************************************/

int ocp_nlp_globalization_fixed_step_find_acceptable_iterate(ocp_nlp_config *nlp_config,
                            ocp_nlp_dims *nlp_dims,
                            ocp_nlp_in *nlp_in,
                            ocp_nlp_out *nlp_out,
                            ocp_nlp_memory *nlp_mem,
                            ocp_nlp_workspace *nlp_work,
                            ocp_nlp_opts *nlp_opts);


void ocp_nlp_globalization_fixed_step_config_initialize_default(ocp_nlp_globalization_config *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_GLOBALIZATION_FIXED_STEP_H_
/// @}
/// @}