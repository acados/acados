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


#include "acados/ocp_nlp/ocp_nlp_globalization_common.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"
// acados
#include "acados/utils/mem.h"

/************************************************
 * config
 ************************************************/

acados_size_t ocp_nlp_globalization_config_calculate_size()
{
    return sizeof(ocp_nlp_globalization_config);
}


ocp_nlp_globalization_config *ocp_nlp_globalization_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_nlp_globalization_config *config = (ocp_nlp_globalization_config *) c_ptr;
    c_ptr += sizeof(ocp_nlp_globalization_config);

    return config;
}

/************************************************
 * dims
 ************************************************/
// NOTE: we use ocp_nlp_dims for the globalization module.

/************************************************
 * options
 ************************************************/
acados_size_t ocp_nlp_globalization_opts_calculate_size(void *config, void *dims)
{
    return sizeof(ocp_nlp_globalization_opts);
}

void ocp_nlp_globalization_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_globalization_opts *opts = opts_;

    opts->use_SOC = 0;
    opts->line_search_use_sufficient_descent = 0;
    opts->full_step_dual = 0;
    opts->alpha_min = 0.05;
    opts->alpha_reduction = 0.7;
    opts->eps_sufficient_descent = 1e-4; // Leineweber1999: MUSCOD-I eps_T = 1e-4 (p.89); Note: eps_T = 0.1 originally proposed by Powell 1978 (Leineweber 1999, p. 53)

    return;
}

void *ocp_nlp_globalization_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    // char *c_ptr = (char *) raw_memory;
    // // align_char_to(8, &c_ptr);

    // ocp_nlp_globalization_opts *opts = (ocp_nlp_globalization_opts *) c_ptr;
    // c_ptr += sizeof(ocp_nlp_globalization_opts);

    // assert((char *) raw_memory + ocp_nlp_globalization_opts_calculate_size(config_, dims_) >= c_ptr);

    // return opts;
    return raw_memory;
}


void ocp_nlp_globalization_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_globalization_opts *opts = (ocp_nlp_globalization_opts *) opts_;

    // nlp_globalization_opts
    if (!strcmp(field, "alpha_reduction"))
    {
        double* alpha_reduction = (double *) value;
        opts->alpha_reduction = *alpha_reduction;
    }
    else if (!strcmp(field, "alpha_min"))
    {
        double* alpha_min = (double *) value;
        opts->alpha_min = *alpha_min;
    }
    else if (!strcmp(field, "eps_sufficient_descent"))
    {
        double* eps_sufficient_descent = (double *) value;
        opts->eps_sufficient_descent = *eps_sufficient_descent;
    }
    else if (!strcmp(field, "full_step_dual"))
    {
        int* full_step_dual = (int *) value;
        opts->full_step_dual = *full_step_dual;
    }
    else if (!strcmp(field, "line_search_use_sufficient_descent"))
    {
        int* line_search_use_sufficient_descent = (int *) value;
        opts->line_search_use_sufficient_descent = *line_search_use_sufficient_descent;
    }
    else if (!strcmp(field, "use_SOC"))
    {
        int* use_SOC = (int *) value;
        opts->use_SOC = *use_SOC;
    }
    else
    {
        printf("\nerror: ocp_nlp_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}
