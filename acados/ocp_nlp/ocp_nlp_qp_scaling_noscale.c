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


#include "acados/ocp_nlp/ocp_nlp_qp_scaling_noscale.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "acados/ocp_nlp/ocp_nlp_qp_scaling_common.h"
#include "acados/utils/math.h"

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_nlp_qp_scaling_noscale_opts_calculate_size(void)
{
    acados_size_t size = 0;

    return size;
}



void *ocp_nlp_qp_scaling_noscale_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_qp_scaling_noscale_opts_initialize_default(void *config_, ocp_nlp_qp_scaling_dims *dims, void *opts_)
{

    return;
}



void ocp_nlp_qp_scaling_noscale_opts_set(void *config_, void *opts_, const char *field, void* value)
{

    printf("\nerror: field %s not available in ocp_nlp_qp_scaling_noscale_opts_set\n", field);
    exit(1);

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_qp_scaling_noscale_memory_calculate_size(void *config_, ocp_nlp_qp_scaling_dims *dims, void *opts_)
{
    acados_size_t size = 0;

    return size;
}



void *ocp_nlp_qp_scaling_noscale_memory_assign(void *config_, ocp_nlp_qp_scaling_dims *dims, void *opts_, void *raw_memory)
{
    return raw_memory;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_qp_scaling_noscale_scale_qp(void *config, ocp_nlp_qp_scaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in)
{
    // printf("ocp_nlp_qp_scaling_noscale_scale_qp: nothing to do\n");
    return;
}


void ocp_nlp_qp_scaling_noscale_rescale_solution(void *config, ocp_nlp_qp_scaling_dims *dims, void *opts_, void *mem_, ocp_qp_in *qp_in, ocp_qp_out *qp_out)
{
    // printf("ocp_nlp_qp_scaling_noscale_rescale_solution: nothing to do\n");
    return;
}


void ocp_nlp_qp_scaling_noscale_config_initialize_default(ocp_nlp_qp_scaling_config *config)
{
    // dims
    config->dims_calculate_size = &ocp_nlp_qp_scaling_dims_calculate_size;
    config->dims_assign = &ocp_nlp_qp_scaling_dims_assign;
    config->dims_set = &ocp_nlp_qp_scaling_dims_set;
    // opts
    config->opts_calculate_size = &ocp_nlp_qp_scaling_noscale_opts_calculate_size;
    config->opts_assign = &ocp_nlp_qp_scaling_noscale_opts_assign;
    config->opts_initialize_default = &ocp_nlp_qp_scaling_noscale_opts_initialize_default;
    config->opts_set = &ocp_nlp_qp_scaling_noscale_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_qp_scaling_noscale_memory_calculate_size;
    config->memory_assign = &ocp_nlp_qp_scaling_noscale_memory_assign;
    // functions
    config->scale_qp = &ocp_nlp_qp_scaling_noscale_scale_qp;
    config->rescale_solution = &ocp_nlp_qp_scaling_noscale_rescale_solution;
}

