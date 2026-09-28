/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include "acados/ocp_nlp/ocp_nlp_reg_noreg.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/utils/math.h"

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_nlp_reg_noreg_opts_calculate_size(void)
{
    acados_size_t size = 0;

    return size;
}



void *ocp_nlp_reg_noreg_opts_assign(void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_reg_noreg_opts_initialize_default(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{

    return;
}



void ocp_nlp_reg_noreg_opts_set(void *config_, void *opts_, const char *field, void* value)
{

    printf("\nerror: field %s not available in ocp_nlp_reg_noreg_opts_set\n", field);
    exit(1);

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_reg_noreg_memory_calculate_size(void *config_, ocp_nlp_reg_dims *dims, void *opts_)
{
    acados_size_t size = 0;

    return size;
}



void *ocp_nlp_reg_noreg_memory_assign(void *config_, ocp_nlp_reg_dims *dims, void *opts_, void *raw_memory)
{
    return raw_memory;
}



void ocp_nlp_reg_noreg_memory_set(void *config_, ocp_nlp_reg_dims *dims, void *memory_, char *field, void *value)
{
    if(!strcmp(field, "RSQrq_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "rq_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "BAbt_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "b_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "idxb_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "DCt_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "ux_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "pi_ptr"))
    {
        // no-op
    }
    else if(!strcmp(field, "lam_ptr"))
    {
        // no-op
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_reg_noreg_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * functions
 ************************************************/

void ocp_nlp_reg_noreg_regularize(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    return;
}


void ocp_nlp_reg_noreg_regularize_lhs(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    return;
}


void ocp_nlp_reg_noreg_regularize_rhs(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    return;
}

void ocp_nlp_reg_noreg_correct_dual_sol(void *config, ocp_nlp_reg_dims *dims, void *opts_, void *mem_)
{
    return;
}



void ocp_nlp_reg_noreg_config_initialize_default(ocp_nlp_reg_config *config)
{
    // dims
    config->dims_calculate_size = &ocp_nlp_reg_dims_calculate_size;
    config->dims_assign = &ocp_nlp_reg_dims_assign;
    config->dims_set = &ocp_nlp_reg_dims_set;
    // opts
    config->opts_calculate_size = &ocp_nlp_reg_noreg_opts_calculate_size;
    config->opts_assign = &ocp_nlp_reg_noreg_opts_assign;
    config->opts_initialize_default = &ocp_nlp_reg_noreg_opts_initialize_default;
    config->opts_set = &ocp_nlp_reg_noreg_opts_set;
    // memory
    config->memory_calculate_size = &ocp_nlp_reg_noreg_memory_calculate_size;
    config->memory_assign = &ocp_nlp_reg_noreg_memory_assign;
    config->memory_set = &ocp_nlp_reg_noreg_memory_set;
    // functions
    config->regularize = &ocp_nlp_reg_noreg_regularize;
    config->regularize_lhs = &ocp_nlp_reg_noreg_regularize_lhs;
    config->regularize_rhs = &ocp_nlp_reg_noreg_regularize_rhs;
    config->correct_dual_sol = &ocp_nlp_reg_noreg_correct_dual_sol;
}

