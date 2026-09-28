/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include "acados/ocp_nlp/ocp_nlp_dynamics_common.h"

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

acados_size_t ocp_nlp_dynamics_config_calculate_size()
{
    acados_size_t size = 0;

    size += sizeof(ocp_nlp_dynamics_config);

    size += sim_config_calculate_size();

    return size;
}



ocp_nlp_dynamics_config *ocp_nlp_dynamics_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    ocp_nlp_dynamics_config *config = (ocp_nlp_dynamics_config *) c_ptr;
    c_ptr += sizeof(ocp_nlp_dynamics_config);

    config->sim_solver = sim_config_assign(c_ptr);
    c_ptr += sim_config_calculate_size();

    return config;
}
