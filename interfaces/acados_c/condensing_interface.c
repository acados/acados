/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include "acados_c/condensing_interface.h"

// external
#include <assert.h>
#include <stdlib.h>
#include <string.h>

// acados
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"
#include "acados/utils/mem.h"

ocp_qp_xcond_config *ocp_qp_condensing_config_create(condensing_plan *plan)
{
    acados_size_t bytes = ocp_qp_condensing_config_calculate_size();
    void *ptr = calloc(1, bytes);
    ocp_qp_xcond_config *config = ocp_qp_condensing_config_assign(ptr);

    switch (plan->condensing_type)
    {
        case PARTIAL_CONDENSING:
            ocp_qp_partial_condensing_config_initialize_default(config);
            break;
        case FULL_CONDENSING:
            ocp_qp_full_condensing_config_initialize_default(config);
            break;
    }
    return config;
}

void *ocp_qp_condensing_opts_create(ocp_qp_xcond_config *config, void *dims_)
{
    acados_size_t bytes = config->opts_calculate_size(dims_);

    void *ptr = calloc(1, bytes);

    void *opts = config->opts_assign(dims_, ptr);

    config->opts_initialize_default(dims_, opts);

    return opts;
}

acados_size_t ocp_qp_condensing_calculate_size(ocp_qp_xcond_config *config, void *dims_, void *opts_)
{
    acados_size_t bytes = sizeof(condensing_module);

    bytes += config->memory_calculate_size(dims_, opts_);
    bytes += config->workspace_calculate_size(dims_, opts_);

    return bytes;
}

condensing_module *ocp_qp_condensing_assign(ocp_qp_xcond_config *config, void *dims_,
                                            void *opts_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    condensing_module *module = (condensing_module *) c_ptr;
    c_ptr += sizeof(condensing_module);

    module->config = config;
    module->dims = dims_;
    module->opts = opts_;

    module->mem = config->memory_assign(dims_, opts_, c_ptr);
    c_ptr += config->memory_calculate_size(dims_, opts_);

    module->work = (void *) c_ptr;
    c_ptr += config->workspace_calculate_size(dims_, opts_);

    assert((char *) raw_memory + ocp_qp_condensing_calculate_size(config, dims_, opts_) >= c_ptr);

    return module;
}

condensing_module *ocp_qp_condensing_create(ocp_qp_xcond_config *config, void *dims_,
                                            void *opts_)
{
    config->opts_update(dims_, opts_);
    acados_size_t bytes = ocp_qp_condensing_calculate_size(config, dims_, opts_);

    void *ptr = calloc(1, bytes);

    condensing_module *module = ocp_qp_condensing_assign(config, dims_, opts_, ptr);

    return module;
}

int ocp_qp_condense(condensing_module *module, void *qp_in, void *qp_out)
{
    return module->config->condensing(qp_in, qp_out, module->opts, module->mem, module->work);
}

int ocp_qp_expand(condensing_module *module, void *qp_in, void *qp_out)
{
    return module->config->expansion(qp_in, qp_out, module->opts, module->mem, module->work);
}
