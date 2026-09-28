/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef INTERFACES_ACADOS_C_CONDENSING_INTERFACE_H_
#define INTERFACES_ACADOS_C_CONDENSING_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/ocp_qp/ocp_qp_partial_condensing.h"

typedef enum {
    PARTIAL_CONDENSING,
    FULL_CONDENSING,
} condensing_t;

typedef struct
{
    condensing_t condensing_type;
} condensing_plan;

typedef struct
{
    ocp_qp_xcond_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} condensing_module;

ocp_qp_xcond_config *ocp_qp_condensing_config_create(condensing_plan *plan);
//
void *ocp_qp_condensing_opts_create(ocp_qp_xcond_config *config, void *dims_);
//
acados_size_t ocp_qp_condensing_calculate_size(ocp_qp_xcond_config *config, void *dims_, void *opts_);
//
condensing_module *ocp_qp_condensing_assign(ocp_qp_xcond_config *config, void *dims_,
                                            void *opts_, void *raw_memory);
//
condensing_module *ocp_qp_condensing_create(ocp_qp_xcond_config *config, void *dims_,
                                            void *opts_);
//
int ocp_qp_condense(condensing_module *module, void *qp_in, void *qp_out);
//
int ocp_qp_expand(condensing_module *module, void *qp_in, void *qp_out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_CONDENSING_INTERFACE_H_
