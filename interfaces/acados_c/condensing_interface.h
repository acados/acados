/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
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
int ocp_qp_condensing_calculate_size(ocp_qp_xcond_config *config, void *dims_, void *opts_);
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
