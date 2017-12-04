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

#ifndef ACADOS_C_OCP_LINEARIZATION_H_
#define ACADOS_C_OCP_LINEARIZATION_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/ocp_linearization/ocp_linearization_common.h>
#include <acados/utils/types.h>

typedef enum {
    ocp_linearization_HPIPM,
    ocp_linearization_QORE,
    ocp_linearization_QPOASES
} ocp_linearization_method_t;

typedef struct {
    ocp_linearization_method_t linearization_method;
} ocp_linearization_config;

//
int ocp_linearization_calculate_args_size(ocp_linearization_config *config, ocp_linearization_dims *dims);
//
void *ocp_linearization_assign_args(ocp_linearization_config *config, ocp_linearization_dims *dims, void *raw_memory);
//
void *ocp_linearization_create_args(ocp_linearization_config *config, ocp_linearization_dims *dims);
//
void ocp_linearization_assign_default_args(ocp_linearization_config *config, void *args_);
//
int ocp_linearization_calculate_memory_size(ocp_linearization_dims *dims, void *args_);
//
void *ocp_linearization_assign_memory(ocp_linearization_dims *dims, void *args_, void *raw_memory);
//
void *ocp_linearization_create_memory(ocp_linearization_dims *dims, void *args_);
//
int ocp_linearization_calculate_workspace_size(ocp_linearization_dims *dims, void *args_);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_OCP_LINEARIZATION_QP_H_
