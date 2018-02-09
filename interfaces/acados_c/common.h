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

#ifndef ACADOS_C_COMMON_H_
#define ACADOS_C_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/utils/timing.h>
#include <acados/utils/types.h>

// TODO(nielsvd): is this struct still needed?
typedef struct {
    int (*fun)(void *);
    int (*calculate_args_size)(void *);
    void *(*assign_args)(void *);
    void *(*copy_args)(void *);
    void (*initialize_default_args)(void *);
    int (*calculate_memory_size)(void *);
    void *(*assign_memory)(void *);
    int (*calculate_workspace_size)(void *);
} module_fcn_ptrs;

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_COMMON_H_