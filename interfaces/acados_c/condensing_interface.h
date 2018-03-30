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

#ifndef ACADOS_C_CONDENSING_INTERFACE_H_
#define ACADOS_C_CONDENSING_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    PARTIAL_CONDENSING,
    FULL_CONDENSING,
} condensing_t;

typedef struct {
    condensing_t condensing_type;
} condensing_plan;

typedef struct {
    condensing_config *config;  // TODO
    void *dims_;
    void *opts;
    void *mem;
    void *work;
} condensing_module;

condensing_config *condensing_config_create(condensing_plan *plan);
//
void *condensing_opts_create(condensing_config *config, void *dims_);
//
int condensing_calculate_size(condensing_config *config, void *dims_, void *opts_);
//
condensing_module *condensing_assign(condensing_config *config, void *dims_, void *opts_, void *raw_memory);
//
condensing_module *condensing_create(condensing_config *config, void *dims_, void *opts_);
//
int condense(condensing_module *module, void *qp_in, void *qp_out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_CONDENSING_INTERFACE_H_
