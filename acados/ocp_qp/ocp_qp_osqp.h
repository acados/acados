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

#ifndef ACADOS_OCP_QP_OCP_QP_OSQP_H_
#define ACADOS_OCP_QP_OCP_QP_OSQP_H_

#ifdef __cplusplus
extern "C" {
#endif

// osqp
#include "osqp/include/types.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/types.h"

typedef struct ocp_qp_osqp_opts_
{
    c_int verbose; // enable or disable printing
    c_int polish; // enable or disable polishing
    OSQPSettings *osqp_opts;
} ocp_qp_osqp_opts;


typedef struct ocp_qp_osqp_memory_
{
    c_int first_run;

    c_float *q;
    c_float *l;
    c_float *u;

    c_int P_nnzmax;
    c_int *P_i;
    c_int *P_p;
    c_float *P_x;

    c_int A_nnzmax;
    c_int *A_i;
    c_int *A_p;
    c_float *A_x;

    OSQPData *osqp_data;
    OSQPWorkspace *osqp_work;

} ocp_qp_osqp_memory;

int ocp_qp_osqp_opts_calculate_size(void *config, void *dims);
//
void *ocp_qp_osqp_opts_assign(void *config, void *dims, void *raw_memory);
//
void ocp_qp_osqp_opts_initialize_default(void *config, void *dims, void *opts_);
//
void ocp_qp_osqp_opts_update(void *config, void *dims, void *opts_);
//
int ocp_qp_osqp_memory_calculate_size(void *config, void *dims, void *opts_);
//
void *ocp_qp_osqp_memory_assign(void *config, void *dims, void *opts_, void *raw_memory);
//
int ocp_qp_osqp_workspace_calculate_size(void *config, void *dims, void *opts_);
//
int ocp_qp_osqp(void *config, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_);
//
void ocp_qp_osqp_config_initialize_default(void *config);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_OSQP_H_
