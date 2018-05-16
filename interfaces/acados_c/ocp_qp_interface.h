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

#ifndef INTERFACES_ACADOS_C_OCP_QP_INTERFACE_H_
#define INTERFACES_ACADOS_C_OCP_QP_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"

typedef enum {
    PARTIAL_CONDENSING_HPIPM,
#ifdef ACADOS_WITH_HPMPC
    PARTIAL_CONDENSING_HPMPC,
#endif
#ifdef ACADOS_WITH_OOQP
    PARTIAL_CONDENSING_OOQP,
#endif
#ifdef ACADOS_WITH_QPDUNES
    PARTIAL_CONDENSING_QPDUNES,
#endif
    FULL_CONDENSING_HPIPM,
#ifdef ACADOS_WITH_QPOASES
    FULL_CONDENSING_QPOASES,
#endif
#ifdef ACADOS_WITH_QORE
    FULL_CONDENSING_QORE
#endif
} ocp_qp_solver_t;

typedef struct
{
    ocp_qp_solver_t qp_solver;
} ocp_qp_solver_plan;

typedef struct
{
    ocp_qp_xcond_solver_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} ocp_qp_solver;

//
ocp_qp_xcond_solver_config *ocp_qp_config_create(ocp_qp_solver_plan plan);
//
ocp_qp_dims *ocp_qp_dims_create(int N);
//
ocp_qp_in *ocp_qp_in_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims);
//
ocp_qp_out *ocp_qp_out_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims);
//
void *ocp_qp_opts_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims);
//
int ocp_qp_calculate_size(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_);
//
ocp_qp_solver *ocp_qp_assign(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_,
                             void *raw_memory);
//
ocp_qp_solver *ocp_qp_create(ocp_qp_xcond_solver_config *config, ocp_qp_dims *dims, void *opts_);
//
int ocp_qp_solve(ocp_qp_solver *solver, ocp_qp_in *qp_in, ocp_qp_out *qp_out);
//
void ocp_qp_inf_norm_residuals(ocp_qp_dims *dims, ocp_qp_in *qp_in, ocp_qp_out *qp_out,
                               double *res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_OCP_QP_INTERFACE_H_
