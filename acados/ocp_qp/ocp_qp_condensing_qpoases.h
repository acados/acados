/*
 *    This file is part of ACADOS.
 *
 *    ACADOS is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADOS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADOS; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_
#define ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_

#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"

typedef struct ocp_qp_condensing_qpoases_args_ {
    real_t dummy;
} ocp_qp_condensing_qpoases_args;

int_t ocp_qp_condensing_qpoases(ocp_qp_in *input, ocp_qp_out *output,
    ocp_qp_condensing_qpoases_args *args, real_t *work);

int_t ocp_qp_condensing_qpoases_workspace_size(ocp_qp_in *input,
    ocp_qp_condensing_qpoases_args *args);

void initialise_qpoases(ocp_qp_in *input);

#endif  // ACADOS_OCP_QP_OCP_QP_CONDENSING_QPOASES_H_
