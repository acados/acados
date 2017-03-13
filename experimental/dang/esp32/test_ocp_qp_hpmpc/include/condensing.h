/*    acados/ocp_qp/condensing.h
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

#ifndef ACADOS_OCP_QP_CONDENSING_H_
#define ACADOS_OCP_QP_CONDENSING_H_

#include "ocp_qp_common.h"
#include "types.h"

typedef struct condensing_in_ {
    ocp_qp_in *qp_input;
} condensing_in;

typedef struct condensing_out_ {
    real_t *H;
    real_t *h;
    real_t *lb;
    real_t *ub;
    real_t *A;
    real_t *lbA;
    real_t *ubA;
} condensing_out;

typedef struct condensing_memory_ {
    real_t dummy;
} condensing_memory;

typedef struct condensing_workspace_ {
    int_t nconvars;
    int_t nconstraints;
    int_t *nstate_bounds;
    real_t ***G;
    real_t **g;
    real_t ***D;
    real_t *W1_x;
    real_t *W2_x;
    real_t *W1_u;
    real_t *W2_u;
    real_t *w1;
    real_t *w2;
    real_t *Sx0;
} condensing_workspace;

void condensing_N2_fixed_initial_state(condensing_in *in, condensing_out *out,
    condensing_workspace *work);

void condensingN2_free_initial_state();

#endif  // ACADOS_OCP_QP_CONDENSING_H_
