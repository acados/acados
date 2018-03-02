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

#include "acados_c/ocp_qp/ocp_qp_hpmpc.h"

#include <string.h>



void *ocp_qp_hpmpc_copy_args(ocp_qp_dims *dims, void *raw_memory, void *source_)
{
    ocp_qp_hpmpc_args *source = (ocp_qp_hpmpc_args *)source_;
    ocp_qp_hpmpc_args *dest;

    dest = (ocp_qp_hpmpc_args *)ocp_qp_hpmpc_assign_args(dims, raw_memory);

    dest->tol = source->tol;
    dest->max_iter = source->max_iter;
    dest->mu0 = source->mu0;
    dest->warm_start = source->warm_start;
    dest->N2 = source->N2;
    dest->out_iter = source->out_iter;
    dest->sigma_mu = source->sigma_mu;
    dest->alpha_min = source->alpha_min;
    dest->N = source->N;
    dest->M = source->M;

    int_t N = dims->N;
    int_t sz;
    for (int_t i = 0; i <= N; i++) {
        sz = (dims->nu[i] + dims->nx[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->ux0));
        memcpy(dest->ux0[i], source->ux0[i], sz);
    }
    for (int_t i = 1; i <= N; i++) {
        sz = dims->nx[i] * sizeof(**(((ocp_qp_hpmpc_args *)0)->pi0));
        memcpy(dest->pi0[i], source->pi0[i], sz);
    }
    for (int_t i = 0; i <= N; i++) {
        sz = (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->lam0));
        memcpy(dest->lam0[i], source->lam0[i], sz);
    }
    for (int_t i = 0; i <= N; i++) {
        sz = (2 * dims->nb[i] + 2 * dims->ng[i]) * sizeof(**(((ocp_qp_hpmpc_args *)0)->t0));
        memcpy(dest->t0[i], source->t0[i], sz);
    }
    sz = 5 * sizeof(*(((ocp_qp_hpmpc_args *)0)->inf_norm_res));
    memcpy(dest->inf_norm_res, source->inf_norm_res, sz);

    return (void *)dest;
}
