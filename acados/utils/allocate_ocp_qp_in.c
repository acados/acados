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

#include <stdio.h>
#include <stdlib.h>

#include "acados/utils/allocate_ocp_qp_in.h"
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_i_aux.h"

static void allocate_QP_basic(int_t N, int_t *nx, int_t *nu, ocp_qp_in *const qp) {
    qp->N = N;
    int_zeros((int_t **) &qp->nx, N+1, 1);
    int_zeros((int_t **) &qp->nu, N, 1);
    int_zeros((int_t **) &qp->nb, N+1, 1);
    int_zeros((int_t **) &qp->nc, N+1, 1);

    // TODO(dimitris): change it to p_zeros once implemented
    qp->A = (const real_t **) malloc(sizeof(*qp->A) * N);
    qp->B = (const real_t **) malloc(sizeof(*qp->B) * N);
    qp->b = (const real_t **) malloc(sizeof(*qp->b) * N);
    qp->Q = (const real_t **) malloc(sizeof(*qp->Q) * (N+1));
    qp->S = (const real_t **) malloc(sizeof(*qp->S) * N);
    qp->R = (const real_t **) malloc(sizeof(*qp->R) * N);
    qp->q = (const real_t **) malloc(sizeof(*qp->q) * (N+1));
    qp->r = (const real_t **) malloc(sizeof(*qp->r) * N);
    for (int_t ii = 0; ii < N; ii++) {
        d_zeros((real_t **) &qp->A[ii], nx[ii], nx[ii]);
        d_zeros((real_t **) &qp->B[ii], nx[ii], nu[ii]);
        d_zeros((real_t **) &qp->b[ii], nx[ii], 1);
        d_zeros((real_t **) &qp->Q[ii], nx[ii], nx[ii]);
        d_zeros((real_t **) &qp->S[ii], nu[ii], nx[ii]);
        d_zeros((real_t **) &qp->R[ii], nu[ii], nu[ii]);
        d_zeros((real_t **) &qp->q[ii], nx[ii], 1);
        d_zeros((real_t **) &qp->r[ii], nu[ii], 1);
    }
    d_zeros((real_t **) &qp->q[N], nx[N], 1);
    d_zeros((real_t **) &qp->Q[N], nx[N], nx[N]);
}


static void free_QP_basic(ocp_qp_in *const qp) {
    int_free((int_t *)qp->nx);
    int_free((int_t *)qp->nu);
    int_free((int_t *)qp->nb);
    int_free((int_t *)qp->nc);

    for (int_t i = 0; i < qp->N; i++) {
        d_free((real_t*)qp->A[i]);
        d_free((real_t*)qp->B[i]);
        d_free((real_t*)qp->b[i]);
        d_free((real_t*)qp->Q[i]);
        d_free((real_t*)qp->S[i]);
        d_free((real_t*)qp->R[i]);
        d_free((real_t*)qp->q[i]);
        d_free((real_t*)qp->r[i]);
    }
    d_free((real_t*)qp->Q[qp->N]);
    d_free((real_t*)qp->q[qp->N]);

    free((real_t**)qp->A);
    free((real_t**)qp->B);
    free((real_t**)qp->b);
    free((real_t**)qp->Q);
    free((real_t**)qp->S);
    free((real_t**)qp->R);
    free((real_t**)qp->q);
    free((real_t**)qp->r);
}


static void allocate_QP_bounds(int_t N, int_t *nb, ocp_qp_in *const qp) {
    qp->lb = (const real_t **) malloc(sizeof(*qp->lb) * (N+1));
    qp->ub = (const real_t **) malloc(sizeof(*qp->ub) * (N+1));
    qp->idxb = (const int_t **) malloc(sizeof(*qp->idxb) * (N+1));

    for (int_t ii = 0; ii < N+1; ii++) {
        int_zeros((int_t **) &qp->idxb[ii], nb[ii], 1);
        d_zeros((real_t **) &qp->lb[ii], nb[ii], 1);
        d_zeros((real_t **) &qp->ub[ii], nb[ii], 1);
    }
}


static void free_QP_bounds(ocp_qp_in *const qp) {
    for (int_t ii = 0; ii < qp->N+1; ii++) {
        d_free((real_t*)qp->lb[ii]);
        d_free((real_t*)qp->ub[ii]);
        int_free((int_t*)qp->idxb[ii]);
    }
    free((real_t**)qp->lb);
    free((real_t**)qp->ub);
    free((int_t**)qp->idxb);
}


static void allocate_QP_polyhedral(int_t N, int_t *nc, ocp_qp_in *const qp) {
    qp->lc = (const real_t **) malloc(sizeof(*qp->lb) * (N+1));
    qp->uc = (const real_t **) malloc(sizeof(*qp->ub) * (N+1));
    qp->Cx = (const real_t **) malloc(sizeof(*qp->Cx) * (N+1));
    qp->Cu = (const real_t **) malloc(sizeof(*qp->Cu) * (N));

    for (int_t ii = 0; ii < N; ii++) {
        d_zeros((real_t **) &qp->lc[ii], nc[ii], 1);
        d_zeros((real_t **) &qp->uc[ii], nc[ii], 1);
        d_zeros((real_t **) &qp->Cx[ii], nc[ii], 1);
        d_zeros((real_t **) &qp->Cu[ii], nc[ii], 1);
    }
    d_zeros((real_t **) &qp->lc[N], nc[N], 1);
    d_zeros((real_t **) &qp->uc[N], nc[N], 1);
    d_zeros((real_t **) &qp->Cx[N], nc[N], 1);
}


static void free_QP_polyhedral(ocp_qp_in *const qp) {
    for (int_t ii = 0; ii < qp->N; ii++) {
        d_free((real_t*)qp->lc[ii]);
        d_free((real_t*)qp->uc[ii]);
        d_free((real_t*)qp->Cx[ii]);
        d_free((real_t*)qp->Cu[ii]);
    }
    d_free((real_t*)qp->lc[qp->N]);
    d_free((real_t*)qp->uc[qp->N]);
    d_free((real_t*)qp->Cx[qp->N]);

    free((real_t**)qp->lc);
    free((real_t**)qp->uc);
    free((real_t**)qp->Cx);
    free((real_t**)qp->Cu);
}


static void allocate_QP_x0(int_t nx0, ocp_qp_in *const qp) {
    qp->lb = (const real_t **) malloc(sizeof(*qp->lb));
    qp->ub = (const real_t **) malloc(sizeof(*qp->ub));
    qp->idxb = (const int_t **) malloc(sizeof(*qp->idxb));
    d_zeros((real_t **) &qp->lb[0], nx0, 1);
    d_zeros((real_t **) &qp->ub[0], nx0, 1);
    int_zeros((int_t **) &qp->idxb[0], nx0, 1);
}


static void free_QP_x0(ocp_qp_in *const qp) {
    d_free((real_t *)qp->lb[0]);
    d_free((real_t *)qp->ub[0]);
    int_free((int_t *)qp->idxb[0]);
    free((real_t **)qp->lb);
    free((real_t **)qp->ub);
    free((int_t **)qp->idxb);
}


void allocate_QP_full(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *nc, ocp_qp_in *const qp) {
    allocate_QP_basic(N, nx, nu, qp);
    allocate_QP_bounds(N, nb, qp);
    allocate_QP_polyhedral(N, nc, qp);
}


void free_QP_full(ocp_qp_in *const qp) {
    free_QP_basic(qp);
    free_QP_bounds(qp);
    free_QP_polyhedral(qp);
}


void allocate_QP_unconstrained(int_t N, int_t *nx, int_t *nu, ocp_qp_in *const qp) {
    allocate_QP_basic(N, nx, nu, qp);
    allocate_QP_x0(nx[0], qp);
}


void free_QP_unconstrained(ocp_qp_in *const qp) {
    free_QP_basic(qp);
    free_QP_x0(qp);
}

// // toy code for validation (and valgrind)
// int main( ) {
//     /* code */
//     int_t N = 3;
//     int_t nx[] = {3, 2, 2, 4};
//     int_t nu[] = {3, 2, 2};
//     int_t nb[] = {1, 1, 1, 1};
//     int_t nc[] = {1, 2, 0, 1};
//
//     ocp_qp_in qp_constrained;
//     ocp_qp_in qp_unconstrained;
//
//     allocate_QP_full(N, nx, nu, nb, nc, &qp_constrained);
//     allocate_QP_unconstrained(N, nx, nu, &qp_unconstrained);
//
//     free_QP_full(&qp_constrained);
//     free_QP_unconstrained(&qp_unconstrained);
//
//     return 0;
// }
