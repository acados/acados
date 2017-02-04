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
#include <string.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_i_aux.h"

#include "acados/utils/allocate_ocp_qp.h"

static void allocate_ocp_qp_in_basic(int_t N, int_t *nx, int_t *nu, ocp_qp_in *const qp) {
    int_zeros((int_t **) &qp->nx, N+1, 1);
    int_zeros((int_t **) &qp->nu, N+1, 1);
    int_zeros((int_t **) &qp->nb, N+1, 1);
    int_zeros((int_t **) &qp->nc, N+1, 1);

    qp->N = N;
    memcpy((void *) qp->nx, (void *) nx, sizeof(*nx)*(N+1));
    memcpy((void *) qp->nu, (void *) nu, sizeof(*nu)*(N));

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
    d_zeros((real_t **) &qp->Q[N], nx[N], nx[N]);
    d_zeros((real_t **) &qp->q[N], nx[N], 1);
}


static void free_ocp_qp_in_basic(ocp_qp_in *const qp) {
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


static void allocate_ocp_qp_in_bounds(int_t N, int_t *nb, ocp_qp_in *const qp) {
    qp->lb = (const real_t **) malloc(sizeof(*qp->lb) * (N+1));
    qp->ub = (const real_t **) malloc(sizeof(*qp->ub) * (N+1));
    qp->idxb = (const int_t **) malloc(sizeof(*qp->idxb) * (N+1));

    qp-> N = N;
    memcpy((void *) qp->nb, (void *) nb, sizeof(*nb)*(N+1));

    for (int_t ii = 0; ii <= N; ii++) {
        int_zeros((int_t **) &qp->idxb[ii], nb[ii], 1);
        d_zeros((real_t **) &qp->lb[ii], nb[ii], 1);
        d_zeros((real_t **) &qp->ub[ii], nb[ii], 1);
    }
}


static void free_ocp_qp_in_bounds(ocp_qp_in *const qp) {
    for (int_t ii = 0; ii < qp->N+1; ii++) {
        d_free((real_t*)qp->lb[ii]);
        d_free((real_t*)qp->ub[ii]);
        int_free((int_t*)qp->idxb[ii]);
    }
    free((real_t**)qp->lb);
    free((real_t**)qp->ub);
    free((int_t**)qp->idxb);
}


static void allocate_ocp_qp_in_polyhedral(int_t N, int_t *nc, ocp_qp_in *const qp) {
    qp->lc = (const real_t **) malloc(sizeof(*qp->lc) * (N+1));
    qp->uc = (const real_t **) malloc(sizeof(*qp->uc) * (N+1));
    qp->Cx = (const real_t **) malloc(sizeof(*qp->Cx) * (N+1));
    qp->Cu = (const real_t **) malloc(sizeof(*qp->Cu) * (N));

    qp-> N = N;
    memcpy((void *) qp->nc, (void *) nc, sizeof(*nc)*(N+1));

    for (int_t ii = 0; ii <= N; ii++) {
        d_zeros((real_t **) &qp->lc[ii], nc[ii], 1);
        d_zeros((real_t **) &qp->uc[ii], nc[ii], 1);
        d_zeros((real_t **) &qp->Cx[ii], nc[ii], qp->nx[ii]);
        if (ii < N) d_zeros((real_t **) &qp->Cu[ii], nc[ii], qp->nu[ii]);
    }
}


static void free_ocp_qp_in_polyhedral(ocp_qp_in *const qp) {
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


void allocate_ocp_qp_in(int_t N, int_t *nx, int_t *nu, int_t *nb, int_t *nc, ocp_qp_in *const qp) {
    allocate_ocp_qp_in_basic(N, nx, nu, qp);
    allocate_ocp_qp_in_bounds(N, nb, qp);
    allocate_ocp_qp_in_polyhedral(N, nc, qp);
}


void free_ocp_qp_in(ocp_qp_in *const qp) {
    free_ocp_qp_in_basic(qp);
    free_ocp_qp_in_bounds(qp);
    free_ocp_qp_in_polyhedral(qp);
}


void allocate_ocp_qp_out(ocp_qp_in *const in, ocp_qp_out *out) {
    int_t kk;
    int_t nPrimalVars = 0;
    int_t nPis = 0;
    int_t nLambdas = 0;
    int_t accPrimal = 0;
    int_t accPis = 0;
    int_t accLamdas = 0;
    int_t N = in->N;
    real_t *primalVars, *pis, *lambdas;

    for (kk = 0; kk < N; kk++) {
        nPrimalVars += in->nx[kk] + in->nu[kk];
        nPis += in->nx[kk+1];
        nLambdas += 2*(in->nb[kk] + in->nc[kk]);
    }
    nLambdas += 2*(in->nb[N] + in->nc[N]);
    nPrimalVars += in->nx[N];
    // printf("\nProblem with:\n");
    // printf("- %d\tprimal variables.\n", nPrimalVars);
    // printf("- %d\tmultipliers of eq. constraints \n", nPis);
    // printf("- %d\tmultipliers of ineq. constraints \n", nLambdas);

    d_zeros(&primalVars, nPrimalVars, 1);
    d_zeros(&pis, nPis, 1);
    d_zeros(&lambdas, nLambdas, 1);

    // TODO(dimitris): Reference writes N+1 pointers for u, is it already updated?
    out->x = (real_t **) malloc(sizeof(*out->x)*(N+1));
    out->u = (real_t **) malloc(sizeof(*out->u)*N);
    out->pi = (real_t **) malloc(sizeof(*out->pi)*N);
    out->lam = (real_t **) malloc(sizeof(*out->lam)*(N+1));

    for (kk = 0; kk < N; kk++) {
        out->x[kk] = &primalVars[accPrimal];
        out->u[kk] = &primalVars[accPrimal+in->nx[kk]];
        out->pi[kk] = &pis[accPis];
        out->lam[kk] = &lambdas[accLamdas];
        accPrimal += in->nx[kk] + in->nu[kk];
        accPis += in->nx[kk+1];
        accLamdas += 2*(in->nb[kk] + in->nc[kk]);
    }
    out->x[N] = &primalVars[accPrimal];
    out->lam[N] = &lambdas[accLamdas];
}


void free_ocp_qp_out(ocp_qp_out *out) {
    d_free((real_t*)out->x[0]);
    d_free((real_t*)out->pi[0]);
    d_free((real_t*)out->lam[0]);
    free(out->lam);
    free(out->pi);
    free(out->u);
    free(out->x);
}
