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
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "acados/utils/print.h"
#include "OOQP/include/cQpGenSparse.h"
// #include "blasfeo/include/blasfeo_i_aux.h"

int_t ocp_qp_ooqp_create_workspace(const ocp_qp_in *in, ocp_qp_ooqp_workspace *work) {
    int_t ii;
    int_t N = in->N;
    work->nx = 0;    // # of primal optimization variables
    work->nnzQ = 0;  // # non-zeros in lower part of Hessian
    work->nnzA = 0;  // TODO(dimitris): # non-zeros in matrix of equality constraints
    work->nnzC = 0;  // TODO(dimitris): # non-zeros in matrix of inequality constraints
    work->my = 0;    // # of equality constraints
    work->mz = 0;    // # of inequality constraints

    for (ii = 0; ii < N; ii++) {
        work->nx += in->nx[ii] + in->nu[ii];
        work->nnzQ += (in->nx[ii]*in->nx[ii] - in->nx[ii])/2 + in->nx[ii];
        work->nnzQ += (in->nu[ii]*in->nu[ii] - in->nu[ii])/2 + in->nu[ii];
        work->nnzQ += in->nx[ii]*in->nu[ii];
        work->my += N*in->nx[ii+1];  // TODO(dimitris): double check this
        work->mz += in->nb[ii];
    }
    work->nx += in->nx[N];
    work->nnzQ += (in->nx[N]*in->nx[N] -in->nx[N])/2 + in->nx[N];
    work->mz += in->nb[N];

    // TODO(dimitris): Check with giaf how x0 constr. is handled

    // memory allocation
    newQpGenSparse(&work->c, work->nx,
        &work->irowQ, work->nnzQ, &work->jcolQ, &work->dQ,
        &work->xlow, &work->ixlow, &work->xupp, &work->ixupp,
        &work->irowA, work->nnzA, &work->jcolA, &work->dA,
        &work-> bA, work->my,
        &work->irowC, work->nnzC, &work->jcolC, &work->dC,
        &work->clow, work->mz, &work->iclow, &work->cupp, &work->icupp, &work->ierr);

    work->x = malloc(sizeof(*work->x)*work->nx);
    work->y = malloc(sizeof(*work->y)*work->my);
    work->lambda = malloc(sizeof(*work->lambda)*work->mz);
    work->pi = malloc(sizeof(*work->pi)*work->mz);
    work->gamma = malloc(sizeof(*work->gamma)*work->nx);
    work->phi = malloc(sizeof(*work->phi)*work->nx);

    return work->ierr;
}

void ocp_qp_ooqp_free_workspace(ocp_qp_ooqp_workspace *work) {
    freeQpGenSparse(&work->c,
        &work->irowQ, &work->jcolQ, &work->dQ,
        &work->xlow, &work->ixlow, &work->xupp, &work->ixupp,
        &work->irowA, &work->jcolA, &work->dA,
        &work-> bA,
        &work->irowC, &work->jcolC, &work->dC,
        &work->clow, &work->iclow, &work->cupp, &work->icupp);
}

static void fill_in_workspace(const ocp_qp_in *in, ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn;
    int_t offset = 0;  // offset on row/column index of Hessian

    nn = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk]; ii++) work->c[nn++] = in->q[kk][ii];
        for (ii = 0; ii < in->nu[kk]; ii++) work->c[nn++] = in->r[kk][ii];
    }
    for (ii = 0; ii < in->nx[in->N]; ii++) work->c[nn++] = in->q[in->N][ii];

    nn = 0;
    // TODO(dimitris): For the moment I assume full matrices Q,R etc (we need to def. sparsities)
    for (kk = 0; kk < in->N; kk++) {
        // printf("===== @stage %d\n", kk);
        for (jj = 0; jj< in->nx[kk]; jj++) {
            // printf("===== @column %d of Q\n", jj);
            for (ii = jj; ii < in->nx[kk]; ii++) {  // we write only the lower triangular part
                // printf("===== @row %d of Q\n", jj);
                work->dQ[nn] = in->Q[kk][jj*in->nx[kk]+ii];
                work->irowQ[nn] = offset + ii;
                work->jcolQ[nn] = offset + jj;
                // printf("===== Q[%d, %d] = %f\n", work->irowQ[nn]+1, work->jcolQ[nn]+1,work->dQ[nn]);
                nn += 1;
            }
        }
        for (jj = 0; jj< in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nu[kk]; ii++) {
                work->dQ[nn] = in->S[kk][jj*in->nu[kk]+ii];
                work->irowQ[nn] = offset + in->nx[kk] + ii;
                work->jcolQ[nn] = offset + jj;
                nn += 1;
            }
        }
        for (jj = 0; jj< in->nu[kk]; jj++) {
            for (ii = jj; ii < in->nu[kk]; ii++) {
                work->dQ[nn] = in->R[kk][jj*in->nu[kk]+ii];
                work->irowQ[nn] = offset + in->nx[kk] + ii;
                work->jcolQ[nn] = offset + in->nx[kk] + jj;
                nn += 1;
            }
        }
        offset += kk*(in->nx[kk]+in->nu[kk]);
    }
    for (jj = 0; jj< in->nx[in->N]; jj++) {
        for (ii = jj; ii < in->nx[in->N]; ii++) {
            work->dQ[nn] = in->Q[in->N][jj*in->nx[in->N]+ii];
            work->irowQ[nn] = offset + ii;
            work->jcolQ[nn] = offset + jj;
            nn += 1;
        }
    }
    doubleLexSort( work->irowQ, work->nnzQ, work->jcolQ, work-> dQ);
    // for (ii = 0; ii < nn; ii++) {
    //     printf("===== Q[%d, %d] = %f\n", work->irowQ[ii]+1, work->jcolQ[ii]+1,work->dQ[ii]);
    // }
    // printf("===== c:\n");
    // for(ii = 0; ii < work->nx; ii++) printf("%f\n", work->c[ii]);

    work->print_level = 20;  // TODO(dimitris): set to zero when no debugging
}

int_t ocp_qp_ooqp(ocp_qp_in *in, ocp_qp_out *out, void *args_, void *work_) {
    ocp_qp_ooqp_workspace *work = (ocp_qp_ooqp_workspace *) work_;
    fill_in_workspace(in, work);
    // printf("WORKSPACE FILLED IN++++++++++++++++++++++++\n");

    // call sparse OOQP
    qpsolvesp(work->c, work->nx,
        work->irowQ, work->nnzQ, work->jcolQ, work->dQ,
        work->xlow, work->ixlow, work->xupp, work->ixupp,
        work->irowA, work->nnzA, work->jcolA, work->dA,
        work->bA, work->my,
        work->irowC, work->nnzC, work->jcolC, work->dC,
        work->clow, work->mz, work->iclow, work->cupp, work->icupp,
        work->x, work->gamma, work->phi, work->y, work->z, work->lambda, work->pi,
        &work->objectiveValue, work->print_level, &work->ierr);

    printf("===================== OOQP FLAG: %d\n", work->ierr);
    printf("===================== OOQP SOLUTION:\n");
    for (int ii=0; ii<work->nx;ii++) printf("%f\n",work->x[ii]);


return work->ierr;
}
