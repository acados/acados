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

int_t ocp_qp_ooqp_create_workspace(const ocp_qp_in *in, ocp_qp_ooqp_workspace *work) {
    int_t kk;
    int_t N = in->N;
    work->nx = 0;    // # of primal optimization variables
    work->nnzQ = 0;  // # non-zeros in lower part of Hessian
    work->nnzA = 0;  // TODO(dimitris): # non-zeros in matrix of equality constraints
    work->nnzC = 0;  // TODO(dimitris): # non-zeros in matrix of inequality constraints
    work->my = 0;    // # of equality constraints
    work->mz = 0;    // # of inequality constraints

    for (kk = 0; kk < N; kk++) {
        work->nx += in->nx[kk] + in->nu[kk];
        work->nnzQ += (in->nx[kk]*in->nx[kk] - in->nx[kk])/2 + in->nx[kk];
        work->nnzQ += (in->nu[kk]*in->nu[kk] - in->nu[kk])/2 + in->nu[kk];
        work->nnzQ += in->nx[kk]*in->nu[kk];
        work->nnzA += in->nx[kk+1]*(in->nx[kk] + in->nu[kk] + 1);
        work->my += in->nx[kk+1];
        work->mz += in->nb[kk];
    }
    work->nx += in->nx[N];
    work->nnzQ += (in->nx[N]*in->nx[N] -in->nx[N])/2 + in->nx[N];
    work->mz += in->nb[N];

    // TODO(dimitris): Discuss with giaf how  to handle x0 constr.

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

        // TODO(dimitris): free rest of workspace! (outputs)
        // TODO(dimitris): ask gianf why all these warnings
        // d_free(work->x);
}

static void fill_in_workspace(const ocp_qp_in *in, ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn;
    int_t offset, offsetRows, offsetCols;

    // TODO(dimitris): For the moment I assume full matrices Q,R,A,B... (we need to def. sparsities)

    // ------- Build objective
    nn = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk]; ii++) work->c[nn++] = in->q[kk][ii];
        for (ii = 0; ii < in->nu[kk]; ii++) work->c[nn++] = in->r[kk][ii];
    }
    for (ii = 0; ii < in->nx[in->N]; ii++) work->c[nn++] = in->q[in->N][ii];

    nn = 0; offset = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (jj = 0; jj< in->nx[kk]; jj++) {
            for (ii = jj; ii < in->nx[kk]; ii++) {  // we write only the lower triangular part
                work->dQ[nn] = in->Q[kk][jj*in->nx[kk]+ii];
                work->irowQ[nn] = offset + ii;
                work->jcolQ[nn] = offset + jj;
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
        offset += in->nx[kk]+in->nu[kk];
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

    // ------- Build equality  constraints
    nn = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk+1]; ii++) work->bA[nn++] = -in->b[kk][ii];
    }

    nn = 0; offsetRows = 0; offsetCols = 0;
    for (kk = 0; kk < in->N; kk++) {
        // write matrix A[kk] (nx[k+1] x nx[k])
        for (jj = 0; jj< in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk+1]; ii++) {
                // printf("writing A_%d[%d,%d]\n", kk, ii, jj);
                work->dA[nn] = in->A[kk][jj*in->nx[kk+1]+ii];
                work->irowA[nn] = offsetRows + ii;
                work->jcolA[nn] = offsetCols + jj;
                nn += 1;
            }
        }
        // write matrix B[kk] (nx[k+1] x nu[k])
        for (jj = 0; jj< in->nu[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk+1]; ii++) {
                work->dA[nn] = in->B[kk][jj*in->nx[kk+1]+ii];
                work->irowA[nn] = offsetRows + ii;
                work->jcolA[nn] = offsetCols + in->nx[kk] + jj;
                nn += 1;
            }
        }
        // write -I (nx[k+1] x nx[k+1])
        for (jj = 0; jj< in->nx[kk+1]; jj++) {
            work->dA[nn] = -1;
            work->irowA[nn] = offsetRows + jj;
            work->jcolA[nn] = offsetCols + in->nx[kk] + in->nu[kk] + jj;
            nn += 1;
        }
        offsetCols += in->nx[kk] + in->nu[kk];
        offsetRows += in->nx[kk+1];
    }
    doubleLexSort( work->irowA, work->nnzA, work->jcolA, work-> dA);
    // for (ii = 0; ii < nn; ii++) {
    //     printf("===== A[%d, %d] = %f\n", work->irowA[ii]+1, work->jcolA[ii]+1,work->dA[ii]);
    // }
    // printf("===== bA:\n");
    // for(ii = 0; ii < work->my; ii++) printf("%f\n", work->bA[ii]);

    // TODO(dimitris): add x0 constr. in bounds!

    work->print_level = 0;  // TODO(dimitris): set to zero when no debugging
}

static void print_ooqp_inputs(ocp_qp_ooqp_workspace *work) {
    printf("\n----------> OOQP INPUTS <----------\n\n");
    printf("NUMBER OF PRIMAL VARIABLES: %d\n", work->nx);
    printf("NUMBER OF NON-ZEROS in HESSIAN: %d\n", work->nnzQ);
    printf("NUMBER OF EQUALITY CONSTRAINTS: %d\n", work->my);
    printf("NUMBER OF NON-ZEROS in EQUALITIES: %d\n", work->nnzA);
    printf("NUMBER OF INEQUALITY CONSTRAINTS: %d\n", work->mz);
    printf("NUMBER OF NON-ZEROS in INEQUALITIES: %d\n", work->nnzC);
    printf("PRINT LEVEL: %d", work->print_level);
    // TODO(dimitris): complete this list
    printf("\n-----------------------------------\n\n");
}

static void print_ooqp_outputs(ocp_qp_ooqp_workspace *work) {
        printf("\n----------> OOQP OUTPUTS <---------\n\n");
        printf("RETURN STATUS: %d\n", work->ierr);
        printf("OBJECTIVE VALUE: %f\n", work->objectiveValue);
        printf("FIRST AND LAST ELEMENT OF SOLUTION:\n");
        printf("x[0] = %f\n", work->x[0]);
        printf("x[%d] = %f\n", work->nx, work->x[work->nx-1]);
        printf("\n----------------------------------\n\n");
}

int_t ocp_qp_ooqp(ocp_qp_in *in, ocp_qp_out *out, void *args_, void *work_) {
    ocp_qp_ooqp_workspace *work = (ocp_qp_ooqp_workspace *) work_;
    fill_in_workspace(in, work);

    printf("!!!!!!!! TEMPORARY REMOVING ALL CONSTRAINTS !!!!!!!!\n");
    //TODO(dimitris): WHERE DO THE BOUND VALUES COME FROM?!
    for (int ii=0; ii < work->nx; ii++) {
        // printf("xlow[%d] before %f\n", ii, work->xlow[ii]);
        work->ixlow[ii] = (char)0;  // TODO(dimitris): cast prob. not needed
        work->ixupp[ii] = (char)0;
        work->xlow[ii] = 0.0;
        work->xupp[ii] = 0.0;
        // printf("xlow[%d] after %f\n", ii, work->xlow[ii]);
    }
    for (int ii=0; ii < work->mz; ii++) {
        work->iclow[ii] = (char)0;
        work->icupp[ii] = (char)0;
        work->clow[ii] = 0.0;
        work->cupp[ii] = 0.0;
    }
    work->mz = 0;

    print_ooqp_inputs(work);

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

    print_ooqp_outputs(work);
    // printf("===================== OOQP SOLUTION:\n");
    // for (int ii=0; ii<work->nx;ii++) printf("%f\n",work->x[ii]);


return work->ierr;
}
