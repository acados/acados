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

static void calculate_problem_size(const ocp_qp_in *in, ocp_qp_ooqp_args *args, int_t *nx,
    int_t *my, int_t *mz, int_t *nnzQ, int_t *nnzA, int_t *nnzC) {

        int_t kk;
        int_t N = in->N;

        // dummy command
        if (args->fixHessianSparsity) kk = 0;

        *nx = 0;    // # of primal optimization variables
        *nnzQ = 0;  // # non-zeros in lower part of Hessian
        *nnzA = 0;  // # non-zeros in matrix of equality constraints
        *nnzC = 0;  // # non-zeros in matrix of inequality constraints
        *my = 0;    // # of equality constraints
        *mz = 0;    // # of inequality constraints

        for (kk = 0; kk < N; kk++) {
            *nx += in->nx[kk] + in->nu[kk];
            *nnzQ += (in->nx[kk]*in->nx[kk] - in->nx[kk])/2 + in->nx[kk];
            *nnzQ += (in->nu[kk]*in->nu[kk] - in->nu[kk])/2 + in->nu[kk];
            *nnzQ += in->nx[kk]*in->nu[kk];
            *nnzA += in->nx[kk+1]*(in->nx[kk] + in->nu[kk] + 1);
            *nnzC += in->nc[kk]*(in->nx[kk] + in->nu[kk]);
            *my += in->nx[kk+1];
            *mz += in->nc[kk];
        }
        *nx += in->nx[N];
        *nnzQ += (in->nx[N]*in->nx[N] -in->nx[N])/2 + in->nx[N];
        *mz += in->nc[N];
        *my += in->nx[0];  // constraint on x0
        *nnzA += in->nx[0];
        *nnzC += in->nx[N]*in->nc[N];
}


static void fill_in_structs(const ocp_qp_in *in,  const ocp_qp_ooqp_args *args,
    ocp_qp_ooqp_memory *mem) {

    int_t ii, jj, kk, nn;
    int_t offset, offsetRows, offsetCols, lim;

    // TODO(dimitris): For the moment I assume full matrices Q,R,A,B... (we need to def. sparsities)

    // ------- Build objective
    nn = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk]; ii++) mem->c[nn++] = in->q[kk][ii];
        for (ii = 0; ii < in->nu[kk]; ii++) mem->c[nn++] = in->r[kk][ii];
    }
    for (ii = 0; ii < in->nx[in->N]; ii++) mem->c[nn++] = in->q[in->N][ii];

    nn = 0; offset = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (jj = 0; jj< in->nx[kk]; jj++) {
            for (ii = jj; ii < in->nx[kk]; ii++) {  // we write only the lower triangular part
                mem->dQ[nn] = in->Q[kk][jj*in->nx[kk]+ii];
                mem->irowQ[nn] = offset + ii;
                mem->jcolQ[nn] = offset + jj;
                nn += 1;
            }
        }
        for (jj = 0; jj< in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nu[kk]; ii++) {
                mem->dQ[nn] = in->S[kk][jj*in->nu[kk]+ii];
                mem->irowQ[nn] = offset + in->nx[kk] + ii;
                mem->jcolQ[nn] = offset + jj;
                nn += 1;
            }
        }
        for (jj = 0; jj< in->nu[kk]; jj++) {
            for (ii = jj; ii < in->nu[kk]; ii++) {
                mem->dQ[nn] = in->R[kk][jj*in->nu[kk]+ii];
                mem->irowQ[nn] = offset + in->nx[kk] + ii;
                mem->jcolQ[nn] = offset + in->nx[kk] + jj;
                nn += 1;
            }
        }
        offset += in->nx[kk]+in->nu[kk];
    }
    for (jj = 0; jj< in->nx[in->N]; jj++) {
        for (ii = jj; ii < in->nx[in->N]; ii++) {
            mem->dQ[nn] = in->Q[in->N][jj*in->nx[in->N]+ii];
            mem->irowQ[nn] = offset + ii;
            mem->jcolQ[nn] = offset + jj;
            nn += 1;
        }
    }
    doubleLexSortC(mem->irowQ, mem->nnzQ, mem->jcolQ, mem-> dQ);

    // ------- Build equality  constraints
    nn = in->nx[0];
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk+1]; ii++) mem->bA[nn++] = -in->b[kk][ii];
    }

    nn = 0;
    // write I on first nx[0] rows/cols
    for (ii = 0; ii < in->nx[0]; ii++) {
        mem->dA[nn] = 1;
        mem->irowA[nn] = ii;
        mem->jcolA[nn] = ii;
        nn += 1;
    }
    offsetRows = in->nx[0]; offsetCols = 0;
    for (kk = 0; kk < in->N; kk++) {
        // write matrix A[kk] (nx[k+1] x nx[k])
        for (jj = 0; jj< in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk+1]; ii++) {
                // printf("writing A_%d[%d,%d]\n", kk, ii, jj);
                mem->dA[nn] = in->A[kk][jj*in->nx[kk+1]+ii];
                mem->irowA[nn] = offsetRows + ii;
                mem->jcolA[nn] = offsetCols + jj;
                nn += 1;
            }
        }
        // write matrix B[kk] (nx[k+1] x nu[k])
        for (jj = 0; jj< in->nu[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk+1]; ii++) {
                mem->dA[nn] = in->B[kk][jj*in->nx[kk+1]+ii];
                mem->irowA[nn] = offsetRows + ii;
                mem->jcolA[nn] = offsetCols + in->nx[kk] + jj;
                nn += 1;
            }
        }
        // write -I (nx[k+1] x nx[k+1])
        for (jj = 0; jj< in->nx[kk+1]; jj++) {
            mem->dA[nn] = -1;
            mem->irowA[nn] = offsetRows + jj;
            mem->jcolA[nn] = offsetCols + in->nx[kk] + in->nu[kk] + jj;
            nn += 1;
        }
        offsetCols += in->nx[kk] + in->nu[kk];
        offsetRows += in->nx[kk+1];
    }
    doubleLexSortC(mem->irowA, mem->nnzA, mem->jcolA, mem-> dA);

    // ------- Build bounds
    offset = 0;
    for (kk = 0; kk < in->N+1; kk++) {
        nn = 0;
        if (kk < in->N) {
            lim = in->nx[kk]+in->nu[kk];
        } else {
            lim = in->nx[kk];
        }
        for (ii = 0; ii < lim; ii++) {
            if (in->nb[kk] > 0) {
                if (in->idxb[kk][nn] == ii) {  // element has bounds
                    mem->ixlow[offset+ii] = (char)1;  // TODO(dimitris): check if cast is redundant
                    mem->xlow[offset+ii] = in->lb[kk][nn];
                    mem->ixupp[offset+ii] = (char)1;
                    mem->xupp[offset+ii] = in->ub[kk][nn];
                    nn += 1;
                }
            } else {
                mem->ixlow[offset+ii] = (char)0;
                mem->xlow[offset+ii] = 0.0;
                mem->ixupp[offset+ii] = (char)0;
                mem->xupp[offset+ii] = 0.0;
            }
        }
        offset += lim;
    }
    // removing bounds on x0 since it is constrained in the equality constraints
    // TODO(dimitris): skip this for MHE
    for (ii = 0; ii < in->nx[0]; ii++) {
        mem->ixlow[ii] = (char)0;
        mem->ixupp[ii] = (char)0;
        mem->xlow[ii] = 0.0;
        mem->xupp[ii] = 0.0;
    }

    // ------- Build inequality constraints
    nn = 0;
    for (kk = 0; kk < in->N+1; kk++) {
        for (ii = 0; ii < in->nc[kk]; ii++) {
            mem->iclow[nn] = (char) 1;
            mem->clow[nn] = in->lc[kk][ii];
            mem->icupp[nn] = (char) 1;
            mem->cupp[nn] = in->uc[kk][ii];
            nn += 1;
        }
    }

    nn = 0;
    offsetRows = 0; offsetCols = 0;
    for (kk = 0; kk < in->N+1; kk++) {
        // write matrix Cx[k] (nc[k] x nx[k])
        for (jj = 0; jj< in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nc[kk]; ii++) {
                // printf("writing C_%d[%d,%d]\n", kk, ii, jj);
                mem->dC[nn] = in->Cx[kk][jj*in->nc[kk]+ii];
                mem->irowC[nn] = offsetRows + ii;
                mem->jcolC[nn] = offsetCols + jj;
                nn += 1;
            }
        }
        if (kk < in->N) {
            // write matrix Cu[k] (nc[k] x nu[k])
            for (jj = 0; jj< in->nu[kk]; jj++) {
                for (ii = 0; ii < in->nc[kk]; ii++) {
                    mem->dC[nn] = in->Cu[kk][jj*in->nc[kk]+ii];
                    mem->irowC[nn] = offsetRows + ii;
                    mem->jcolC[nn] = offsetCols + in->nx[kk] + jj;
                    nn += 1;
                }
            }
        }
        offsetCols += in->nx[kk] + in->nu[kk];
        offsetRows += in->nc[kk];
    }
    doubleLexSortC(mem->irowC, mem->nnzC, mem->jcolC, mem-> dC);

    mem->print_level = args->printLevel;

    if (mem->firstRun == 0) mem->firstRun = 1;
}


static void print_inputs(ocp_qp_ooqp_memory *mem) {
    // int ii;
    printf("\n----------> OOQP INPUTS <----------\n\n");
    printf("NUMBER OF PRIMAL VARIABLES: %d\n", mem->nx);
    printf("NUMBER OF NON-ZEROS in HESSIAN: %d\n", mem->nnzQ);
    printf("NUMBER OF EQUALITY CONSTRAINTS: %d\n", mem->my);
    printf("NUMBER OF NON-ZEROS in EQUALITIES: %d\n", mem->nnzA);
    printf("NUMBER OF INEQUALITY CONSTRAINTS: %d\n", mem->mz);
    printf("NUMBER OF NON-ZEROS in INEQUALITIES: %d\n", mem->nnzC);
    printf("PRINT LEVEL: %d", mem->print_level);
    printf("\n-----------------------------------\n\n");

    // printf("\nEQUALITY CONSTRAINTS:\n");
    // for (ii = 0; ii < mem->nnzA; ii++) {
    //     printf("=====> A[%d, %d] = %f\n", mem->irowA[ii]+1, mem->jcolA[ii]+1, mem->dA[ii]);
    // }
    // for (ii = 0; ii < mem->my; ii++) printf("===> bA[%d] = %f\n", ii+1, mem->bA[ii]);
    //
    // printf("\nBOUNDS:\n");
    // for (ii = 0; ii < mem->nx; ii++) {
    //     printf("ixlow[%d] = %d \t xlow[%d] = %4.2f \t ixupp[%d] = %d \t xupp[%d] = %4.2f\n",
    //         ii, mem->ixlow[ii], ii, mem->xlow[ii], ii, mem->ixupp[ii], ii, mem->xupp[ii]);
    // }
    //
    // printf("\nINEQUALITY CONSTRAINTS:\n");
    // for (ii = 0; ii < mem->nnzC; ii++) {
    //     printf("=====> C[%d, %d] = %f\n", mem->irowC[ii]+1, mem->jcolC[ii]+1, mem->dC[ii]);
    // }
    // for (ii = 0; ii < mem->mz; ii++) {
    //     printf("===> clow[%d] = %4.2f \t cupp[%d] = %4.2f\n",
    //     ii+1, mem->clow[ii], ii+1, mem->cupp[ii]);
    // }
}


static void print_outputs(ocp_qp_ooqp_memory *mem, ocp_qp_ooqp_workspace *work, int return_value) {
        printf("\n----------> OOQP OUTPUTS <---------\n\n");
        printf("RETURN STATUS: %d\n", return_value);
        printf("OBJECTIVE VALUE: %f\n", work->objectiveValue);
        printf("FIRST AND LAST ELEMENT OF SOLUTION:\n");
        printf("x[0] = %f\n", work->x[0]);
        printf("x[%d] = %f\n", mem->nx, work->x[mem->nx-1]);
        printf("\n----------------------------------\n\n");
}


// TODO(dimitris): change order of functions
static void fill_in_qp_out(ocp_qp_in *in, ocp_qp_out *out, ocp_qp_ooqp_workspace *work) {
    int kk, ii, nn;

    nn = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk]; ii++) out->x[kk][ii] = work->x[nn++];
        for (ii = 0; ii < in->nu[kk]; ii++) out->u[kk][ii] = work->x[nn++];
    }
    for (ii = 0; ii < in->nx[in->N]; ii++) out->x[in->N][ii] = work->x[nn++];

    // TODO(dimitris): fill-in all qp_out fields
}


static void embed_x0(const real_t *x0, int_t nx, ocp_qp_ooqp_memory *mem) {
    int ii;
    for (ii = 0; ii < nx; ii++) mem->bA[ii] = x0[ii];
}


int_t ocp_qp_ooqp_create_memory(const ocp_qp_in *in, void *args_, void *mem_) {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args*) args_;
    ocp_qp_ooqp_memory *mem = (ocp_qp_ooqp_memory *) mem_;

    // TODO(dimitris): only perform actions if firstRun
    int_t return_value;

    calculate_problem_size(in, args, &mem->nx, &mem->my, &mem->mz,
        &mem->nnzQ, &mem->nnzA, &mem->nnzC);

    newQpGenSparse(&mem->c, mem->nx,
        &mem->irowQ, mem->nnzQ, &mem->jcolQ, &mem->dQ,
        &mem->xlow, &mem->ixlow, &mem->xupp, &mem->ixupp,
        &mem->irowA, mem->nnzA, &mem->jcolA, &mem->dA,
        &mem->bA, mem->my,
        &mem->irowC, mem->nnzC, &mem->jcolC, &mem->dC,
        &mem->clow, mem->mz, &mem->iclow, &mem->cupp, &mem->icupp, &return_value);

    return return_value;
}


int_t ocp_qp_ooqp_create_workspace(const ocp_qp_in *in, void *args_, void *work_) {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args*) args_;
    ocp_qp_ooqp_workspace *work = (ocp_qp_ooqp_workspace *) work_;

    int nx, my, mz, nnzQ, nnzA, nnzC;

    // TODO(dimitris): do not call the function twice if memory already initialized before
    calculate_problem_size(in, args, &nx, &my, &mz, &nnzQ, &nnzA, &nnzC);

    work->x = (real_t*)malloc(sizeof(*work->x)*nx);
    work->gamma = (real_t*)malloc(sizeof(*work->gamma)*nx);
    work->phi = (real_t*)malloc(sizeof(*work->phi)*nx);
    work->y = (real_t*)malloc(sizeof(*work->y)*my);
    work->z = (real_t*)malloc(sizeof(*work->z)*mz);
    work->lambda = (real_t*)malloc(sizeof(*work->lambda)*mz);
    work->pi = (real_t*)malloc(sizeof(*work->pi)*mz);

    // TODO(dimitris): implement this
    return 0;
}


void ocp_qp_ooqp_free_workspace(void *work_) {
    ocp_qp_ooqp_workspace *work = (ocp_qp_ooqp_workspace *) work_;
        free(work->x);
        free(work->gamma);
        free(work->phi);
        free(work->y);
        free(work->z);
        free(work->lambda);
        free(work->pi);
}


void ocp_qp_ooqp_free_memory(void *mem_) {
    ocp_qp_ooqp_memory *mem = (ocp_qp_ooqp_memory *) mem_;

    freeQpGenSparse(&mem->c,
        &mem->irowQ, &mem->jcolQ, &mem->dQ,
        &mem->xlow, &mem->ixlow, &mem->xupp, &mem->ixupp,
        &mem->irowA, &mem->jcolA, &mem->dA,
        &mem-> bA,
        &mem->irowC, &mem->jcolC, &mem->dC,
        &mem->clow, &mem->iclow, &mem->cupp, &mem->icupp);
}


int_t ocp_qp_ooqp(ocp_qp_in *in, ocp_qp_out *out, void *args_, void *memory_, void *work_) {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args*) args_;
    ocp_qp_ooqp_memory *mem = (ocp_qp_ooqp_memory *) memory_;
    ocp_qp_ooqp_workspace *work = (ocp_qp_ooqp_workspace *) work_;

    int return_value;

    fill_in_structs(in, args, mem);

    // here I assume that x0 is already encoded in ub/lb
    embed_x0(in->lb[0], in->nx[0], mem);

    if (0) print_inputs(mem);

    // call sparse OOQP
    // TODO(dimitris): implement dense OOQP
    // TODO(dimitris): Discuss with giaf how  to handle x0 constr.

    qpsolvesp(mem->c, mem->nx,
        mem->irowQ, mem->nnzQ, mem->jcolQ, mem->dQ,
        mem->xlow, mem->ixlow, mem->xupp, mem->ixupp,
        mem->irowA, mem->nnzA, mem->jcolA, mem->dA,
        mem->bA, mem->my,
        mem->irowC, mem->nnzC, mem->jcolC, mem->dC,
        mem->clow, mem->mz, mem->iclow, mem->cupp, mem->icupp,
        work->x, work->gamma, work->phi, work->y, work->z, work->lambda, work->pi,
        &work->objectiveValue, mem->print_level, &return_value);

    if (0) print_outputs(mem, work, return_value);
    fill_in_qp_out(in, out, work);

return return_value;
}
