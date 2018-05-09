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

#include "acados/ocp_qp/ocp_qp_ooqp.h"

#include <stdio.h>
#include <stdlib.h>

#include "ooqp/cQpGenSparse.h"

#include "acados/utils/timing.h"

#define TIMINGS \
    0  // 0: do not print any timings inside here
       // 1: print only time to solve QP
       // 2: print detailed timings

int_t *rows;
int_t *cols;
int_t lda;

static int_t max_of_three(int_t a, int_t b, int_t c) {
    int_t ans = a;
    (void)((ans < b) && (ans = b));
    (void)((ans < c) && (ans = c));
    return ans;
}

// comparator for qsort
static int_t comparator(const void *p1, const void *p2) {
    int_t ans1, ans2;
    int_t ind1 = *((int *)p1);
    int_t ind2 = *((int *)p2);

    ans1 = rows[ind1] * lda + cols[ind1];
    ans2 = rows[ind2] * lda + cols[ind2];

    return ans1 - ans2;
}

static void sort_matrix_structure_row_major(int_t *order, int_t *irow, int_t nnz, int_t *jcol,
                                            int_t *tmp) {
    int_t ii;

    for (ii = 0; ii < nnz; ii++) {
        tmp[ii] = irow[order[ii]];
    }
    for (ii = 0; ii < nnz; ii++) {
        irow[ii] = tmp[ii];
    }

    for (ii = 0; ii < nnz; ii++) {
        tmp[ii] = jcol[order[ii]];
    }
    for (ii = 0; ii < nnz; ii++) {
        jcol[ii] = tmp[ii];
    }
}

static void sort_matrix_data_row_major(int_t *order, int_t nnz, real_t *d, real_t *tmp) {
    int_t ii;

    for (ii = 0; ii < nnz; ii++) {
        tmp[ii] = d[order[ii]];
    }
    for (ii = 0; ii < nnz; ii++) {
        d[ii] = tmp[ii];
    }
}

static int_t get_number_of_primal_vars(const ocp_qp_in *in) {
    int_t nx = 0;
    int_t kk;
    for (kk = 0; kk <= in->N; kk++) {
        nx += in->nx[kk] + in->nu[kk];
    }
    return nx;
}

static int_t get_number_of_equalities(const ocp_qp_in *in) {
    int_t my = 0;
    int_t kk;
    for (kk = 0; kk < in->N; kk++) {
        my += in->nx[kk + 1];
    }
    return my;
}

static int_t get_number_of_inequalities(const ocp_qp_in *in) {
    int_t mz = 0;
    int_t kk;
    for (kk = 0; kk < in->N + 1; kk++) {
        mz += in->nc[kk];
    }
    return mz;
}

static int_t get_nnzQ(const ocp_qp_in *in, const ocp_qp_ooqp_args *args) {
    int_t kk;
    int_t nnzQ = 0;

    // dummy command
    if (args->printLevel) kk = 0;

    for (kk = 0; kk <= in->N; kk++) {
        nnzQ += (in->nx[kk] * in->nx[kk] - in->nx[kk]) / 2 + in->nx[kk];
        nnzQ += (in->nu[kk] * in->nu[kk] - in->nu[kk]) / 2 + in->nu[kk];
        nnzQ += in->nx[kk] * in->nu[kk];
    }
    return nnzQ;
}

static int_t get_nnzA(const ocp_qp_in *in, const ocp_qp_ooqp_args *args) {
    int_t kk;
    int_t nnzA = 0;

    // dummy command
    if (args->printLevel) kk = 0;

    for (kk = 0; kk < in->N; kk++) {
        nnzA += in->nx[kk + 1] * (in->nx[kk] + in->nu[kk] + 1);
    }
    return nnzA;
}

static int_t get_nnzC(const ocp_qp_in *in, const ocp_qp_ooqp_args *args) {
    int_t kk;
    int_t nnzC = 0;

    // dummy command
    if (args->printLevel) kk = 0;

    for (kk = 0; kk <= in->N; kk++) {
        nnzC += in->nc[kk] * (in->nx[kk] + in->nu[kk]);
    }
    return nnzC;
}

static void update_gradient(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem) {
    int_t ii, kk, nn;

    nn = 0;
    for (kk = 0; kk <= in->N; kk++) {
        for (ii = 0; ii < in->nx[kk]; ii++) mem->c[nn++] = in->q[kk][ii];
        for (ii = 0; ii < in->nu[kk]; ii++) mem->c[nn++] = in->r[kk][ii];
    }
}

static void update_hessian_structure(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                     ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn, offset;

    // TODO(dimitris): For the moment I assume full matrices Q,R,A,B... (we need
    // to def. sparsities) printf("------------> updating Hessian sparsity\n");
    nn = 0;
    offset = 0;
    for (kk = 0; kk <= in->N; kk++) {
        // writing Q[kk]
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = jj; ii < in->nx[kk]; ii++) {  // we write only the lower triangular part
                mem->irowQ[nn] = offset + ii;
                mem->jcolQ[nn] = offset + jj;
                nn += 1;
            }
        }
        // writing S[kk]
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nu[kk]; ii++) {
                mem->irowQ[nn] = offset + in->nx[kk] + ii;
                mem->jcolQ[nn] = offset + jj;
                nn += 1;
            }
        }
        // writing R[kk]
        for (jj = 0; jj < in->nu[kk]; jj++) {
            for (ii = jj; ii < in->nu[kk]; ii++) {
                mem->irowQ[nn] = offset + in->nx[kk] + ii;
                mem->jcolQ[nn] = offset + in->nx[kk] + jj;
                nn += 1;
            }
        }
        offset += in->nx[kk] + in->nu[kk];
    }
    rows = mem->irowQ;
    cols = mem->jcolQ;
    lda = mem->nx;
    qsort(mem->orderQ, mem->nnzQ, sizeof(*mem->orderQ), comparator);
    sort_matrix_structure_row_major(mem->orderQ, mem->irowQ, mem->nnzQ, mem->jcolQ, work->tmpInt);
}

static void update_hessian_data(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn, offset;

    // printf("------------> updating Hessian data\n");
    nn = 0;
    offset = 0;
    for (kk = 0; kk <= in->N; kk++) {
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = jj; ii < in->nx[kk]; ii++) {  // we write only the lower triangular part
                mem->dQ[nn++] = in->Q[kk][jj * in->nx[kk] + ii];
            }
        }
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nu[kk]; ii++) {
                mem->dQ[nn++] = in->S[kk][jj * in->nu[kk] + ii];
            }
        }
        for (jj = 0; jj < in->nu[kk]; jj++) {
            for (ii = jj; ii < in->nu[kk]; ii++) {
                mem->dQ[nn++] = in->R[kk][jj * in->nu[kk] + ii];
            }
        }
        offset += in->nx[kk] + in->nu[kk];
    }
    sort_matrix_data_row_major(mem->orderQ, mem->nnzQ, mem->dQ, work->tmpReal);
}

static void update_b_vector(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem) {
    int_t ii, kk;
    int_t nn = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk + 1]; ii++) mem->bA[nn++] = -in->b[kk][ii];
    }
}

static void update_dynamics_structure(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                      ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn, offsetRows, offsetCols;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk < in->N; kk++) {
        // writing A[kk] (nx[k+1] x nx[k])
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk + 1]; ii++) {
                // printf("writing A_%d[%d,%d]\n", kk, ii, jj);
                mem->irowA[nn] = offsetRows + ii;
                mem->jcolA[nn] = offsetCols + jj;
                nn += 1;
            }
        }
        // writing B[kk] (nx[k+1] x nu[k])
        for (jj = 0; jj < in->nu[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk + 1]; ii++) {
                mem->irowA[nn] = offsetRows + ii;
                mem->jcolA[nn] = offsetCols + in->nx[kk] + jj;
                nn += 1;
            }
        }
        // writing -I (nx[k+1] x nx[k+1])
        for (jj = 0; jj < in->nx[kk + 1]; jj++) {
            mem->irowA[nn] = offsetRows + jj;
            mem->jcolA[nn] = offsetCols + in->nx[kk] + in->nu[kk] + jj;
            nn += 1;
        }
        offsetCols += in->nx[kk] + in->nu[kk];
        offsetRows += in->nx[kk + 1];
    }
    rows = mem->irowA;
    cols = mem->jcolA;
    lda = mem->nx;
    qsort(mem->orderA, mem->nnzA, sizeof(*mem->orderA), comparator);
    sort_matrix_structure_row_major(mem->orderA, mem->irowA, mem->nnzA, mem->jcolA, work->tmpInt);
}

static void update_dynamics_data(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                 ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn, offsetRows, offsetCols;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk + 1]; ii++) {
                mem->dA[nn++] = in->A[kk][jj * in->nx[kk + 1] + ii];
            }
        }
        for (jj = 0; jj < in->nu[kk]; jj++) {
            for (ii = 0; ii < in->nx[kk + 1]; ii++) {
                mem->dA[nn++] = in->B[kk][jj * in->nx[kk + 1] + ii];
            }
        }
        for (jj = 0; jj < in->nx[kk + 1]; jj++) {
            mem->dA[nn++] = -1;
        }
        offsetCols += in->nx[kk] + in->nu[kk];
        offsetRows += in->nx[kk + 1];
    }
    sort_matrix_data_row_major(mem->orderA, mem->nnzA, mem->dA, work->tmpReal);
}

static void update_bounds(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem) {
    int_t ii, kk;
    int_t offset = 0;
    int_t idx;

    for (kk = 0; kk <= in->N; kk++) {
        for (ii = 0; ii < in->nx[kk] + in->nu[kk]; ii++) {
            mem->ixlow[offset + ii] = (char)0;
            mem->ixupp[offset + ii] = (char)0;
            mem->xlow[offset + ii] = 0.0;
            mem->xupp[offset + ii] = 0.0;
        }
        for (ii = 0; ii < in->nb[kk]; ii++) {
#ifdef FLIP_BOUNDS
            if (in->idxb[kk][ii] < in->nu[kk]) {
                idx = in->idxb[kk][ii] + in->nx[kk];
            } else {
                idx = in->idxb[kk][ii] - in->nu[kk];
            }
// printf("OOQP with flipped bounds\n"); exit(1);
#else
            idx = in->idxb[kk][ii];
// printf("OOQP with normal bounds\n"); exit(1);
#endif
            // TODO(dimitris): check if cast is redundant
            // NOTE(dimitris): OOQP can give wrong results if there are 1e12 bounds
            if (in->lb[kk][ii] > -1e10) {  // TODO(dimitris): use acados inf
                mem->ixlow[offset + idx] = (char)1;
                mem->xlow[offset + idx] = in->lb[kk][ii];
            }
            if (in->ub[kk][ii] < 1e10) {  // TODO(dimitris): same here
                mem->ixupp[offset + idx] = (char)1;
                mem->xupp[offset + idx] = in->ub[kk][ii];
            }
        }
        offset += in->nx[kk] + in->nu[kk];
    }
}

static void update_ineq_bounds(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem) {
    int_t ii, kk;
    int_t nn = 0;

    for (kk = 0; kk < in->N + 1; kk++) {
        for (ii = 0; ii < in->nc[kk]; ii++) {
            mem->iclow[nn] = (char)1;
            mem->clow[nn] = in->lc[kk][ii];
            mem->icupp[nn] = (char)1;
            mem->cupp[nn] = in->uc[kk][ii];
            nn += 1;
        }
    }
}

static void update_inequalities_structure(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                          ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn, offsetRows, offsetCols;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk <= in->N; kk++) {
        // writing Cx[k] (nc[k] x nx[k])
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nc[kk]; ii++) {
                // printf("writing C_%d[%d,%d]\n", kk, ii, jj);
                mem->irowC[nn] = offsetRows + ii;
                mem->jcolC[nn] = offsetCols + jj;
                nn += 1;
            }
        }
        // writing Cu[k] (nc[k] x nu[k])
        for (jj = 0; jj < in->nu[kk]; jj++) {
            for (ii = 0; ii < in->nc[kk]; ii++) {
                mem->irowC[nn] = offsetRows + ii;
                mem->jcolC[nn] = offsetCols + in->nx[kk] + jj;
                nn += 1;
            }
        }
        offsetCols += in->nx[kk] + in->nu[kk];
        offsetRows += in->nc[kk];
    }
    rows = mem->irowC;
    cols = mem->jcolC;
    lda = mem->nx;
    qsort(mem->orderC, mem->nnzC, sizeof(*mem->orderC), comparator);
    sort_matrix_structure_row_major(mem->orderC, mem->irowC, mem->nnzC, mem->jcolC, work->tmpInt);
}

static void update_inequalities_data(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                     ocp_qp_ooqp_workspace *work) {
    int_t ii, jj, kk, nn, offsetRows, offsetCols;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk <= in->N; kk++) {
        for (jj = 0; jj < in->nx[kk]; jj++) {
            for (ii = 0; ii < in->nc[kk]; ii++) {
                mem->dC[nn++] = in->Cx[kk][jj * in->nc[kk] + ii];
            }
        }
        for (jj = 0; jj < in->nu[kk]; jj++) {
            for (ii = 0; ii < in->nc[kk]; ii++) {
                mem->dC[nn++] = in->Cu[kk][jj * in->nc[kk] + ii];
            }
        }
        offsetCols += in->nx[kk] + in->nu[kk];
        offsetRows += in->nc[kk];
    }
    sort_matrix_data_row_major(mem->orderC, mem->nnzC, mem->dC, work->tmpReal);
}

static void ocp_qp_ooqp_update_memory(const ocp_qp_in *in, const ocp_qp_ooqp_args *args,
                                      ocp_qp_ooqp_memory *mem, ocp_qp_ooqp_workspace *work) {
    int_t ii;

    if (mem->firstRun == 1) {
        for (ii = 0; ii < mem->nnzQ; ii++) mem->orderQ[ii] = ii;
        for (ii = 0; ii < mem->nnzA; ii++) mem->orderA[ii] = ii;
        for (ii = 0; ii < mem->nnzC; ii++) mem->orderC[ii] = ii;
    }

    // ------- Update objective
    update_gradient(in, mem);

    if (mem->firstRun == 1 || (args->fixHessianSparsity == 0 && args->fixHessian == 0)) {
        update_hessian_structure(in, mem, work);
    }
    if (mem->firstRun == 1 || args->fixHessian == 0) {
        update_hessian_data(in, mem, work);
    }

    // ------- Update equality constraints
    update_b_vector(in, mem);

    if (mem->firstRun == 1 || (args->fixDynamicsSparsity == 0 && args->fixDynamics == 0)) {
        update_dynamics_structure(in, mem, work);
    }
    if (mem->firstRun == 1 || args->fixDynamics == 0) {
        update_dynamics_data(in, mem, work);
    }

    // ------- Update bounds
    update_bounds(in, mem);

    // ------- Update inequality constraints
    update_ineq_bounds(in, mem);

    if (mem->firstRun == 1 || (args->fixInequalitiesSparsity == 0 && args->fixInequalities == 0)) {
        update_inequalities_structure(in, mem, work);
    }
    if (mem->firstRun == 1 || args->fixInequalities == 0) {
        update_inequalities_data(in, mem, work);
    }

    mem->firstRun = 0;
}

static void print_inputs(ocp_qp_ooqp_memory *mem) {
    printf("\n----------> OOQP INPUTS <----------\n\n");
    printf("NUMBER OF PRIMAL VARIABLES: %d\n", mem->nx);
    printf("NUMBER OF NON-ZEROS in HESSIAN: %d\n", mem->nnzQ);
    printf("NUMBER OF EQUALITY CONSTRAINTS: %d\n", mem->my);
    printf("NUMBER OF NON-ZEROS in EQUALITIES: %d\n", mem->nnzA);
    printf("NUMBER OF INEQUALITY CONSTRAINTS: %d\n", mem->mz);
    printf("NUMBER OF NON-ZEROS in INEQUALITIES: %d\n", mem->nnzC);
    printf("\n-----------------------------------\n\n");

    int_t ii;
    printf("\nOBJECTIVE FUNCTION:\n");
    for (ii = 0; ii < mem->nnzQ; ii++) {
        printf("=====> Q[%d, %d] = %f\n", mem->irowQ[ii] + 1, mem->jcolQ[ii] + 1, mem->dQ[ii]);
    }
    for (ii = 0; ii < mem->nx; ii++) printf("===> c[%d] = %f\n", ii + 1, mem->c[ii]);

    printf("\nEQUALITY CONSTRAINTS:\n");
    for (ii = 0; ii < mem->nnzA; ii++) {
        printf("=====> A[%d, %d] = %f\n", mem->irowA[ii] + 1, mem->jcolA[ii] + 1, mem->dA[ii]);
    }
    for (ii = 0; ii < mem->my; ii++) printf("===> bA[%d] = %f\n", ii + 1, mem->bA[ii]);

    printf("\nBOUNDS:\n");
    for (ii = 0; ii < mem->nx; ii++) {
        printf(
            "ixlow[%d] = %d \t xlow[%d] = %4.2f \t ixupp[%d] = %d \t xupp[%d] "
            "= %4.2f\n",
            ii + 1, mem->ixlow[ii], ii + 1, mem->xlow[ii], ii + 1, mem->ixupp[ii], ii + 1,
            mem->xupp[ii]);
    }

    printf("\nINEQUALITY CONSTRAINTS:\n");
    for (ii = 0; ii < mem->nnzC; ii++) {
        printf("=====> C[%d, %d] = %f\n", mem->irowC[ii] + 1, mem->jcolC[ii] + 1, mem->dC[ii]);
    }
    for (ii = 0; ii < mem->mz; ii++) {
        printf("===> clow[%d] = %4.2f \t cupp[%d] = %4.2f\n", ii + 1, mem->clow[ii], ii + 1,
               mem->cupp[ii]);
    }
}

static void print_outputs(ocp_qp_ooqp_memory *mem, ocp_qp_ooqp_workspace *work, int return_value) {
    int_t ii;

    printf("\n----------> OOQP OUTPUTS <---------\n\n");
    printf("RETURN STATUS: %d\n", return_value);
    printf("OBJECTIVE VALUE: %f\n", work->objectiveValue);
    printf("FIRST AND LAST ELEMENT OF SOLUTION:\n");
    printf("x[0] = %f\n", work->x[0]);
    printf("x[%d] = %f\n", mem->nx, work->x[mem->nx - 1]);
    printf("\n----------------------------------\n\n");

    printf("\nPRIMAL SOLUTION:\n");
    for (ii = 0; ii < mem->nx; ii++) {
        printf("=====> x[%d] = %f\n", ii + 1, work->x[ii]);
    }
}

static void fill_in_qp_out(const ocp_qp_in *in, ocp_qp_out *out, ocp_qp_ooqp_workspace *work) {
    int_t kk, ii, nn;

    nn = 0;
    for (kk = 0; kk <= in->N; kk++) {
        for (ii = 0; ii < in->nx[kk]; ii++) out->x[kk][ii] = work->x[nn++];
        for (ii = 0; ii < in->nu[kk]; ii++) out->u[kk][ii] = work->x[nn++];
    }
    nn = 0;
    for (kk = 0; kk < in->N; kk++) {
        for (ii = 0; ii < in->nx[kk + 1]; ii++) out->pi[kk][ii] = -work->y[nn++];
    }
    // TODO(dimitris): fill-in multipliers of inequalities
}

ocp_qp_ooqp_args *ocp_qp_ooqp_create_arguments() {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args *)malloc(sizeof(ocp_qp_ooqp_args));

    args->printLevel = 0;
    args->fixHessianSparsity = 1;
    args->fixDynamicsSparsity = 1;
    args->fixInequalitiesSparsity = 1;
    args->fixHessian = 0;
    args->fixDynamics = 0;
    args->fixInequalities = 0;

    return args;
}

static void ocp_qp_ooqp_cast_workspace(ocp_qp_ooqp_workspace *work, ocp_qp_ooqp_memory *mem) {
    char *ptr = (char *)work;

    ptr += sizeof(ocp_qp_ooqp_workspace);
    work->x = (real_t *)ptr;
    ptr += (mem->nx) * sizeof(real_t);
    work->gamma = (real_t *)ptr;
    ptr += (mem->nx) * sizeof(real_t);
    work->phi = (real_t *)ptr;
    ptr += (mem->nx) * sizeof(real_t);
    work->y = (real_t *)ptr;
    ptr += (mem->my) * sizeof(real_t);
    work->z = (real_t *)ptr;
    ptr += (mem->mz) * sizeof(real_t);
    work->lambda = (real_t *)ptr;
    ptr += (mem->mz) * sizeof(real_t);
    work->pi = (real_t *)ptr;
    ptr += (mem->mz) * sizeof(real_t);
    work->tmpInt = (int_t *)ptr;
    ptr += (mem->nnz) * sizeof(int_t);
    work->tmpReal = (real_t *)ptr;
    // ptr += (mem->nnz)*sizeof(real_t);
}

ocp_qp_ooqp_memory *ocp_qp_ooqp_create_memory(const ocp_qp_in *in, void *args_) {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args *)args_;
    ocp_qp_ooqp_memory *mem = (ocp_qp_ooqp_memory *)malloc(sizeof(ocp_qp_ooqp_memory));

    mem->firstRun = 1;
    mem->nx = get_number_of_primal_vars(in);
    mem->my = get_number_of_equalities(in);
    mem->mz = get_number_of_inequalities(in);
    mem->nnzQ = get_nnzQ(in, args);
    mem->nnzA = get_nnzA(in, args);
    mem->nnzC = get_nnzC(in, args);
    mem->nnz = max_of_three(mem->nnzQ, mem->nnzA, mem->nnzC);

    int ooqp_failed;
    newQpGenSparse(&mem->c, mem->nx, &mem->irowQ, mem->nnzQ, &mem->jcolQ, &mem->dQ, &mem->xlow,
                   &mem->ixlow, &mem->xupp, &mem->ixupp, &mem->irowA, mem->nnzA, &mem->jcolA,
                   &mem->dA, &mem->bA, mem->my, &mem->irowC, mem->nnzC, &mem->jcolC, &mem->dC,
                   &mem->clow, mem->mz, &mem->iclow, &mem->cupp, &mem->icupp, &ooqp_failed);

    mem->orderQ = (int_t *)calloc(mem->nnzQ, sizeof(*mem->orderQ));
    mem->orderA = (int_t *)calloc(mem->nnzA, sizeof(*mem->orderA));
    mem->orderC = (int_t *)calloc(mem->nnzC, sizeof(*mem->orderC));

    if (ooqp_failed) return NULL;
    return mem;
}

int_t ocp_qp_ooqp_calculate_workspace_size(const ocp_qp_in *in, void *args_) {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args *)args_;

    int_t size = 0;
    int_t nx, my, mz, nnzQ, nnzA, nnzC, nnz;

    nx = get_number_of_primal_vars(in);
    my = get_number_of_equalities(in);
    mz = get_number_of_inequalities(in);
    nnzQ = get_nnzQ(in, args);
    nnzA = get_nnzA(in, args);
    nnzC = get_nnzC(in, args);
    nnz = max_of_three(nnzQ, nnzA, nnzC);

    size += sizeof(ocp_qp_ooqp_workspace);
    size += sizeof(real_t) * (3 * nx + my + 3 * mz);
    size += sizeof(int_t) * nnz;
    size += sizeof(real_t) * nnz;

    return size;
}

void ocp_qp_ooqp_free_memory(void *mem_) {
    ocp_qp_ooqp_memory *mem = (ocp_qp_ooqp_memory *)mem_;

    freeQpGenSparse(&mem->c, &mem->irowQ, &mem->jcolQ, &mem->dQ, &mem->xlow, &mem->ixlow,
                    &mem->xupp, &mem->ixupp, &mem->irowA, &mem->jcolA, &mem->dA, &mem->bA,
                    &mem->irowC, &mem->jcolC, &mem->dC, &mem->clow, &mem->iclow, &mem->cupp,
                    &mem->icupp);

    free(mem->orderQ);
    free(mem->orderA);
    free(mem->orderC);
}

int_t ocp_qp_ooqp(const ocp_qp_in *in, ocp_qp_out *out, void *args_, void *memory_, void *work_) {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args *)args_;
    ocp_qp_ooqp_memory *mem = (ocp_qp_ooqp_memory *)memory_;
    ocp_qp_ooqp_workspace *work = (ocp_qp_ooqp_workspace *)work_;

    if (0) print_inputs(mem);

    // NOTE: has to be called after setting up the memory which contains the problem dimensions
    ocp_qp_ooqp_cast_workspace(work, mem);
    ocp_qp_ooqp_update_memory(in, args, mem, work);

    // TODO(dimitris): implement dense OOQP
    // call sparse OOQP
    int ooqp_status;
    qpsolvesp(mem->c, mem->nx, mem->irowQ, mem->nnzQ, mem->jcolQ, mem->dQ, mem->xlow, mem->ixlow,
              mem->xupp, mem->ixupp, mem->irowA, mem->nnzA, mem->jcolA, mem->dA, mem->bA, mem->my,
              mem->irowC, mem->nnzC, mem->jcolC, mem->dC, mem->clow, mem->mz, mem->iclow, mem->cupp,
              mem->icupp, work->x, work->gamma, work->phi, work->y, work->z, work->lambda, work->pi,
              &work->objectiveValue, args->printLevel, &ooqp_status);

    if (0) print_outputs(mem, work, return_value);
    fill_in_qp_out(in, out, work);

    int acados_status = ooqp_status;
    if (ooqp_status == SUCCESSFUL_TERMINATION) acados_status = ACADOS_SUCCESS;
    if (ooqp_status == MAX_ITS_EXCEEDED) acados_status = ACADOS_MAXITER;
    return acados_status;
}

void ocp_qp_ooqp_initialize(const ocp_qp_in *qp_in, void *args_, void **mem, void **work) {
    ocp_qp_ooqp_args *args = (ocp_qp_ooqp_args *)args_;

    *mem = ocp_qp_ooqp_create_memory(qp_in, args);
    int_t work_space_size = ocp_qp_ooqp_calculate_workspace_size(qp_in, args);
    *work = (void *)malloc(work_space_size);
}

void ocp_qp_ooqp_destroy(void *mem_, void *work) {
    free(work);
    ocp_qp_ooqp_free_memory(mem_);
}
