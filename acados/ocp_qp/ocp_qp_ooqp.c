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
#include <assert.h>

#include "ooqp/cQpGenSparse.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// acados
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "acados/utils/timing.h"
#include "acados/utils/print.h"
#include "acados/utils/mem.h"

#define FLIP_BOUNDS

int *rows;
int *cols;
int lda;

static int max_of_three(int a, int b, int c)
{
    int ans = a;

    (void)((ans < b) && (ans = b));
    (void)((ans < c) && (ans = c));

    return ans;
}

// comparator for qsort
static int comparator(const void *p1, const void *p2)
{
    int ans1, ans2;
    int ind1 = *((int *)p1);
    int ind2 = *((int *)p2);

    ans1 = rows[ind1] * lda + cols[ind1];
    ans2 = rows[ind2] * lda + cols[ind2];

    return ans1 - ans2;
}

static void sort_matrix_structure_row_major(int *order, int *irow, int nnz, int *jcol,
                                            int *tmp)
{
    int ii;

    for (ii = 0; ii < nnz; ii++)
    {
        tmp[ii] = irow[order[ii]];
    }

    for (ii = 0; ii < nnz; ii++)
    {
        irow[ii] = tmp[ii];
    }

    for (ii = 0; ii < nnz; ii++)
    {
        tmp[ii] = jcol[order[ii]];
    }

    for (ii = 0; ii < nnz; ii++)
    {
        jcol[ii] = tmp[ii];
    }
}

static void sort_matrix_data_row_major(int *order, int nnz, real_t *d, real_t *tmp)
{
    int ii;

    for (ii = 0; ii < nnz; ii++)
    {
        tmp[ii] = d[order[ii]];
    }

    for (ii = 0; ii < nnz; ii++)
    {
        d[ii] = tmp[ii];
    }
}

static int get_number_of_primal_vars(const ocp_qp_dims *dims)
{
    int nx = 0;
    int kk;

    for (kk = 0; kk <= dims->N; kk++)
    {
        nx += dims->nx[kk] + dims->nu[kk];
    }

    return nx;
}

static int get_number_of_equalities(const ocp_qp_dims *dims)
{
    int my = 0;
    int kk;

    for (kk = 0; kk < dims->N; kk++)
    {
        my += dims->nx[kk + 1];
    }

    return my;
}

static int get_number_of_inequalities(const ocp_qp_dims *dims)
{
    int mz = 0;
    int kk;

    for (kk = 0; kk < dims->N + 1; kk++)
    {
        mz += dims->ng[kk];
    }

    return mz;
}

static int get_nnzQ(const ocp_qp_dims *dims)
{
    int kk;
    int nnzQ = 0;

    for (kk = 0; kk <= dims->N; kk++)
    {
        nnzQ += (dims->nx[kk] * dims->nx[kk] - dims->nx[kk]) / 2 + dims->nx[kk];
        nnzQ += (dims->nu[kk] * dims->nu[kk] - dims->nu[kk]) / 2 + dims->nu[kk];
        nnzQ += dims->nx[kk] * dims->nu[kk];
    }

    return nnzQ;
}

static int get_nnzA(const ocp_qp_dims *dims)
{
    int kk;
    int nnzA = 0;

    for (kk = 0; kk < dims->N; kk++)
    {
        nnzA += dims->nx[kk + 1] * (dims->nx[kk] + dims->nu[kk] + 1);
    }

    return nnzA;
}

static int get_nnzC(const ocp_qp_dims *dims)
{
    int kk;
    int nnzC = 0;

    for (kk = 0; kk <= dims->N; kk++)
    {
        nnzC += dims->ng[kk] * (dims->nx[kk] + dims->nu[kk]);
    }

    return nnzC;
}

static void update_gradient(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem)
{
    int kk, nn;
    ocp_qp_dims *dims = in->dim;

    nn = 0;
    for (kk = 0; kk <= dims->N; kk++)
    {
        blasfeo_unpack_dvec(dims->nx[kk], in->rqz+kk, dims->nu[kk], &mem->c[nn]);
        nn += dims->nx[kk];
        blasfeo_unpack_dvec(dims->nu[kk], in->rqz+kk, 0, &mem->c[nn]);
        nn += dims->nu[kk];
    }
}

static void update_hessian_structure(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                     ocp_qp_ooqp_workspace *work)
{
    int ii, jj, kk, nn, offset;
    ocp_qp_dims *dims = in->dim;

    // TODO(dimitris): For the moment I assume full matrices Q,R,A,B... (we need
    // to def. sparsities) printf("------------> updating Hessian sparsity\n");
    nn = 0;
    offset = 0;
    for (kk = 0; kk <= dims->N; kk++)
    {
        // writing Q[kk]
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = jj; ii < dims->nx[kk]; ii++)
            {  // we write only the lower triangular part
                mem->irowQ[nn] = offset + ii;
                mem->jcolQ[nn] = offset + jj;
                nn++;
            }
        }

        // writing S[kk]
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = 0; ii < dims->nu[kk]; ii++)
            {
                mem->irowQ[nn] = offset + dims->nx[kk] + ii;
                mem->jcolQ[nn] = offset + jj;
                nn++;
            }
        }

        // writing R[kk]
        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            for (ii = jj; ii < dims->nu[kk]; ii++)
            {
                mem->irowQ[nn] = offset + dims->nx[kk] + ii;
                mem->jcolQ[nn] = offset + dims->nx[kk] + jj;
                nn++;
            }
        }

        offset += dims->nx[kk] + dims->nu[kk];
    }

    rows = mem->irowQ;
    cols = mem->jcolQ;
    lda = mem->nx;
    qsort(mem->orderQ, mem->nnzQ, sizeof(*mem->orderQ), comparator);
    sort_matrix_structure_row_major(mem->orderQ, mem->irowQ, mem->nnzQ, mem->jcolQ, work->tmpInt);
}

static void update_hessian_data(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                ocp_qp_ooqp_workspace *work)
{
    int ii, jj, kk, nn, offset;
    ocp_qp_dims *dims = in->dim;

    // printf("------------> updating Hessian data\n");
    nn = 0;
    offset = 0;
    for (kk = 0; kk <= dims->N; kk++)
    {
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = jj; ii < dims->nx[kk]; ii++)
            {  // we write only the lower triangular part
                mem->dQ[nn++] = BLASFEO_DMATEL(&in->RSQrq[kk], ii+dims->nu[kk], jj+dims->nu[kk]);
            }
        }
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = 0; ii < dims->nu[kk]; ii++)
            {
                mem->dQ[nn++] = BLASFEO_DMATEL(&in->RSQrq[kk], jj+dims->nu[kk], ii);
            }
        }
        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            for (ii = jj; ii < dims->nu[kk]; ii++)
            {
                mem->dQ[nn++] = BLASFEO_DMATEL(&in->RSQrq[kk], ii, jj);
            }
        }
        offset += dims->nx[kk] + dims->nu[kk];
    }
    sort_matrix_data_row_major(mem->orderQ, mem->nnzQ, mem->dQ, work->tmpReal);
}

static void update_b_vector(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem)
{
    int ii, kk;
    int nn = 0;
    ocp_qp_dims *dims = in->dim;

    for (kk = 0; kk < dims->N; kk++)
    {
        for (ii = 0; ii < dims->nx[kk + 1]; ii++)
        {
            mem->bA[nn++] = -BLASFEO_DVECEL(&in->b[kk], ii);
        }
    }
}

static void update_dynamics_structure(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                      ocp_qp_ooqp_workspace *work)
{
    int ii, jj, kk, nn, offsetRows, offsetCols;
    ocp_qp_dims *dims = in->dim;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk < dims->N; kk++)
    {
        // writing A[kk] (nx[k+1] x nx[k])
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = 0; ii < dims->nx[kk + 1]; ii++)
            {
                // printf("writing A_%d[%d,%d]\n", kk, ii, jj);
                mem->irowA[nn] = offsetRows + ii;
                mem->jcolA[nn] = offsetCols + jj;
                nn++;
            }
        }
        // writing B[kk] (nx[k+1] x nu[k])
        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            for (ii = 0; ii < dims->nx[kk + 1]; ii++)
            {
                mem->irowA[nn] = offsetRows + ii;
                mem->jcolA[nn] = offsetCols + dims->nx[kk] + jj;
                nn++;
            }
        }
        // writing -I (nx[k+1] x nx[k+1])
        for (jj = 0; jj < dims->nx[kk + 1]; jj++)
        {
            mem->irowA[nn] = offsetRows + jj;
            mem->jcolA[nn] = offsetCols + dims->nx[kk] + dims->nu[kk] + jj;
            nn++;
        }
        offsetCols += dims->nx[kk] + dims->nu[kk];
        offsetRows += dims->nx[kk + 1];
    }
    rows = mem->irowA;
    cols = mem->jcolA;
    lda = mem->nx;
    qsort(mem->orderA, mem->nnzA, sizeof(*mem->orderA), comparator);
    sort_matrix_structure_row_major(mem->orderA, mem->irowA, mem->nnzA, mem->jcolA, work->tmpInt);
}

static void update_dynamics_data(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                 ocp_qp_ooqp_workspace *work)
{
    int ii, jj, kk, nn, offsetRows, offsetCols;
    ocp_qp_dims *dims = in->dim;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk < dims->N; kk++)
    {
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = 0; ii < dims->nx[kk + 1]; ii++)
            {
                mem->dA[nn++] = BLASFEO_DMATEL(&in->BAbt[kk], jj+dims->nu[kk], ii);
            }
        }

        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            for (ii = 0; ii < dims->nx[kk + 1]; ii++)
            {
                mem->dA[nn++] = BLASFEO_DMATEL(&in->BAbt[kk], jj, ii);
            }
        }

        for (jj = 0; jj < dims->nx[kk + 1]; jj++)
        {
            mem->dA[nn++] = -1;
        }
        offsetCols += dims->nx[kk] + dims->nu[kk];
        offsetRows += dims->nx[kk + 1];
    }
    sort_matrix_data_row_major(mem->orderA, mem->nnzA, mem->dA, work->tmpReal);
}

static void update_bounds(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem)
{
    int ii, kk;
    int offset = 0;
    int idx;
    ocp_qp_dims *dims = in->dim;

    for (kk = 0; kk <= dims->N; kk++)
    {
        for (ii = 0; ii < dims->nx[kk] + dims->nu[kk]; ii++)
        {
            mem->ixlow[offset + ii] = (char)0;
            mem->ixupp[offset + ii] = (char)0;
            mem->xlow[offset + ii] = 0.0;
            mem->xupp[offset + ii] = 0.0;
        }

        for (ii = 0; ii < dims->nb[kk]; ii++)
        {
#ifdef FLIP_BOUNDS
            if (in->idxb[kk][ii] < dims->nu[kk])
            {
                idx = in->idxb[kk][ii] + dims->nx[kk];
            }
            else
            {
                idx = in->idxb[kk][ii] - dims->nu[kk];
            }
// printf("OOQP with flipped bounds\n"); exit(1);
#else
            idx = in->idxb[kk][ii];
// printf("OOQP with normal bounds\n"); exit(1);
#endif
            // TODO(dimitris): check if cast is redundant
            // NOTE(dimitris): OOQP can give wrong results if there are 1e12 bounds
            if (BLASFEO_DVECEL(&in->d[kk], ii) > -1e10)
            {  // TODO(dimitris): use acados inf
                mem->ixlow[offset + idx] = (char)1;
                mem->xlow[offset + idx] = BLASFEO_DVECEL(&in->d[kk], ii);
            }
            if (-BLASFEO_DVECEL(&in->d[kk], ii+dims->nb[kk]+dims->ng[kk])  < 1e10)
            {  // TODO(dimitris): same here
                mem->ixupp[offset + idx] = (char)1;
                mem->xupp[offset + idx] = -BLASFEO_DVECEL(&in->d[kk], ii+dims->nb[kk]+dims->ng[kk]);
            }
        }
        offset += dims->nx[kk] + dims->nu[kk];
    }
}

static void update_ineq_bounds(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem)
{
    int ii, kk;
    int nn = 0;
    ocp_qp_dims *dims = in->dim;

    for (kk = 0; kk < dims->N + 1; kk++)
    {
        for (ii = 0; ii < dims->ng[kk]; ii++)
        {
            mem->iclow[nn] = (char)1;
            mem->clow[nn] = BLASFEO_DVECEL(&in->d[kk], ii+dims->nb[kk]);  // lg
            mem->icupp[nn] = (char)1;
            mem->cupp[nn] = -BLASFEO_DVECEL(&in->d[kk], ii+2*dims->nb[kk]+dims->ng[kk]);  // ug
            nn += 1;
        }
    }
}

static void update_inequalities_structure(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                          ocp_qp_ooqp_workspace *work)
{
    int ii, jj, kk, nn, offsetRows, offsetCols;
    ocp_qp_dims *dims = in->dim;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk <= dims->N; kk++)
    {
        // writing Cx[k] (nc[k] x nx[k])
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                // printf("writing C_%d[%d,%d]\n", kk, ii, jj);
                mem->irowC[nn] = offsetRows + ii;
                mem->jcolC[nn] = offsetCols + jj;
                nn++;
            }
        }
        // writing Cu[k] (nc[k] x nu[k])
        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                mem->irowC[nn] = offsetRows + ii;
                mem->jcolC[nn] = offsetCols + dims->nx[kk] + jj;
                nn++;
            }
        }
        offsetCols += dims->nx[kk] + dims->nu[kk];
        offsetRows += dims->ng[kk];
    }
    rows = mem->irowC;
    cols = mem->jcolC;
    lda = mem->nx;
    qsort(mem->orderC, mem->nnzC, sizeof(*mem->orderC), comparator);
    sort_matrix_structure_row_major(mem->orderC, mem->irowC, mem->nnzC, mem->jcolC, work->tmpInt);
}

static void update_inequalities_data(const ocp_qp_in *in, ocp_qp_ooqp_memory *mem,
                                     ocp_qp_ooqp_workspace *work)
{
    int ii, jj, kk, nn, offsetRows, offsetCols;
    ocp_qp_dims *dims = in->dim;

    nn = 0;
    offsetRows = 0;
    offsetCols = 0;
    for (kk = 0; kk <= dims->N; kk++)
    {
        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                mem->dC[nn++] = BLASFEO_DMATEL(&in->DCt[kk], jj+dims->nu[kk], ii);
            }
        }
        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                mem->dC[nn++] = BLASFEO_DMATEL(&in->DCt[kk], jj, ii);
            }
        }
        offsetCols += dims->nx[kk] + dims->nu[kk];
        offsetRows += dims->ng[kk];
    }
    sort_matrix_data_row_major(mem->orderC, mem->nnzC, mem->dC, work->tmpReal);
}

static void ocp_qp_ooqp_update_memory(const ocp_qp_in *in, const ocp_qp_ooqp_opts *opts,
                                      ocp_qp_ooqp_memory *mem, ocp_qp_ooqp_workspace *work)
{
    int ii;

    if (mem->firstRun == 1)
    {
        for (ii = 0; ii < mem->nnzQ; ii++) mem->orderQ[ii] = ii;
        for (ii = 0; ii < mem->nnzA; ii++) mem->orderA[ii] = ii;
        for (ii = 0; ii < mem->nnzC; ii++) mem->orderC[ii] = ii;
    }

    // ------- Update objective
    update_gradient(in, mem);

    if (mem->firstRun == 1 || (opts->fixHessianSparsity == 0 && opts->fixHessian == 0))
    {
        update_hessian_structure(in, mem, work);
    }
    if (mem->firstRun == 1 || opts->fixHessian == 0)
    {
        update_hessian_data(in, mem, work);
    }

    // ------- Update equality constraints
    update_b_vector(in, mem);

    if (mem->firstRun == 1 || (opts->fixDynamicsSparsity == 0 && opts->fixDynamics == 0))
    {
        update_dynamics_structure(in, mem, work);
    }
    if (mem->firstRun == 1 || opts->fixDynamics == 0)
    {
        update_dynamics_data(in, mem, work);
    }

    // ------- Update bounds
    update_bounds(in, mem);

    // ------- Update inequality constraints
    update_ineq_bounds(in, mem);

    if (mem->firstRun == 1 || (opts->fixInequalitiesSparsity == 0 && opts->fixInequalities == 0))
    {
        update_inequalities_structure(in, mem, work);
    }
    if (mem->firstRun == 1 || opts->fixInequalities == 0)
    {
        update_inequalities_data(in, mem, work);
    }

    mem->firstRun = 0;
}

static void print_inputs(ocp_qp_ooqp_memory *mem)
{
    printf("\n----------> OOQP INPUTS <----------\n\n");
    printf("NUMBER OF PRIMAL VARIABLES: %d\n", mem->nx);
    printf("NUMBER OF NON-ZEROS in HESSIAN: %d\n", mem->nnzQ);
    printf("NUMBER OF EQUALITY CONSTRAINTS: %d\n", mem->my);
    printf("NUMBER OF NON-ZEROS in EQUALITIES: %d\n", mem->nnzA);
    printf("NUMBER OF INEQUALITY CONSTRAINTS: %d\n", mem->mz);
    printf("NUMBER OF NON-ZEROS in INEQUALITIES: %d\n", mem->nnzC);
    printf("\n-----------------------------------\n\n");

    int ii;
    printf("\nOBJECTIVE FUNCTION:\n");
    for (ii = 0; ii < mem->nnzQ; ii++)
    {
        printf("=====> Q[%d, %d] = %f\n", mem->irowQ[ii] + 1, mem->jcolQ[ii] + 1, mem->dQ[ii]);
    }
    for (ii = 0; ii < mem->nx; ii++)
        printf("===> c[%d] = %f\n", ii + 1, mem->c[ii]);

    printf("\nEQUALITY CONSTRAINTS:\n");
    for (ii = 0; ii < mem->nnzA; ii++)
    {
        printf("=====> A[%d, %d] = %f\n", mem->irowA[ii] + 1, mem->jcolA[ii] + 1, mem->dA[ii]);
    }
    for (ii = 0; ii < mem->my; ii++)
        printf("===> bA[%d] = %f\n", ii + 1, mem->bA[ii]);

    printf("\nBOUNDS:\n");
    for (ii = 0; ii < mem->nx; ii++)
    {
        printf(
            "ixlow[%d] = %d \t xlow[%d] = %4.2f \t ixupp[%d] = %d \t xupp[%d] "
            "= %4.2f\n",
            ii + 1, mem->ixlow[ii], ii + 1, mem->xlow[ii], ii + 1, mem->ixupp[ii], ii + 1,
            mem->xupp[ii]);
    }

    printf("\nINEQUALITY CONSTRAINTS:\n");
    for (ii = 0; ii < mem->nnzC; ii++)
    {
        printf("=====> C[%d, %d] = %f\n", mem->irowC[ii] + 1, mem->jcolC[ii] + 1, mem->dC[ii]);
    }
    for (ii = 0; ii < mem->mz; ii++)
    {
        printf("===> clow[%d] = %4.2f \t cupp[%d] = %4.2f\n", ii + 1, mem->clow[ii], ii + 1,
               mem->cupp[ii]);
    }
}

static void print_outputs(ocp_qp_ooqp_memory *mem, ocp_qp_ooqp_workspace *work, int return_value)
{
    int ii;

    printf("\n----------> OOQP OUTPUTS <---------\n\n");
    printf("RETURN STATUS: %d\n", return_value);
    printf("OBJECTIVE VALUE: %f\n", work->objectiveValue);
    printf("FIRST AND LAST ELEMENT OF SOLUTION:\n");
    printf("x[0] = %f\n", work->x[0]);
    printf("x[%d] = %f\n", mem->nx, work->x[mem->nx - 1]);
    printf("\n----------------------------------\n\n");

    printf("\nPRIMAL SOLUTION:\n");
    for (ii = 0; ii < mem->nx; ii++)
    {
        printf("=====> x[%d] = %f\n", ii + 1, work->x[ii]);
    }

    // printf("\nGAMMA:\n");
    // for (ii = 0; ii < mem->nx; ii++)
    // {
    //     printf("=====> gamma[%d] = %f\n", ii + 1, work->gamma[ii]);
    // }

    // printf("\nPHI:\n");
    // for (ii = 0; ii < mem->nx; ii++)
    // {
    //     printf("=====> phi[%d] = %f\n", ii + 1, work->phi[ii]);
    // }

    printf("\nLAMBDA:\n");
    for (ii = 0; ii < mem->mz; ii++)
    {
        printf("=====> lambda[%d] = %f\n", ii + 1, work->lambda[ii]);
    }

    printf("\nPI:\n");
    for (ii = 0; ii < mem->mz; ii++)
    {
        printf("=====> pi[%d] = %f\n", ii + 1, work->pi[ii]);
    }
}

static void fill_in_qp_out(const ocp_qp_in *in, ocp_qp_out *out, ocp_qp_ooqp_workspace *work)
{
    int kk, ii, nn, mm;
    int N = in->dim->N;
    int *nb = in->dim->nb;
    int *ng = in->dim->ng;
    int *ns = in->dim->ns;
    int *nx = in->dim->nx;
    int *nu = in->dim->nu;
    int **idxb = in->idxb;

    nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        blasfeo_pack_dvec(nx[kk], &work->x[nn], out->ux+kk, nu[kk]);
        nn += nx[kk];
        blasfeo_pack_dvec(nu[kk], &work->x[nn], out->ux+kk, 0);
        nn += nu[kk];
    }

    nn = 0;
    for (kk = 0; kk < N; kk++)
    {
        blasfeo_pack_dvec(nx[kk + 1], &work->y[nn], out->pi+kk, 0);
        blasfeo_dvecsc(nx[kk + 1], -1.0, out->pi+kk, 0);
        nn += nx[kk + 1];
    }

    nn = 0;
    mm = 0;
    for (kk = 0; kk <= N; kk++)
    {
        for (ii = 0; ii < 2 * nb[kk] + 2 * ng[kk] + 2 * ns[kk]; ii++) out->lam[kk].pa[ii] = 0.0;

        for (ii = 0; ii < nb[kk]; ii++)
        {
            double delta = work->gamma[idxb[kk][ii]+nn] - work->phi[idxb[kk][ii]+nn];
            if (ii < nx[kk])
            {
                if (delta >= 0)
                    out->lam[kk].pa[nu[kk] + ii] = delta;
                else
                    out->lam[kk].pa[nb[kk] + ng[kk] + nu[kk] + ii] = -delta;
            }
            else
            {
                if (delta >= 0)
                    out->lam[kk].pa[ii - nx[kk]] = delta;
                else
                    out->lam[kk].pa[nb[kk] + ng[kk] + ii - nx[kk]] = -delta;
            }
        }

        nn += nx[kk]+nu[kk];

        for (ii = 0; ii < ng[kk]; ii++)
        {
            double delta = work->lambda[ii+mm] - work->pi[ii+mm];
            if (delta >= 0)
                out->lam[kk].pa[nb[kk]+ii] = delta;
            else
                out->lam[kk].pa[2*nb[kk]+ng[kk]+ii] = -delta;
        }

        mm += ng[kk];
    }
}

int ocp_qp_ooqp_opts_calculate_size(void *config_, ocp_qp_dims *dims)
{
    int size = 0;
    size += sizeof(ocp_qp_ooqp_opts);
    return size;
}

void *ocp_qp_ooqp_opts_assign(void *config_, ocp_qp_dims *dims, void *raw_memory)
{
    ocp_qp_ooqp_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (ocp_qp_ooqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_ooqp_opts);

    assert((char *) raw_memory + ocp_qp_ooqp_opts_calculate_size(config_, dims) >= c_ptr);

    return (void *) opts;
}

void ocp_qp_ooqp_opts_initialize_default(void *config_, ocp_qp_dims *dims, void *opts_)
{
    ocp_qp_ooqp_opts *opts = (ocp_qp_ooqp_opts *) opts_;

    opts->printLevel = 0;
    opts->fixHessianSparsity = 1;
    opts->fixDynamicsSparsity = 1;
    opts->fixInequalitiesSparsity = 1;
    opts->fixHessian = 0;
    opts->fixDynamics = 0;
    opts->fixInequalities = 0;

    return;
}

void ocp_qp_ooqp_opts_update(void *config_, ocp_qp_dims *dims, void *opts_)
{
    return;
}

int ocp_qp_ooqp_memory_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    int nx = get_number_of_primal_vars(dims);
    int my = get_number_of_equalities(dims);
    int mz = get_number_of_inequalities(dims);
    int nnzQ = get_nnzQ(dims);
    int nnzA = get_nnzA(dims);
    int nnzC = get_nnzC(dims);
    // int nnz = max_of_three(nnzQ, nnzA, nnzC);

    int size = 0;
    size += sizeof(ocp_qp_ooqp_memory);

    size += 1 *   nx * sizeof(double);  // c
    size += 1 * nnzQ * sizeof(double);  // dQ
    size += 3 * nnzQ * sizeof(int);     // irowQ, jcolQ, orderQ
    size += 2 *   nx * sizeof(double);  // xlow, xupp
    size += 2 *   nx * sizeof(char);    // ixlow, ixupp
    size += 1 * nnzA * sizeof(double);  // dA
    size += 3 * nnzA * sizeof(int);     // irowA, jcolA, orderA
    size += 1 *   my * sizeof(double);  // bA
    size += 1 * nnzC * sizeof(double);  // dC
    size += 3 * nnzC * sizeof(int);     // irowC, jcolC, orderC
    size += 2 *   mz * sizeof(double);  // clow, cupp
    size += 2 *   mz * sizeof(char);    // iclow, icupp

    make_int_multiple_of(8, &size);

    return size;
}

void *ocp_qp_ooqp_memory_assign(void *config_, ocp_qp_dims *dims, void *opts_, void *raw_memory)
{
    // ocp_qp_ooqp_opts *opts = (ocp_qp_ooqp_opts *)opts_;
    ocp_qp_ooqp_memory *mem;

    int nx = get_number_of_primal_vars(dims);
    int my = get_number_of_equalities(dims);
    int mz = get_number_of_inequalities(dims);
    int nnzQ = get_nnzQ(dims);
    int nnzA = get_nnzA(dims);
    int nnzC = get_nnzC(dims);
    // int nnz = max_of_three(nnzQ, nnzA, nnzC);

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (ocp_qp_ooqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_ooqp_memory);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assign_and_advance_double(nx, &mem->c, &c_ptr);
    assign_and_advance_double(nnzQ, &mem->dQ, &c_ptr);
    assign_and_advance_double(nx, &mem->xlow, &c_ptr);
    assign_and_advance_double(nx, &mem->xupp, &c_ptr);
    assign_and_advance_double(nnzA, &mem->dA, &c_ptr);
    assign_and_advance_double(my, &mem->bA, &c_ptr);
    assign_and_advance_double(nnzC, &mem->dC, &c_ptr);
    assign_and_advance_double(mz, &mem->clow, &c_ptr);
    assign_and_advance_double(mz, &mem->cupp, &c_ptr);

    assert((size_t) c_ptr % 8 == 0 && "double not 8-byte aligned!");

    assign_and_advance_int(nnzQ, &mem->irowQ, &c_ptr);
    assign_and_advance_int(nnzQ, &mem->jcolQ, &c_ptr);
    assign_and_advance_int(nnzQ, &mem->orderQ, &c_ptr);
    assign_and_advance_int(nnzA, &mem->irowA, &c_ptr);
    assign_and_advance_int(nnzA, &mem->jcolA, &c_ptr);
    assign_and_advance_int(nnzA, &mem->orderA, &c_ptr);
    assign_and_advance_int(nnzC, &mem->irowC, &c_ptr);
    assign_and_advance_int(nnzC, &mem->jcolC, &c_ptr);
    assign_and_advance_int(nnzC, &mem->orderC, &c_ptr);

    assign_and_advance_char(nx, &mem->ixlow, &c_ptr);
    assign_and_advance_char(nx, &mem->ixupp, &c_ptr);
    assign_and_advance_char(mz, &mem->iclow, &c_ptr);
    assign_and_advance_char(mz, &mem->icupp, &c_ptr);

    // initialize memory
    mem->firstRun = 1;
    mem->nx = nx;
    mem->my = my;
    mem->mz = mz;
    mem->nnzQ = nnzQ;
    mem->nnzA = nnzA;
    mem->nnzC = nnzC;
    mem->nnz = max_of_three(nnzQ, nnzA, nnzC);

    for (int ii = 0; ii < nnzQ; ii++) mem->orderQ[ii] = 0;
    for (int ii = 0; ii < nnzA; ii++) mem->orderA[ii] = 0;
    for (int ii = 0; ii < nnzC; ii++) mem->orderC[ii] = 0;

    assert((char *) raw_memory + ocp_qp_ooqp_memory_calculate_size(config_, dims, opts_) >=
           c_ptr);

    return mem;
}

int ocp_qp_ooqp_workspace_calculate_size(void *config_, ocp_qp_dims *dims, void *opts_)
{
    UNUSED(opts_);

    int size = 0;
    int nx, my, mz, nnzQ, nnzA, nnzC, nnz;

    nx = get_number_of_primal_vars(dims);
    my = get_number_of_equalities(dims);
    mz = get_number_of_inequalities(dims);
    nnzQ = get_nnzQ(dims);
    nnzA = get_nnzA(dims);
    nnzC = get_nnzC(dims);
    nnz = max_of_three(nnzQ, nnzA, nnzC);

    size += sizeof(ocp_qp_ooqp_workspace);
    size += sizeof(double) * (3 * nx + my + 3 * mz);
    size += sizeof(int) * nnz;
    size += sizeof(double) * nnz;

    return size;
}

static void ocp_qp_ooqp_cast_workspace(ocp_qp_ooqp_workspace *work, ocp_qp_ooqp_memory *mem)
{
    char *ptr = (char *)work;

    ptr += sizeof(ocp_qp_ooqp_workspace);
    work->x = (double *)ptr;
    ptr += (mem->nx) * sizeof(double);
    work->gamma = (double *)ptr;
    ptr += (mem->nx) * sizeof(double);
    work->phi = (double *)ptr;
    ptr += (mem->nx) * sizeof(double);
    work->y = (double *)ptr;
    ptr += (mem->my) * sizeof(double);
    work->z = (double *)ptr;
    ptr += (mem->mz) * sizeof(double);
    work->lambda = (double *)ptr;
    ptr += (mem->mz) * sizeof(double);
    work->pi = (double *)ptr;
    ptr += (mem->mz) * sizeof(double);
    work->tmpInt = (int *)ptr;
    ptr += (mem->nnz) * sizeof(int);
    work->tmpReal = (double *)ptr;
    // ptr += (mem->nnz)*sizeof(double);
}

int ocp_qp_ooqp(void *config_, ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *opts_, void *memory_,
                  void *work_)
{
    int N = qp_in->dim->N;
    int *ns = qp_in->dim->ns;

    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
        {
            printf("\nOOQP interface can not handle ns>0 yet: what about implementing it? :)\n");
            return ACADOS_FAILURE;
        }
    }

    acados_timer tot_timer, qp_timer, interface_timer;
    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_tic(&tot_timer);

    ocp_qp_ooqp_opts *opts = (ocp_qp_ooqp_opts *)opts_;
    ocp_qp_ooqp_memory *mem = (ocp_qp_ooqp_memory *)memory_;
    ocp_qp_ooqp_workspace *work = (ocp_qp_ooqp_workspace *)work_;

    acados_tic(&interface_timer);
    // NOTE: has to be called after setting up the memory which contains the problem dimensions
    ocp_qp_ooqp_cast_workspace(work, mem);
    ocp_qp_ooqp_update_memory(qp_in, opts, mem, work);
    info->interface_time = acados_toc(&interface_timer);

    if (0) print_inputs(mem);

    // TODO(dimitris): implement dense OOQP
    // call sparse OOQP
    acados_tic(&qp_timer);
    int ooqp_status;
    qpsolvesp(mem->c, mem->nx, mem->irowQ, mem->nnzQ, mem->jcolQ, mem->dQ, mem->xlow, mem->ixlow,
              mem->xupp, mem->ixupp, mem->irowA, mem->nnzA, mem->jcolA, mem->dA, mem->bA, mem->my,
              mem->irowC, mem->nnzC, mem->jcolC, mem->dC, mem->clow, mem->mz, mem->iclow, mem->cupp,
              mem->icupp, work->x, work->gamma, work->phi, work->y, work->z, work->lambda, work->pi,
              &work->objectiveValue, opts->printLevel, &ooqp_status);
    info->solve_QP_time = acados_toc(&qp_timer);

    if (0) print_outputs(mem, work, ooqp_status);
    acados_tic(&interface_timer);
    fill_in_qp_out(qp_in, qp_out, work);
    ocp_qp_compute_t(qp_in, qp_out);
    info->interface_time += acados_toc(&interface_timer);

    info->total_time = acados_toc(&tot_timer);
    info->num_iter = -1;
    info->t_computed = 1;

    int acados_status = ooqp_status;
    if (ooqp_status == SPARSE_SUCCESSFUL_TERMINATION) acados_status = ACADOS_SUCCESS;
    if (ooqp_status == SPARSE_MAX_ITS_EXCEEDED) acados_status = ACADOS_MAXITER;
    return acados_status;
}

void ocp_qp_ooqp_destroy(void *mem_, void *work)
{
    free(work);
    free(mem_);
}

void ocp_qp_ooqp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->dims_set = &ocp_qp_dims_set;
    config->opts_calculate_size = (int (*)(void *, void *)) & ocp_qp_ooqp_opts_calculate_size;
    config->opts_assign = (void *(*) (void *, void *, void *) ) & ocp_qp_ooqp_opts_assign;
    config->opts_initialize_default =
        (void (*)(void *, void *, void *)) & ocp_qp_ooqp_opts_initialize_default;
    config->opts_update = (void (*)(void *, void *, void *)) & ocp_qp_ooqp_opts_update;
    config->memory_calculate_size =
        (int (*)(void *, void *, void *)) & ocp_qp_ooqp_memory_calculate_size;
    config->memory_assign =
        (void *(*) (void *, void *, void *, void *) ) & ocp_qp_ooqp_memory_assign;
    config->workspace_calculate_size =
        (int (*)(void *, void *, void *)) & ocp_qp_ooqp_workspace_calculate_size;
    config->evaluate = (int (*)(void *, void *, void *, void *, void *, void *)) & ocp_qp_ooqp;
}
