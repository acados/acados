/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include <assert.h>
#include <string.h>

// clarabel
#include "Clarabel.cpp/include/clarabel.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_clarabel.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"


/************************************************
 * helper functions
 ************************************************/

static void print_csc_as_dns(ClarabelCscMatrix *M)
{
    int i, j = 0; // Predefine row index and column index
    int idx;

    // Initialize matrix of zeros
    double *A = (double *) calloc(M->m * M->n, sizeof(double));
    for (int ii=0; ii<M->m*M->n; ii++)
        A[ii] = 1e30;

    // Allocate elements
    for (idx = 0; idx < M->colptr[M->n]; idx++)
    {
        // Get row index i (starting from 1)
        i = M->rowval[idx];

        // Get column index j (increase if necessary) (starting from 1)
        while (M->colptr[j + 1] <= idx) j++;

        // Assign values to A
        A[j * (M->m) + i] = M->nzval[idx];
    }

    for (i = 0; i < M->m; i++)
    {
        for (j = 0; j < M->n; j++)
        {
            if (A[j * (M->m) + i]==1e30)
                printf("  *      ");
            else
                printf("%8.4f ", A[j * (M->m) + i]);
        }
        printf("\n");
    }

    free(A);
}


void print_csc_matrix(ClarabelCscMatrix *M, const char *name)
{
    int j, i, row_start, row_stop;
    int k = 0;

    // Print name
    printf("%s :\n", name);

    for (j = 0; j < M->n; j++) {
        row_start = M->colptr[j];
        row_stop  = M->colptr[j + 1];

        if (row_start == row_stop) continue;
        else {
        for (i = row_start; i < row_stop; i++) {
            printf("\t[%3u,%3u] = %.3g\n", (int)M->rowval[i], (int)j, M->nzval[k++]);
        }
        }
    }
}


static int acados_clarabel_num_vars(ocp_qp_dims *dims)
{
    int n = 0;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        n += dims->nx[ii] + dims->nu[ii] + 2*dims->ns[ii];
    }

    return n;
}



static int acados_clarabel_num_constr(ocp_qp_dims *dims)
{
    int m = 0;

    //printf("\ndims N %d\n", dims->N);

    for (int ii = 0; ii <= dims->N; ii++)
    {
    //printf("dims[%d] nx %d nu %d nb %d ng %d ns %d\n", ii, dims->nx[ii], dims->nu[ii], dims->nb[ii], dims->ng[ii], dims->ns[ii]);

        m += 2 * dims->nb[ii];
        m += 2 * dims->ng[ii];
        m += 2 * dims->ns[ii];

        // equalities
        if (ii < dims->N)
        {
            m += dims->nx[ii+1];
        }
    }

    return m;
}



static int acados_clarabel_nnzmax_P(const ocp_qp_dims *dims)
{
    int nnz = 0;

    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        nnz += (nu[ii]+nx[ii])*(nu[ii]+nx[ii]+1)/2; // triu(RSQ)
        nnz += 2*ns[ii]; // Z
    }

    return nnz;
}



static int acados_clarabel_nnzmax_A(const ocp_qp_dims *dims)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int nnz = 0;

    for (int ii = 0; ii <= N; ii++)
    {
        // inequality constraints
        nnz += 2*nb[ii];         // eye of box constraints
        nnz += 2*ng[ii]*nx[ii];  // C
        nnz += 2*ng[ii]*nu[ii];  // D
        nnz += 2*(nb[ii]+ng[ii])*2*ns[ii]; // soft constraints at worst case, when idxs_rev encoding is used. Typically just 2*ns
        nnz += 2*ns[ii]; // eye of slacks nonnegativity constraints

        // dynamics equality constraints
        if (ii < dims->N)
        {
            nnz += nx[ii+1] * nx[ii];  // A
            nnz += nx[ii+1] * nu[ii];  // B
            nnz += nx[ii+1];           // eye
        }
    }

    return nnz;
}



static void update_gradient(const ocp_qp_in *in, ocp_qp_clarabel_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int kk, nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        blasfeo_unpack_dvec(nu[kk]+nx[kk]+2*ns[kk], in->rqz + kk, 0, &mem->q[nn], 1);
        nn += nu[kk]+nx[kk]+2*ns[kk];
    }

    // actual number of nonzeros
    mem->q_nnz = nn;
}



static void update_hessian_structure(const ocp_qp_in *in, ocp_qp_clarabel_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int ii, jj, kk;

    // CSC format: P_rowval are row indices and P_col_ptr are column pointers
    int nn = 0, offset = 0, col = 0;
    for (kk = 0; kk <= N; kk++)
    {
        // write triu(RSQ[kk])
        for (jj = 0; jj < nx[kk] + nu[kk]; jj++)
        {
            mem->P_col_ptr[col] = nn;
            col++;

            for (ii = 0; ii <= jj; ii++)
            {
                // we write only the upper triangular part
                mem->P_rowval[nn] = offset + ii;
                nn++;
            }
        }
        offset += nx[kk] + nu[kk];

        // write Z[kk]
        for (jj = 0; jj < 2*ns[kk]; jj++)
        {
            mem->P_col_ptr[col] = nn;
            col++;

            // diagonal
            mem->P_rowval[nn] = offset + jj;
            nn++;
        }

        offset += 2*ns[kk];
    }

    mem->P_col_ptr[col] = nn;
    // actual number of nonzeros
    mem->P_nnz = nn;
}



static void update_hessian_data(const ocp_qp_in *in, ocp_qp_clarabel_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int ii, kk;

    // Traversing the matrix in column-major order
    int nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        // writing RSQ[kk]
        // we write the lower triangular part in row-major order
        // that's the same as writing the upper triangular part in
        // column-major order
        for (ii = 0; ii < nx[kk] + nu[kk]; ii++)
        {
            blasfeo_unpack_dmat(1, ii+1, in->RSQrq+kk, ii, 0, mem->P_nzval+nn, 1);
            nn += ii+1;
        }

        // write Z[kk]
        blasfeo_unpack_dvec(2*ns[kk], in->Z+kk, 0, mem->P_nzval+nn, 1);
        nn += 2*ns[kk];
    }

    // TODO ? check that nn==mem->P_nnz
}



static void update_constraints_matrix_structure(const ocp_qp_in *in, ocp_qp_clarabel_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int ii, jj, kk;

    int row_offset_dyn = 0;
    int row_offset_con = 0;
    int row_offset_slk = 0;

    int con_start = 0;
    int slk_start = 0;
    for (kk = 0; kk <= N; kk++)
    {
        con_start += kk < N ? nx[kk+1] : 0;
        slk_start += 2*nb[kk]+2*ng[kk];
    }

    slk_start += con_start;

    // setup cones type
    int m = acados_clarabel_num_constr(dims);
    int m_eq = con_start;
    int m_ineq = m-m_eq;
    mem->cones[0] = ClarabelZeroConeT(m_eq);
    mem->cones[1] = ClarabelNonnegativeConeT(m_ineq);

    // CSC format: A_i are row indices and A_p are column pointers
    int nn = 0, col = 0;
    for (kk = 0; kk <= N; kk++)
    {

        // control variables
        for (jj = 0; jj < nu[kk]; jj++)
        {
            mem->A_col_ptr[col] = nn;
            col++;

            if (kk < dims->N)
            {
                // write column from B
                for (ii = 0; ii < nx[kk+1]; ii++)
                {
                    mem->A_rowval[nn+ii] = row_offset_dyn + ii;
                }
                nn += nx[kk+1];
            }

            // write bound on u (upper and lower)
            for (ii = 0; ii < nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj)
                {
                    mem->A_rowval[nn] = con_start + row_offset_con + ii;
                    nn++;
                    mem->A_rowval[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + ii;
                    nn++;
                    break;
                }
            }

            // write column from D (upper and lower)
            for (ii = 0; ii < ng[kk]; ii++)
            {
                mem->A_rowval[nn+ii] = con_start + row_offset_con + nb[kk] + ii;
                mem->A_rowval[nn+ng[kk]+ii] = con_start + row_offset_con + 2*nb[kk] + ng[kk] + ii;
            }
            nn += 2*ng[kk];
        }

        // state variables
        for (jj = 0; jj < nx[kk]; jj++)
        {
            mem->A_col_ptr[col] = nn;
            col++;

            if (kk > 0)
            {
                // write column from -I
                mem->A_rowval[nn] = row_offset_dyn - nx[kk] + jj;
                nn++;
            }

            if (kk < N)
            {
                // write column from A
                for (ii = 0; ii < nx[kk + 1]; ii++)
                {
                    mem->A_rowval[nn+ii] = row_offset_dyn + ii;
                }
                nn += nx[kk+1];
            }

            // write bound on x (upper and lower)
            for (ii = 0; ii < nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == nu[kk] + jj)
                {
                    mem->A_rowval[nn] = con_start + row_offset_con + ii;
                    nn++;
                    mem->A_rowval[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + ii;
                    nn++;
                    break;
                }
            }

            // write column from C (upper and lower)
            for (ii = 0; ii < ng[kk]; ii++)
            {
                mem->A_rowval[nn+ii] = con_start + row_offset_con + nb[kk] + ii;
                mem->A_rowval[nn+ng[kk]+ii] = con_start + row_offset_con + 2*nb[kk] + ng[kk] + ii;
            }
            nn += 2*ng[kk];
        }

        // slack variables on lower inequalities
        for (jj = 0; jj < ns[kk]; jj++)
        {
            mem->A_col_ptr[col] = nn;
            col++;

            // soft constraint
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    mem->A_rowval[nn] = con_start + row_offset_con + ii;
                    nn++;
                    // no break, there could possibly be multiple
                }
            }

            // nonnegativity constraint
            mem->A_rowval[nn] = slk_start + row_offset_slk + jj;
            nn++;
        }

        // slack variables on upper inequalities
        for (jj = 0; jj < ns[kk]; jj++)
        {
            mem->A_col_ptr[col] = nn;
            col++;

            // soft constraint
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    mem->A_rowval[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + ii;
                    nn++;
                    // no break, there could possibly be multiple
                }
            }

            // nonnegativity constraint
            mem->A_rowval[nn] = slk_start + row_offset_slk + ns[kk] + jj;
            nn++;
        }

        row_offset_con += 2*nb[kk]+2*ng[kk];
        row_offset_dyn += kk < N ? nx[kk + 1] : 0;
        row_offset_slk += 2*ns[kk];
    }

    // end of matrix
    mem->A_col_ptr[col] = nn;
    // actual number of nonzeros
    mem->A_nnz = nn;
}



// TODO move constant stuff like I to structure routine
static void update_constraints_matrix_data(const ocp_qp_in *in, ocp_qp_clarabel_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int ii, jj, kk;


    // Traverse matrix in column-major order
    int nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        // control variables
        for (jj = 0; jj < nu[kk]; jj++)
        {
            if (kk < dims->N)
            {
                // write column from B
                blasfeo_unpack_dmat(1, nx[kk+1], in->BAbt+kk, jj, 0, mem->A_nzval+nn, 1);
                nn += nx[kk+1];
            }

            // write bound on u
            for (ii = 0; ii < dims->nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj)
                {
                    mem->A_nzval[nn] = -1.0; // lower bound
                    nn++;
                    mem->A_nzval[nn] = 1.0; // upper bound
                    nn++;
                    break;
                }
            }

            // write column from D
            blasfeo_unpack_dmat(1, ng[kk], in->DCt+kk, jj, 0, mem->A_nzval+nn+ng[kk], 1);
            for (ii=0; ii<ng[kk]; ii++)
            {
                mem->A_nzval[nn+ii] = - mem->A_nzval[nn+ng[kk]+ii];
            }
            nn += 2*ng[kk];
        }

        // state variables
        for (jj = 0; jj < nx[kk]; jj++)
        {
            if (kk > 0)
            {
                // write column from -I
                mem->A_nzval[nn] = -1.0;
                nn++;
            }

            if (kk < N)
            {
                // write column from A
                blasfeo_unpack_dmat(1, nx[kk+1], in->BAbt+kk, nu[kk]+jj, 0, mem->A_nzval+nn, 1);
                nn += nx[kk+1];
            }

            // write bound on x
            for (ii = 0; ii < nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == nu[kk] + jj)
                {
                    mem->A_nzval[nn] = -1.0; // lower bound
                    nn++;
                    mem->A_nzval[nn] = 1.0; // upper bound
                    nn++;
                    break;
                }
            }

            // write column from C
            blasfeo_unpack_dmat(1, ng[kk], in->DCt+kk, nu[kk]+jj, 0, mem->A_nzval+nn+ng[kk], 1);
            for (ii=0; ii<ng[kk]; ii++)
            {
                mem->A_nzval[nn+ii] = - mem->A_nzval[nn+ng[kk]+ii];
            }
            nn += 2*ng[kk];

        }

        // slack variables on lower inequalities
        for (jj = 0; jj < ns[kk]; jj++)
        {
            // soft constraint
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    //mem->A_nzval[nn] = 1.0;
                    mem->A_nzval[nn] = -1.0;
                    nn++;
                    // no break, there could possibly be multiple
                }
            }

            // nonnegativity constraint
            mem->A_nzval[nn] = -1.0; //1.0;
            nn++;
        }

        // slack variables on upper inequalities
        for (jj = 0; jj < ns[kk]; jj++)
        {
            // soft constraint
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    //mem->A_nzval[nn] = 1.0; //-1.0;
                    mem->A_nzval[nn] = -1.0; //-1.0;
                    nn++;
                    // no break, there could possibly be multiple
                }
            }

            // nonnegativity constraint
            mem->A_nzval[nn] = -1.0; //1.0;
            nn++;
        }

    }

    // TODO ? check that nn==mem->A_nnz
}



static void update_bounds(const ocp_qp_in *in, ocp_qp_clarabel_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    //int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int ii, kk, nn = 0;

    // write b to b
    for (kk = 0; kk < N; kk++)
    {
        // unpack b to b
        blasfeo_unpack_dvec(nx[kk + 1], in->b + kk, 0, &mem->b[nn], 1);
        for (ii = 0; ii < 2*nx[kk+1]; ii++)
        {
            mem->b[nn+ii] = - mem->b[nn+ii];
        }
        nn += nx[kk + 1];
    }

    // write lb lg and ub ug
    for (kk = 0; kk <= N; kk++)
    {
        // unpack lb lg to l and flip sign because in Clarabel lower bounds are casted as upper bounds
        blasfeo_unpack_dvec(nb[kk]+ng[kk], in->d+kk, 0, &mem->b[nn], 1);
        // unpack ub ug to u and flip signs because in HPIPM the signs are flipped for upper bounds
        blasfeo_unpack_dvec(nb[kk]+ng[kk], in->d+kk, nb[kk]+ng[kk], &mem->b[nn+nb[kk]+ng[kk]], 1);
        for (ii = 0; ii < 2*nb[kk]+2*ng[kk]; ii++)
        {
            mem->b[nn+ii] = - mem->b[nn+ii];
        }
        nn += 2*nb[kk] + 2*ng[kk];
    }

    // write ls and us
    for (kk = 0; kk <= N; kk++)
    {
        // unpack ls and us to b because in Clarabel lower bounds are casted as upper bounds
        blasfeo_unpack_dvec(2*ns[kk], in->d+kk, 2*nb[kk]+2*ng[kk], &mem->b[nn], 1);
        for (ii = 0; ii < 2*ns[kk]; ii++)
        {
            mem->b[nn+ii] = - mem->b[nn+ii];
        }
        nn += 2*ns[kk];
    }

    // actual number of nonzeros
    mem->b_nnz = nn;
}



static void clarabel_init_data(ocp_qp_clarabel_memory* mem, ocp_qp_in *qp_in)
{
    //update_constraints_matrix_data(qp_in, mem);
    //update_hessian_data(qp_in, mem);

    //update_bounds(qp_in, mem);
    //update_gradient(qp_in, mem);

    int n = acados_clarabel_num_vars(qp_in->dim);
    int m = acados_clarabel_num_constr(qp_in->dim);

    //printf("\nn %d m %d\n", n, m);

    // allocates and initializes a csc matrix
    clarabel_CscMatrix_init(&mem->A, m, n, mem->A_col_ptr, mem->A_rowval, mem->A_nzval);
    //printf("ocp_qp_clarabel: created A\n");
    //print_csc_matrix(&mem->A, "A_mat");
    //print_csc_as_dns(&mem->A);
    //for (int ii=0; ii<=n; ii++)
    //{
    //    printf("col ptr %d: %d\n", ii, mem->A_col_ptr[ii]);
    //}

    // allocates and initializes a csc matrix
    clarabel_CscMatrix_init(&mem->P, n, n, mem->P_col_ptr, mem->P_rowval, mem->P_nzval);
    //printf("ocp_qp_clarabel: created P\n");
    //print_csc_matrix(&mem->P, "P_mat");
    //print_csc_as_dns(&mem->P);
    //for (int ii=0; ii<=n; ii++)
    //{
    //    printf("col ptr %d: %d\n", ii, mem->P_col_ptr[ii]);
    //}

    //printf("\ndone\n");
    //exit(0);

}



static void ocp_qp_clarabel_update_memory(const ocp_qp_in *in, const ocp_qp_clarabel_opts *opts,
                                      ocp_qp_clarabel_memory *mem)
{
    if (opts->first_run)
    {
        update_hessian_structure(in, mem);
        update_constraints_matrix_structure(in, mem);
    }

    update_hessian_data(in, mem);
    update_constraints_matrix_data(in, mem);
    update_bounds(in, mem);
    update_gradient(in, mem);
}


/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_clarabel_opts_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;
    size += sizeof(ocp_qp_clarabel_opts);
    size += sizeof(ClarabelDefaultSettings);

    return size;
}



void *ocp_qp_clarabel_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_qp_clarabel_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (ocp_qp_clarabel_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_clarabel_opts);

     opts->clarabel_opts = (ClarabelDefaultSettings *) c_ptr;
     c_ptr += sizeof(ClarabelDefaultSettings);

    assert((char *) raw_memory + ocp_qp_clarabel_opts_calculate_size(config_, dims_) == c_ptr);

    return (void *) opts;
}



void ocp_qp_clarabel_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_clarabel_opts *opts = opts_;
    opts->print_level = 0;

    *opts->clarabel_opts = clarabel_DefaultSettings_default();
    opts->clarabel_opts->verbose = false;
    opts->clarabel_opts->presolve_enable = false;

    opts->first_run = 1;

    return;
}



void ocp_qp_clarabel_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_qp_clarabel_opts *opts = (ocp_qp_clarabel_opts *)opts_;

    return;
}

void ocp_qp_clarabel_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_clarabel_opts *opts = opts_;

    // Updating options through this function does not work, only before the first call!
    if (!opts->first_run)
    {
#ifndef ACADOS_SILENT
        printf("\nWARNING: ocp_qp_clarabel_opts_set: attempting to set field: %s. However, options cannot be changed after first run, option is NOT updated. \n", field);
#endif
        return;
    }

    if (!strcmp(field, "iter_max"))
    {
        int *tmp_ptr = value;
        opts->clarabel_opts->max_iter = *tmp_ptr;
    }
    else if (!strcmp(field, "print_level"))
    {
        int* print_level = (int *) value;
        opts->print_level = *print_level;
    }
    else if (!strcmp(field, "tol_stat"))
    {
        double *tol = value;
        // printf("in ocp_qp_clarabel_opts_set, tol_stat %e\n", *tol);
        opts->clarabel_opts->tol_gap_abs = *tol;
    }
    else if (!strcmp(field, "tol_eq"))
    {
        double *tol = value;
        opts->clarabel_opts->tol_infeas_abs = *tol;
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        double *tol = value;
        opts->clarabel_opts->tol_infeas_abs = *tol;
    }
    else if (!strcmp(field, "tol_comp"))
    {
        // do nothing
    }
    else if (!strcmp(field, "warm_start"))
    {
        // do nothing
    }
    else
    {
        printf("\nWARNING: ocp_qp_clarabel_opts_set: field: %s not interfaced yet. Ignoring option and \n", field);
        exit(1);
    }

    return;
}

void ocp_qp_clarabel_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    // ocp_qp_clarabel_opts *opts = opts_;
    printf("\nerror: ocp_qp_clarabel_opts_get: not implemented for field %s\n", field);
    exit(1);
}




/************************************************
 * memory
 ************************************************/

acados_size_t ocp_qp_clarabel_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;

    size_t n = acados_clarabel_num_vars(dims);
    size_t m = acados_clarabel_num_constr(dims);

    size_t A_nnzmax = acados_clarabel_nnzmax_A(dims);
    size_t P_nnzmax = acados_clarabel_nnzmax_P(dims);

    acados_size_t size = 0;
    size += sizeof(ocp_qp_clarabel_memory);

    size += A_nnzmax * sizeof(ClarabelFloat);  // A_nzval
    size += A_nnzmax * sizeof(uintptr_t);      // A_rowval
    size += (n + 1) * sizeof(uintptr_t);       // A_col_ptr

    size += P_nnzmax * sizeof(ClarabelFloat);  // A_nzval
    size += P_nnzmax * sizeof(uintptr_t);      // P_rowval
    size += (n + 1) * sizeof(uintptr_t);       // P_col_ptr

    size += n * sizeof(ClarabelFloat);  // q
    size += m * sizeof(ClarabelFloat);  // b

    size += 1 * 8;

    return size;
}




void *ocp_qp_clarabel_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_qp_dims *dims = dims_;
    ocp_qp_clarabel_memory *mem;

    int n = acados_clarabel_num_vars(dims);
    int m = acados_clarabel_num_constr(dims);
    int P_nnzmax = acados_clarabel_nnzmax_P(dims);
    int A_nnzmax = acados_clarabel_nnzmax_A(dims);

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (ocp_qp_clarabel_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_clarabel_memory);

    mem->P_nnzmax = P_nnzmax;
    mem->A_nnzmax = A_nnzmax;

    mem->solver = NULL;

    align_char_to(8, &c_ptr);

    // doubles
    mem->q = (ClarabelFloat *) c_ptr;
    c_ptr += n * sizeof(ClarabelFloat);

    mem->b = (ClarabelFloat *) c_ptr;
    c_ptr += m * sizeof(ClarabelFloat);

    mem->P_nzval = (ClarabelFloat *) c_ptr;
    c_ptr += (mem->P_nnzmax) * sizeof(ClarabelFloat);

    mem->A_nzval = (ClarabelFloat *) c_ptr;
    c_ptr += (mem->A_nnzmax) * sizeof(ClarabelFloat);

    // ints
    mem->P_rowval = (uintptr_t *) c_ptr;
    c_ptr += (mem->P_nnzmax) * sizeof(uintptr_t);

    mem->P_col_ptr = (uintptr_t *) c_ptr;
    c_ptr += (n + 1) * sizeof(uintptr_t);

    mem->A_rowval = (uintptr_t *) c_ptr;
    c_ptr += (mem->A_nnzmax) * sizeof(uintptr_t);

    mem->A_col_ptr = (uintptr_t *) c_ptr;
    c_ptr += (n + 1) * sizeof(uintptr_t);

    assert((char *) raw_memory + ocp_qp_clarabel_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_qp_clarabel_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    // qp_solver_config *config = config_;
    ocp_qp_clarabel_memory *mem = mem_;

    if (!strcmp(field, "time_qp_solver_call"))
    {
        double *tmp_ptr = value;
        *tmp_ptr = mem->time_qp_solver_call;
    }
    else if (!strcmp(field, "iter"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->iter;
    }
    else if (!strcmp(field, "status"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->status;
    }
    else
    {
        printf("\nerror: ocp_qp_clarabel_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}


void ocp_qp_clarabel_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    // ocp_qp_in *qp_in = qp_in_;
    // reset memory
    printf("acados: reset clarabel_mem not implemented.\n");
    exit(1);
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_clarabel_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    return 0;
}



/************************************************
 * functions
 ************************************************/

static void fill_in_qp_out(const ocp_qp_in *in, ocp_qp_out *out, ocp_qp_clarabel_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int kk, nn;

    //int con_start = 0;
    //int slk_start = 0;
    //for (kk = 0; kk <= N; kk++)
    //{
    //    con_start += kk < N ? nx[kk + 1] : 0;
    //    slk_start += 2*nb[kk] + 2*ng[kk];
    //}

    //slk_start += con_start;

    ClarabelDefaultSolution *sol = &mem->solution;

    // primal variables
    nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        blasfeo_pack_dvec(nx[kk]+nu[kk]+2*ns[kk], &sol->x[nn], 1, out->ux + kk, 0);
        nn += nx[kk] + nu[kk] + 2*ns[kk];
    }

    // dual variables
    nn = 0;
    for (kk = 0; kk < N; kk++)
    {
        blasfeo_pack_dvec(nx[kk + 1], &sol->z[nn], 1, out->pi + kk, 0);
        nn += nx[kk + 1];
    }

    //nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        blasfeo_pack_dvec(2*nb[kk]+2*ng[kk], &sol->z[nn], 1, out->lam+kk, 0);
        nn += 2*nb[kk]+2*ng[kk];
    }

    //nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        blasfeo_pack_dvec(2*ns[kk], &sol->z[nn], 1, out->lam+kk, 2*nb[kk]+2*ng[kk]);
        nn += 2*ns[kk];
    }
}



// clarabel f64 printing stuff
static void print_array_double(double *array, size_t n)
{
    printf("[");
    for (size_t i = 0; i < n; i++)
    {
        printf("%.10f", array[i]);
        if (i < n - 1)
        {
            printf(", ");
        }
    }
    printf("]\n");
}

static void print_solution(ClarabelDefaultSolution_f64 *solution)
{
    printf("Solution (x)\t = ");
    print_array_double(solution->x, solution->x_length);
    printf("Multipliers (z)\t = ");
    print_array_double(solution->z, solution->z_length);
    printf("Slacks (s)\t = ");
    print_array_double(solution->s, solution->s_length);
}

int ocp_qp_clarabel(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;

    //int N = qp_in->dim->N;
    //int *ns = qp_in->dim->ns;

    //print_ocp_qp_dims(qp_in->dim);

    //for (int ii = 0; ii <= N; ii++)
    //{
    //    if (ns[ii] > 0)
    //    {
    //        printf("\nClarabel interface can not handle ns>0 yet.\n");
    //        exit(1);
    //    }
    //}

    // print_ocp_qp_in(qp_in);

    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer, solver_call_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_clarabel_opts *opts = (ocp_qp_clarabel_opts *) opts_;
    ocp_qp_clarabel_memory *mem = (ocp_qp_clarabel_memory *) mem_;

    acados_tic(&interface_timer);
    ocp_qp_clarabel_update_memory(qp_in, opts, mem);
    info->interface_time = acados_toc(&interface_timer);

    acados_tic(&qp_timer);

    if (!opts->first_run)
    {
        // ClarabelDefaultInfo tmp_info;
        clarabel_DefaultSolver_update_P(mem->solver, mem->P_nzval, mem->P_nnz);
        clarabel_DefaultSolver_update_A(mem->solver, mem->A_nzval, mem->A_nnz);
        clarabel_DefaultSolver_update_q(mem->solver, mem->q, mem->q_nnz);
        clarabel_DefaultSolver_update_b(mem->solver, mem->b, mem->b_nnz);
    }
    else
    {
        //printf("\nbefore build solver\n");
        // TODO for now, need to free previous solver (and csc matrices) ...
        if (mem->solver!=NULL)
        {
            //printf("\nfreeing solver before new\n");
            clarabel_DefaultSolver_free(mem->solver);
        }
        // init new csc matrices
        clarabel_init_data(mem, qp_in);
        // Build solver
        mem->solver = clarabel_DefaultSolver_new(&mem->P, mem->q, &mem->A, mem->b, 2, mem->cones, opts->clarabel_opts);
        opts->first_run = 0;
        // opts->first_run = 1;
        // uncomment above to force first_run, and investigate solver updates.
    }

    // solve Clarabel
    acados_tic(&solver_call_timer);
    // Solve
    //clarabel_DefaultSolver_print_to_file(mem->solver, "clarabel_data.txt");
    //clarabel_DefaultSolver_save_to_file(mem->solver, "clarabel_data.txt");
    clarabel_DefaultSolver_solve(mem->solver);
    mem->time_qp_solver_call = acados_toc(&solver_call_timer);


    // Get solution
    mem->solution = clarabel_DefaultSolver_solution(mem->solver);
    if (opts->print_level > 1)
    {
        print_solution(&mem->solution);
    }
    //clarabel_DefaultSolver_free(mem->solver);
    //print_solution(&mem->solution);

    /* fill qp_out */
    fill_in_qp_out(qp_in, qp_out, mem);
    ocp_qp_compute_t(qp_in, qp_out);

    //d_ocp_qp_sol_print(qp_in->dim, qp_out);

    // info
    ClarabelDefaultInfo clarabel_info = clarabel_DefaultSolver_info(mem->solver);
    info->solve_QP_time = acados_toc(&qp_timer);
    //info->solve_QP_time = clarabel_info.solve_time;
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = clarabel_info.iterations;
    mem->iter = clarabel_info.iterations;

    // status
    int acados_status = ACADOS_QP_FAILURE; // generic QP failure
    ClarabelSolverStatus clarabel_status = mem->solution.status;
    if (clarabel_status==ClarabelSolved)
        acados_status = ACADOS_SUCCESS;
    else if (clarabel_status==ClarabelMaxIterations)
        acados_status = ACADOS_MAXITER;

    return acados_status;
}



void ocp_qp_clarabel_eval_adj_sens(void *config_, void *param_qp_in_, void *seed, void *sens_qp_out_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_clarabel_eval_adj_sens: not implemented yet\n");
    exit(1);
}

void ocp_qp_clarabel_eval_forw_sens(void *config_, void *param_qp_in_, void *seed, void *sens_qp_out_, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_clarabel_eval_forw_sens: not implemented yet\n");
    exit(1);
}

void ocp_qp_clarabel_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: ocp_qp_clarabel_solver_get: not implemented yet\n");
    exit(1);
}


void ocp_qp_clarabel_terminate(void *config_, void *mem_, void *work_)
{
    ocp_qp_clarabel_memory *mem = (ocp_qp_clarabel_memory *) mem_;
    // Free the matrices and the solver
    clarabel_DefaultSolver_free(mem->solver);
}




void ocp_qp_clarabel_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = &ocp_qp_clarabel_opts_calculate_size;
    config->opts_assign = &ocp_qp_clarabel_opts_assign;
    config->opts_initialize_default = &ocp_qp_clarabel_opts_initialize_default;
    config->opts_update = &ocp_qp_clarabel_opts_update;
    config->opts_set = &ocp_qp_clarabel_opts_set;
    config->opts_get = &ocp_qp_clarabel_opts_get;
    config->memory_calculate_size = &ocp_qp_clarabel_memory_calculate_size;
    config->memory_assign = &ocp_qp_clarabel_memory_assign;
    config->memory_get = &ocp_qp_clarabel_memory_get;
    config->workspace_calculate_size = &ocp_qp_clarabel_workspace_calculate_size;
    config->evaluate = &ocp_qp_clarabel;
    config->terminate = &ocp_qp_clarabel_terminate;
    config->eval_forw_sens = &ocp_qp_clarabel_eval_forw_sens;
    config->eval_adj_sens = &ocp_qp_clarabel_eval_adj_sens;
    config->memory_reset = &ocp_qp_clarabel_memory_reset;
    config->solver_get = &ocp_qp_clarabel_solver_get;

    return;
}
