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

// blasfeo
#include "blasfeo/include/blasfeo_d_blasfeo_api.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_osqp.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

// osqp
#include "osqp/include/auxil.h"
#include "osqp/include/constants.h"
#include "osqp/include/glob_opts.h"
#include "osqp/include/osqp.h"
#include "osqp/include/scaling.h"
#include "osqp/include/types.h"
#include "osqp/include/util.h"
#include "osqp/include/lin_sys.h"



#define OSQP_INF 1e8;



/************************************************
 * helper functions
 ************************************************/

#if 0
static void print_csc_as_dns(csc *M)
{
    c_int i, j = 0; // Predefine row index and column index
    c_int idx;

    // Initialize matrix of zeros
    c_float *A = (c_float *)c_calloc(M->m * M->n, sizeof(c_float));
    for (c_int ii=0; ii<M->m*M->n; ii++)
        A[ii] = 1e30;

    // Allocate elements
    for (idx = 0; idx < M->p[M->n]; idx++)
    {
        // Get row index i (starting from 1)
        i = M->i[idx];

        // Get column index j (increase if necessary) (starting from 1)
        while (M->p[j + 1] <= idx) j++;

        // Assign values to A
        A[j * (M->m) + i] = M->x[idx];
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
#endif



static void cpy_vec(c_int n, c_float *from, c_float *to)
{
    for (c_int ii = 0; ii < n; ii++) to[ii] = from[ii];
}



static void cpy_int_vec(c_int n, c_int *from, c_int *to)
{
    for (c_int ii = 0; ii < n; ii++) to[ii] = from[ii];
}



static void init_csc_matrix(c_int m, c_int n, c_int nzmax, c_float *x, c_int *i, c_int *p, csc *M)
{
    M->m = m;
    M->n = n;
    M->nz = -1;
    M->nzmax = nzmax;
    M->x = x;
    M->i = i;
    M->p = p;
}



static void cpy_csc_matrix(csc *from, csc *to)
{
    to->m = from->m;
    to->n = from->n;
    to->nz = -1;
    to->nzmax = from->nzmax;
    cpy_vec(from->nzmax, from->x, to->x);
    cpy_int_vec(from->nzmax, from->i, to->i);
    cpy_int_vec(from->n + 1, from->p, to->p);
}



static void set_vec(c_int n, c_float val, c_float *vec)
{
    for (c_int ii = 0; ii < n; ii++) vec[ii] = val;
}



static void set_int_vec(c_int n, c_int val, c_int *vec)
{
    for (c_int ii = 0; ii < n; ii++) vec[ii] = val;
}



static void cpy_osqp_settings(OSQPSettings *from, OSQPSettings *to)
{
    to->rho = from->rho;
    to->sigma = from->sigma;
    to->scaling = from->scaling;
    to->adaptive_rho = from->adaptive_rho;
    to->adaptive_rho_interval = from->adaptive_rho_interval;
    to->adaptive_rho_tolerance = from->adaptive_rho_tolerance;
    to->max_iter = from->max_iter;
    to->eps_abs = from->eps_abs;
    to->eps_rel = from->eps_rel;
    to->eps_prim_inf = from->eps_prim_inf;
    to->eps_dual_inf = from->eps_dual_inf;
    to->alpha = from->alpha;
    to->linsys_solver = from->linsys_solver;
    to->delta = from->delta;
    to->polish = from->polish;
    to->polish_refine_iter = from->polish_refine_iter;
    to->verbose = from->verbose;
    to->scaled_termination = from->scaled_termination;
    to->check_termination = from->check_termination;
    to->warm_start = from->warm_start;
}



// static void print_inputs(OSQPData *data)
// {
    // printf("\n----------> OSQP INPUTS <----------\n\n");
    // printf("NUMBER OF VARIABLES: %d\n", data->n);
    // printf("NUMBER OF CONSTRAINTS: %d\n", data->m);
    // printf("NUMBER OF NON-ZEROS in HESSIAN: %d\n", data->P->nzmax);
    // printf("NUMBER OF NON-ZEROS in CONSTRAINTS: %d\n", data->A->nzmax);
    // printf("\n-----------------------------------\n\n");

    // int ii;

    // printf("\nOBJECTIVE FUNCTION:\n");
    // for (ii = 0; ii < data->P->nzmax; ii++)
    //     printf("=====> P_x[%d] = %f, P_i[%d] = %d\n", ii + 1, data->P->x[ii], ii + 1,
    //            data->P->i[ii]);

    // for (ii = 0; ii < mem->A_nnzmax; ii++)
    // printf("=====> A_x[%d] = %f, A_i[%d] = %d\n", ii + 1, mem->A_x[ii], ii+1, mem->A_i[ii]);
    // print_csc_matrix(data->P, "Matrix P");
    // for (ii = 0; ii < mem->osqp_data->n; ii++)
    //     printf("=====> q[%d] = %f\n", ii + 1, mem->q[ii]);
    // for (ii = 0; ii < mem->osqp_data->n+1; ii++)
    //     printf("=====> P_p[%d] = %d\n", ii + 1, mem->P_p[ii]);

    // print_csc_as_dns(mem->osqp_data->P);

    // printf("\nBOUNDS:\n");
    // for (ii = 0; ii < data->m; ii++)
    //     printf("=====> l[%d] = %f, u[%d] = %f\n", ii + 1, data->l[ii], ii + 1, data->u[ii]);

    // printf("\nCONSTRAINTS MATRIX:\n");
    // print_csc_matrix(data->A, "Matrix A");
    // print_csc_as_dns(mem->osqp_data->A);
    // for (ii = 0; ii < data->A->nzmax; ii++)
    //     printf("=====> A_x[%d] = %f, A_i[%d] = %d\n", ii + 1, data->A->x[ii], ii + 1,
    //            data->A->i[ii]);
    // for (ii = 0; ii < mem->osqp_data->n+1; ii++)
    //     printf("=====> A_p[%d] = %d\n", ii + 1, mem->A_p[ii]);
// }



static int acados_osqp_num_vars(ocp_qp_dims *dims)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int n = 0;

    for (int ii = 0; ii <= N; ii++)
    {
        n += nx[ii];   // states
        n += nu[ii];   // controls
        n += 2*ns[ii]; // slacks
    }

    return n;
}



static int acados_osqp_num_constr(ocp_qp_dims *dims)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int m = 0;

    // inequality constraints
    for (int ii = 0; ii <= N; ii++)
    {
        m += nb[ii];   // box constraints
        m += ng[ii];   // general constraints
        m += ns[ii];   // replicated box/general softed constraint
        m += 2*ns[ii]; // slacks nonnegativity constraints
    }

    // dynamics equality constraints
    for (int ii = 0; ii < N; ii++)
    {
        m += nx[ii + 1];
    }

    return m;
}



static int acados_osqp_nnzmax_P(const ocp_qp_dims *dims)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int nnz = 0;

    for (int ii = 0; ii <= N; ii++)
    {
        nnz += nx[ii] * nx[ii];      // Q
        nnz += nu[ii] * nu[ii];      // R
        nnz += 2 * nx[ii] * nu[ii];  // S // TODO 1*, likely only the L or U are needed
        nnz += 2 * ns[ii];           // Z
    }

    return nnz;
}



static int acados_osqp_nnzmax_A(const ocp_qp_dims *dims)
{
    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int nnz = 0;

    // inequality constraints
    for (int ii = 0; ii <= N; ii++)
    {
        nnz += nb[ii];           // eye of box constraints
        nnz += ng[ii] * nx[ii];  // C
        nnz += ng[ii] * nu[ii];  // D
        nnz += ns[ii] * (nu[ii] + nx[ii]);  // replicated box/general softed constraint at worst case
        nnz += (nb[ii] + ng[ii]) * 2 * ns[ii]; // soft constraints at worst case, when idxs_rev encoding is used. Typically just 2*ns
        nnz += 2 * ns[ii];       // eye of slacks nonnegativity constraints
    }

    // dynamics equality constraints
    for (int ii = 0; ii < N; ii++)
    {
        nnz += nx[ii + 1] * nx[ii];  // A
        nnz += nx[ii + 1] * nu[ii];  // B
        nnz += nx[ii + 1];           // eye
    }

    return nnz;
}



static void update_gradient(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
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
}



static void update_hessian_structure(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int ii, jj, kk;

    // CSC format: P_i are row indices and P_p are column pointers
    c_int nn = 0, offset = 0, col = 0;
    for (kk = 0; kk <= N; kk++)
    {
        // write RSQ[kk]
        for (jj = 0; jj < nx[kk] + nu[kk]; jj++)
        {
            mem->P_p[col] = nn;
            col++;

            for (ii = 0; ii <= jj; ii++)
            {
                // we write only the upper triangular part
                mem->P_i[nn] = offset + ii;
                nn++;
            }
        }
        offset += nx[kk] + nu[kk];

        // write Z[kk]
        for (jj = 0; jj < 2*ns[kk]; jj++)
        {
            mem->P_p[col] = nn;
            col++;

            // diagonal
            mem->P_i[nn] = offset + jj;
            nn++;
        }

        offset += 2*ns[kk];
    }

    mem->P_p[col] = nn;
}



static void update_hessian_data(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ns = dims->ns;

    int ii, kk;

    // Traversing the matrix in column-major order
    c_int nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        // writing RSQ[kk]
        // we write the lower triangular part in row-major order
        // that's the same as writing the upper triangular part in
        // column-major order
        for (ii = 0; ii < nx[kk] + nu[kk]; ii++)
        {
            blasfeo_unpack_dmat(1, ii+1, in->RSQrq+kk, ii, 0, mem->P_x+nn, 1);
            nn += ii+1;
        }

        // write Z[kk]
        blasfeo_unpack_dvec(2*ns[kk], in->Z+kk, 0, mem->P_x+nn, 1);
        nn += 2*ns[kk];
    }
}



static void update_constraints_matrix_structure(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int ii, jj, kk;

    c_int row_offset_dyn = 0;
    c_int row_offset_con = 0;
    c_int row_offset_slk = 0;

    c_int con_start = 0;
    c_int slk_start = 0;
    for (kk = 0; kk <= N; kk++)
    {
        con_start += kk < N ? nx[kk + 1] : 0;
        slk_start += nb[kk]+ng[kk]+ns[kk];
    }

    slk_start += con_start;

    // CSC format: A_i are row indices and A_p are column pointers
    c_int nn = 0, col = 0;
    for (kk = 0; kk <= N; kk++)
    {

        // compute number of softed box constraints
        int nsb = 0;
        for (ii=0; ii<nb[kk]; ii++)
        {
            if (in->idxs_rev[kk][ii]>=0)
            {
                nsb++;
            }
        }

        // control variables
        for (jj = 0; jj < nu[kk]; jj++)
        {
            mem->A_p[col] = nn;
            col++;

            if (kk < dims->N)
            {
                // write column from B
                for (ii = 0; ii < nx[kk + 1]; ii++)
                {
                    mem->A_i[nn] = row_offset_dyn + ii;
                    nn++;
                }
            }

            // write bound on u
            for (ii = 0; ii < nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj)
                {
                    mem->A_i[nn] = con_start + row_offset_con + ii;
                    nn++;
                    break;
                }
            }
            int idxbu = ii;

            // write column from D
            for (ii = 0; ii < ng[kk]; ii++)
            {
                mem->A_i[nn] = con_start + row_offset_con + nb[kk] + ii;
                nn++;
            }

            // replicated softed bound on u
            if (idxbu<nb[kk]) // bounded input
            {
                if (in->idxs_rev[kk][idxbu]>=0) // softed bounded input
                {
                    // compute position in "packed" soft constraints, i.e. it is the itmp-th one
                    int itmp = 0;
                    for (ii=0; ii<idxbu; ii++)
                    {
                        if (in->idxs_rev[kk][ii]>=0)
                        {
                            itmp++;
                        }
                    }
                    mem->A_i[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + itmp;
                    nn++;
                }
            }

            // replicated softed D
            {
                int itmp = 0;
                for (ii = 0; ii < ng[kk]; ii++)
                {
                    if (in->idxs_rev[kk][nb[kk]+ii]>=0) // softed
                    {
                        mem->A_i[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + nsb + itmp;
                        nn++;
                        itmp++;
                    }
                }
            }

        }

        // state variables
        for (jj = 0; jj < nx[kk]; jj++)
        {
            mem->A_p[col] = nn;
            col++;

            if (kk > 0)
            {
                // write column from -I
                mem->A_i[nn] = row_offset_dyn - nx[kk] + jj;
                nn++;
            }

            if (kk < N)
            {
                // write column from A
                for (ii = 0; ii < nx[kk + 1]; ii++)
                {
                    mem->A_i[nn] = row_offset_dyn + ii;
                    nn++;
                }
            }

            // write bound on x
            for (ii = 0; ii < nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == nu[kk] + jj)
                {
                    mem->A_i[nn] = con_start + row_offset_con + ii;
                    nn++;
                    break;
                }
            }
            int idxbx = ii;

            // write column from C
            for (ii = 0; ii < ng[kk]; ii++)
            {
                mem->A_i[nn] = con_start + row_offset_con + nb[kk] + ii;
                nn++;
            }

            // replicated softed bound on x
            if (idxbx<nb[kk]) // bounded input
            {
                if (in->idxs_rev[kk][idxbx]>=0) // softed bounded input
                {
                    // compute position in "packed" soft constraints, i.e. it is the itmp-th one
                    int itmp = 0;
                    for (ii=0; ii<idxbx; ii++)
                    {
                        if (in->idxs_rev[kk][ii]>=0)
                        {
                            itmp++;
                        }
                    }
                    mem->A_i[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + itmp;
                    nn++;
                }
            }

            // replicated softed C
            {
                int itmp = 0;
                for (ii = 0; ii < ng[kk]; ii++)
                {
                    if (in->idxs_rev[kk][nb[kk]+ii]>=0) // softed
                    {
                        mem->A_i[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + nsb + itmp;
                        nn++;
                        itmp++;
                    }
                }
            }

        }

        // slack variables on lower inequalities (original)
        for (jj = 0; jj < ns[kk]; jj++)
        {
            mem->A_p[col] = nn;
            col++;

            // soft constraint
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    mem->A_i[nn] = con_start + row_offset_con + ii;
                    nn++;
                    // no break, there could possibly be multiple
                }
            }

            // nonnegativity constraint
            mem->A_i[nn] = slk_start + row_offset_slk + jj;
            nn++;
        }

        // slack variables on upper inequalities (replicated)
        for (jj = 0; jj < ns[kk]; jj++)
        {
            mem->A_p[col] = nn;
            col++;

            // soft constraint
            int itmp = 0;
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    mem->A_i[nn] = con_start + row_offset_con + nb[kk] + ng[kk] + itmp;
                    nn++;
                    // no break, there could possibly be multiple
                }
                if (in->idxs_rev[kk][ii]>=0)
                {
                    itmp++;
                }
            }

            // nonnegativity constraint
            mem->A_i[nn] = slk_start + row_offset_slk + ns[kk] + jj;
            nn++;
        }

        row_offset_con += nb[kk]+ng[kk]+ns[kk];
        row_offset_dyn += kk < N ? nx[kk + 1] : 0;
        row_offset_slk += 2*ns[kk];
    }

    // end of matrix
    mem->A_p[col] = nn;
}



// TODO move constant stuff like I to structure routine
static void update_constraints_matrix_data(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
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
    c_int nn = 0;
    for (kk = 0; kk <= N; kk++)
    {

        // control variables
        for (jj = 0; jj < nu[kk]; jj++)
        {
            if (kk < dims->N)
            {
                // write column from B
                blasfeo_unpack_dmat(1, nx[kk+1], in->BAbt+kk, jj, 0, mem->A_x+nn, 1);
                nn += nx[kk+1];
            }

            // write bound on u
            for (ii = 0; ii < dims->nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj)
                {
                    mem->A_x[nn] = 1.0;
                    nn++;
                    break;
                }
            }
            int idxbu = ii;

            // write column from D
            blasfeo_unpack_dmat(1, ng[kk], in->DCt+kk, jj, 0, mem->A_x+nn, 1);
            nn += ng[kk];

            // replicated softed bound on u
            if (idxbu<nb[kk]) // bounded input
            {
                if (in->idxs_rev[kk][idxbu]>=0) // softed bounded input
                {
                    mem->A_x[nn] = 1.0;
                    nn++;
                }
            }

            // replicated softed D
            for (ii = 0; ii < ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][nb[kk]+ii]>=0) // softed
                {
                    mem->A_x[nn] = BLASFEO_DMATEL(in->DCt+kk, jj, ii);
                    nn++;
                }
            }

        }

        // state variables
        for (jj = 0; jj < nx[kk]; jj++)
        {
            if (kk > 0)
            {
                // write column from -I
                mem->A_x[nn] = -1.0;
                nn++;
            }

            if (kk < N)
            {
                // write column from A
                blasfeo_unpack_dmat(1, nx[kk+1], in->BAbt+kk, nu[kk]+jj, 0, mem->A_x+nn, 1);
                nn += nx[kk+1];
            }

            // write bound on x
            for (ii = 0; ii < dims->nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == dims->nu[kk] + jj)
                {
                    mem->A_x[nn] = 1.0;
                    nn++;
                    break;
                }
            }
            int idxbx = ii;

            // write column from C
            blasfeo_unpack_dmat(1, ng[kk], in->DCt+kk, nu[kk]+jj, 0, mem->A_x+nn, 1);
            nn += ng[kk];

            // replicated softed bound on x
            if (idxbx<nb[kk]) // bounded input
            {
                if (in->idxs_rev[kk][idxbx]>=0) // softed bounded input
                {
                    mem->A_x[nn] = 1.0;
                    nn++;
                }
            }

            // replicated softed C
            for (ii = 0; ii < ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][nb[kk]+ii]>=0) // softed
                {
                    mem->A_x[nn] = BLASFEO_DMATEL(in->DCt+kk, nu[kk]+jj, ii);
                    nn++;
                }
            }

        }

        // slack variables on lower inequalities (original)
        for (jj = 0; jj < ns[kk]; jj++)
        {

            // soft constraint
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    mem->A_x[nn] = 1.0;
                    nn++;
                    // no break, there could possibly be multiple
                }
            }

            // nonnegativity constraint
            mem->A_x[nn] = 1.0;
            nn++;
        }

        // slack variables on upper inequalities (replicated)
        for (jj = 0; jj < ns[kk]; jj++)
        {

            // soft constraint
            for (ii=0; ii<nb[kk]+ng[kk]; ii++)
            {
                if (in->idxs_rev[kk][ii]==jj)
                {
                    mem->A_x[nn] = -1.0;
                    nn++;
                    // no break, there could possibly be multiple
                }
            }

            // nonnegativity constraint
            mem->A_x[nn] = 1.0;
            nn++;
        }

    }

}



static void update_bounds(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    //int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int ii, kk, nn = 0;

    // write -b to l and u
    for (kk = 0; kk < N; kk++)
    {
        // unpack b to l
        blasfeo_unpack_dvec(nx[kk + 1], in->b + kk, 0, &mem->l[nn], 1);

        // change sign of l (to get -b) and copy to u
        for (ii = 0; ii < nx[kk + 1]; ii++)
        {
            mem->l[nn + ii] = -mem->l[nn + ii];
            mem->u[nn + ii] = mem->l[nn + ii];
        }

        nn += nx[kk + 1];
    }

    // write lb lg and ub ug
    for (kk = 0; kk <= N; kk++)
    {
        // unpack lb lg to l
        blasfeo_unpack_dvec(nb[kk]+ng[kk], in->d + kk, 0, &mem->l[nn], 1);
        // set replicated to -inf
        for (ii=0; ii<ns[kk]; ii++)
        {
            mem->l[nn+nb[kk]+ng[kk]+ii] = -OSQP_INF;
        }

        // unpack ub ug to u and flip signs because in HPIPM the signs are flipped for upper bounds
        // keep in original place if not softed; set to inf and replicated under if softed
        int itmp = 0;
        for (ii = 0; ii < nb[kk] + ng[kk]; ii++)
        {
            if (in->idxs_rev[kk][ii]==-1) // not softed
            {
                mem->u[nn + ii] = -BLASFEO_DVECEL(&in->d[kk], nb[kk]+ng[kk]+ii);
            }
            else
            {
                mem->u[nn + ii] = OSQP_INF;
                mem->u[nn + nb[kk]+ng[kk]+itmp] = -BLASFEO_DVECEL(&in->d[kk], nb[kk]+ng[kk]+ii);
                itmp++;
            }
        }

        nn += nb[kk] + ng[kk] + ns[kk];
    }

    // write ls and us
    for (kk = 0; kk <= N; kk++)
    {
        // unpack ls and us to l
        blasfeo_unpack_dvec(2*ns[kk], in->d + kk, 2*nb[kk]+2*ng[kk], &mem->l[nn], 1);

        // OSQP_INF at upper bound
        for (ii = 0; ii < 2*ns[kk]; ii++)
        {
            mem->u[nn + ii] = OSQP_INF;
        }

        nn += 2*ns[kk];
    }

}



static void ocp_qp_osqp_update_memory(const ocp_qp_in *in, const ocp_qp_osqp_opts *opts,
                                      ocp_qp_osqp_memory *mem)
{
    if (mem->first_run)
    {
        update_hessian_structure(in, mem);
        update_constraints_matrix_structure(in, mem);
    }

    update_bounds(in, mem);
    update_gradient(in, mem);
    update_hessian_data(in, mem);
    update_constraints_matrix_data(in, mem);

    //printf("\nP\n");
    //print_csc_as_dns(mem->osqp_data->P);
    //printf("\nA\n");
    //print_csc_as_dns(mem->osqp_data->A);
    //exit(1);
}


/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_osqp_opts_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;
    size += sizeof(ocp_qp_osqp_opts);
    size += sizeof(OSQPSettings);

    return size;
}



void *ocp_qp_osqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_qp_osqp_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (ocp_qp_osqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_osqp_opts);

    opts->osqp_opts = (OSQPSettings *) c_ptr;
    c_ptr += sizeof(OSQPSettings);

    assert((char *) raw_memory + ocp_qp_osqp_opts_calculate_size(config_, dims_) == c_ptr);

    return (void *) opts;
}



void ocp_qp_osqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_osqp_opts *opts = opts_;

    osqp_set_default_settings(opts->osqp_opts);
    opts->osqp_opts->verbose = 0;
    opts->osqp_opts->polish = 1;
    opts->osqp_opts->check_termination = 5;
    opts->osqp_opts->warm_start = 1;

    return;
}



void ocp_qp_osqp_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_qp_osqp_opts *opts = (ocp_qp_osqp_opts *)opts_;

    return;
}

void ocp_qp_osqp_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_osqp_opts *opts = opts_;

    // NOTE/TODO(oj): options are copied into OSQP at first call.
    // Updating options through this function does not work, only before the first call!
    if (!strcmp(field, "iter_max"))
    {
        int *tmp_ptr = value;
        opts->osqp_opts->max_iter = *tmp_ptr;
    }
    else if (!strcmp(field, "tol_stat"))
    {
        double *tol = value;
        // printf("in ocp_qp_osqp_opts_set, tol_stat %e\n", *tol);

        // opts->osqp_opts->eps_rel = *tol;
        // opts->osqp_opts->eps_dual_inf = *tol;

        opts->osqp_opts->eps_rel = fmax(*tol, 1e-5);
        opts->osqp_opts->eps_dual_inf = fmax(*tol, 1e-5);

        if (*tol <= 1e-3)
        {
            opts->osqp_opts->polish = 1;
            opts->osqp_opts->polish_refine_iter = 5;
        }
    }
    else if (!strcmp(field, "tol_eq"))
    {
        double *tol = value;
        opts->osqp_opts->eps_prim_inf = *tol;
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        double *tol = value;
        opts->osqp_opts->eps_prim_inf = *tol;
    }
    else if (!strcmp(field, "tol_comp"))
    {
        // "OSQP always satisfies complementary slackness conditions
        //  with machine precision by construction." - Strellato2020
    }
    else if (!strcmp(field, "warm_start"))
    {
        // XXX after the first call to the solver, this doesn't work any more, as in osqp the settings are copied in the work !!!!!
        // XXX i.e. as it is, it gets permanently set to zero if warm start is disabled at the fist iteration !!!!!
        int *tmp_ptr = value;
        // int tmp_ptr[] = {1};
        opts->osqp_opts->warm_start = *tmp_ptr;
        // printf("\nwarm start %d\n", opts->osqp_opts->warm_start);
    }
    else
    {
        printf("\nerror: ocp_qp_osqp_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * memory
 ************************************************/

static acados_size_t osqp_workspace_calculate_size(int n, int m, int P_nnzmax, int A_nnzmax)
{
    acados_size_t size = 0;

    size += sizeof(OSQPWorkspace);
    size += sizeof(OSQPData);
    size += 2 * sizeof(csc);

    size += 1 * n * sizeof(c_float);  // q
    size += 2 * m * sizeof(c_float);  // l, u

    size += P_nnzmax * sizeof(c_float);  // P_x
    size += P_nnzmax * sizeof(c_int);    // P_i
    size += (n + 1) * sizeof(c_int);     // P_p

    size += A_nnzmax * sizeof(c_float);  // A_x
    size += A_nnzmax * sizeof(c_int);    // A_i
    size += (n + 1) * sizeof(c_int);     // A_p

    size += 2 * m * sizeof(c_float);  // rho_vec, rho_inv_vec
    size += m * sizeof(c_int);        // constr_type

    size += n * sizeof(c_float);        // x
    size += m * sizeof(c_float);        // z
    size += (n + m) * sizeof(c_float);  // xz_tilde
    size += n * sizeof(c_float);        // x_prev
    size += m * sizeof(c_float);        // z_prev
    size += m * sizeof(c_float);        // y

    size += m * sizeof(c_float);  // Ax
    size += n * sizeof(c_float);  // Px
    size += n * sizeof(c_float);  // Aty

    size += m * sizeof(c_float);  // delta_y
    size += n * sizeof(c_float);  // Atdelta_y

    size += n * sizeof(c_float);  // delta_x
    size += n * sizeof(c_float);  // Pdelta_x
    size += m * sizeof(c_float);  // Adelta_x

    size += sizeof(OSQPSettings);  // settings

    size += sizeof(OSQPScaling);  // scaling
    size += n * sizeof(c_float);  // scaling->D
    size += n * sizeof(c_float);  // scaling->Dinv
    size += m * sizeof(c_float);  // scaling->E
    size += m * sizeof(c_float);  // scaling->Einv

    size += n * sizeof(c_float);  // D_temp
    size += n * sizeof(c_float);  // D_temp_A
    size += m * sizeof(c_float);  // E_temp

    size += sizeof(OSQPPolish);   // pol
    size += m * sizeof(c_int);    // pol->Alow_to_A
    size += m * sizeof(c_int);    // pol->Aupp_to_A
    size += m * sizeof(c_int);    // pol->A_to_Alow
    size += m * sizeof(c_int);    // pol->A_to_Aupp
    size += n * sizeof(c_float);  // pol->x
    size += m * sizeof(c_float);  // pol->z
    size += m * sizeof(c_float);  // pol->y

    size += sizeof(OSQPSolution);  // solution
    size += n * sizeof(c_float);   // solution->x
    size += m * sizeof(c_float);   // solution->y

    size += sizeof(OSQPInfo);  // info

    size += 1 * 8;

    return size;
}


acados_size_t ocp_qp_osqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;

    size_t n = acados_osqp_num_vars(dims);
    size_t m = acados_osqp_num_constr(dims);

    size_t P_nnzmax = acados_osqp_nnzmax_P(dims);
    size_t A_nnzmax = acados_osqp_nnzmax_A(dims);

    acados_size_t size = 0;
    size += sizeof(ocp_qp_osqp_memory);

    size += 1 * n * sizeof(c_float);  // q
    size += 2 * m * sizeof(c_float);  // l, u

    size += P_nnzmax * sizeof(c_float);  // P_x
    size += P_nnzmax * sizeof(c_int);    // P_i
    size += (n + 1) * sizeof(c_int);     // P_p

    size += A_nnzmax * sizeof(c_float);  // A_x
    size += A_nnzmax * sizeof(c_int);    // A_i
    size += (n + 1) * sizeof(c_int);     // A_p

    size += sizeof(OSQPData);
    size += 2 * sizeof(csc);  // matrices P and A
    size += osqp_workspace_calculate_size(n, m, P_nnzmax, A_nnzmax);

    size += 1 * 8;

    return size;
}



static void *osqp_workspace_assign(int n, int m, int P_nnzmax, int A_nnzmax, void *raw_memory)
{
    OSQPWorkspace *work;

    // char pointer
    char *c_ptr = (char *) raw_memory;

    work = (OSQPWorkspace *) c_ptr;
    c_ptr += sizeof(OSQPWorkspace);

    work->data = (OSQPData *) c_ptr;
    c_ptr += sizeof(OSQPData);

    work->data->P = (csc *) c_ptr;
    c_ptr += sizeof(csc);

    work->data->A = (csc *) c_ptr;
    c_ptr += sizeof(csc);

    work->settings = (OSQPSettings *) c_ptr;
    c_ptr += sizeof(OSQPSettings);

    work->scaling = (OSQPScaling *) c_ptr;
    c_ptr += sizeof(OSQPScaling);

    work->pol = (OSQPPolish *) c_ptr;
    c_ptr += sizeof(OSQPPolish);

    work->solution = (OSQPSolution *) c_ptr;
    c_ptr += sizeof(OSQPSolution);

    work->info = (OSQPInfo *) c_ptr;
    c_ptr += sizeof(OSQPInfo);

    align_char_to(8, &c_ptr);

    // doubles
    work->data->q = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->data->l = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->data->u = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->data->P->x = (c_float *) c_ptr;
    c_ptr += P_nnzmax * sizeof(c_float);

    work->data->A->x = (c_float *) c_ptr;
    c_ptr += A_nnzmax * sizeof(c_float);

    work->rho_vec = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->rho_inv_vec = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->x = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->z = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->xz_tilde = (c_float *) c_ptr;
    c_ptr += (n + m) * sizeof(c_float);

    work->x_prev = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->z_prev = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->y = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->Ax = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->Px = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->Aty = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->delta_y = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->Atdelta_y = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->delta_x = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->Pdelta_x = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->Adelta_x = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->scaling->D = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->scaling->Dinv = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->scaling->E = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->scaling->Einv = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->D_temp = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->D_temp_A = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->E_temp = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->pol->x = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->pol->z = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->pol->y = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    work->solution->x = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    work->solution->y = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    // integers

    work->data->P->i = (c_int *) c_ptr;
    c_ptr += P_nnzmax * sizeof(c_int);

    work->data->P->p = (c_int *) c_ptr;
    c_ptr += (n + 1) * sizeof(c_int);

    work->data->A->i = (c_int *) c_ptr;
    c_ptr += A_nnzmax * sizeof(c_int);

    work->data->A->p = (c_int *) c_ptr;
    c_ptr += (n + 1) * sizeof(c_int);

    work->constr_type = (c_int *) c_ptr;
    c_ptr += m * sizeof(c_int);

    work->pol->Alow_to_A = (c_int *) c_ptr;
    c_ptr += m * sizeof(c_int);

    work->pol->Aupp_to_A = (c_int *) c_ptr;
    c_ptr += m * sizeof(c_int);

    work->pol->A_to_Alow = (c_int *) c_ptr;
    c_ptr += m * sizeof(c_int);

    work->pol->A_to_Aupp = (c_int *) c_ptr;
    c_ptr += m * sizeof(c_int);

    return work;
}



static int osqp_init_data(OSQPData *data, OSQPSettings *settings, OSQPWorkspace *work)
{
    c_int n = data->n;
    c_int m = data->m;

    // Copy problem data into workspace
    work->data->n = n;
    work->data->m = m;
    cpy_csc_matrix(data->P, work->data->P);
    cpy_csc_matrix(data->A, work->data->A);
    cpy_vec(m, data->l, work->data->l);
    cpy_vec(m, data->u, work->data->u);
    cpy_vec(n, data->q, work->data->q);

    // Vectorized rho parameter
    set_vec(m, 0.0, work->rho_vec);
    set_vec(m, 0.0, work->rho_inv_vec);

    // Type of constraints
    set_int_vec(m, 0, work->constr_type);

    // Allocate internal solver variables (ADMM steps)
    set_vec(n, 0.0, work->x);
    set_vec(m, 0.0, work->z);
    set_vec(n + m, 0.0, work->xz_tilde);
    set_vec(n, 0.0, work->x_prev);
    set_vec(m, 0.0, work->z_prev);
    set_vec(m, 0.0, work->y);

    // Initialize variables x, y, z to 0
    cold_start(work);

    // Primal and dual residuals variables
    set_vec(m, 0.0, work->Ax);
    set_vec(n, 0.0, work->Px);
    set_vec(n, 0.0, work->Aty);

    // Primal infeasibility variables
    set_vec(m, 0.0, work->delta_y);
    set_vec(n, 0.0, work->Atdelta_y);

    // Dual infeasibility variables
    set_vec(n, 0.0, work->delta_x);
    set_vec(n, 0.0, work->Pdelta_x);
    set_vec(m, 0.0, work->Adelta_x);

    // Copy settings
    cpy_osqp_settings(settings, work->settings);

    if (settings->scaling)
        scale_data(work);  // Scale data
    else
        work->scaling = OSQP_NULL;

    // Set type of constraints
    set_rho_vec(work);

    // Load linear system solver
    if (load_linsys_solver(work->settings->linsys_solver))
    {
        c_eprint("%s linear system solver not available.\nTried to obtain it from shared library",
        LINSYS_SOLVER_NAME[work->settings->linsys_solver]);
        return 0;
    }

    // Initialize linear system solver structure
    // NOTE: mallocs, memory is freed in ocp_qp_osqp_terminate
    if (init_linsys_solver(&(work->linsys_solver), work->data->P, work->data->A, work->settings->sigma,
                                             work->rho_vec, work->settings->linsys_solver, 0)){
        c_eprint("Failed to initialize %s linear system solver",
                 LINSYS_SOLVER_NAME[work->settings->linsys_solver]);
        return 0;
    }


    // Initialize solution to 0
    set_vec(n, 0.0, work->solution->x);
    set_vec(m, 0.0, work->solution->y);

    // Initialize info
    update_status(work->info, OSQP_UNSOLVED);
    work->info->iter = 0;
    work->info->status_polish = 0;
    work->info->obj_val = 0.0;
    work->info->pri_res = 0.0;
    work->info->dua_res = 0.0;
    work->info->rho_updates = 0;
    work->info->rho_estimate = work->settings->rho;

    if (work->settings->verbose) print_setup_header(work);
    work->summary_printed = 0; // Initialize last summary  to not printed

    if (work->settings->adaptive_rho && !work->settings->adaptive_rho_interval)
    {
        if (work->settings->check_termination)
        {
            // If check_termination is enabled, we set it to a multiple of the check
            // termination interval
            work->settings->adaptive_rho_interval =
                ADAPTIVE_RHO_MULTIPLE_TERMINATION * work->settings->check_termination;
        }
        else
        {
            // If check_termination is disabled we set it to a predefined fix number
            work->settings->adaptive_rho_interval = ADAPTIVE_RHO_FIXED;
        }
    }

    return 1;
}



void *ocp_qp_osqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    UNUSED(opts_);

    ocp_qp_dims *dims = dims_;
    ocp_qp_osqp_memory *mem;

    int n = acados_osqp_num_vars(dims);
    int m = acados_osqp_num_constr(dims);
    int P_nnzmax = acados_osqp_nnzmax_P(dims);
    int A_nnzmax = acados_osqp_nnzmax_A(dims);

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (ocp_qp_osqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_osqp_memory);

    mem->P_nnzmax = P_nnzmax;
    mem->A_nnzmax = A_nnzmax;
    mem->first_run = 1;

    align_char_to(8, &c_ptr);

    // doubles
    mem->q = (c_float *) c_ptr;
    c_ptr += n * sizeof(c_float);

    mem->l = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    mem->u = (c_float *) c_ptr;
    c_ptr += m * sizeof(c_float);

    mem->P_x = (c_float *) c_ptr;
    c_ptr += (mem->P_nnzmax) * sizeof(c_float);

    mem->A_x = (c_float *) c_ptr;
    c_ptr += (mem->A_nnzmax) * sizeof(c_float);

    // ints
    mem->P_i = (c_int *) c_ptr;
    c_ptr += (mem->P_nnzmax) * sizeof(c_int);

    mem->P_p = (c_int *) c_ptr;
    c_ptr += (n + 1) * sizeof(c_int);

    mem->A_i = (c_int *) c_ptr;
    c_ptr += (mem->A_nnzmax) * sizeof(c_int);

    mem->A_p = (c_int *) c_ptr;
    c_ptr += (n + 1) * sizeof(c_int);

    mem->osqp_data = (OSQPData *) c_ptr;
    c_ptr += sizeof(OSQPData);

    mem->osqp_data->P = (csc *) c_ptr;
    c_ptr += sizeof(csc);

    mem->osqp_data->A = (csc *) c_ptr;
    c_ptr += sizeof(csc);

    // mem->osqp_work = (OSQPWorkspace *) c_ptr;
    // c_ptr += sizeof(OSQPWorkspace);
    mem->osqp_work = osqp_workspace_assign(n, m, P_nnzmax, A_nnzmax, c_ptr);
    c_ptr += osqp_workspace_calculate_size(n, m, P_nnzmax, A_nnzmax);

    // initialize data pointers
    OSQPData *data = mem->osqp_data;
    data->n = n;
    data->m = m;

    data->q = mem->q;
    data->l = mem->l;
    data->u = mem->u;

    init_csc_matrix(n, n, P_nnzmax, mem->P_x, mem->P_i, mem->P_p, data->P);
    init_csc_matrix(m, n, A_nnzmax, mem->A_x, mem->A_i, mem->A_p, data->A);

    assert((char *) raw_memory + ocp_qp_osqp_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_qp_osqp_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    // qp_solver_config *config = config_;
    ocp_qp_osqp_memory *mem = mem_;

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
        printf("\nerror: ocp_qp_osqp_memory_get: field %s not available\n", field);
        exit(1);
    }

    return;

}


void ocp_qp_osqp_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    // ocp_qp_in *qp_in = qp_in_;
    // reset memory
    printf("acados: reset osqp_mem not implemented.\n");
    exit(1);
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_osqp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    return 0;
}



/************************************************
 * functions
 ************************************************/

static void fill_in_qp_out(const ocp_qp_in *in, ocp_qp_out *out, ocp_qp_osqp_memory *mem)
{
    ocp_qp_dims *dims = in->dim;

    int N = dims->N;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *nb = dims->nb;
    int *ng = dims->ng;
    int *ns = dims->ns;

    int ii, kk, nn;

    c_int con_start = 0;
    c_int slk_start = 0;
    for (kk = 0; kk <= N; kk++)
    {
        con_start += kk < N ? nx[kk + 1] : 0;
        slk_start += nb[kk] + ng[kk] + ns[kk];
    }

    slk_start += con_start;

    OSQPSolution *sol = mem->osqp_work->solution;

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
        blasfeo_pack_dvec(nx[kk + 1], &sol->y[nn], 1, out->pi + kk, 0);
        nn += nx[kk + 1];
    }

    nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        //blasfeo_dvecse(2*nb[kk]+2*ng[kk]+2*ns[kk], 0.0, out->lam+kk, 0);
        blasfeo_dvecse(2*nb[kk]+2*ng[kk], 0.0, out->lam+kk, 0);

        int itmp = 0;
        for (ii = 0; ii < nb[kk]+ng[kk]; ii++)
        {
            if (in->idxs_rev[kk][ii]==-1) // not softed: two sided
            {
                double lam = sol->y[con_start+nn+ii];
                if (lam <= 0)
                {
                    BLASFEO_DVECEL(out->lam+kk, ii) = -lam;
                }
                else
                {
                    BLASFEO_DVECEL(out->lam+kk, nb[kk]+ng[kk] + ii) = lam;
                }
            }
            else // softed: replicated one sided
            {
                BLASFEO_DVECEL(out->lam+kk, ii) = -sol->y[con_start+nn+ii];
                BLASFEO_DVECEL(out->lam+kk, nb[kk]+ng[kk] + ii) = sol->y[con_start+nn+nb[kk]+ng[kk]+itmp];
                itmp++;
            }
        }
        nn += nb[kk]+ng[kk]+ns[kk];
    }

    nn = 0;
    for (kk = 0; kk <= N; kk++)
    {
        blasfeo_pack_dvec(2*ns[kk], &sol->y[slk_start+nn], 1, out->lam+kk, 2*nb[kk]+2*ng[kk]);
        blasfeo_dvecsc(2*ns[kk], -1.0, out->lam+kk, 2*nb[kk]+2*ng[kk]);
        nn += 2*ns[kk];
    }
}


void ocp_qp_osqp_terminate(void *config_, void *mem_, void *work_)
{
    ocp_qp_osqp_memory *mem = (ocp_qp_osqp_memory *) mem_;
    OSQPWorkspace *work = mem->osqp_work;
    work->linsys_solver->free(work->linsys_solver);
}




int ocp_qp_osqp(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;

    // print_ocp_qp_dims(qp_in->dim);
    // print_ocp_qp_in(qp_in);

    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer, solver_call_timer;
    acados_tic(&tot_timer);

    // cast data structures
    ocp_qp_osqp_opts *opts = (ocp_qp_osqp_opts *) opts_;
    ocp_qp_osqp_memory *mem = (ocp_qp_osqp_memory *) mem_;

    acados_tic(&interface_timer);
    ocp_qp_osqp_update_memory(qp_in, opts, mem);
    info->interface_time = acados_toc(&interface_timer);

    acados_tic(&qp_timer);

    // update osqp workspace with new data
    if (!mem->first_run)
    {
        osqp_update_lin_cost(mem->osqp_work, mem->q);
        osqp_update_P_A(mem->osqp_work, mem->P_x, NULL, mem->P_nnzmax, mem->A_x, NULL,
                        mem->A_nnzmax);
        osqp_update_bounds(mem->osqp_work, mem->l, mem->u);
        // TODO(oj): update OSQP options here if they were updated?
    }
    else
    {
        // mem->osqp_work = osqp_setup(mem->osqp_data, opts->osqp_opts);
        osqp_init_data(mem->osqp_data, opts->osqp_opts, mem->osqp_work);
        mem->first_run = 0;
    }

    // check settings:
    // OSQPSettings *settings = mem->osqp_work->settings;
    // printf("OSQP settings: warm_start %d\n", settings->warm_start);
    // printf("polish %d, polish_refine_iter %d, delta: %e\n", settings->polish, settings->polish_refine_iter, settings->delta);
    // printf("eps_abs %e, eps_rel %e, eps_prim_inf: %e, eps_dual_inf: %e\n", settings->eps_abs, settings->eps_rel, settings->eps_prim_inf, settings->eps_dual_inf);

    // solve OSQP
    acados_tic(&solver_call_timer);
    osqp_solve(mem->osqp_work);
    mem->time_qp_solver_call = acados_toc(&solver_call_timer);
    mem->iter = mem->osqp_work->info->iter;

    // fill qp_out
    fill_in_qp_out(qp_in, qp_out, mem);
    ocp_qp_compute_t(qp_in, qp_out);

    // info
    info->solve_QP_time = acados_toc(&qp_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = mem->osqp_work->info->iter;
    info->t_computed = 1;

    c_int osqp_status = mem->osqp_work->info->status_val;
    int acados_status = osqp_status;

    // check exit conditions
    if (osqp_status == OSQP_SOLVED)
    {
        //printf("\nOSQP solved\n");
        acados_status = ACADOS_SUCCESS;
    }
    if (osqp_status == OSQP_MAX_ITER_REACHED)
    {
        //printf("\nOSQP max iter reached\n");
        acados_status = ACADOS_MAXITER;
    }
    mem->status = acados_status;

    return acados_status;
}



void ocp_qp_osqp_eval_sens(void *config_, void *qp_in, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_osqp_eval_sens: not implemented yet\n");
    exit(1);
}


void ocp_qp_osqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: ocp_qp_osqp_solver_get: not implemented yet\n");
    exit(1);
}


void ocp_qp_osqp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = &ocp_qp_osqp_opts_calculate_size;
    config->opts_assign = &ocp_qp_osqp_opts_assign;
    config->opts_initialize_default = &ocp_qp_osqp_opts_initialize_default;
    config->opts_update = &ocp_qp_osqp_opts_update;
    config->opts_set = &ocp_qp_osqp_opts_set;
    config->memory_calculate_size = &ocp_qp_osqp_memory_calculate_size;
    config->memory_assign = &ocp_qp_osqp_memory_assign;
    config->memory_get = &ocp_qp_osqp_memory_get;
    config->workspace_calculate_size = &ocp_qp_osqp_workspace_calculate_size;
    config->evaluate = &ocp_qp_osqp;
    config->terminate = &ocp_qp_osqp_terminate;
    config->eval_sens = &ocp_qp_osqp_eval_sens;
    config->memory_reset = &ocp_qp_osqp_memory_reset;
    config->solver_get = &ocp_qp_osqp_solver_get;

    return;
}
