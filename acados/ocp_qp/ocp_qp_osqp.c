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

// blasfeo
#include "blasfeo_d_blasfeo_api.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_osqp.h"
#include "acados/utils/mem.h"
#include "acados/utils/math.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

// osqp
#include "osqp/include/public/osqp.h"




/************************************************
 * helper functions
 ************************************************/



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
    OSQPInt nn = 0, offset = 0, col = 0;
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
    OSQPInt nn = 0;
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

    OSQPInt row_offset_dyn = 0;
    OSQPInt row_offset_con = 0;
    OSQPInt row_offset_slk = 0;

    OSQPInt con_start = 0;
    OSQPInt slk_start = 0;
    for (kk = 0; kk <= N; kk++)
    {
        con_start += kk < N ? nx[kk + 1] : 0;
        slk_start += nb[kk]+ng[kk]+ns[kk];
    }

    slk_start += con_start;

    // CSC format: A_i are row indices and A_p are column pointers
    OSQPInt nn = 0, col = 0;
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
    OSQPInt nn = 0;
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
            mem->l[nn+nb[kk]+ng[kk]+ii] = -OSQP_INFTY;
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
                mem->u[nn + ii] = OSQP_INFTY;
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

        // OSQP_INFTY at upper bound
        for (ii = 0; ii < 2*ns[kk]; ii++)
        {
            mem->u[nn + ii] = OSQP_INFTY;
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

    assert((char *) raw_memory + ocp_qp_osqp_opts_calculate_size(config_, dims_) >= c_ptr);

    return (void *) opts;
}



void ocp_qp_osqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_osqp_opts *opts = opts_;

    osqp_set_default_settings(opts->osqp_opts);
    opts->osqp_opts->verbose = 0;
    opts->osqp_opts->polishing = 1;
    opts->osqp_opts->check_termination = 5;
    opts->osqp_opts->warm_starting = 1;

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

    // NOTE: settings are passed to OSQP at every call, via osqp_setup on the first call and
    // osqp_update_settings afterwards; some settings can only be set before the first call,
    // see the documentation of osqp_update_settings.
    if (!strcmp(field, "iter_max"))
    {
        int *tmp_ptr = value;
        opts->osqp_opts->max_iter = *tmp_ptr;
    }
    else if (!strcmp(field, "print_level"))
    {
        int *tmp_ptr = value;
        opts->print_level = *tmp_ptr;
        if (opts->print_level > 0)
        {
            opts->osqp_opts->verbose = 1;
        }
    }
    else if (!strcmp(field, "tol_stat"))
    {
        double *tol = value;
        // printf("in ocp_qp_osqp_opts_set, tol_stat %e\n", *tol);

        // opts->osqp_opts->eps_rel = *tol;
        // opts->osqp_opts->eps_dual_inf = *tol;

        opts->osqp_opts->eps_rel = MAX(*tol, 1e-5);
        opts->osqp_opts->eps_dual_inf = MAX(*tol, 1e-5);

        if (*tol <= 1e-3)
        {
            opts->osqp_opts->polishing = 1;
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
        int *tmp_ptr = value;
        opts->osqp_opts->warm_starting = *tmp_ptr;
    }
    else
    {
        printf("\nerror: ocp_qp_osqp_opts_set: wrong field: %s\n", field);
        exit(1);
    }

    return;
}


void ocp_qp_osqp_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    // ocp_qp_osqp_opts *opts = opts_;
    printf("\nerror: ocp_qp_osqp_opts_get: not implemented for field %s\n", field);
    exit(1);
}



/************************************************
 * memory
 ************************************************/



acados_size_t ocp_qp_osqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;

    size_t n = acados_osqp_num_vars(dims);
    size_t m = acados_osqp_num_constr(dims);

    size_t P_nnzmax = acados_osqp_nnzmax_P(dims);
    size_t A_nnzmax = acados_osqp_nnzmax_A(dims);

    acados_size_t size = 0;
    size += sizeof(ocp_qp_osqp_memory);

    size += 1 * n * sizeof(OSQPFloat);  // q
    size += 2 * m * sizeof(OSQPFloat);  // l, u

    size += P_nnzmax * sizeof(OSQPFloat);  // P_x
    size += P_nnzmax * sizeof(OSQPInt);    // P_i
    size += (n + 1) * sizeof(OSQPInt);     // P_p

    size += A_nnzmax * sizeof(OSQPFloat);  // A_x
    size += A_nnzmax * sizeof(OSQPInt);    // A_i
    size += (n + 1) * sizeof(OSQPInt);     // A_p

    size += 2 * sizeof(OSQPCscMatrix);  // matrices P and A

    size += 1 * 8;

    return size;
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
    mem->q = (OSQPFloat *) c_ptr;
    c_ptr += n * sizeof(OSQPFloat);

    mem->l = (OSQPFloat *) c_ptr;
    c_ptr += m * sizeof(OSQPFloat);

    mem->u = (OSQPFloat *) c_ptr;
    c_ptr += m * sizeof(OSQPFloat);

    mem->P_x = (OSQPFloat *) c_ptr;
    c_ptr += (mem->P_nnzmax) * sizeof(OSQPFloat);

    mem->A_x = (OSQPFloat *) c_ptr;
    c_ptr += (mem->A_nnzmax) * sizeof(OSQPFloat);

    // ints
    mem->P_i = (OSQPInt *) c_ptr;
    c_ptr += (mem->P_nnzmax) * sizeof(OSQPInt);

    mem->P_p = (OSQPInt *) c_ptr;
    c_ptr += (n + 1) * sizeof(OSQPInt);

    mem->A_i = (OSQPInt *) c_ptr;
    c_ptr += (mem->A_nnzmax) * sizeof(OSQPInt);

    mem->A_p = (OSQPInt *) c_ptr;
    c_ptr += (n + 1) * sizeof(OSQPInt);

    mem->P = (OSQPCscMatrix *) c_ptr;
    c_ptr += sizeof(OSQPCscMatrix);

    mem->A = (OSQPCscMatrix *) c_ptr;
    c_ptr += sizeof(OSQPCscMatrix);

    // initialize matrix structs; the array pointers remain acados-owned (owned = 0)
    OSQPCscMatrix_set_data(mem->P, n, n, P_nnzmax, mem->P_x, mem->P_i, mem->P_p);
    OSQPCscMatrix_set_data(mem->A, m, n, A_nnzmax, mem->A_x, mem->A_i, mem->A_p);

    // the OSQPSolver is opaque and cannot be placed in acados-managed memory;
    // it is allocated by osqp_setup at the first call and freed in ocp_qp_osqp_terminate
    mem->osqp_solver = NULL;

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

    OSQPInt con_start = 0;
    OSQPInt slk_start = 0;
    for (kk = 0; kk <= N; kk++)
    {
        con_start += kk < N ? nx[kk + 1] : 0;
        slk_start += nb[kk] + ng[kk] + ns[kk];
    }

    slk_start += con_start;

    OSQPSolution *sol = mem->osqp_solver->solution;

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
    if (mem->osqp_solver != NULL)
        osqp_cleanup(mem->osqp_solver);
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

    // update osqp solver with new data
    if (!mem->first_run)
    {
        osqp_update_data_vec(mem->osqp_solver, mem->q, mem->l, mem->u);
        osqp_update_data_mat(mem->osqp_solver, mem->P_x, NULL, mem->P->p[mem->P->n],
                             mem->A_x, NULL, mem->A->p[mem->A->n]);
        osqp_update_settings(mem->osqp_solver, opts->osqp_opts);
    }
    else
    {
        if (osqp_setup(&mem->osqp_solver, mem->P, mem->q, mem->A, mem->l, mem->u,
                       mem->A->m, mem->P->n, opts->osqp_opts) != 0)
        {
            printf("\nerror: ocp_qp_osqp: osqp_setup failed\n");
            exit(1);
        }
        mem->first_run = 0;
    }

    // solve OSQP
    acados_tic(&solver_call_timer);
    osqp_solve(mem->osqp_solver);
    mem->time_qp_solver_call = acados_toc(&solver_call_timer);
    mem->iter = mem->osqp_solver->info->iter;

    // fill qp_out
    fill_in_qp_out(qp_in, qp_out, mem);
    ocp_qp_compute_t(qp_in, qp_out);

    // info
    info->solve_QP_time = acados_toc(&qp_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = mem->osqp_solver->info->iter;
    info->t_computed = 1;

    OSQPInt osqp_status = mem->osqp_solver->info->status_val;
    int acados_status = osqp_status;

    // check exit conditions
    if (osqp_status == OSQP_SOLVED)
    {
        //printf("\nOSQP solved\n");
        acados_status = ACADOS_SUCCESS;
    }
    else if (osqp_status == OSQP_MAX_ITER_REACHED)
    {
        //printf("\nOSQP max iter reached\n");
        acados_status = ACADOS_MAXITER;
    }
    else if (osqp_status == OSQP_PRIMAL_INFEASIBLE)
    {
        //printf("\nOSQP primal infeasible\n");
        acados_status = ACADOS_INFEASIBLE;
    }
    else
    {
        acados_status = ACADOS_UNKNOWN;
    }
    mem->status = acados_status;

    return acados_status;
}



void ocp_qp_osqp_eval_forw_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_osqp_eval_forw_sens: not implemented yet\n");
    exit(1);
}

void ocp_qp_osqp_eval_adj_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_osqp_eval_adj_sens: not implemented yet\n");
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
    config->opts_get = &ocp_qp_osqp_opts_get;
    config->memory_calculate_size = &ocp_qp_osqp_memory_calculate_size;
    config->memory_assign = &ocp_qp_osqp_memory_assign;
    config->memory_get = &ocp_qp_osqp_memory_get;
    config->workspace_calculate_size = &ocp_qp_osqp_workspace_calculate_size;
    config->evaluate = &ocp_qp_osqp;
    config->terminate = &ocp_qp_osqp_terminate;
    config->eval_forw_sens = &ocp_qp_osqp_eval_forw_sens;
    config->eval_adj_sens = &ocp_qp_osqp_eval_adj_sens;
    config->memory_reset = &ocp_qp_osqp_memory_reset;
    config->solver_get = &ocp_qp_osqp_solver_get;

    return;
}
