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

#include <assert.h>

// osqp
#include "osqp/include/auxil.h"
#include "osqp/include/constants.h"
#include "osqp/include/glob_opts.h"
#include "osqp/include/osqp.h"
#include "osqp/include/scaling.h"
#include "osqp/include/types.h"
#include "osqp/include/util.h"
#include "osqp/include/lin_sys.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_osqp.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"


/************************************************
 * helper functions
 ************************************************/

// static void print_csc_as_dns(csc *M)
// {
//     c_int i, j = 0; // Predefine row index and column index
//     c_int idx;

//     // Initialize matrix of zeros
//     c_float *A = (c_float *)c_calloc(M->m * M->n, sizeof(c_float));

//     // Allocate elements
//     for (idx = 0; idx < M->p[M->n]; idx++)
//     {
//         // Get row index i (starting from 1)
//         i = M->i[idx];

//         // Get column index j (increase if necessary) (starting from 1)
//         while (M->p[j + 1] <= idx) j++;

//         // Assign values to A
//         A[j * (M->m) + i] = M->x[idx];
//     }

//     for (i = 0; i < M->m; i++)
//     {
//         for (j = 0; j < M->n; j++)
//         {
//             printf("%f ", A[j * (M->m) + i]);
//         }
//         printf("\n");
//     }

//     free(A);
// }



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
    int n = 0;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        n += dims->nx[ii] + dims->nu[ii];
    }

    return n;
}



static int acados_osqp_num_constr(ocp_qp_dims *dims)
{
    int m = 0;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        m += dims->nb[ii];
        m += dims->ng[ii];

        if (ii < dims->N)
        {
            m += dims->nx[ii + 1];
        }
    }

    return m;
}



static int acados_osqp_nnzmax_P(const ocp_qp_dims *dims)
{
    int nnz = 0;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        nnz += dims->nx[ii] * dims->nx[ii];      // Q
        nnz += dims->nu[ii] * dims->nu[ii];      // R
        nnz += 2 * dims->nx[ii] * dims->nu[ii];  // S
    }

    return nnz;
}



static int acados_osqp_nnzmax_A(const ocp_qp_dims *dims)
{
    int nnz = 0;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        // inequality constraints
        nnz += dims->nb[ii];                 // eye
        nnz += dims->ng[ii] * dims->nx[ii];  // C
        nnz += dims->ng[ii] * dims->nu[ii];  // D

        // equality constraints
        if (ii < dims->N)
        {
            nnz += dims->nx[ii + 1] * dims->nx[ii];  // A
            nnz += dims->nx[ii + 1] * dims->nu[ii];  // B
            nnz += dims->nx[ii + 1];                 // eye
        }
    }

    return nnz;
}



static void update_gradient(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    int kk, nn = 0;
    ocp_qp_dims *dims = in->dim;

    for (kk = 0; kk <= dims->N; kk++)
    {
        blasfeo_unpack_dvec(dims->nu[kk] + dims->nx[kk], in->rqz + kk, 0, &mem->q[nn]);
        nn += dims->nu[kk] + dims->nx[kk];
    }
}



static void update_hessian_structure(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    c_int ii, jj, kk, nn = 0, offset = 0, col = 0;
    ocp_qp_dims *dims = in->dim;

    // CSC format: P_i are row indices and P_p are column pointers
    for (kk = 0; kk <= dims->N; kk++)
    {
        // writing RSQ[kk]
        for (jj = 0; jj < dims->nx[kk] + dims->nu[kk]; jj++)
        {
            mem->P_p[col++] = nn;

            for (ii = 0; ii <= jj; ii++)
            {
                // we write only the upper triangular part
                mem->P_i[nn++] = offset + ii;
            }
        }

        offset += dims->nx[kk] + dims->nu[kk];
    }

    mem->P_p[col] = nn;
}



static void update_hessian_data(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    c_int ii, jj, kk, nn = 0;
    ocp_qp_dims *dims = in->dim;

    // Traversing the matrix in column-major order
    for (kk = 0; kk <= dims->N; kk++)
    {
        // writing RSQ[kk]
        for (ii = 0; ii < dims->nx[kk] + dims->nu[kk]; ii++)
        {
            for (jj = 0; jj <= ii; jj++)
            {
                // we write the lower triangular part in row-major order
                // that's the same as writing the upper triangular part in
                // column-major order
                mem->P_x[nn++] = BLASFEO_DMATEL(&in->RSQrq[kk], ii, jj);
            }
        }
    }
}



static void update_constraints_matrix_structure(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    c_int ii, jj, kk, nn = 0, col = 0;
    c_int con_start = 0, bnd_start = 0;
    c_int row_offset_dyn = 0, row_offset_con = 0, row_offset_bnd = 0;
    ocp_qp_dims *dims = in->dim;

    for (kk = 0; kk <= dims->N; kk++)
    {
        con_start += kk < dims->N ? dims->nx[kk + 1] : 0;
        bnd_start += dims->ng[kk];
    }

    bnd_start += con_start;

    // CSC format: A_i are row indices and A_p are column pointers
    for (kk = 0; kk <= dims->N; kk++)
    {
        int nbu = 0;

        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            mem->A_p[col++] = nn;

            if (kk < dims->N)
            {
                // write column from B
                for (ii = 0; ii < dims->nx[kk + 1]; ii++)
                {
                    mem->A_i[nn++] = ii + row_offset_dyn;
                }
            }

            // write column from D
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                mem->A_i[nn++] = ii + con_start + row_offset_con;
            }

            // write bound on u
            for (ii = 0; ii < dims->nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj)
                {
                    mem->A_i[nn++] = ii + bnd_start + row_offset_bnd;
                    nbu++;
                    break;
                }
            }
        }

        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            mem->A_p[col++] = nn;

            if (kk > 0)
            {
                // write column from -I
                mem->A_i[nn++] = jj + row_offset_dyn - dims->nx[kk];
            }

            if (kk < dims->N)
            {
                // write column from A
                for (ii = 0; ii < dims->nx[kk + 1]; ii++)
                {
                    mem->A_i[nn++] = ii + row_offset_dyn;
                }
            }

            // write column from C
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                mem->A_i[nn++] = ii + con_start + row_offset_con;
            }

            // write bound on x
            for (ii = 0; ii < dims->nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj + dims->nu[kk])
                {
                    mem->A_i[nn++] = ii + bnd_start + row_offset_bnd;
                    break;
                }
            }
        }

        row_offset_bnd += dims->nb[kk];
        row_offset_con += dims->ng[kk];
        row_offset_dyn += kk < dims->N ? dims->nx[kk + 1] : 0;
    }

    mem->A_p[col] = nn;
}



static void update_constraints_matrix_data(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    c_int ii, jj, kk, nn = 0;
    ocp_qp_dims *dims = in->dim;

    // Traverse matrix in column-major order
    for (kk = 0; kk <= dims->N; kk++)
    {
        int nbu = 0;

        for (jj = 0; jj < dims->nu[kk]; jj++)
        {
            if (kk < dims->N)
            {
                // write column from B
                for (ii = 0; ii < dims->nx[kk + 1]; ii++)
                {
                    mem->A_x[nn++] = BLASFEO_DMATEL(&in->BAbt[kk], jj, ii);
                }
            }

            // write column from D
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                mem->A_x[nn++] = BLASFEO_DMATEL(&in->DCt[kk], jj, ii);
            }

            // write bound on u
            for (ii = 0; ii < dims->nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj)
                {
                    mem->A_x[nn++] = 1.0;
                    nbu++;
                    break;
                }
            }
        }

        for (jj = 0; jj < dims->nx[kk]; jj++)
        {
            if (kk > 0)
            {
                // write column from -I
                mem->A_x[nn++] = -1.0;
            }

            if (kk < dims->N)
            {
                // write column from A
                for (ii = 0; ii < dims->nx[kk + 1]; ii++)
                {
                    mem->A_x[nn++] = BLASFEO_DMATEL(&in->BAbt[kk], jj + dims->nu[kk], ii);
                }
            }

            // write column from C
            for (ii = 0; ii < dims->ng[kk]; ii++)
            {
                mem->A_x[nn++] = BLASFEO_DMATEL(&in->DCt[kk], jj + dims->nu[kk], ii);
            }

            // write bound on x
            for (ii = 0; ii < dims->nb[kk]; ii++)
            {
                if (in->idxb[kk][ii] == jj + dims->nu[kk])
                {
                    mem->A_x[nn++] = 1.0;
                }
            }
        }
    }
}



static void update_bounds(const ocp_qp_in *in, ocp_qp_osqp_memory *mem)
{
    int ii, kk, nn = 0;
    ocp_qp_dims *dims = in->dim;

    // write -b to l and u
    for (kk = 0; kk < dims->N; kk++)
    {
        // unpack b to l
        blasfeo_unpack_dvec(dims->nx[kk + 1], in->b + kk, 0, &mem->l[nn]);

        // change sign of l (to get -b) and copy to u
        for (ii = 0; ii < dims->nx[kk + 1]; ii++)
        {
            mem->l[nn + ii] = -mem->l[nn + ii];
            mem->u[nn + ii] = mem->l[nn + ii];
        }

        nn += dims->nx[kk + 1];
    }

    // write lg and ug
    for (kk = 0; kk <= dims->N; kk++)
    {
        // unpack lg to l
        blasfeo_unpack_dvec(dims->ng[kk], in->d + kk, dims->nb[kk], &mem->l[nn]);

        // unpack ug to u and flip signs because in HPIPM the signs are flipped for upper bounds
        for (ii = 0; ii < dims->ng[kk]; ii++)
        {
            mem->u[nn + ii] = -BLASFEO_DVECEL(&in->d[kk], ii + 2 * dims->nb[kk] + dims->ng[kk]);
        }

        nn += dims->ng[kk];
    }

    // write lb and ub
    for (kk = 0; kk <= dims->N; kk++)
    {
        // unpack lb to l
        blasfeo_unpack_dvec(dims->nb[kk], in->d + kk, 0, &mem->l[nn]);

        // unpack ub to u and flip signs because in HPIPM the signs are flipped for upper bounds
        for (ii = 0; ii < dims->nb[kk]; ii++)
        {
            mem->u[nn + ii] = -BLASFEO_DVECEL(&in->d[kk], ii + dims->nb[kk] + dims->ng[kk]);
        }

        nn += dims->nb[kk];
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

int ocp_qp_osqp_opts_calculate_size(void *config_, void *dims_)
{
    int size = 0;
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

    opts->verbose = 0;  // default value, disable printing
    opts->polish = 1; // default value, enable polishing

    assert((char *) raw_memory + ocp_qp_osqp_opts_calculate_size(config_, dims_) == c_ptr);

    return (void *) opts;
}



void ocp_qp_osqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_osqp_opts *opts = opts_;

    osqp_set_default_settings(opts->osqp_opts);
    opts->osqp_opts->verbose = opts->verbose;
    opts->osqp_opts->polish = opts->polish;

    return;
}



void ocp_qp_osqp_opts_update(void *config_, void *dims_, void *opts_)
{
    // ocp_qp_osqp_opts *opts = (ocp_qp_osqp_opts *)opts_;

    return;
}

/************************************************
 * memory
 ************************************************/

static int osqp_workspace_calculate_size(int n, int m, int P_nnzmax, int A_nnzmax)
{
    int size = 0;

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

int ocp_qp_osqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;

    int n = acados_osqp_num_vars(dims);
    int m = acados_osqp_num_constr(dims);

    int P_nnzmax = acados_osqp_nnzmax_P(dims);
    int A_nnzmax = acados_osqp_nnzmax_A(dims);

    int size = 0;
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
    work->linsys_solver = init_linsys_solver(work->data->P, work->data->A, work->settings->sigma,
                                             work->rho_vec, work->settings->linsys_solver, 0);

    if (!work->linsys_solver) return 0;

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

/************************************************
 * workspace
 ************************************************/

int ocp_qp_osqp_workspace_calculate_size(void *config_, void *dims_, void *opts_) { return 0; }

/************************************************
 * functions
 ************************************************/

static void fill_in_qp_out(const ocp_qp_in *in, ocp_qp_out *out, ocp_qp_osqp_memory *mem)
{
    int ii, kk, nn = 0, mm, con_start = 0, bnd_start = 0;
    ocp_qp_dims *dims = in->dim;
    OSQPSolution *sol = mem->osqp_work->solution;

    for (kk = 0; kk <= dims->N; kk++)
    {
        blasfeo_pack_dvec(dims->nx[kk] + dims->nu[kk], &sol->x[nn], out->ux + kk, 0);
        nn += dims->nx[kk] + dims->nu[kk];

        con_start += kk < dims->N ? dims->nx[kk + 1] : 0;
        bnd_start += dims->ng[kk];
    }

    bnd_start += con_start;

    nn = 0;
    for (kk = 0; kk < dims->N; kk++)
    {
        blasfeo_pack_dvec(dims->nx[kk + 1], &sol->y[nn], out->pi + kk, 0);
        nn += dims->nx[kk + 1];
    }

    nn = 0;
    mm = 0;
    for (kk = 0; kk <= dims->N; kk++)
    {
        for (ii = 0; ii < 2 * dims->nb[kk] + 2 * dims->ng[kk] + 2 * dims->ns[kk]; ii++)
            out->lam[kk].pa[ii] = 0.0;

        for (ii = 0; ii < dims->nb[kk]; ii++)
        {
            double lam = sol->y[bnd_start + nn + ii];
            if (lam <= 0)
                out->lam[kk].pa[ii] = -lam;
            else
                out->lam[kk].pa[dims->nb[kk] + dims->ng[kk] + ii] = lam;
        }

        nn += dims->nb[kk];

        for (ii = 0; ii < dims->ng[kk]; ii++)
        {
            double lam = sol->y[con_start + mm + ii];
            if (lam <= 0)
                out->lam[kk].pa[dims->nb[kk] + ii] = -lam;
            else
                out->lam[kk].pa[2 * dims->nb[kk] + dims->ng[kk] + ii] = lam;
        }

        mm += dims->ng[kk];
    }
}



int ocp_qp_osqp(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;

    int N = qp_in->dim->N;
    int *ns = qp_in->dim->ns;

    // print_ocp_qp_dims(qp_in->dim);

    for (int ii = 0; ii <= N; ii++)
    {
        if (ns[ii] > 0)
        {
            printf("\nOSQP interface can not handle ns>0 yet: what about implementing it? :)\n");
            return ACADOS_FAILURE;
        }
    }

    // print_ocp_qp_in(qp_in);

    ocp_qp_info *info = (ocp_qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;

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
    }
    else
    {
        // mem->osqp_work = osqp_setup(mem->osqp_data, opts->osqp_opts);
        osqp_init_data(mem->osqp_data, opts->osqp_opts, mem->osqp_work);
        mem->first_run = 0;
    }

    // solve OSQP
    osqp_solve(mem->osqp_work);
    fill_in_qp_out(qp_in, qp_out, mem);
    ocp_qp_compute_t(qp_in, qp_out);

    info->solve_QP_time = acados_toc(&qp_timer);
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = mem->osqp_work->info->iter;
    info->t_computed = 1;

    c_int osqp_status = mem->osqp_work->info->status_val;
    int acados_status = osqp_status;

    // check exit conditions
    if (osqp_status == OSQP_SOLVED) acados_status = ACADOS_SUCCESS;
    if (osqp_status == OSQP_MAX_ITER_REACHED) acados_status = ACADOS_MAXITER;
    return acados_status;
}



void ocp_qp_osqp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = &ocp_qp_osqp_opts_calculate_size;
    config->opts_assign = &ocp_qp_osqp_opts_assign;
    config->opts_initialize_default = &ocp_qp_osqp_opts_initialize_default;
    config->opts_update = &ocp_qp_osqp_opts_update;
    config->memory_calculate_size = &ocp_qp_osqp_memory_calculate_size;
    config->memory_assign = &ocp_qp_osqp_memory_assign;
    config->workspace_calculate_size = &ocp_qp_osqp_workspace_calculate_size;
    config->evaluate = &ocp_qp_osqp;

    return;
}
