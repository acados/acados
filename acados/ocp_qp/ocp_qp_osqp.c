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
#include "osqp/include/types.h"
#include "osqp/include/osqp.h"
#include "osqp/include/util.h"
#include "osqp/include/constants.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_osqp.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/print.h"


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



static void print_inputs(ocp_qp_osqp_memory *mem)
{
    printf("\n----------> OSQP INPUTS <----------\n\n");
    printf("NUMBER OF VARIABLES: %d\n", mem->osqp_data->n);
    printf("NUMBER OF CONSTRAINTS: %d\n", mem->osqp_data->m);
    printf("NUMBER OF NON-ZEROS in HESSIAN: %d\n", mem->P_nnzmax);
    printf("NUMBER OF NON-ZEROS in CONSTRAINTS: %d\n", mem->A_nnzmax);
    printf("\n-----------------------------------\n\n");

    int ii;

    printf("\nOBJECTIVE FUNCTION:\n");
    // for (ii = 0; ii < mem->P_nnzmax; ii++)
    //     printf("=====> P_x[%d] = %f, P_i[%d] = %d\n", ii + 1, mem->P_x[ii], ii+1, mem->P_i[ii]);

    // for (ii = 0; ii < mem->A_nnzmax; ii++)
        // printf("=====> A_x[%d] = %f, A_i[%d] = %d\n", ii + 1, mem->A_x[ii], ii+1, mem->A_i[ii]);
    print_csc_matrix(mem->osqp_data->P, "Matrix P");
    // for (ii = 0; ii < mem->osqp_data->n; ii++)
    //     printf("=====> q[%d] = %f\n", ii + 1, mem->q[ii]);
    // for (ii = 0; ii < mem->osqp_data->n+1; ii++)
    //     printf("=====> P_p[%d] = %d\n", ii + 1, mem->P_p[ii]);

    // print_csc_as_dns(mem->osqp_data->P);

    printf("\nBOUNDS:\n");
    for (ii = 0; ii < mem->osqp_data->m; ii++)
        printf("=====> l[%d] = %f, u[%d] = %f\n", ii + 1, mem->l[ii], ii + 1, mem->u[ii]);

    printf("\nCONSTRAINTS MATRIX:\n");
    print_csc_matrix(mem->osqp_data->A, "Matrix A");
    // print_csc_as_dns(mem->osqp_data->A);
    // for (ii = 0; ii < mem->A_nnzmax; ii++)
    //     printf("=====> A_x[%d] = %f, A_i[%d] = %d\n", ii + 1, mem->A_x[ii], ii+1, mem->A_i[ii]);
    // for (ii = 0; ii < mem->osqp_data->n+1; ii++)
    //     printf("=====> A_p[%d] = %d\n", ii + 1, mem->A_p[ii]);

}



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
            m += dims->nx[ii+1];
        }
    }

    return m;
}



static int acados_osqp_nnzmax_P(const ocp_qp_dims *dims)
{
    int nnz = 0;

    for (int ii = 0; ii <= dims->N; ii++)
    {
        nnz +=     dims->nx[ii] * dims->nx[ii];  // Q
        nnz +=     dims->nu[ii] * dims->nu[ii];  // R
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
        nnz += dims->nb[ii];  // eye
        nnz += dims->ng[ii]*dims->nx[ii]; // C
        nnz += dims->ng[ii]*dims->nu[ii]; // D

        // equality constraints
        if (ii < dims->N)
        {
            nnz += dims->nx[ii+1] * dims->nx[ii];  // A
            nnz += dims->nx[ii+1] * dims->nu[ii];  // B
            nnz += dims->nx[ii+1];  // eye
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
        blasfeo_unpack_dvec(dims->nu[kk]+dims->nx[kk], in->rqz+kk, 0, &mem->q[nn]);
        nn += dims->nu[kk]+dims->nx[kk];
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
        con_start += kk < dims->N ? dims->nx[kk+1] : 0;
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
                for (ii = 0; ii < dims->nx[kk+1]; ii++)
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
                for (ii = 0; ii < dims->nx[kk+1]; ii++)
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
        row_offset_dyn += kk < dims->N ? dims->nx[kk+1] : 0;
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
                for (ii = 0; ii < dims->nx[kk+1]; ii++)
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
                for (ii = 0; ii < dims->nx[kk+1]; ii++)
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
    int ii, jj, kk, nn = 0;
    ocp_qp_dims *dims = in->dim;

    // write -b to l and u
    for (kk = 0; kk < dims->N; kk++)
    {
        // unpack b to l
        blasfeo_unpack_dvec(dims->nx[kk+1], in->b+kk, 0, &mem->l[nn]);

        // change sign of l (to get -b) and copy to u
        for (ii = 0; ii < dims->nx[kk+1]; ii++)
        {
            mem->l[nn+ii] = -mem->l[nn+ii];
            mem->u[nn+ii] = mem->l[nn+ii];
        }

        nn += dims->nx[kk+1];
    }

    // write lg and ug
    for (kk = 0; kk <= dims->N; kk++)
    {
        // unpack lg to l
        blasfeo_unpack_dvec(dims->ng[kk], in->d+kk, dims->nb[kk], &mem->l[nn]);

        // unpack ug to u and flip signs because in HPIPM the signs are flipped for upper bounds
        for (ii = 0; ii < dims->ng[kk]; ii++)
        {
            mem->u[nn+ii] = -BLASFEO_DVECEL(&in->d[kk], ii+2*dims->nb[kk]+dims->ng[kk]);
        }

        nn += dims->ng[kk];
    }

    // write lb and ub
    for (kk = 0; kk <= dims->N; kk++)
    {
        // unpack lb to l
        blasfeo_unpack_dvec(dims->nb[kk], in->d+kk, 0, &mem->l[nn]);

        // unpack ub to u and flip signs because in HPIPM the signs are flipped for upper bounds
        for (ii = 0; ii < dims->nb[kk]; ii++)
        {
            mem->u[nn+ii] = -BLASFEO_DVECEL(&in->d[kk], ii+dims->nb[kk]+dims->ng[kk]);
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

        OSQPData *data = mem->osqp_data;

        data->n = acados_osqp_num_vars(in->dim);
        data->m = acados_osqp_num_constr(in->dim);

        mem->osqp_data->q = mem->q;
        mem->osqp_data->l = mem->l;
        mem->osqp_data->u = mem->u;

        mem->osqp_data->P = csc_matrix(data->n, data->n, mem->P_nnzmax,
                                       mem->P_x, mem->P_i, mem->P_p);
        mem->osqp_data->A = csc_matrix(data->m, data->n, mem->A_nnzmax,
                                       mem->A_x, mem->A_i, mem->A_p);

        update_bounds(in, mem);
        update_gradient(in, mem);
        update_hessian_data(in, mem);
        update_constraints_matrix_data(in, mem);
    }
    else
    {
        update_bounds(in, mem);
        update_gradient(in, mem);
        update_hessian_data(in, mem);
        update_constraints_matrix_data(in, mem);
    }
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

    opts->verbose = 0; // default value, disable printing

    assert((char *) raw_memory + ocp_qp_osqp_opts_calculate_size(config_, dims_) == c_ptr);

    return (void *) opts;
}



void ocp_qp_osqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_osqp_opts *opts = opts_;

    osqp_set_default_settings(opts->osqp_opts);
    opts->osqp_opts->verbose = opts->verbose;

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

int ocp_qp_osqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;

    int n = acados_osqp_num_vars(dims);
    int m = acados_osqp_num_constr(dims);

    int P_nnzmax = acados_osqp_nnzmax_P(dims);
    int A_nnzmax = acados_osqp_nnzmax_A(dims);

    int size = 0;
    size += sizeof(ocp_qp_osqp_memory);

    size += 1*n*sizeof(c_float);  // q
    size += 2*m*sizeof(c_float);  // l, u

    size += P_nnzmax*sizeof(c_float);  // P_x
    size += P_nnzmax*sizeof(c_int);  // P_i
    size += (n+1)*sizeof(c_int);  // P_p

    size += A_nnzmax*sizeof(c_float);  // A_x
    size += A_nnzmax*sizeof(c_int);  // A_i
    size += (n+1)*sizeof(c_int);  // A_p

    size += sizeof(OSQPData);
    size += sizeof(OSQPWorkspace);

    size += 3 * 8;

    return size;
}



void *ocp_qp_osqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_qp_dims *dims = dims_;
    ocp_qp_osqp_opts *opts = opts_;
    ocp_qp_osqp_memory *mem;

    int n = acados_osqp_num_vars(dims);
    int m = acados_osqp_num_constr(dims);

    // char pointer
    char *c_ptr = (char *) raw_memory;

    mem = (ocp_qp_osqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_osqp_memory);

    mem->P_nnzmax = acados_osqp_nnzmax_P(dims);
    mem->A_nnzmax = acados_osqp_nnzmax_A(dims);
    mem->first_run = 1;

    align_char_to(8, &c_ptr);

    // doubles
    mem->q = (c_float *) c_ptr;
    c_ptr += n*sizeof(c_float);

    mem->l = (c_float *) c_ptr;
    c_ptr += m*sizeof(c_float);

    mem->u = (c_float *) c_ptr;
    c_ptr += m*sizeof(c_float);

    mem->P_x = (c_float *) c_ptr;
    c_ptr += (mem->P_nnzmax)*sizeof(c_float);

    mem->A_x = (c_float *) c_ptr;
    c_ptr += (mem->A_nnzmax)*sizeof(c_float);

    // ints
    mem->P_i = (c_int *) c_ptr;
    c_ptr += (mem->P_nnzmax)*sizeof(c_int);

    mem->P_p = (c_int *) c_ptr;
    c_ptr += (n+1)*sizeof(c_int);

    mem->A_i = (c_int *) c_ptr;
    c_ptr += (mem->A_nnzmax)*sizeof(c_int);

    mem->A_p = (c_int *) c_ptr;
    c_ptr += (n+1)*sizeof(c_int);

    align_char_to(8, &c_ptr);

    mem->osqp_data = (OSQPData *) c_ptr;
    c_ptr += sizeof(OSQPData);

    align_char_to(8, &c_ptr);

    mem->osqp_work = (OSQPWorkspace *) c_ptr;
    c_ptr += sizeof(OSQPWorkspace);

    assert((char *)raw_memory + ocp_qp_osqp_memory_calculate_size(config_, dims, opts_) >= c_ptr);

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
    int ii, jj, kk, nn = 0, mm, con_start = 0, bnd_start = 0;
    ocp_qp_dims *dims = in->dim;
    OSQPSolution *sol = mem->osqp_work->solution;

    for (kk = 0; kk <= dims->N; kk++)
    {
        blasfeo_pack_dvec(dims->nx[kk]+dims->nu[kk], &sol->x[nn], out->ux+kk, 0);
        nn += dims->nx[kk]+dims->nu[kk];

        con_start += kk < dims->N ? dims->nx[kk+1] : 0;
        bnd_start += dims->ng[kk];
    }

    bnd_start += con_start;

    nn = 0;
    for (kk = 0; kk < dims->N; kk++)
    {
        blasfeo_pack_dvec(dims->nx[kk + 1], &sol->y[nn], out->pi+kk, 0);
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

    // print_inputs(mem);

    // update osqp workspace with new data
    if (!mem->first_run)
    {
        osqp_update_lin_cost(mem->osqp_work, mem->q);
        osqp_update_P_A(mem->osqp_work, mem->P_x, NULL, mem->P_nnzmax,
                        mem->A_x, NULL, mem->A_nnzmax);
        osqp_update_bounds(mem->osqp_work, mem->l, mem->u);
    }
    else
    {
        mem->osqp_work = osqp_setup(mem->osqp_data, opts->osqp_opts);
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
