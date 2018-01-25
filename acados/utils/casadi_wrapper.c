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

#include "acados/utils/casadi_wrapper.h"

#include <assert.h>
#include <stdlib.h>



int nnz_output(const int *sparsity)
{
    int nnz = 0;
    if (sparsity != NULL) {
        const int nrow = sparsity[0];
        const int ncol = sparsity[1];
        const int dense = sparsity[2];
        if (dense) {
            nnz = nrow * ncol;
        } else {
            const int *colind = sparsity + 2;
            for (int i = 0; i < ncol; ++i) {
                nnz += colind[i + 1] - colind[i];
            }
        }
    }

    return nnz;
}



void densify(const double *sparse_in, double *dense_out, const int *sparsity)
{
    const int nrow = sparsity[0];
    const int ncol = sparsity[1];
    const int dense = sparsity[2];

    if (dense) {
        for (int i = 0; i < ncol * nrow; i++) dense_out[i] = sparse_in[i];
    } else {
        // Fill with zeros
        for (int i = 0; i < ncol; i++) {
            for (int j = 0; j < nrow; j++) dense_out[i * nrow + j] = 0;
        }
        // Additional data
        const double *x = sparse_in;
        const int *colind = sparsity + 2, *row = sparsity + ncol + 3;
        // Copy nonzeros
        for (int i = 0; i < ncol; ++i) {
            for (int el = colind[i]; el != colind[i + 1]; ++el) {
                dense_out[row[el] + i * nrow] = *x++;
            }
        }
    }
}



int casadi_wrapper_calculate_args_size(casadi_wrapper_dims *dims)
{
    int size = sizeof(casadi_wrapper_args);

    return size;
}



void *casadi_wrapper_assign_args(casadi_wrapper_dims *dims, void *raw_memory)
{
    casadi_wrapper_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (casadi_wrapper_args *) c_ptr;
    c_ptr += sizeof(casadi_wrapper_args);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + casadi_wrapper_calculate_args_size(dims) == c_ptr);

    return (void *)args;
}



void casadi_wrapper_initialize_default_args(casadi_wrapper_args *args)
{
    args->fun = NULL;
    args->dims = NULL;
    args->sparsity = NULL;
}



int casadi_wrapper_calculate_memory_size(casadi_wrapper_dims *dims, casadi_wrapper_args *args)
{
    int size = sizeof(casadi_wrapper_memory);

    return size;
}



void *casadi_wrapper_assign_memory(casadi_wrapper_dims *dims, casadi_wrapper_args *args, void *raw_memory)
{
    casadi_wrapper_memory *mem;

    char *c_ptr = (char *) raw_memory;

    args = (casadi_wrapper_memory *) c_ptr;
    c_ptr += sizeof(casadi_wrapper_memory);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + casadi_wrapper_calculate_memory_size(dims, args) == c_ptr);

    return (void *)mem;
}



int casadi_wrapper_calculate_workspace_size(casadi_wrapper_dims *dims, casadi_wrapper_args *args)
{
    int size = sizeof(casadi_wrapper_workspace);

    int sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

    // Allocate for three inputs and outputs
    if (sz_arg < 5) sz_arg = 5;
    if (sz_res < 4) sz_res = 4;

    // Fixed input dimension: x, u, xdot, p, mul
    size += sz_arg * sizeof(double *);
    // Fixed output dimension: f, df/d(x,u,z,xdot), (d^2 f)/d(x,u)^2
    size += sz_res * sizeof(double *);
    size += sz_iw * sizeof(int);
    size += sz_w * sizeof(double);

    // Workspace for temporary sparse matrix storage
    size += 4 * sizeof(double *);
    for (int i=0; i<4; i++) {
        size += nnz_output(args->sparsity(i)) * sizeof(double);
    }

    return size;
}



static void cast_workspace(casadi_wrapper_dims *dims, casadi_wrapper_args *args, casadi_wrapper_memory *mem, casadi_wrapper_workspace *work)
{
    int sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

    // Allocate for three inputs and outputs
    if (sz_arg < 5) sz_arg = 5;
    if (sz_res < 4) sz_res = 4;

    char *c_ptr = (char *) work;
    c_ptr += sizeof(casadi_wrapper_workspace);

    work->arg = (const double **)c_ptr;
    c_ptr += sz_arg * sizeof(double *);

    work->res = (double **)c_ptr;
    c_ptr += sz_res * sizeof(double *);

    work->iw = (int *)c_ptr;
    c_ptr += sz_iw * sizeof(int);

    work->w = (double *)c_ptr;
    c_ptr += sz_w * sizeof(double);

    work->sparse_res = (double **)c_ptr;
    c_ptr += 4 * sizeof(double *);

    for (int i=0; i<4; i++) {
        work->sparse_res[i] = (double *)c_ptr;
        c_ptr += nnz_output(args->sparsity(i)) * sizeof(double);
    }
}



int casadi_wrapper(casadi_wrapper_in *cw_in, casadi_wrapper_out *cw_out, casadi_wrapper_args *args, casadi_wrapper_memory *mem, casadi_wrapper_workspace *work)
{
    work->arg[0] = cw_in->x;
    work->arg[1] = cw_in->u;
    work->arg[2] = cw_in->z;
    work->arg[3] = cw_in->xdot;
    work->arg[4] = cw_in->p;
    work->arg[5] = cw_in->mul;

    work->res[0] = cw_in->compute_y ? work->sparse_res[0] : NULL;
    work->res[1] = cw_in->compute_jac_y ? work->sparse_res[1] : NULL;
    work->res[2] = cw_in->compute_grad_mul_y ? work->sparse_res[2] : NULL;
    work->res[3] = cw_in->compute_hess_mul_y ? work->sparse_res[3] : NULL;

    int output = args->fun(work->arg, work->res, work->iw, work->w, 0);

    // Densify
    if (cw_in->compute_y)
        densify(work->sparse_res[0], cw_out->y, args->sparsity(0));
    if (cw_in->compute_jac_y)
        densify(work->sparse_res[1], cw_out->jac_y, args->sparsity(1));
    if (cw_in->compute_grad_mul_y)
        densify(work->sparse_res[2], cw_out->grad_mul_y, args->sparsity(2));
    if (cw_in->compute_hess_mul_y)
        densify(work->sparse_res[3], cw_out->hess_mul_y, args->sparsity(3));

    return output;
}
