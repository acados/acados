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

int_t nnz_output(const int_t *sparsity) {
    int_t nnz = 0;
    if (sparsity != NULL) {
        const int_t nrow = sparsity[0];
        const int_t ncol = sparsity[1];
        const int_t dense = sparsity[2];
        if (dense) {
            nnz = nrow * ncol;
        } else {
            const int_t *colind = sparsity + 2;
            for (int_t i = 0; i < ncol; ++i) {
                nnz += colind[i + 1] - colind[i];
            }
        }
    }

    return nnz;
}

void densify(const real_t *sparse_in, real_t *dense_out,
             const int_t *sparsity) {
    const int_t nrow = sparsity[0];
    const int_t ncol = sparsity[1];
    const int_t dense = sparsity[2];

    if (dense) {
        for (int_t i = 0; i < ncol * nrow; i++) dense_out[i] = sparse_in[i];
    } else {
        // Fill with zeros
        for (int_t i = 0; i < ncol; i++) {
            for (int_t j = 0; j < nrow; j++) dense_out[i * nrow + j] = 0;
        }
        // Additional data
        const real_t *x = sparse_in;
        const int_t *colind = sparsity + 2, *row = sparsity + ncol + 3;
        // Copy nonzeros
        for (int_t i = 0; i < ncol; ++i) {
            for (int_t el = colind[i]; el != colind[i + 1]; ++el) {
                dense_out[row[el] + i * nrow] = *x++;
            }
        }
    }
}

casadi_wrapper_args *casadi_wrapper_create_arguments() {
    casadi_wrapper_args *args =
        (casadi_wrapper_args *)malloc(sizeof(casadi_wrapper_args));
    args->fun = NULL;
    args->dims = NULL;

    return args;
}

int_t casadi_wrapper_calculate_workspace_size(const casadi_wrapper_in *cw_in,
                                              casadi_wrapper_args *args) {
    int_t size = sizeof(casadi_wrapper_workspace);

    int_t sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

    // Allocate for three inputs and outputs
    if (sz_arg < 3) sz_arg = 3;
    if (sz_res < 3) sz_res = 3;

    // Fixed input dimension: x, u, p
    size += sz_arg * sizeof(real_t *);
    // Fixed output dimension: f, df/d(x,u), (d^2 f)/d(x,u)^2
    size += sz_res * sizeof(real_t *);
    size += sz_iw * sizeof(int_t);
    size += sz_w * sizeof(real_t);

    // Workspace for temporary sparse matrix storage
    size += 3 * sizeof(real_t *);
    size += nnz_output(args->sparsity(0)) * sizeof(real_t);
    size += nnz_output(args->sparsity(1)) * sizeof(real_t);
    size += nnz_output(args->sparsity(2)) * sizeof(real_t);

    return size;
}

char *casadi_wrapper_assign_workspace(const casadi_wrapper_in *cw_in,
                                      casadi_wrapper_args *args,
                                      casadi_wrapper_workspace **work,
                                      void *raw_memory) {
    int_t sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

    // Allocate for three inputs and outputs
    if (sz_arg < 3) sz_arg = 3;
    if (sz_res < 3) sz_res = 3;

    char *c_ptr = (char *)raw_memory;

    *work = (casadi_wrapper_workspace *)c_ptr;
    c_ptr += sizeof(casadi_wrapper_workspace);

    (*work)->arg = (const real_t **)c_ptr;
    c_ptr += sz_arg * sizeof(real_t *);

    (*work)->res = (real_t **)c_ptr;
    c_ptr += sz_res * sizeof(real_t *);

    (*work)->iw = (int_t *)c_ptr;
    c_ptr += sz_iw * sizeof(int_t);

    (*work)->w = (real_t *)c_ptr;
    c_ptr += sz_w * sizeof(real_t);

    (*work)->sparse_res = (real_t **)c_ptr;
    c_ptr += 3 * sizeof(real_t *);

    (*work)->sparse_res[0] = (real_t *)c_ptr;
    c_ptr += nnz_output(args->sparsity(0)) * sizeof(real_t);

    (*work)->sparse_res[1] = (real_t *)c_ptr;
    c_ptr += nnz_output(args->sparsity(1)) * sizeof(real_t);

    (*work)->sparse_res[2] = (real_t *)c_ptr;
    c_ptr += nnz_output(args->sparsity(2)) * sizeof(real_t);

    return c_ptr;
}

casadi_wrapper_workspace *casadi_wrapper_create_workspace(
    const casadi_wrapper_in *cw_in, casadi_wrapper_args *args) {
    casadi_wrapper_workspace *work;

    int_t workspace_size = casadi_wrapper_calculate_workspace_size(cw_in, args);
    void *raw_memory_ptr = malloc(workspace_size);

    char *ptr_end =
        casadi_wrapper_assign_workspace(cw_in, args, &work, raw_memory_ptr);
    assert((char *)raw_memory_ptr + workspace_size >= ptr_end);
    (void)ptr_end;

    return work;
}

int_t casadi_wrapper(const casadi_wrapper_in *cw_in, casadi_wrapper_out *cw_out,
                     casadi_wrapper_args *args,
                     casadi_wrapper_workspace *work) {
    work->arg[0] = cw_in->x;
    work->arg[1] = cw_in->u;
    work->arg[2] = cw_in->p;

    work->res[0] = work->sparse_res[0];
    work->res[1] = cw_in->compute_jac ? work->sparse_res[1] : NULL;
    work->res[2] = cw_in->compute_hess ? work->sparse_res[2] : NULL;

    int_t output = args->fun(work->arg, work->res, work->iw, work->w, 0);

    // Densify
    densify(work->sparse_res[0], cw_out->y, args->sparsity(0));
    if (cw_in->compute_jac)
        densify(work->sparse_res[1], cw_out->jac_y, args->sparsity(1));
    if (cw_in->compute_hess)
        densify(work->sparse_res[2], cw_out->hess_y, args->sparsity(2));

    return output;
}

void casadi_wrapper_initialize(const casadi_wrapper_in *cw_in,
                               casadi_wrapper_args *args,
                               casadi_wrapper_workspace **work) {
    *work = casadi_wrapper_create_workspace(cw_in, args);
}

void casadi_wrapper_destroy(casadi_wrapper_workspace *work) { free(work); }
