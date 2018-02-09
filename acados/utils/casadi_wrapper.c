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

#include <acados/utils/mem.h>



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



int casadi_wrapper_calculate_args_size(external_function_dims *dims, void *submodules_)
{
    int size = sizeof(casadi_wrapper_args);

    return size;
}



void *casadi_wrapper_assign_args(external_function_dims *dims, void **submodules_, void *raw_memory)
{
    casadi_wrapper_args *args;

    char *c_ptr = (char *) raw_memory;

    args = (casadi_wrapper_args *) c_ptr;
    c_ptr += sizeof(casadi_wrapper_args);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + casadi_wrapper_calculate_args_size(dims, *submodules_) == c_ptr);

    // Update submodules pointer
    *submodules_ = NULL;

    return (void *)args;
}



void *casadi_wrapper_copy_args(external_function_dims *dims, void *raw_memory, void *source_)
{
    casadi_wrapper_args *source = (casadi_wrapper_args *)source_;
    casadi_wrapper_args *dest;

    void *submodules = NULL;
    dest = casadi_wrapper_assign_args(dims, &submodules, raw_memory);

    dest->fun = source->fun;

    dest->dims = source->dims;

    dest->sparsity = source->sparsity;
    
    return (void *)dest;
}



void casadi_wrapper_initialize_default_args(void *args_)
{
    casadi_wrapper_args *args = (casadi_wrapper_args *)args_;

    args->fun = NULL;
    args->dims = NULL;
    args->sparsity = NULL;
}



int casadi_wrapper_calculate_memory_size(external_function_dims *dims, void *args_)
{
    int size = sizeof(casadi_wrapper_memory);

    size += 1*8;

    return size;
}



void *casadi_wrapper_assign_memory(external_function_dims *dims, void *args_, void *raw_memory)
{
    casadi_wrapper_args *args = (casadi_wrapper_args *)args_;

    casadi_wrapper_memory *mem;

    char *c_ptr = (char *) raw_memory;

    mem = (casadi_wrapper_memory *) c_ptr;
    c_ptr += sizeof(casadi_wrapper_memory);

    align_char_to(8, &c_ptr);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + casadi_wrapper_calculate_memory_size(dims, args) >= c_ptr);

    return (void *)mem;
}



int casadi_wrapper_calculate_workspace_size(external_function_dims *dims, void *args_)
{
    casadi_wrapper_args *args = (casadi_wrapper_args *)args_;
    int size = sizeof(casadi_wrapper_workspace);

    int sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

    size += sz_arg * sizeof(double *);
    size += sz_res * sizeof(double *);
    size += sz_iw * sizeof(int);
    size += sz_w * sizeof(double);

    // Workspace for temporary sparse matrix storage
    size += sz_res * sizeof(double *);
    for (int i=0; i<sz_res; i++) {
        size += nnz_output(args->sparsity(i)) * sizeof(double);
    }

    return size;
}



static void cast_workspace(casadi_wrapper_args *args, casadi_wrapper_memory *mem, casadi_wrapper_workspace *work)
{
    int sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

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
    c_ptr += sz_res * sizeof(double *);

    for (int i=0; i<sz_res; i++) {
        work->sparse_res[i] = (double *)c_ptr;
        c_ptr += nnz_output(args->sparsity(i)) * sizeof(double);
    }
}



int casadi_wrapper(external_function_in *ef_in, external_function_out *ef_out, void *args_, void *mem_, void *work_)
{
    casadi_wrapper_args *args = (casadi_wrapper_args *) args_;
    casadi_wrapper_memory *mem = (casadi_wrapper_memory *) mem_;
    casadi_wrapper_workspace *work = (casadi_wrapper_workspace *) work_;

    cast_workspace(args, mem, work);

    int sz_arg, sz_res;
    args->dims(&sz_arg, &sz_res, NULL, NULL);

    for (int i=0; i<sz_arg; i++) {
        work->arg[i] = ef_in->inputs[i];
    }

    for (int i=0; i<sz_res; i++) {
        work->res[i] = ef_in->compute_output[i] ? work->sparse_res[i] : NULL;
    }

    int output = args->fun(work->arg, work->res, work->iw, work->w, 0);

    // Densify
    for (int i=0; i<sz_res; i++) {
        if (ef_in->compute_output[i])
            densify(work->sparse_res[i], ef_out->outputs[i], args->sparsity(i));
    }

    return output;
}
