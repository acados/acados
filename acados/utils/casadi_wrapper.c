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

// TODO(nielsvd): only perform assert in debug mode?
#include <assert.h>

#include <stdlib.h>

casadi_wrapper_args *casadi_wrapper_create_arguments() {
    
    casadi_wrapper_args *args = (casadi_wrapper_args *)malloc(sizeof(casadi_wrapper_args));
    args->nx = 0;
    args->nu = 0;
    args->np = 0;

    args->fun = NULL;
    args->dims = NULL;

    return args;
}

int_t casadi_wrapper_calculate_workspace_size(const casadi_wrapper_in *cw_in, void *args_) {
    casadi_wrapper_args *args = (casadi_wrapper_args *) args_;
    
    int_t size = sizeof(casadi_wrapper_workspace);

    int_t sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

    size += 3*sizeof(real_t *);
    size += 3*sizeof(real_t *);
    size += sz_iw*sizeof(int_t);
    size += sz_w*sizeof(real_t);

    assert(sz_arg <= 3);
    // nielsvd: apparently it is possible casadi adds additional real_t-pointers in the
    //          res-array, which is used in the body of the code. This wrapper is not 
    //          compatible with this case. Therefore, check whether #outputs <= 3.
    assert(sz_res <= 3);

    return size;
}

char *casadi_wrapper_assign_workspace(const casadi_wrapper_in *cw_in, void *args_, void **work_, void *raw_memory) {
    casadi_wrapper_args *args = (casadi_wrapper_args *)args_;

    int_t sz_arg, sz_res, sz_iw, sz_w;
    args->dims(&sz_arg, &sz_res, &sz_iw, &sz_w);

    casadi_wrapper_workspace **work = (casadi_wrapper_workspace **)work_;
    char *c_ptr = (char *)raw_memory;

    *work = (casadi_wrapper_workspace *)c_ptr;
    c_ptr += sizeof(casadi_wrapper_workspace);

    (*work)->arg = (const real_t **) c_ptr;
    c_ptr += sz_arg*sizeof(real_t *);

    (*work)->res = (real_t **) c_ptr;
    c_ptr += sz_res*sizeof(real_t *);

    (*work)->iw = (int_t *) c_ptr;
    c_ptr += sz_iw*sizeof(int_t);

    (*work)->w = (real_t *) c_ptr;
    c_ptr += sz_w*sizeof(real_t);

    return c_ptr;
}

casadi_wrapper_workspace *casadi_wrapper_create_workspace(const casadi_wrapper_in *cw_in, void *args_) {
    casadi_wrapper_workspace *work;

    int_t workspace_size = casadi_wrapper_calculate_workspace_size(cw_in, args_);
    void *raw_memory_ptr = malloc(workspace_size);

    char *ptr_end = casadi_wrapper_assign_workspace(cw_in, args_, (void **)&work, raw_memory_ptr);
    assert((char *)raw_memory_ptr + workspace_size >= ptr_end); (void)ptr_end;

    return work;
}

int_t casadi_wrapper(const casadi_wrapper_in *cw_in, casadi_wrapper_out *cw_out, void *args_, void *workspace_) {
    casadi_wrapper_args *args = (casadi_wrapper_args *)args_;
    casadi_wrapper_workspace *work = (casadi_wrapper_workspace *)workspace_;

    work->arg[0] = cw_in->x;
    work->arg[1] = cw_in->u;
    work->arg[2] = cw_in->p;

    work->res[0] = cw_out->y;
    work->res[1] = cw_in->compute_jac ? cw_out->jac_y : NULL;
    work->res[2] = cw_in->compute_hess ? cw_out->hess_y: NULL;

    int_t output = args->fun(work->arg, work->res, work->iw, work->w, 0);

    return output;
}

void casadi_wrapper_initialize(const casadi_wrapper_in *cw_in, void *args_, void **work) {
    *work = casadi_wrapper_create_workspace(cw_in, args_);
}

void casadi_wrapper_destroy(void *work) {
    // TODO(nielsvd): replace dummy commands once interface completed
    (void)work;
}