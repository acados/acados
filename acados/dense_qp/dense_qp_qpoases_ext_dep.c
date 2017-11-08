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

// external
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
// acados
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/dense_qp/dense_qp_qpoases_ext_dep.h"
#include "acados/dense_qp/dense_qp_common.h"



dense_qp_qpoases_args *dense_qp_qpoases_create_arguments(dense_qp_in *qp_in) {

    int memory_size = dense_qp_qpoases_calculate_args_size(qp_in);
    void *ptr = malloc(memory_size);
    dense_qp_qpoases_args *args = dense_qp_qpoases_assign_args(qp_in, ptr);
    dense_qp_qpoases_initialize_default_args(args);

    return args;
}



dense_qp_qpoases_memory *dense_qp_qpoases_create_memory(dense_qp_in *qp_in, void *args_) {
    dense_qp_qpoases_args *args = (dense_qp_qpoases_args *) args_;

    dense_qp_qpoases_memory *mem;

    int memory_size = dense_qp_qpoases_calculate_memory_size(qp_in, args);
    void *ptr = malloc(memory_size);
    char *ptr_end = dense_qp_qpoases_assign_memory(qp_in, args, (void **) &mem, ptr);
    assert((char *) ptr + memory_size >= ptr_end); (void) ptr_end;
    // TODO(dimitris): ASSERT DOES NOT WORK (REDEFINED IN QPOASES)
    assert(1==0);
    return mem;
}
