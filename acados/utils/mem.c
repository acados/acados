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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
// acados
#include "acados/utils/mem.h"

// #define _USE_VALGRIND_  // uncomment to bypass assignment and do new memory allocation instead

#define _USE_MALLOC_  // acados_malloc = malloc / acados_malloc = calloc

void make_int_multiple_of(int num, int *size) {
    *size = (*size + num - 1) / num * num;
}



int align_char_to(int num, char **c_ptr) {
    size_t s_ptr = (size_t)*c_ptr;
    s_ptr = (s_ptr + num - 1) / num * num;
    int offset = num - (int)(s_ptr - (size_t)(*c_ptr));
    *c_ptr = (char *)s_ptr;
    return offset;
}



// print warning when by-passing pointer and allocating new memory (for debugging)
static void print_warning ()
{
    printf(" -- using dynamically allocated memory for debugging --\n");
}



void *acados_malloc(size_t nitems, size_t size)
{
#if defined(_USE_MALLOC_)
    void *ptr = malloc(nitems*size);
#else
    void *ptr = calloc(nitems, size);
#endif
    return ptr;
}



void assign_double_ptrs(int n, double ***v, char **ptr)
{
    assert((size_t)*ptr % 8 == 0 && "pointer not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    *v = (double **)acados_malloc(n, sizeof(double *));
#else
    *v = (double **)*ptr;
    *ptr += sizeof(double *) * n;
#endif
}



void assign_int_ptrs(int n, int ***v, char **ptr)
{
    assert((size_t)*ptr % 8 == 0 && "pointer not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    *v = (int **)acados_malloc(n, sizeof(int *));
#else
    *v = (int **)*ptr;
    *ptr += sizeof(int *) * n;
#endif
}



void assign_strvec_ptrs(int n, struct d_strvec **sv, char **ptr)
{
    assert((size_t)*ptr % 8 == 0 && "pointer not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    *sv = (struct d_strvec *)acados_malloc(n, sizeof(struct d_strvec));
#else
    *sv = (struct d_strvec *) *ptr;
    *ptr += sizeof(struct d_strvec) * n;
#endif
}



void assign_strmat_ptrs(int n, struct d_strmat **sm, char **ptr)
{
    assert((size_t)*ptr % 8 == 0 && "pointer not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    *sm = (struct d_strmat *)acados_malloc(n, sizeof(struct d_strmat));
#else
    *sm = (struct d_strmat *) *ptr;
    *ptr += sizeof(struct d_strmat) * n;
#endif
}


void assign_strmat_ptrs_to_ptrs(int n, struct d_strmat ***sm, char **ptr)
{
    assert((size_t)*ptr % 8 == 0 && "pointer not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    *sm = (struct d_strmat **) acados_malloc(n, sizeof(struct d_strmat *));
#else
    *sm = (struct d_strmat **) *ptr;
    *ptr += sizeof(struct d_strmat *) * n;
#endif
}



void assign_int(int n, int **v, char **ptr)
{
#ifdef _USE_VALGRIND_
    *v = (int *)acados_malloc(n, sizeof(int));
    print_warning();
#else
    *v = (int *)*ptr;
    *ptr += sizeof(int) * n;
#endif
}



void assign_double(int n, double **v, char **ptr)
{
    assert((size_t)*ptr % 8 == 0 && "double not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    *v = (double *)acados_malloc(n, sizeof(double));
    print_warning();
#else
    *v = (double *)*ptr;
    *ptr += sizeof(double) * n;
#endif
}



void assign_strvec(int n, struct d_strvec *sv, char **ptr)
{
    assert((size_t)*ptr % 8 == 0 && "strvec not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    d_allocate_strvec(n, sv);
    print_warning();
#else
    d_create_strvec(n, sv, *ptr);
    *ptr += sv->memory_size;
#endif
}



void assign_strmat(int m, int n, struct d_strmat *sA, char **ptr)
{
#ifdef LA_HIGH_PERFORMANCE
    assert((size_t)*ptr % 64 == 0 && "strmat not 64-byte aligned!");
#else
    assert((size_t)*ptr % 8 == 0 && "strmat not 8-byte aligned!");
#endif

#ifdef _USE_VALGRIND_
    d_allocate_strmat(m, n, sA);
    print_warning();
#else
    d_create_strmat(m, n, sA, *ptr);
    *ptr += sA->memory_size;
#endif
}



// TODO(dimitris): probably does not belong here
void copy_module_pointers_to_args(void *solver_in_args_, void *solver_)
{
    module_solver *solver_in_args = solver_in_args_;
    module_solver *solver = solver_;

    solver_in_args->calculate_args_size = solver->calculate_args_size;
    solver_in_args->assign_args = solver->assign_args;
    solver_in_args->initialize_default_args = solver->initialize_default_args;
    solver_in_args->calculate_memory_size = solver->calculate_memory_size;
    solver_in_args->assign_memory = solver->assign_memory;
    solver_in_args->calculate_workspace_size = solver->calculate_workspace_size;
    solver_in_args->fun = solver->fun;
}
