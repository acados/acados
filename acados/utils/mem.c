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
#include <stdio.h>
#include <stdlib.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
// acados
#include "acados/utils/mem.h"

// #define WINDOWS_SKIP_PTR_ALIGNMENT_CHECK

// #define _USE_VALGRIND_  // uncomment to bypass assignment and do new memory allocation instead

#define _USE_MALLOC_  // acados_malloc = malloc / acados_malloc = calloc

void make_int_multiple_of(int num, int *size) { *size = (*size + num - 1) / num * num; }
int align_char_to(int num, char **c_ptr)
{
    size_t s_ptr = (size_t) *c_ptr;
    s_ptr = (s_ptr + num - 1) / num * num;
    int offset = num - (int) (s_ptr - (size_t)(*c_ptr));
    *c_ptr = (char *) s_ptr;
    return offset;
}

#ifdef _USE_VALGRIND_
// print warning when by-passing pointer and allocating new memory (for debugging)
static void print_warning() { printf(" -- using dynamically allocated memory for debugging --\n"); }
#endif

void *acados_malloc(size_t nitems, size_t size)
{
#if defined(_USE_MALLOC_)
    void *ptr = malloc(nitems * size);
#else
    void *ptr = calloc(nitems, size);
#endif
    return ptr;
}

void *acados_calloc(size_t nitems, size_t size)
{
    void *ptr = calloc(nitems, size);
    return ptr;
}

void assign_and_advance_double_ptrs(int n, double ***v, char **ptr)
{
#ifndef WINDOWS_SKIP_PTR_ALIGNMENT_CHECK
    assert((size_t) *ptr % 8 == 0 && "pointer not 8-byte aligned!");
#endif
#ifdef _USE_VALGRIND_
    *v = (double **) acados_malloc(n, sizeof(double *));
#else
    *v = (double **) *ptr;
    *ptr += sizeof(double *) * n;
#endif
}

void assign_and_advance_int_ptrs(int n, int ***v, char **ptr)
{
#ifndef WINDOWS_SKIP_PTR_ALIGNMENT_CHECK
    assert((size_t) *ptr % 8 == 0 && "pointer not 8-byte aligned!");
#endif
#ifdef _USE_VALGRIND_
    *v = (int **) acados_malloc(n, sizeof(int *));
#else
    *v = (int **) *ptr;
    *ptr += sizeof(int *) * n;
#endif
}

void assign_and_advance_blasfeo_dvec_structs(int n, struct blasfeo_dvec **sv, char **ptr)
{
#ifndef WINDOWS_SKIP_PTR_ALIGNMENT_CHECK
    assert((size_t) *ptr % 8 == 0 && "pointer not 8-byte aligned!");
#endif
#ifdef _USE_VALGRIND_
    *sv = (struct blasfeo_dvec *) acados_malloc(n, sizeof(struct blasfeo_dvec));
#else
    *sv = (struct blasfeo_dvec *) *ptr;
    *ptr += sizeof(struct blasfeo_dvec) * n;
#endif
}

void assign_and_advance_blasfeo_dmat_structs(int n, struct blasfeo_dmat **sm, char **ptr)
{
#ifndef WINDOWS_SKIP_PTR_ALIGNMENT_CHECK
    assert((size_t) *ptr % 8 == 0 && "pointer not 8-byte aligned!");
#endif
#ifdef _USE_VALGRIND_
    *sm = (struct blasfeo_dmat *) acados_malloc(n, sizeof(struct blasfeo_dmat));
#else
    *sm = (struct blasfeo_dmat *) *ptr;
    *ptr += sizeof(struct blasfeo_dmat) * n;
#endif
}

void assign_and_advance_blasfeo_dmat_ptrs(int n, struct blasfeo_dmat ***sm, char **ptr)
{
#ifndef WINDOWS_SKIP_PTR_ALIGNMENT_CHECK
    assert((size_t) *ptr % 8 == 0 && "pointer not 8-byte aligned!");
#endif
#ifdef _USE_VALGRIND_
    *sm = (struct blasfeo_dmat **) acados_malloc(n, sizeof(struct blasfeo_dmat *));
#else
    *sm = (struct blasfeo_dmat **) *ptr;
    *ptr += sizeof(struct blasfeo_dmat *) * n;
#endif
}

void assign_and_advance_char(int n, char **v, char **ptr)
{
#ifdef _USE_VALGRIND_
    *v = (char *) acados_malloc(n, sizeof(char));
    print_warning();
#else
    *v = (char *) *ptr;
    *ptr += sizeof(char) * n;
#endif
}

void assign_and_advance_int(int n, int **v, char **ptr)
{
#ifdef _USE_VALGRIND_
    *v = (int *) acados_malloc(n, sizeof(int));
    print_warning();
#else
    *v = (int *) *ptr;
    *ptr += sizeof(int) * n;
#endif
}

void assign_and_advance_double(int n, double **v, char **ptr)
{
    assert((size_t) *ptr % 8 == 0 && "double not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    *v = (double *) acados_malloc(n, sizeof(double));
    print_warning();
#else
    *v = (double *) *ptr;
    *ptr += sizeof(double) * n;
#endif
}

void assign_and_advance_blasfeo_dvec_mem(int n, struct blasfeo_dvec *sv, char **ptr)
{
    assert((size_t) *ptr % 8 == 0 && "strvec not 8-byte aligned!");

#ifdef _USE_VALGRIND_
    blasfeo_allocate_dvec(n, sv);
    print_warning();
#else
    blasfeo_create_dvec(n, sv, *ptr);
    *ptr += sv->memsize;
#endif
}

void assign_and_advance_blasfeo_dmat_mem(int m, int n, struct blasfeo_dmat *sA, char **ptr)
{
#ifdef LA_HIGH_PERFORMANCE
    assert((size_t) *ptr % 64 == 0 && "strmat not 64-byte aligned!");
#else
    assert((size_t) *ptr % 8 == 0 && "strmat not 8-byte aligned!");
#endif

#ifdef _USE_VALGRIND_
    blasfeo_allocate_dmat(m, n, sA);
    print_warning();
#else
    blasfeo_create_dmat(m, n, sA, *ptr);
    *ptr += sA->memsize;
#endif
}
