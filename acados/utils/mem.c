/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


// external
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
// acados
#include "acados/utils/mem.h"


// #define _USE_VALGRIND_  // uncomment to bypass assignment and do new memory allocation instead

// #define _USE_MALLOC_  // acados_malloc = malloc / acados_malloc = calloc

int acados_get_pointer_size()
{
    // returns pointer size in bytes
    int pointer_size = sizeof(void *);
    if (pointer_size != 4 && pointer_size != 8)
    {
        printf("\nerror: acados expected pointers to be 4 or 8 bytes.\n");
        exit(1);
    }
    return pointer_size;
}

void make_int_multiple_of(acados_size_t num, acados_size_t *size)
{
    // avoid changes for num < 2
    if(num>1)
    // round integer size up to next multiple of num;
	*size = (*size + num - 1) / num * num;
}



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

void *acados_malloc(size_t nitems, acados_size_t size)
{
#if defined(_USE_MALLOC_)
    void *ptr = malloc(nitems * size);
#else
    void *ptr = calloc(nitems, size);
#endif
    return ptr;
}

void *acados_calloc(size_t nitems, acados_size_t size)
{
    void *ptr = calloc(nitems, size);
    return ptr;
}

void assign_and_advance_double_ptrs(int n, double ***v, char **ptr)
{
    assert((size_t) *ptr % acados_get_pointer_size() == 0 && "pointer not aligned properly!");
#ifdef _USE_VALGRIND_
    *v = (double **) acados_malloc(n, sizeof(double *));
#else
    *v = (double **) *ptr;
    *ptr += sizeof(double *) * n;
#endif
}

void assign_and_advance_int_ptrs(int n, int ***v, char **ptr)
{
    assert((size_t) *ptr % acados_get_pointer_size() == 0 && "pointer not aligned properly!");
#ifdef _USE_VALGRIND_
    *v = (int **) acados_malloc(n, sizeof(int *));
#else
    *v = (int **) *ptr;
    *ptr += sizeof(int *) * n;
#endif
}

void assign_and_advance_blasfeo_dvec_structs(int n, struct blasfeo_dvec **sv, char **ptr)
{
    assert(((size_t) *ptr % acados_get_pointer_size()) == 0 && "pointer not aligned properly!");
#ifdef _USE_VALGRIND_
    *sv = (struct blasfeo_dvec *) acados_malloc(n, sizeof(struct blasfeo_dvec));
#else
    *sv = (struct blasfeo_dvec *) *ptr;
    *ptr += sizeof(struct blasfeo_dvec) * n;
#endif
}

void assign_and_advance_blasfeo_dmat_structs(int n, struct blasfeo_dmat **sm, char **ptr)
{
    assert((size_t) *ptr % acados_get_pointer_size() == 0 && "pointer not aligned properly!");
#ifdef _USE_VALGRIND_
    *sm = (struct blasfeo_dmat *) acados_malloc(n, sizeof(struct blasfeo_dmat));
#else
    *sm = (struct blasfeo_dmat *) *ptr;
    *ptr += sizeof(struct blasfeo_dmat) * n;
#endif
}

void assign_and_advance_blasfeo_dmat_ptrs(int n, struct blasfeo_dmat ***sm, char **ptr)
{
    assert((size_t) *ptr % acados_get_pointer_size() == 0 && "pointer not aligned properly!");
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

void assign_and_advance_bool(int n, bool **v, char **ptr)
{
#ifdef _USE_VALGRIND_
    *v = (bool *) acados_malloc(n, sizeof(bool));
    print_warning();
#else
    *v = (bool *) *ptr;
    *ptr += sizeof(bool) * n;
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
