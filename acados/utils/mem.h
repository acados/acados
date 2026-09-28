/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef ACADOS_UTILS_MEM_H_
#define ACADOS_UTILS_MEM_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>

#include "types.h"

// blasfeo
#include "blasfeo_d_aux.h"
#include "blasfeo_d_aux_ext_dep.h"

// make int counter of memory multiple of a number (typically 8 or 64)
void make_int_multiple_of(acados_size_t num, acados_size_t *size);

// align char pointer to number (typically 8 for pointers and doubles,
// 64 for blasfeo structs) and return offset
int align_char_to(int num, char **c_ptr);

// switch between malloc and calloc (for valgrinding)
void *acados_malloc(size_t nitems, acados_size_t size);

// uses always calloc
void *acados_calloc(size_t nitems, acados_size_t size);

// allocate vector of pointers to vectors of doubles and advance pointer
void assign_and_advance_double_ptrs(int n, double ***v, char **ptr);

// allocate vector of pointers to vectors of ints and advance pointer
void assign_and_advance_int_ptrs(int n, int ***v, char **ptr);

// allocate vector of pointers to strvecs and advance pointer
void assign_and_advance_blasfeo_dvec_structs(int n, struct blasfeo_dvec **sv, char **ptr);

// allocate vector of pointers to strmats and advance pointer
void assign_and_advance_blasfeo_dmat_structs(int n, struct blasfeo_dmat **sm, char **ptr);

// allocate vector of pointers to vector of pointers to strmats and advance pointer
void assign_and_advance_blasfeo_dmat_ptrs(int n, struct blasfeo_dmat ***sm, char **ptr);

// allocate vector of chars and advance pointer
void assign_and_advance_char(int n, char **v, char **ptr);

// allocate vector of ints and advance pointer
void assign_and_advance_int(int n, int **v, char **ptr);

// allocate vector of bools and advance pointer
void assign_and_advance_bool(int n, bool **v, char **ptr);

// allocate vector of doubles and advance pointer
void assign_and_advance_double(int n, double **v, char **ptr);

// allocate strvec and advance pointer
void assign_and_advance_blasfeo_dvec_mem(int n, struct blasfeo_dvec *sv, char **ptr);

// allocate strmat and advance pointer
void assign_and_advance_blasfeo_dmat_mem(int m, int n, struct blasfeo_dmat *sA, char **ptr);

// print pointer alignment
void print_pointer_alignment(char **ptr);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_MEM_H_
