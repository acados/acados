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

#ifndef ACADOS_UTILS_MEM_H_
#define ACADOS_UTILS_MEM_H_

#ifdef __cplusplus
extern "C" {
#endif

// external
#include <stdio.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"


// TODO(dimitris): probably does not belong here
typedef struct {
    int (*fun)(void *);
    int (*calculate_args_size)(void *);
    void *(*assign_args)(void *);
    void (*initialize_default_args)(void *);
    int (*calculate_memory_size)(void *);
    void *(*assign_memory)(void *);
    int (*calculate_workspace_size)(void *);
} module_solver;


void copy_module_pointers_to_args(void *solver_in_args_, void *solver_);


// make int counter of memory multiple of a number (typically 8 or 64)
void make_int_multiple_of(int num, int *size);


// align char pointer to number (typically 8 for pointers and doubles,
// 64 for blasfeo structs) and return offset
int align_char_to(int num, char **c_ptr);


// switch between malloc and calloc (for valgrinding)
void *acados_malloc(size_t nitems, size_t size);


// allocate vector of pointers to vectors of doubles and advance pointer
void assign_double_ptrs(int n, double ***v, char **ptr);


// allocate vector of pointers to vectors of ints and advance pointer
void assign_int_ptrs(int n, int ***v, char **ptr);


// allocate vector of pointers to strvecs and advance pointer
void assign_strvec_ptrs(int n, struct d_strvec **sv, char **ptr);


// allocate vector of pointers to strmats and advance pointer
void assign_strmat_ptrs(int n, struct d_strmat **sm, char **ptr);


// allocate vector of ints and advance pointer
void assign_int(int n, int **v, char **ptr);


// allocate vector of doubles and advance pointer
void assign_double(int n, double **v, char **ptr);


// allocate strvec and advance pointer
void assign_strvec(int n, struct d_strvec *sv, char **ptr);


// allocate strmat and advance pointer
void assign_strmat(int m, int n, struct d_strmat *sA, char **ptr);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_MEM_H_
