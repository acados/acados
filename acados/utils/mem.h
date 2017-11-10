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


// make int counter of memory multiple of a number (typically 8 or 64)
void make_int_multiple_of(int num, int *size);


// align char pointer to number (typically 8 for pointers and doubles, 64 for blasfeo structs)
void align_char_to(int num, char **c_ptr);


// switch between malloc and calloc (for valgrinding)
void *acados_malloc(size_t nitems, size_t size);


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