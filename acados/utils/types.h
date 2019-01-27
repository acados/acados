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

#ifndef ACADOS_UTILS_TYPES_H_
#define ACADOS_UTILS_TYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#define MAX_STR_LEN 256
#define ACADOS_EPS 1e-12
#define ACADOS_NEG_INFTY -1.0e9
#define ACADOS_POS_INFTY +1.0e9
#define UNUSED(x) ((void)(x))



typedef double real_t;
typedef int int_t;



typedef int (*casadi_function_t)(const double** arg, double** res, int* iw, double* w, void* mem);



// enum of return values
enum return_values
{
    ACADOS_SUCCESS,
    ACADOS_FAILURE,
    ACADOS_MAXITER,
    ACADOS_MINSTEP,
    ACADOS_QP_FAILURE,
    ACADOS_READY,
};



// opts values ( please keep in alphabetical order ! )
// TODO outdated, proably remove !!!!!!!!!
enum acados_opts
{
    COMPUTE_ADJ,
//    COMPUTE_DUAL_SOL,
};



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_TYPES_H_
