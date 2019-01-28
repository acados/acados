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

#ifndef ACADOS_UTILS_EXTERNAL_FUNCTION_GENERIC_H_
#define ACADOS_UTILS_EXTERNAL_FUNCTION_GENERIC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

/************************************************
 * generic external function
 ************************************************/

// type of arguments
typedef enum {
    COLMAJ,
    BLASFEO_DMAT,
    BLASFEO_DVEC,
    COLMAJ_ARGS,
    BLASFEO_DMAT_ARGS,
    BLASFEO_DVEC_ARGS,
    IGNORE_ARGUMENT
} ext_fun_arg_t;

struct colmaj_args
{
    double *A;
    int lda;
};

struct blasfeo_dmat_args
{
    struct blasfeo_dmat *A;
    int ai;
    int aj;
};

struct blasfeo_dvec_args
{
    struct blasfeo_dvec *x;
    int xi;
};

// prototype of an external function
typedef struct
{
    // public members (have to be before private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    // private members
    // .....
} external_function_generic;

/************************************************
 * casadi external function
 ************************************************/

typedef struct
{
    // public members (have to be the same as in the prototype, and before the private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    // private members
    void *ptr_ext_mem;  // pointer to external memory
    int (*casadi_fun)(const double **, double **, int *, double *, void *);
    int (*casadi_work)(int *, int *, int *, int *);
    const int *(*casadi_sparsity_in)(int);
    const int *(*casadi_sparsity_out)(int);
    int (*casadi_n_in)();
    int (*casadi_n_out)();
    double **args;
    double **res;
    double *w;
    int *iw;
    int *args_size;     // size of args[i]
    int *res_size;      // size of res[i]
    int args_num;       // number of args arrays
    int args_size_tot;  // total size of args arrays
    int res_num;        // number of res arrays
    int res_size_tot;   // total size of res arrays
    int in_num;         // number of input arrays
    int out_num;        // number of output arrays
    int iw_size;        // number of ints for worksapce
    int w_size;         // number of dobules for workspace
} external_function_casadi;

//
int external_function_casadi_struct_size();
//
void external_function_casadi_set_fun(external_function_casadi *fun, void *value);
//
void external_function_casadi_set_work(external_function_casadi *fun, void *value);
//
void external_function_casadi_set_sparsity_in(external_function_casadi *fun, void *value);
//
void external_function_casadi_set_sparsity_out(external_function_casadi *fun, void *value);
//
void external_function_casadi_set_n_in(external_function_casadi *fun, void *value);
//
void external_function_casadi_set_n_out(external_function_casadi *fun, void *value);
//
int external_function_casadi_calculate_size(external_function_casadi *fun);
//
void external_function_casadi_assign(external_function_casadi *fun, void *mem);
//
void external_function_casadi_wrapper(void *self, ext_fun_arg_t *type_in, void **in,
                                      ext_fun_arg_t *type_out, void **out);

/************************************************
 * casadi external parametric function
 ************************************************/

typedef struct
{
    // public members (have to be the same as in the prototype, and before the private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    // private members
    void (*set_param)(void *, double *);
    void *ptr_ext_mem;  // pointer to external memory
    int (*casadi_fun)(const double **, double **, int *, double *, void *);
    int (*casadi_work)(int *, int *, int *, int *);
    const int *(*casadi_sparsity_in)(int);
    const int *(*casadi_sparsity_out)(int);
    int (*casadi_n_in)();
    int (*casadi_n_out)();
    double **args;
    double **res;
    double *w;
    double *p;  // parameters
    int *iw;
    int *args_size;     // size of args[i]
    int *res_size;      // size of res[i]
    int args_num;       // number of args arrays
    int args_size_tot;  // total size of args arrays
    int res_num;        // number of res arrays
    int res_size_tot;   // total size of res arrays
    int in_num;         // number of input arrays
    int out_num;        // number of output arrays
    int iw_size;        // number of ints for worksapce
    int w_size;         // number of dobules for workspace
    int np;             // number of parameters
} external_function_param_casadi;

//
int external_function_param_casadi_struct_size();
//
void external_function_param_casadi_set_fun(external_function_param_casadi *fun, void *value);
//
void external_function_param_casadi_set_work(external_function_param_casadi *fun, void *value);
//
void external_function_param_casadi_set_sparsity_in(external_function_param_casadi *fun, void *value);
//
void external_function_param_casadi_set_sparsity_out(external_function_param_casadi *fun, void *value);
//
void external_function_param_casadi_set_n_in(external_function_param_casadi *fun, void *value);
//
void external_function_param_casadi_set_n_out(external_function_param_casadi *fun, void *value);
//
int external_function_param_casadi_calculate_size(external_function_param_casadi *fun, int np);
//
void external_function_param_casadi_assign(external_function_param_casadi *fun, void *mem);
//
void external_function_param_casadi_wrapper(void *self, ext_fun_arg_t *type_in, void **in,
                                            ext_fun_arg_t *type_out, void **out);
//
void external_function_param_casadi_set_param(void *self, double *p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_EXTERNAL_FUNCTION_GENERIC_H_
