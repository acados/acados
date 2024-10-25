/*
 * Copyright (c) The acados authors.
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


// external_function_opts
typedef struct
{
    bool external_workspace;
    bool with_global_data;
} external_function_opts;



// prototype of an external function
typedef struct
{
    // public members (have to be before private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    size_t (*get_external_workspace_requirement)(void *);
    void (*set_external_workspace)(void *, void *);
    // private members
    // .....
} external_function_generic;


size_t external_function_get_workspace_requirement_if_defined(external_function_generic *fun);

void external_function_set_fun_workspace_if_defined(external_function_generic *fun, void *work_);

void external_function_opts_set_to_default(external_function_opts *opts);




/************************************************
 * generic external parametric function
 ************************************************/

// prototype of a parametric external function
typedef struct
{
    // public members for core (have to be before private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    size_t (*get_external_workspace_requirement)(void *);
    void (*set_external_workspace)(void *, void *);
    // public members for interfaces
    void (*get_nparam)(void *, int *);
    void (*set_param)(void *, double *);
    void (*set_param_sparse)(void *, int n_update, int *idx, double *);
    // private members
    void *ptr_ext_mem;  // pointer to external memory
    int (*fun)(void **, void **, void *);
    double *p;  // parameters
    int np;     // number of parameters
    external_function_opts opts;
    // .....
} external_function_param_generic;

//
acados_size_t external_function_param_generic_struct_size();
//
acados_size_t external_function_param_generic_calculate_size(external_function_param_generic *fun, int np, external_function_opts *opts_);
//
void external_function_param_generic_assign(external_function_param_generic *fun, void *mem);
//
void external_function_param_generic_wrapper(void *self, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out);
//
void external_function_param_generic_get_nparam(void *self, int *np);
//
void external_function_param_generic_set_param(void *self, double *p);
//
size_t external_function_param_generic_get_external_workspace_requirement(void *self);
//
void external_function_param_generic_set_external_workspace(void *self, void *workspace);

/************************************************
 * casadi external function
 ************************************************/

typedef struct
{
    // public members (have to be the same as in the prototype, and before the private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    size_t (*get_external_workspace_requirement)(void *);
    void (*set_external_workspace)(void *, void *);
    // private members
    void *ptr_ext_mem;  // pointer to external memory
    int (*casadi_fun)(const double **, double **, int *, double *, void *);
    int (*casadi_work)(int *, int *, int *, int *);
    const int *(*casadi_sparsity_in)(int);
    const int *(*casadi_sparsity_out)(int);
    int (*casadi_n_in)(void);
    int (*casadi_n_out)(void);
    double **args;
    double **res;
    double *float_work;
    int *int_work;
    int *args_size;     // size of args[i]
    int *res_size;      // size of res[i]
    int *args_dense;    // indicates if args[i] is dense
    int *res_dense;     // indicates if res[i] is dense
    int args_num;       // number of args arrays
    int args_size_tot;  // total size of args arrays
    int res_num;        // number of res arrays
    int res_size_tot;   // total size of res arrays
    int in_num;         // number of input arrays
    int out_num;        // number of output arrays
    int int_work_size;        // number of ints for worksapce
    int float_work_size;         // number of doubles for workspace
    external_function_opts opts;
} external_function_casadi;

//
acados_size_t external_function_casadi_struct_size();
//
acados_size_t external_function_casadi_calculate_size(external_function_casadi *fun, external_function_opts *opts_);
//
void external_function_casadi_assign(external_function_casadi *fun, void *mem);
//
void external_function_casadi_wrapper(void *self, ext_fun_arg_t *type_in, void **in,
                                      ext_fun_arg_t *type_out, void **out);
//
size_t external_function_casadi_get_external_workspace_requirement(void *self);
//
void external_function_casadi_set_external_workspace(void *self, void *workspace);

/************************************************
 * casadi external parametric function
 ************************************************/

typedef struct
{
    // public members for core (have to be the same as in the prototype, and before the private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    size_t (*get_external_workspace_requirement)(void *);
    void (*set_external_workspace)(void *, void *);
    // public members for interfaces
    void (*get_nparam)(void *, int *);
    void (*set_param)(void *, double *);
    void (*set_param_sparse)(void *, int n_update, int *idx, double *);
    // private members
    void *ptr_ext_mem;  // pointer to external memory
    int (*casadi_fun)(const double **, double **, int *, double *, void *);
    int (*casadi_work)(int *, int *, int *, int *);
    const int *(*casadi_sparsity_in)(int);
    const int *(*casadi_sparsity_out)(int);
    int (*casadi_n_in)(void);
    int (*casadi_n_out)(void);
    double **args;
    double **res;
    double *float_work;
    int *int_work;
    int *args_size;     // size of args[i]
    int *res_size;      // size of res[i]
    int *args_dense;    // indicates if args[i] is dense
    int *res_dense;     // indicates if res[i] is dense
    int args_num;       // number of args arrays
    int args_size_tot;  // total size of args arrays
    int res_num;        // number of res arrays
    int res_size_tot;   // total size of res arrays
    int in_num;         // number of input arrays
    int out_num;        // number of output arrays
    int int_work_size;        // number of ints for worksapce
    int float_work_size;         // number of doubles for workspace
    int np;             // number of parameters
    int idx_in_p;
    external_function_opts opts;
} external_function_param_casadi;

//
acados_size_t external_function_param_casadi_struct_size();
//
acados_size_t external_function_param_casadi_calculate_size(external_function_param_casadi *fun, int np, external_function_opts *opts_);
//
void external_function_param_casadi_assign(external_function_param_casadi *fun, void *mem);
//
void external_function_param_casadi_wrapper(void *self, ext_fun_arg_t *type_in, void **in,
                                            ext_fun_arg_t *type_out, void **out);
//
void external_function_param_casadi_get_nparam(void *self, int *np);
//
size_t external_function_param_casadi_get_external_workspace_requirement(void *self);
//
void external_function_param_casadi_set_external_workspace(void *self, void *workspace);


/************************************************
 * external_function_external_param_casadi
 ************************************************/

typedef struct
{
    // public members for core (have to be the same as in the prototype, and before the private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    size_t (*get_external_workspace_requirement)(void *);
    void (*set_external_workspace)(void *, void *);
    void (*set_global_data_pointer)(void *, double *);
    // public members for interfaces
    void (*set_param_pointer)(void *, double *);
    // private members
    void *ptr_ext_mem;  // pointer to external memory
    int (*casadi_fun)(const double **, double **, int *, double *, void *);
    int (*casadi_work)(int *, int *, int *, int *);
    const int *(*casadi_sparsity_in)(int);
    const int *(*casadi_sparsity_out)(int);
    int (*casadi_n_in)(void);
    int (*casadi_n_out)(void);
    double **args;
    double **res;
    double *float_work;
    int *int_work;
    int *args_size;     // size of args[i]
    int *res_size;      // size of res[i]
    int *args_dense;    // indicates if args[i] is dense
    int *res_dense;     // indicates if res[i] is dense
    int args_num;       // number of args arrays
    int args_size_tot;  // total size of args arrays
    int res_num;        // number of res arrays
    int res_size_tot;   // total size of res arrays
    int in_num;         // number of input arrays
    int out_num;        // number of output arrays
    int int_work_size;        // number of ints for worksapce
    int float_work_size;         // number of doubles for workspace

    bool param_mem_is_set;  // indicates if param memory is set;
    bool global_data_ptr_is_set;  // indicates if global data pointer is set;
    int idx_in_p;
    int idx_in_global_data;

    external_function_opts opts;
} external_function_external_param_casadi;

//
acados_size_t external_function_external_param_casadi_struct_size();
//
acados_size_t external_function_external_param_casadi_calculate_size(external_function_external_param_casadi *fun, external_function_opts *opts_);
//
void external_function_external_param_casadi_assign(external_function_external_param_casadi *fun, void *mem);
//
void external_function_external_param_casadi_wrapper(void *self, ext_fun_arg_t *type_in, void **in,
                                            ext_fun_arg_t *type_out, void **out);
//
size_t external_function_external_param_casadi_get_external_workspace_requirement(void *self);
//
void external_function_external_param_casadi_set_external_workspace(void *self, void *workspace);


/************************************************
 * external_function_external_param_generic
 ************************************************/

// prototype of a parametric external function
typedef struct
{
    // public members for core (have to be before private ones)
    void (*evaluate)(void *, ext_fun_arg_t *, void **, ext_fun_arg_t *, void **);
    size_t (*get_external_workspace_requirement)(void *);
    void (*set_external_workspace)(void *, void *);
    void (*set_global_data_pointer)(void *, double *);
    // public members for interfaces
    void (*set_param_pointer)(void *, double *);

    // private members
    void *ptr_ext_mem;  // pointer to external memory
    int (*fun)(void **, void **, void *);
    double *p;  // parameters
    bool param_mem_is_set;
    external_function_opts opts;

} external_function_external_param_generic;

//
acados_size_t external_function_external_param_generic_struct_size();
//
acados_size_t external_function_external_param_generic_calculate_size(external_function_external_param_generic *fun, external_function_opts *opts_);
//
void external_function_external_param_generic_assign(external_function_external_param_generic *fun, void *mem);
//
void external_function_external_param_generic_wrapper(void *self, ext_fun_arg_t *type_in, void **in, ext_fun_arg_t *type_out, void **out);
//
void external_function_external_param_generic_set_param_ptr(void *self, double *p);
//
size_t external_function_external_param_generic_get_external_workspace_requirement(void *self);
//
void external_function_external_param_generic_set_external_workspace(void *self, void *workspace);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_EXTERNAL_FUNCTION_GENERIC_H_
