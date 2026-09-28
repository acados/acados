/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#ifndef INTERFACES_ACADOS_C_EXTERNAL_FUNCTION_INTERFACE_H_
#define INTERFACES_ACADOS_C_EXTERNAL_FUNCTION_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/external_function_generic.h"



/************************************************
 * generic external parametric function
 ************************************************/

//
void external_function_param_generic_create(external_function_param_generic *fun, int np, external_function_opts *opts_);
//
void external_function_param_generic_free(external_function_param_generic *fun);



/************************************************
 * casadi external function
 ************************************************/

//
void external_function_casadi_create(external_function_casadi *fun, external_function_opts *opts_);
//
void external_function_casadi_free(external_function_casadi *fun);
//
void external_function_casadi_create_array(int size, external_function_casadi *funs, external_function_opts *opts_);
//
void external_function_casadi_free_array(int size, external_function_casadi *funs);



/************************************************
 * casadi external parametric function
 ************************************************/

//
void external_function_param_casadi_create(external_function_param_casadi *fun, int np, external_function_opts *opts_);
//
void external_function_param_casadi_free(external_function_param_casadi *fun);
//
void external_function_param_casadi_create_array(int size, external_function_param_casadi *funs,
                                                 int np, external_function_opts *opts_);
//
void external_function_param_casadi_free_array(int size, external_function_param_casadi *funs);



/************************************************
 * external_function_external_param_casadi
 ************************************************/

//
void external_function_external_param_casadi_create(external_function_external_param_casadi *fun, external_function_opts *opts_);
//
void external_function_external_param_casadi_free(external_function_external_param_casadi *fun);

/************************************************
 * external_function_external_param_generic
 ************************************************/

//
void external_function_external_param_generic_create(external_function_external_param_generic *fun, external_function_opts *opts_);
//
void external_function_external_param_generic_free(external_function_external_param_generic *fun);





#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_EXTERNAL_FUNCTION_INTERFACE_H_
