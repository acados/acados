/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "acados_c/external_function_interface.h"

#include "acados/utils/external_function_generic.h"

#include "acados/utils/mem.h"



/************************************************
 * generic external parametric function
 ************************************************/

void external_function_param_generic_create(external_function_param_generic *fun, int np, external_function_opts *opts_)
{
    acados_size_t fun_size = external_function_param_generic_calculate_size(fun, np, opts_);
    void *fun_mem = acados_malloc(1, fun_size);
    assert(fun_mem != 0);
    external_function_param_generic_assign(fun, fun_mem);

    return;
}



void external_function_param_generic_free(external_function_param_generic *fun)
{
    free(fun->ptr_ext_mem);

    return;
}



/************************************************
 * casadi external function
 ************************************************/

void external_function_casadi_create(external_function_casadi *fun, external_function_opts *opts_)
{
    acados_size_t fun_size = external_function_casadi_calculate_size(fun, opts_);
    void *fun_mem = acados_malloc(1, fun_size);
    assert(fun_mem != 0);
    external_function_casadi_assign(fun, fun_mem);

    return;
}



void external_function_casadi_create_array(int size, external_function_casadi *funs, external_function_opts *opts_)
{
    // loop index
    int ii;

    char *c_ptr;

    // create size array
    acados_size_t *funs_size = (acados_size_t *) acados_malloc(1, size * sizeof(acados_size_t));
    assert(funs_size != 0);
    // acados_size_t *funs_size = malloc(size * sizeof(acados_size_t));
    acados_size_t funs_size_tot = 0;

    // compute sizes
    for (ii = 0; ii < size; ii++)
    {
        funs_size[ii] = external_function_casadi_calculate_size(funs + ii, opts_);
        funs_size_tot += funs_size[ii];
    }

    // allocate memory
    void *funs_mem = acados_malloc(1, funs_size_tot);
    assert(funs_mem != 0);

    // assign
    c_ptr = funs_mem;
    for (ii = 0; ii < size; ii++)
    {
        external_function_casadi_assign(funs + ii, c_ptr);
        c_ptr += funs_size[ii];
    }

    // free size array
    free(funs_size);

    return;
}



void external_function_casadi_free(external_function_casadi *fun)
{
    free(fun->ptr_ext_mem);

    return;
}



void external_function_casadi_free_array(int size, external_function_casadi *funs)
{
    free(funs[0].ptr_ext_mem);

    return;
}



/************************************************
 * casadi external parametric function
 ************************************************/

void external_function_param_casadi_create(external_function_param_casadi *fun, int np, external_function_opts *opts_)
{
    acados_size_t fun_size = external_function_param_casadi_calculate_size(fun, np, opts_);
    void *fun_mem = acados_malloc(1, fun_size);
    assert(fun_mem != 0);
    external_function_param_casadi_assign(fun, fun_mem);

    return;
}



void external_function_param_casadi_create_array(int size, external_function_param_casadi *funs, int np, external_function_opts *opts_)
{
    // loop index
    int ii;

    char *c_ptr;

    // create size array
    acados_size_t *funs_size = (acados_size_t *) acados_malloc(1, size * sizeof(acados_size_t));
    assert(funs_size != 0);
    // acados_size_t *funs_size = malloc(size * sizeof(acados_size_t));
    acados_size_t funs_size_tot = 0;

    // compute sizes
    for (ii = 0; ii < size; ii++)
    {
        funs_size[ii] = external_function_param_casadi_calculate_size(funs + ii, np, opts_);
        funs_size_tot += funs_size[ii];
    }

    // allocate memory
    void *funs_mem = acados_malloc(1, funs_size_tot);
    assert(funs_mem != 0);

    // assign
    c_ptr = funs_mem;
    for (ii = 0; ii < size; ii++)
    {
        external_function_param_casadi_assign(funs + ii, c_ptr);
        c_ptr += funs_size[ii];
    }

    // free size array
    free(funs_size);

    return;
}



void external_function_param_casadi_free(external_function_param_casadi *fun)
{
    free(fun->ptr_ext_mem);

    return;
}



void external_function_param_casadi_free_array(int size, external_function_param_casadi *funs)
{
    free(funs[0].ptr_ext_mem);

    return;
}


// external_function_external_param

void external_function_external_param_casadi_create(external_function_external_param_casadi *fun, external_function_opts *opts_)
{
    acados_size_t fun_size = external_function_external_param_casadi_calculate_size(fun, opts_);
    void *fun_mem = acados_malloc(1, fun_size);
    assert(fun_mem != 0);
    external_function_external_param_casadi_assign(fun, fun_mem);

    return;
}


void external_function_external_param_casadi_free(external_function_external_param_casadi *fun)
{
    free(fun->ptr_ext_mem);

    return;
}


// external_function_external_param

void external_function_external_param_generic_create(external_function_external_param_generic *fun, external_function_opts *opts_)
{
    acados_size_t fun_size = external_function_external_param_generic_calculate_size(fun, opts_);
    void *fun_mem = acados_malloc(1, fun_size);
    assert(fun_mem != 0);
    external_function_external_param_generic_assign(fun, fun_mem);

    return;
}


void external_function_external_param_generic_free(external_function_external_param_generic *fun)
{
    free(fun->ptr_ext_mem);

    return;
}