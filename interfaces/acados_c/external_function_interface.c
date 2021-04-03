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

void external_function_param_generic_create(external_function_param_generic *fun, int np)
{
    acados_size_t fun_size = external_function_param_generic_calculate_size(fun, np);
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

void external_function_casadi_create(external_function_casadi *fun)
{
    acados_size_t fun_size = external_function_casadi_calculate_size(fun);
    void *fun_mem = acados_malloc(1, fun_size);
    assert(fun_mem != 0);
    external_function_casadi_assign(fun, fun_mem);

    return;
}



void external_function_casadi_create_array(int size, external_function_casadi *funs)
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
        funs_size[ii] = external_function_casadi_calculate_size(funs + ii);
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

void external_function_param_casadi_create(external_function_param_casadi *fun, int np)
{
    acados_size_t fun_size = external_function_param_casadi_calculate_size(fun, np);
    void *fun_mem = acados_malloc(1, fun_size);
    assert(fun_mem != 0);
    external_function_param_casadi_assign(fun, fun_mem);

    return;
}



void external_function_param_casadi_create_array(int size, external_function_param_casadi *funs, int np)
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
        funs_size[ii] = external_function_param_casadi_calculate_size(funs + ii, np);
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
