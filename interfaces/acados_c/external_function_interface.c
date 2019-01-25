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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados_c/external_function_interface.h"

#include "acados/utils/external_function_generic.h"

/************************************************
 * casadi external function
 ************************************************/

void external_function_casadi_create(external_function_casadi *fun)
{
    int fun_size = external_function_casadi_calculate_size(fun);
    void *fun_mem = malloc(fun_size);
    external_function_casadi_assign(fun, fun_mem);

    return;
}



void external_function_casadi_create_array(int size, external_function_casadi *funs)
{
    // loop index
    int ii;

    char *c_ptr;

    // create size array
    int *funs_size = malloc(size * sizeof(int));
    int funs_size_tot = 0;

    // compute sizes
    for (ii = 0; ii < size; ii++)
    {
        funs_size[ii] = external_function_casadi_calculate_size(funs + ii);
        funs_size_tot += funs_size[ii];
    }

    // allocate memory
    void *funs_mem = malloc(funs_size_tot);

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
    int fun_size = external_function_param_casadi_calculate_size(fun, np);
    void *fun_mem = malloc(fun_size);
    external_function_param_casadi_assign(fun, fun_mem);

    return;
}



void external_function_param_casadi_create_array(int size, external_function_param_casadi *funs,
                                                 int np)
{
    // loop index
    int ii;

    char *c_ptr;

    // create size array
    int *funs_size = malloc(size * sizeof(int));
    int funs_size_tot = 0;

    // compute sizes
    for (ii = 0; ii < size; ii++)
    {
        funs_size[ii] = external_function_param_casadi_calculate_size(funs + ii, np);
        funs_size_tot += funs_size[ii];
    }

    // allocate memory
    void *funs_mem = malloc(funs_size_tot);

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
