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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "acados/utils/math.h"

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"



/************************************************
 * config
 ************************************************/

int ocp_nlp_reg_config_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_config);
}



void *ocp_nlp_reg_config_assign(void *raw_memory)
{
    return raw_memory;
}



/************************************************
 * dims
 ************************************************/

int ocp_nlp_reg_dims_calculate_size(int N)
{
    int size = sizeof(ocp_nlp_reg_dims);

    size += 5*(N+1)*sizeof(int); // nx nu nbu nbx ng

    return size;
}



ocp_nlp_reg_dims *ocp_nlp_reg_dims_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // dims
    ocp_nlp_reg_dims *dims = (ocp_nlp_reg_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_reg_dims);
    // nx
    dims->nx = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // nu
    dims->nu = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // nbu
    dims->nbu = (int *) c_ptr;
    // nbx
    dims->nbx = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // ng
    dims->ng = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);

    dims->N = N;

    assert((char *) raw_memory + ocp_nlp_reg_dims_calculate_size(N) >= c_ptr);

    return dims;
}



void ocp_nlp_reg_dims_set(void *config_, ocp_nlp_reg_dims *dims, int stage, char *field, int* value)
{

    if (!strcmp(field, "nx"))
    {
        dims->nx[stage] = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu[stage] = *value;
    }
    else if (!strcmp(field, "nbu"))
    {
        dims->nbu[stage] = *value;
    }
    else if (!strcmp(field, "nbx"))
    {
        dims->nbx[stage] = *value;
    }
    else if (!strcmp(field, "ng"))
    {
        dims->ng[stage] = *value;
    }
    else
    {
        printf("\nerror: field %s not available in module ocp_nlp_reg_dims_set\n", field);
        exit(1);
    }

    return;
}



/************************************************
 * regularization help functions
 ************************************************/

// reconstruct A = V * d * V'
void acados_reconstruct_A(int dim, double *A, double *V, double *d)
{
    int i, j, k;

    for (i=0; i<dim; i++)
    {
        for (j=0; j<=i; j++)
        {
            A[i*dim+j] = 0.0;
            for (k=0; k<dim; k++)
                A[i*dim+j] += V[i*dim+k] * d[k] * V[j*dim+k];
            A[j*dim+i] = A[i*dim+j];
        }
    }
}



// mirroring regularization
void acados_mirror(int dim, double *A, double *V, double *d, double *e, double epsilon)
{
    int i;

    acados_eigen_decomposition(dim, A, V, d, e);

    for (i = 0; i < dim; i++)
    {
        // project
        if (d[i] >= -epsilon && d[i] <= epsilon)
            d[i] = epsilon;
        // mirror
        else if (d[i] < 0)
            d[i] = -d[i];
    }

    acados_reconstruct_A(dim, A, V, d);
}



// projecting regularization
void acados_project(int dim, double *A, double *V, double *d, double *e, double epsilon)
{
    int i;

    acados_eigen_decomposition(dim, A, V, d, e);

    // project
    for (i = 0; i < dim; i++)
    {
        if (d[i] < epsilon)
            d[i] = epsilon;
    }

    acados_reconstruct_A(dim, A, V, d);
}




