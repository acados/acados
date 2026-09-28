/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "acados/utils/math.h"

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"



/************************************************
 * config
 ************************************************/

acados_size_t ocp_nlp_reg_config_calculate_size(void)
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

acados_size_t ocp_nlp_reg_dims_calculate_size(int N)
{
    acados_size_t size = sizeof(ocp_nlp_reg_dims);

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
    c_ptr += (N+1)*sizeof(int);
    // nbx
    dims->nbx = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // ng
    dims->ng = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);

    dims->N = N;

    // initialize to zero by default
    int ii;
    // nx
    for(ii=0; ii<=N; ii++)
        dims->nx[ii] = 0;
    // nu
    for(ii=0; ii<=N; ii++)
        dims->nu[ii] = 0;
    // nbx
    for(ii=0; ii<=N; ii++)
        dims->nbx[ii] = 0;
    // nbu
    for(ii=0; ii<=N; ii++)
        dims->nbu[ii] = 0;
    // ng
    for(ii=0; ii<=N; ii++)
        dims->ng[ii] = 0;

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

    // mirror
    for (i = 0; i < dim; i++)
    {
        if (d[i] >= -epsilon && d[i] <= epsilon)
            d[i] = epsilon;
        else if (d[i] < 0)
            d[i] = -d[i];
    }

    acados_reconstruct_A(dim, A, V, d);
}

void acados_mirror_adaptive_eps(int dim, double *A, double *V, double *d, double *e, double max_cond_block, double min_eps)
{
    int i;
    acados_eigen_decomposition(dim, A, V, d, e);
    double max_eig = 0.0;
    double eps;

    // compute max and min eigenvalues
    for (i=0; i < dim; i++)
    {
        max_eig = MAX(max_eig, fabs(d[i]));
    }
    eps = MAX(max_eig/max_cond_block, min_eps);

    // mirror
    for (i = 0; i < dim; i++)
    {
        if (d[i] >= -eps && d[i] <= eps)
            d[i] = eps;
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


void acados_project_adaptive_eps(int dim, double *A, double *V, double *d, double *e, double max_cond_block, double min_eps)
{
    int i;
    acados_eigen_decomposition(dim, A, V, d, e);
    double max_eig = 0.0;
    double eps;

    // compute max and min eigenvalues
    for (i=0; i < dim; i++)
    {
        max_eig = MAX(max_eig, d[i]);
    }
    eps = MAX(max_eig/max_cond_block, min_eps);

    // project
    for (i = 0; i < dim; i++)
    {
        if (d[i] < eps)
            d[i] = eps;
    }

    acados_reconstruct_A(dim, A, V, d);
}
