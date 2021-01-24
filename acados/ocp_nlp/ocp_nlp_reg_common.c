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


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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




