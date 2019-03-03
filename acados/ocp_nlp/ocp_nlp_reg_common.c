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

#include "acados/utils/math.h"

#include "acados/ocp_nlp/ocp_nlp_reg_common.h"



int ocp_nlp_reg_opts_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_opts);
}



void *ocp_nlp_reg_opts_assign(void *raw_memory)
{
    return raw_memory;
}



int ocp_nlp_reg_config_calculate_size(void)
{
    return sizeof(ocp_nlp_reg_config);
}



void *ocp_nlp_reg_config_assign(void *raw_memory)
{
    return raw_memory;
}



/* regularization help functions */

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




