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


#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "acados/utils/math.h"

#include "acados/ocp_nlp/ocp_nlp_qpscaling_common.h"

#include "blasfeo_d_aux.h"
#include "blasfeo_d_blas.h"


/************************************************
 * config
 ************************************************/

acados_size_t ocp_nlp_qpscaling_config_calculate_size(void)
{
    return sizeof(ocp_nlp_qpscaling_config);
}



void *ocp_nlp_qpscaling_config_assign(void *raw_memory)
{
    return raw_memory;
}



/************************************************
 * dims
 ************************************************/

acados_size_t ocp_nlp_qpscaling_dims_calculate_size(int N)
{
    acados_size_t size = sizeof(ocp_nlp_qpscaling_dims);

    size += 3*(N+1)*sizeof(int); // nx nu ng

    return size;
}



ocp_nlp_qpscaling_dims *ocp_nlp_qpscaling_dims_assign(int N, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    // dims
    ocp_nlp_qpscaling_dims *dims = (ocp_nlp_qpscaling_dims *) c_ptr;
    c_ptr += sizeof(ocp_nlp_qpscaling_dims);
    // nx
    dims->nx = (int *) c_ptr;
    c_ptr += (N+1)*sizeof(int);
    // nu
    dims->nu = (int *) c_ptr;
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

    assert((char *) raw_memory + ocp_nlp_qpscaling_dims_calculate_size(N) >= c_ptr);

    return dims;
}



void ocp_nlp_qpscaling_dims_set(void *config_, ocp_nlp_qpscaling_dims *dims, int stage, char *field, int* value)
{
    if (!strcmp(field, "nx"))
    {
        dims->nx[stage] = *value;
    }
    else if (!strcmp(field, "nu"))
    {
        dims->nu[stage] = *value;
    }
    else if (!strcmp(field, "ng"))
    {
        dims->ng[stage] = *value;
    }
    else
    {
        printf("\nerror: field %s not available in module ocp_nlp_qpscaling_dims_set\n", field);
        exit(1);
    }

    return;
}


/************************************************
 * functions
 ************************************************/

void ocp_qp_scale_objective(ocp_qp_in *qp_in, double factor)
{
    int *nx = qp_in->dim->nx;
    int *nu = qp_in->dim->nu;
    int *ns = qp_in->dim->ns;

    struct blasfeo_dvec *rqz = qp_in->rqz;
    struct blasfeo_dmat *RSQrq = qp_in->RSQrq;
    struct blasfeo_dvec *Z = qp_in->Z;

    for (int stage = 0; stage <= qp_in->dim->N; stage++)
    {
        // scale cost
        blasfeo_dvecsc(nx[stage]+nu[stage]+2*ns[stage], factor, rqz+stage, 0);
        blasfeo_dvecsc(2*ns[stage], factor, Z+stage, 0);
        blasfeo_dgecpsc(nx[stage]+nu[stage], nx[stage]+nu[stage], factor, RSQrq+stage, 0, 0, RSQrq+stage, 0, 0);
    }
}



void ocp_qp_out_scale_duals(ocp_qp_out *qp_out, double factor)
{
    struct d_ocp_qp_dim *qp_dim = qp_out->dim;
    struct blasfeo_dvec *pi = qp_out->pi;
    struct blasfeo_dvec *lam = qp_out->lam;
    // print_ocp_qp_dims(qp_dim);
    for (int stage = 0; stage <= qp_dim->N; stage++)
    {
        blasfeo_dvecsc(2*qp_dim->nb[stage]+2*qp_dim->ng[stage]+2*qp_dim->ns[stage], factor, lam+stage, 0);
        if (stage < qp_dim->N)
            blasfeo_dvecsc(qp_dim->nx[stage+1], factor, pi+stage, 0);
    }
}
