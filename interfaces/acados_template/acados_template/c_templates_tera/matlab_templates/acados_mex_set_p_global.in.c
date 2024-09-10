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


// standard
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_solver_{{ name }}.h"
{% if dims.np_global > 0 %}
#include "{{ model.name }}_model/{{ model.name }}_model.h"
#include "{{ name }}_p_global_precompute_fun.h"
{% endif %}
// mex
#include "mex.h"
#include "mex_macros.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    long long *ptr;
    char fun_name[20] = "ocp_set_p_global";
    char buffer [500]; // for error messages

    /* RHS */
    int nrhs_ref = 2;

    // C ocp
    const mxArray *C_ocp = prhs[0];
    // capsule
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "capsule" ) );
    {{ name }}_solver_capsule *capsule = ({{ name }}_solver_capsule *) ptr[0];
    // dims
    ptr = (long long *) mxGetData( mxGetField( C_ocp, 0, "dims" ) );
    ocp_nlp_dims *dims = (ocp_nlp_dims *) ptr[0];

    // value
    double *value = mxGetPr( prhs[1] );

    // for checks
    int matlab_size = (int) mxGetNumberOfElements( prhs[1] );
    int nrow = (int) mxGetM( prhs[1] );
    int ncol = (int) mxGetN( prhs[1] );

    if (nrhs == nrhs_ref)
    {
    {% if dims.np_global > 0 %}
        external_function_casadi* fun = &capsule->p_global_precompute_fun;
        fun->args[0] = value;
        int np_global = {{ dims.np_global }};

        if (matlab_size != np_global)
        {
            sprintf(buffer, "acados_mex_set_p_global: np_global = %d should match data_len = %d. Exiting.\n", np_global, matlab_size);
            mexErrMsgTxt(buffer);
        }

    {% set n_pools = casadi_pool_names | length %}
    {% for ip in range(end=n_pools) %}
    {% set pool_name = casadi_pool_names[ip] %}
    {% set fun_name_split = pool_name | split(pat='|') %}
        fun->res[{{ ip }}] = {{ fun_name_split[0] }}_get_pool_double("{{ pool_name }}");
    {%- endfor %}

        fun->casadi_fun((const double **) fun->args, fun->res, fun->iw, fun->w, NULL);

    {%- else %}
        sprintf(buffer, "acados_mex_set_p_global: p_global is not defined, {{ name }}_acados_set_p_global does nothing.\n");
        mexWarnMsgTxt(buffer);
    {%- endif %}
    }
    else
    {
        sprintf(buffer, "acados_mex_set_p_global: wrong nrhs: %d\n", nrhs);
        mexErrMsgTxt(buffer);
    }

    return;
}

