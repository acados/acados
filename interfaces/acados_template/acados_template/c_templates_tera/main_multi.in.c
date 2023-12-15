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

{%- if not solver_options.custom_update_filename %}
    {%- set custom_update_filename = "" %}
{% else %}
    {%- set custom_update_filename = solver_options.custom_update_filename %}
{%- endif %}

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_{{ name }}.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

{%- set nx_values = [] -%}
{%- for jj in range(end=n_phases) %}
    {%- set_global nx_values = nx_values | concat(with=(phases_dims[jj].nx)) %}
{%- endfor %}
{%- set nx_max = nx_values | sort | last %}

{%- set nu_values = [] -%}
{%- for jj in range(end=n_phases) %}
    {%- set_global nu_values = nu_values | concat(with=(phases_dims[jj].nu)) %}
{%- endfor %}
{%- set nu_max = nu_values | sort | last %}

#define {{ name | upper }}_N      {{ N_horizon }}


int main()
{

    int nx_max = {{ nx_max }};
    int nu_max = {{ nu_max }};

    {{ name }}_solver_capsule *acados_ocp_capsule = {{ name }}_acados_create_capsule();
    // there is an opportunity to change the number of shooting intervals in C without new code generation
    int N = {{ name | upper }}_N;
    // allocate the array and fill it accordingly
    double* new_time_steps = NULL;
    int status = {{ name }}_acados_create_with_discretization(acados_ocp_capsule, N, new_time_steps);

    if (status)
    {
        printf("{{ name }}_acados_create() returned status %d. Exiting.\n", status);
        exit(1);
    }

    ocp_nlp_config *nlp_config = {{ name }}_acados_get_nlp_config(acados_ocp_capsule);
    ocp_nlp_dims *nlp_dims = {{ name }}_acados_get_nlp_dims(acados_ocp_capsule);
    ocp_nlp_in *nlp_in = {{ name }}_acados_get_nlp_in(acados_ocp_capsule);
    ocp_nlp_out *nlp_out = {{ name }}_acados_get_nlp_out(acados_ocp_capsule);
    ocp_nlp_solver *nlp_solver = {{ name }}_acados_get_nlp_solver(acados_ocp_capsule);
    void *nlp_opts = {{ name }}_acados_get_nlp_opts(acados_ocp_capsule);

    // initialization for state values
    double* x_init = calloc(nx_max, sizeof(double));
    double* u_init = calloc(nu_max, sizeof(double));
    double* x_val = calloc(nx_max, sizeof(double));
    double* u_val = calloc(nu_max, sizeof(double));

    // prepare evaluation
    int NTIMINGS = 1;
    double min_time = 1e12;
    double kkt_norm_inf;
    double elapsed_time;
    int sqp_iter;

    // solve ocp in loop
    int rti_phase = 0;

    for (int ii = 0; ii < NTIMINGS; ii++)
    {
        // initialize solution
        for (int i = 0; i < N; i++)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x_init);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u_init);
        }
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, N, "x", x_init);
        ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_phase", &rti_phase);
        status = {{ name }}_acados_solve(acados_ocp_capsule);
        ocp_nlp_get(nlp_config, nlp_solver, "time_tot", &elapsed_time);

        min_time = MIN(elapsed_time, min_time);
    }

    int nx;
    int nu;
    /* print solution and statistics */
    printf("\n--- x trajectory ---\n");
    for (int ii = 0; ii <= nlp_dims->N; ii++)
    {
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "x", x_val);
        nx = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, ii, "x");
        d_print_exp_tran_mat(nx, 1, x_val, 1);
    }
    printf("\n--- u trajectory ---\n");
    for (int ii = 0; ii < nlp_dims->N; ii++)
    {
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "u", u_val);
        nu = ocp_nlp_dims_get_from_attr(nlp_config, nlp_dims, nlp_out, ii, "u");
        d_print_exp_tran_mat(nu, 1, u_val, 1);    }
    printf("\nsolved ocp %d times, solution printed above\n\n", NTIMINGS);

    if (status == ACADOS_SUCCESS)
    {
        printf("{{ name }}_acados_solve(): SUCCESS!\n");
    }
    else
    {
        printf("{{ name }}_acados_solve() failed with status %d.\n", status);
    }


{%- if custom_update_filename != "" %}
    {{ name }}_acados_custom_update(acados_ocp_capsule, xtraj, NX*(N+1));
{%- endif %}

    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "kkt_norm_inf", &kkt_norm_inf);
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", &sqp_iter);

    {{ name }}_acados_print_stats(acados_ocp_capsule);

    printf("\nSolver info:\n");
    printf("SQP iterations %2d\n minimum time for %d solve %f [ms]\n KKT %e\n",
           sqp_iter, NTIMINGS, min_time*1000, kkt_norm_inf);

    // free solver
    status = {{ name }}_acados_free(acados_ocp_capsule);
    if (status) {
        printf("{{ name }}_acados_free() returned status %d. \n", status);
    }
    // free solver capsule
    status = {{ name }}_acados_free_capsule(acados_ocp_capsule);
    if (status) {
        printf("{{ name }}_acados_free_capsule() returned status %d. \n", status);
    }

    return status;
}
