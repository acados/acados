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


// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_{{ model.name }}.h"


int main()
{

    int status = 0;
    status = acados_create();

    if (status)
    {
        printf("acados_create() returned status %d. Exiting.\n", status);
        exit(1);
    }

    // initial condition
    int idxbx0[{{ dims.nbx_0 }}];
    {% for i in range(end=dims.nbx_0) %}
    idxbx0[{{ i }}] = {{ constraints.idxbx_0[i] }};
    {%- endfor %}

    double lbx0[{{ dims.nbx_0 }}];
    double ubx0[{{ dims.nbx_0 }}];
    {% for i in range(end=dims.nbx_0) %}
    lbx0[{{ i }}] = {{ constraints.lbx_0[i] }};
    ubx0[{{ i }}] = {{ constraints.ubx_0[i] }};
    {%- endfor %}

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);

    // initialization for state values
    double x_init[{{ dims.nx }}];
    {%- for i in range(end=dims.nx) %}
    x_init[{{ i }}] = 0.0;
    {%- endfor %}

    // initial value for control input
    double u0[{{ dims.nu }}];
    {%- for i in range(end=dims.nu) %}
    u0[{{ i }}] = 0.0;
    {%- endfor %}


  {%- if dims.np > 0 %}
    // set parameters
    double p[{{ dims.np }}];
    {% for item in parameter_values %}
    p[{{ loop.index0 }}] = {{ item }};
    {% endfor %}


    {%- if solver_options.integrator_type == "IRK" -%}
    for (int ii = 0; ii < {{ dims.N }}; ii++)
    {
        impl_dae_fun[ii].set_param(impl_dae_fun+ii, p);
        impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, p);
        impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, p);
    }
    {% elif solver_options.integrator_type == "ERK" %}
    for (int ii = 0; ii < {{ dims.N }}; ii++)
    {
        forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, p);
    }
    {%- endif %}
    for (int ii = 0; ii < {{ dims.N }}; ii++) {
        {%- if constraints.constr_type == "BGP" %}
        // r_constraint[ii].set_param(r_constraint+ii, p);
        phi_constraint[ii].set_param(phi_constraint+ii, p);
        {%- endif %}
        {%- if dims.nh > 0 %}
        h_constraint[ii].set_param(h_constraint+ii, p);
        {% endif %}
    }
    {%- if constraints.constr_type_e == "BGP" %}
    // r_e_constraint.set_param(&r_e_constraint, p);
    phi_e_constraint.set_param(&phi_e_constraint, p);
    {% endif %}
    {%- if dims.nh_e > 0 %}
    h_e_constraint.set_param(&h_e_constraint, p);
    {% endif %}
  {% endif %}{# if np > 0 #}

    // prepare evaluation
    int NTIMINGS = 1;
    double min_time = 1e12;
    double kkt_norm_inf;
    double elapsed_time;
    int sqp_iter;

    double xtraj[{{ dims.nx }} * ({{ dims.N }}+1)];
    double utraj[{{ dims.nu }} * ({{ dims.N }})];

    
    // solve ocp in loop

    int rti_phase = 0; 

    for (int ii = 0; ii < NTIMINGS; ii++)
    {
        // initialize solution
        for (int i = 0; i <= nlp_dims->N; i++)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x_init);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
        }
        ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_phase", rti_phase);
        status = acados_solve();
        ocp_nlp_get(nlp_config, nlp_solver, "time_tot", &elapsed_time);
        min_time = MIN(elapsed_time, min_time);
    }

    /* print solution and statistics */
    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "x", &xtraj[ii*{{ dims.nx }}]);
    for (int ii = 0; ii < nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "u", &utraj[ii*{{ dims.nu }}]);

    printf("\n--- xtraj ---\n");
    d_print_exp_tran_mat( {{ dims.nx }}, {{ dims.N }}+1, xtraj, {{ dims.nx }} );
    printf("\n--- utraj ---\n");
    d_print_exp_tran_mat( {{ dims.nu }}, {{ dims.N }}, utraj, {{ dims.nu }} );
    // ocp_nlp_out_print(nlp_solver->dims, nlp_out);

    printf("\nsolved ocp %d times, solution printed above\n\n", NTIMINGS);

    if (status == ACADOS_SUCCESS)
    {
        printf("acados_solve(): SUCCESS!\n");
    }
    else
    {
        printf("acados_solve() failed with status %d.\n", status);
    }

    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "kkt_norm_inf", &kkt_norm_inf);
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", &sqp_iter);

    printf("\nSolver info:\n");
    printf(" SQP iterations %2d\n minimum time for 1 solve %f [ms]\n KKT %e\n",
           sqp_iter, min_time*1000, kkt_norm_inf);

    // free solver
    status = acados_free();
    if (status) {
        printf("acados_free() returned status %d. \n", status);
    }

    return status;
}
