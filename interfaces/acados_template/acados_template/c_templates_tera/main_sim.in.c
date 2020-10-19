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
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_sim_solver_{{ model.name }}.h"


int main()
{
    int status = 0;
    status = {{ model.name }}_acados_sim_create();

    if (status)
    {
        printf("acados_create() returned status %d. Exiting.\n", status);
        exit(1);
    }

    // initial condition
    double x0[{{ dims.nx }}];
    {% if constraints.lbx_0 %}
    {%- for i in range(end=dims.nx) %}
    x0[{{ i }}] = {{ constraints.lbx_0[i] }};
    {%- endfor %}
    {%- else %}
    {%- for i in range(end=dims.nx) %}
    x0[{{ i }}] = 0.0;
    {%- endfor %}
    {%- endif %}

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

    {{ model.name }}_acados_sim_update_params(p, {{ dims.np }});
  {% endif %}{# if np > 0 #}

    int n_sim_steps = 3;
    // solve ocp in loop
    for (int ii = 0; ii < n_sim_steps; ii++)
    {
        sim_in_set({{ model.name }}_sim_config, {{ model.name }}_sim_dims,
            {{ model.name }}_sim_in, "x", x0);
        status = {{ model.name }}_acados_sim_solve();

        if (status != ACADOS_SUCCESS)
        {
            printf("acados_solve() failed with status %d.\n", status);
        }

        sim_out_get({{ model.name }}_sim_config, {{ model.name }}_sim_dims,
               {{ model.name }}_sim_out, "x", x0);
        
        printf("\nx0, %d\n", ii);
        for (int jj = 0; jj < {{ dims.nx }}; jj++)
        {
            printf("%e\n", x0[jj]);
        }
    }

    printf("\nPerformed %d simulation steps with acados integrator successfully.\n\n", n_sim_steps);

    // free solver
    status = {{ model.name }}_acados_sim_free();
    if (status) {
        printf("{{ model.name }}_acados_sim_free() returned status %d. \n", status);
    }

    return status;
}
