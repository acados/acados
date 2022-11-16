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

#include <stdlib.h>
#include <stdio.h>

#include "acados_solver_pendulum_ode.h"
#include "acados_c/ocp_nlp_interface.h"

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

static void print_x_trajectory(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    ocp_nlp_config *nlp_config = solver->config;
    ocp_nlp_dims *nlp_dims = solver->dims;

    int N = nlp_dims->N;
    int nx = nlp_dims->nx[0];

    double *buffer = malloc(nx*(N+1) * sizeof(double)); // TODO: buffer has to be big enough; malloc is slow;
    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "x", &buffer[ii*nx]);
    printf("xtraj in custom function:\t");
    d_print_exp_tran_mat( nx, N+1, buffer, nx);
    free(buffer);
}

static void print_u_trajectory(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out)
{
    ocp_nlp_config *nlp_config = solver->config;
    ocp_nlp_dims *nlp_dims = solver->dims;

    int N = nlp_dims->N;
    int nu = nlp_dims->nu[0];

    double *buffer = malloc(nu*(N+1) * sizeof(double)); // TODO: buffer has to be big enough; malloc is slow;
    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "u", &buffer[ii*nu]);
    printf("utraj in custom function:\t");
    d_print_exp_mat( nu, N+1, buffer, nu);
    free(buffer);
}



int pendulum_ode_custom_update_function(pendulum_ode_solver_capsule* capsule)
{
    // printf("\nInside the actual custom_function\n");
    ocp_nlp_config *nlp_config = pendulum_ode_acados_get_nlp_config(capsule);
    ocp_nlp_dims *nlp_dims = pendulum_ode_acados_get_nlp_dims(capsule);
    ocp_nlp_in *nlp_in = pendulum_ode_acados_get_nlp_in(capsule);
    ocp_nlp_out *nlp_out = pendulum_ode_acados_get_nlp_out(capsule);
    ocp_nlp_solver *nlp_solver = pendulum_ode_acados_get_nlp_solver(capsule);
    void *nlp_opts = pendulum_ode_acados_get_nlp_opts(capsule);

    // // EXAMPLE 1: print x trajectory
    // print_x_trajectory(nlp_solver, nlp_in, nlp_out);
    // EXAMPLE 2: print u trajectory
    print_u_trajectory(nlp_solver, nlp_in, nlp_out);

}

