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
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_pendulum_ode.h"
#include "acados_sim_solver_pendulum_ode.h"


int main() {

    // test integrator first
    int sim_status = 0;
    sim_status = pendulum_ode_acados_sim_create();

    // Set sim input
    double x_sim[4];
    
    x_sim[0] = 0.0;
    
    x_sim[1] = 3.14;
    
    x_sim[2] = 0.0;
    
    x_sim[3] = 0.0;
    

    // Set initial condition
    double u_sim[1];
    
    u_sim[0] = 0;
    

    // Set forward seeds
    double S_forw[4 * (4 + 1)];
    for (int ii = 0; ii < 4 * (4 + 1); ii++)
        S_forw[ii] = 0.0;
    for (int ii = 0; ii < 4; ii++)
        S_forw[ii * (4 + 1)] = 1.0;

    // Set discretization time
    double Td = 2.0/ 50;

    sim_in_set(pendulum_ode_sim_config, pendulum_ode_sim_dims, pendulum_ode_sim_in, "T", &Td);
    // Set initial state
    sim_in_set(pendulum_ode_sim_config, pendulum_ode_sim_dims, pendulum_ode_sim_in, "x", x_sim);
    // Set initial control
    sim_in_set(pendulum_ode_sim_config, pendulum_ode_sim_dims, pendulum_ode_sim_in, "u", u_sim);
    sim_in_set(pendulum_ode_sim_config, pendulum_ode_sim_dims, pendulum_ode_sim_in, "S_forw", &S_forw);



    sim_status = pendulum_ode_acados_sim_solve();
    // get and print output
    double *xn_out = calloc( 4, sizeof(double));
    sim_out_get(pendulum_ode_sim_config, pendulum_ode_sim_dims, pendulum_ode_sim_out, "xn", xn_out);
    printf("\nxn: \n");
    d_print_exp_mat(1, 4, xn_out, 1);

    double *S_forw_out = calloc(4*(4+1), sizeof(double));
    if (pendulum_ode_sim_opts->sens_forw){
        sim_out_get(pendulum_ode_sim_config, pendulum_ode_sim_dims, pendulum_ode_sim_out, "S_forw", S_forw_out);
        printf("\nS_forw_out: \n");
        d_print_exp_mat(4, 4 * (4 + 1), S_forw_out, 4);
    }

    int status = 0;
    status = acados_create();

    if (status) {
        printf("acados_create() returned status %d. Exiting.\n", status);
        exit(1); }

    // set initial condition
    double x0[4];
    
    x0[0] = 0.0;
    
    x0[1] = 3.14;
    
    x0[2] = 0.0;
    
    x0[3] = 0.0;
    

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", x0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", x0);

    

    

    double kkt_norm_inf = 1e12, elapsed_time;

#if 1
    int NTIMINGS = 100;
    double min_time = 1e12;
    for (int ii = 0; ii < NTIMINGS; ii ++) {
        for (int i = 0; i <= nlp_dims->N; ++i) {
            blasfeo_dvecse(nlp_dims->nu[i]+nlp_dims->nx[i], 0.0, nlp_out->ux+i, 0);
        }
        status = acados_solve();
        elapsed_time = nlp_out->total_time;
        if (elapsed_time < min_time) min_time = elapsed_time;
    }
    elapsed_time = min_time;
#else
    status = acados_solve();
#endif
    kkt_norm_inf = nlp_out->inf_norm_res;
    elapsed_time = nlp_out->total_time;
    printf(" iterations %2d | time  %f |  KKT %e\n", nlp_out->sqp_iter, elapsed_time, kkt_norm_inf);

    printf("\n--- solution ---\n");
    ocp_nlp_out_print(nlp_solver->dims, nlp_out);
    if (status) {
        printf("acados_solve() returned status %d.\n", status);
    }

    status = acados_free();

    if (status) {
        printf("acados_free() returned status %d. \n", status);
    }

    return status;
}