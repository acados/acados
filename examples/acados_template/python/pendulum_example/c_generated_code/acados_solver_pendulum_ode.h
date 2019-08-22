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

#ifndef ACADOS_SOLVER_pendulum_ode_H_
#define ACADOS_SOLVER_pendulum_ode_H_

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

int acados_create();
int acados_solve();
int acados_free();

ocp_nlp_in * acados_get_nlp_in();
ocp_nlp_out * acados_get_nlp_out();
ocp_nlp_solver * acados_get_nlp_solver();
ocp_nlp_config * acados_get_nlp_config();
void * acados_get_nlp_opts();
ocp_nlp_dims * acados_get_nlp_dims();

#ifdef __cplusplus
}
#endif


// ** global data **
extern ocp_nlp_in * nlp_in;
extern ocp_nlp_out * nlp_out;
extern ocp_nlp_solver * nlp_solver;
extern void * nlp_opts;
extern ocp_nlp_plan * nlp_solver_plan;
extern ocp_nlp_config * nlp_config;
extern ocp_nlp_dims * nlp_dims;


extern external_function_casadi * forw_vde_casadi;








#endif  // ACADOS_SOLVER_pendulum_ode_H_