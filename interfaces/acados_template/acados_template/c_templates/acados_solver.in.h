#ifndef ACADOS_SOLVER_{{ra.model_name.upper()}}_H_
#define ACADOS_SOLVER_{{ra.model_name.upper()}}_H_

#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

ocp_nlp_in * nlp_in;
ocp_nlp_out * nlp_out;
ocp_nlp_solver * nlp_solver;
void * nlp_opts;
ocp_nlp_solver_plan * nlp_solver_plan;
ocp_nlp_solver_config * nlp_config;
ocp_nlp_dims * nlp_dims;

int acados_create();
int acados_solve();
int acados_free();

void * acados_get_nlp_in();
void * acados_get_nlp_out();
void * acados_get_nlp_solver();
void * acados_get_nlp_config();
void * acados_get_nlp_opts();
void * acados_get_nlp_dims();

#endif  // ACADOS_SOLVER_{{ra.model_name.upper()}}_H_
