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

int acados_create(ocp_nlp_solver ** _nlp_solver, ocp_nlp_in ** _nlp_in, 
        ocp_nlp_out ** _nlp_out, void ** _nlp_opts, ocp_nlp_solver_config ** _nlp_config, 
        ocp_nlp_solver_plan ** _solver_plan, ocp_nlp_dims ** _nlp_dims);

int acados_solve(ocp_nlp_solver * nlp_solver, ocp_nlp_in * nlp_in, 
        ocp_nlp_out * nlp_out, void * nlp_opts, ocp_nlp_solver_config * nlp_config, 
        ocp_nlp_solver_plan * nlp_solver_plan, ocp_nlp_dims * nlp_dims);

int acados_free(ocp_nlp_solver * nlp_solver, ocp_nlp_in * nlp_in, 
        ocp_nlp_out * nlp_out, void * nlp_opts, ocp_nlp_solver_config * nlp_config, 
        ocp_nlp_solver_plan * nlp_solver_plan, ocp_nlp_dims * nlp_dims);

#endif  // ACADOS_SOLVER_{{ra.model_name.upper()}}_H_
