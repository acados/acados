/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_
#define INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
// acados_c
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/sim_interface.h"

typedef enum
{
    SQP_GN,
} ocp_nlp_solver_t;



typedef enum
{
    LINEAR_LS,
    NONLINEAR_LS,
    EXTERNALLY_PROVIDED,
} ocp_nlp_cost_t;



typedef enum
{
    CONTINUOUS_MODEL,
    DISCRETE_MODEL,
} ocp_nlp_dynamics_t;



typedef enum
{
    BGH,
    BGHP,
} ocp_nlp_constraints_t;



typedef struct
{
    ocp_qp_solver_plan ocp_qp_solver_plan;
    sim_solver_plan *sim_solver_plan;
    ocp_nlp_solver_t nlp_solver;
    ocp_nlp_cost_t *nlp_cost;
    ocp_nlp_dynamics_t *nlp_dynamics;
    ocp_nlp_constraints_t *nlp_constraints;
} ocp_nlp_solver_plan;



typedef struct
{
    ocp_nlp_solver_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} ocp_nlp_solver;



//
ocp_nlp_solver_plan *ocp_nlp_plan_create(int N);
//
ocp_nlp_solver_config *ocp_nlp_config_create(ocp_nlp_solver_plan plan, int N);
//
ocp_nlp_dims *ocp_nlp_dims_create(void *config);
//
ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);
//
int nlp_set_model_in_stage(ocp_nlp_solver_config *config, ocp_nlp_in *in, int stage,
                           const char *fun_type, void *fun_ptr);
//
ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);
//
void *ocp_nlp_opts_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_solver *ocp_nlp_create(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *opts_);
//
int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *qp_in, ocp_nlp_out *qp_out);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_
