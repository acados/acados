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
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/sim/sim_gnsf.h"
// acados_c
#include "acados_c/ocp_qp_interface.h"
#include "acados_c/sim_interface.h"

typedef enum
{
    SQP,
    SQP_RTI,
    INVALID_SCHEME,
} ocp_nlp_solver_t;



typedef enum
{
    LINEAR_LS,
    NONLINEAR_LS,
    EXTERNALLY_PROVIDED,
    INVALID_COST,
} ocp_nlp_cost_t;



typedef enum
{
    CONTINUOUS_MODEL,
    DISCRETE_MODEL,
    INVALID_MODEL,
} ocp_nlp_dynamics_t;



typedef enum
{
    BGH,
    BGHP,
    INVALID_CONSTRAINT,
} ocp_nlp_constraints_t;



typedef enum
{
    NO_REGULARIZATION,
    MIRROR,
    CONVEXIFICATION,
} ocp_nlp_reg_t;


typedef struct
{
    ocp_qp_solver_plan ocp_qp_solver_plan;
    sim_solver_plan *sim_solver_plan;
    ocp_nlp_solver_t nlp_solver;
    ocp_nlp_reg_t regularization;
    ocp_nlp_cost_t *nlp_cost;
    ocp_nlp_dynamics_t *nlp_dynamics;
    ocp_nlp_constraints_t *nlp_constraints;
    int N;
} ocp_nlp_plan;



typedef struct
{
    ocp_nlp_config *config;
    void *dims;
    void *opts;
    void *mem;
    void *work;
} ocp_nlp_solver;


/* plan */
ocp_nlp_plan *ocp_nlp_plan_create(int N);
//
void ocp_nlp_plan_destroy(void* plan_);

/* config */
ocp_nlp_config *ocp_nlp_config_create(ocp_nlp_plan plan);
//
void ocp_nlp_config_destroy(void *config_);

/* dims */
ocp_nlp_dims *ocp_nlp_dims_create(void *config_);
//
void ocp_nlp_dims_destroy(void *dims_);
//

/* in */
ocp_nlp_in *ocp_nlp_in_create(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void ocp_nlp_in_destroy(void *in);
//
void ocp_nlp_in_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in, int stage,
        const char *field, void *value);
//
int ocp_nlp_dynamics_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                           int stage, const char *fun_type, void *fun_ptr);
//
// TODO remove and use ocp_nlp_dynamics_model_set instead !!!
int nlp_set_discrete_model_in_stage(ocp_nlp_config *config, ocp_nlp_in *in, int stage,
                                    void *fun_ptr);
//
int ocp_nlp_cost_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_in *in,
                           int stage, const char *field, void *value);
//
int ocp_nlp_constraints_model_set(ocp_nlp_config *config, ocp_nlp_dims *dims,
                                     ocp_nlp_in *in, int stage, const char *field, void *value);

/* out */
//
ocp_nlp_out *ocp_nlp_out_create(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void ocp_nlp_out_destroy(void *out);
//
void ocp_nlp_out_set(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
                     int stage, const char *field, void *value);
//
void ocp_nlp_out_get(ocp_nlp_config *config, ocp_nlp_dims *dims, ocp_nlp_out *out,
                     int stage, const char *field, void *value);

/* opts */
//
void *ocp_nlp_opts_create(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void ocp_nlp_opts_destroy(void *opts);
//
void ocp_nlp_opts_set(ocp_nlp_config *config, void *opts_,
                      const char *field, const void* value);
//
int ocp_nlp_dynamics_opts_set(ocp_nlp_config *config, void *opts_, int stage,
                                         const char *field, void *value);
//
void ocp_nlp_opts_update(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_);


/* solver */
//
ocp_nlp_solver *ocp_nlp_solver_create(ocp_nlp_config *config, ocp_nlp_dims *dims, void *opts_);
//
void ocp_nlp_solver_destroy(void *solver);
//
int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out);
//
int ocp_nlp_precompute(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out);

/* get */
void ocp_nlp_get(ocp_nlp_config *config, ocp_nlp_solver *solver,
                 const char *field, void *return_value_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // INTERFACES_ACADOS_C_OCP_NLP_INTERFACE_H_
