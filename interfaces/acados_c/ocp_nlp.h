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

#ifndef ACADOS_C_OCP_NLP_H_
#define ACADOS_C_OCP_NLP_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/ocp_nlp/ocp_nlp_common.h>
#include <acados/utils/types.h>
// acados_c
#include <acados_c/sim.h>
#include <acados_c/ocp_qp.h>

typedef enum {
    SQP_GN
} ocp_nlp_solver_t;

typedef struct {
    ocp_qp_config *ocp_qp_config;
    sim_config **sim_config;
} ocp_nlp_config;

typedef struct {
    ocp_nlp_solver_fcn_ptrs *fcn_ptrs;
    void *dims;
    void *args;
    void *mem;
    void *work;
} ocp_nlp_solver;

// INPUT AND OUTPUT
//
ocp_nlp_in *create_ocp_nlp_in(ocp_nlp_dims *dims);
//
ocp_nlp_out *create_ocp_nlp_out(ocp_nlp_dims *dims);

// BASIC INTERFACE
//
int ocp_nlp_calculate_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
ocp_nlp_solver *ocp_nlp_assign(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory);
//
ocp_nlp_solver *ocp_nlp_create(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
int ocp_nlp_solve(ocp_nlp_solver *solver, ocp_nlp_in *qp_in, ocp_nlp_out *qp_out);
//
void ocp_nlp_initialize_default_args(ocp_nlp_solver *solver);

// EXPERT INTERFACE
//
//
int ocp_nlp_calculate_args_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void *ocp_nlp_assign_args(ocp_nlp_config  *config, ocp_nlp_dims *dims, void *raw_memory);
//
void *ocp_nlp_create_args(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
int ocp_nlp_calculate_memory_size(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
void *ocp_nlp_assign_memory(ocp_nlp_config *config, ocp_nlp_dims *dims, void *raw_memory);
//
void *ocp_nlp_create_memory(ocp_nlp_config *config, ocp_nlp_dims *dims);
//
int ocp_nlp_calculate_workspace_size(ocp_nlp_config *config, ocp_nlp_dims *dims);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_OCP_NLP_H_