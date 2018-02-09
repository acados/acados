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

#ifndef ACADOS_C_SIM_H_
#define ACADOS_C_SIM_H_

#ifdef __cplusplus
extern "C" {
#endif

// acados
#include <acados/sim/sim_common.h>
#include <acados/utils/types.h>
// acados_c
#include "acados_c/common.h"
#include "acados_c/external_function.h"



typedef enum {
    PREVIOUS,
    ERK,
    IRK,
    LIFTED_IRK
} sim_solver_t;



typedef struct {
    sim_solver_t sim_solver;
    external_function_config extfun;
} sim_solver_config;



typedef struct {
    sim_solver_fcn_ptrs *fcn_ptrs;
    void *dims;
    void *args;
    void *mem;
    void *work;
} sim_solver;



// INPUT, OUTPUT AND OPTIONS
//
sim_dims *create_sim_dims();
//
sim_in *create_sim_in(sim_dims *dims);
//
sim_out *create_sim_out(sim_dims *dims);
//
int sim_calculate_args_size(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims);
//
void *sim_assign_args(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *raw_memory);
//
void *sim_create_args(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims);
//
void *sim_copy_args(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *raw_memory, void *source);

// BASIC INTERFACE
//
int sim_calculate_size(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *args_);
//
sim_solver *sim_assign(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *args_, void *raw_memory);
//
sim_solver *sim_create(sim_solver_fcn_ptrs *fcn_ptrs, sim_dims *dims, void *args_);
//
int sim_solve(sim_solver *solver, sim_in *qp_in, sim_out *qp_out);

// OPTIONS BASED CONFIGURATION STRATEGY
//
int sim_calculate_submodules_size(sim_solver_config *config, sim_dims *dims);
//
void *sim_assign_submodules(sim_solver_config *config, sim_dims *dims, void *raw_memory);
//
int calculate_sim_solver_fcn_ptrs_size(sim_solver_config *config, sim_dims *dims);
//
void *assign_sim_solver_fcn_ptrs(sim_solver_config *config, sim_dims *dims, void *raw_memory);
//
void *create_sim_solver_fcn_ptrs(sim_solver_config *config, sim_dims *dims);

// EXPERT INTERFACE
//
int set_sim_solver_fcn_ptrs(sim_solver_config *config, sim_solver_fcn_ptrs *fcn_ptrs);




#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_C_SIM_H_