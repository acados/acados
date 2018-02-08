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

#ifndef ACADOS_SIM_SIM_COMMON_H_
#define ACADOS_SIM_SIM_COMMON_H_

#include <stdbool.h>

#include "acados/sim/sim_collocation.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"



typedef struct {
    int num_stages;
    int nx;
    int nu;
    int np;
} sim_dims;



typedef struct {
    int nx;   // NX
    int nu;   // NU
    int np;   // NP

    // int nz;   // ALGEBRAIC VARIABLES: currently only internal, similar to ACADO code generation
    double *x;  // x[NX]
    double *u;  // u[NU]
    double *p;  // p[NU]

    double *S_forw;  // forward seed
    double *S_adj;   // backward seed

    double step;

} sim_in;



typedef struct {
    double CPUtime;
    double LAtime;
    double ADtime;
} sim_info;



typedef struct {
    double *xn;      // xn[NX]
    double *S_forw;  // S_forw[NX*(NX+NU)]
    double *S_adj;   //
    double *S_hess;  //

    double *grad;  // gradient correction

    sim_info *info;
} sim_out;



typedef struct {
    int (*fun)(sim_in *in, sim_out *out, void *args, void *mem, void *work);
    int (*calculate_args_size)(sim_dims *dims, void *submodules);
    void *(*assign_args)(sim_dims *dims, void **submodules, void *raw_memory);
    void *(*copy_args)(sim_dims *dims, void *raw_memory, void *source);
    void (*initialize_default_args)(sim_dims *dims, void *args);
    int (*calculate_memory_size)(sim_dims *dims, void *args);
    void *(*assign_memory)(sim_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(sim_dims *dims, void *args);
    void *submodules;
} sim_solver_fcn_ptrs;



//
int sim_dims_calculate_size();
//
sim_dims *assign_sim_dims(void *raw_memory);
//
int sim_in_calculate_size(sim_dims *dims);
//
sim_in *assign_sim_in(sim_dims *dims, void *raw_memory);
//
int sim_out_calculate_size(sim_dims *dims);
//
sim_out *assign_sim_out(sim_dims *dims, void *raw_memory);



#endif  // ACADOS_SIM_SIM_COMMON_H_
