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

#include "acados/sim/sim_collocation_utils.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "acados/utils/external_function_generic.h"



typedef struct {
    int num_stages;
    int nx;
    int nu;
} sim_dims;



typedef struct {
    int nx;   // NX
    int nu;   // NU

    // int nz;   // ALGEBRAIC VARIABLES: currently only internal, similar to ACADO code generation
    double *x;  // x[NX]
    double *u;  // u[NU]

    double *S_forw;  // forward seed
    double *S_adj;   // backward seed

	/* external functions */

	// explicit ode
	external_function_generic *ode_expl;
	// jacobian explicit ode
	external_function_generic *jac_ode_expl;
	// hessian explicit ode
	external_function_generic *hess_ode_expl;
	// forward explicit vde
	external_function_generic *forw_vde_expl;
	// adjoint explicit vde
	external_function_generic *adj_vde_expl;

	// implicit ode
	external_function_generic *ode_impl;
	// jac_x implicit ode
	external_function_generic *jac_x_ode_impl;
	// jac_xdot implicit ode
	external_function_generic *jac_xdot_ode_impl;
	// jac_u implicit ode
	external_function_generic *jac_u_ode_impl;

    double step;

} sim_in;



typedef struct
{
    double CPUtime; // in seconds
    double LAtime; // in seconds
    double ADtime; // in seconds
} sim_info;



typedef struct
{
    double *xn;      // xn[NX]
    double *S_forw;  // S_forw[NX*(NX+NU)]
    double *S_adj;   //
    double *S_hess;  //

    double *grad;  // gradient correction

    sim_info *info;
} sim_out;



typedef struct
{

    double interval;
    int num_stages;

    int num_steps;
    int num_forw_sens;

    double *A_mat;
    double *c_vec;
    double *b_vec;

    bool sens_forw;
    bool sens_adj;
    bool sens_hess;

    // for explicit integrators: newton_iter == 0 && scheme == NULL 
    // && jac_reuse=false
    int newton_iter;
    bool jac_reuse;
    Newton_scheme *scheme;

	// work space
	void *work;

} sim_rk_opts;



typedef struct {
    int (*fun)(sim_in *in, sim_out *out, void *args, void *mem, void *work);
    int (*calculate_args_size)(sim_dims *dims);
    void *(*assign_args)(sim_dims *dims, void *raw_memory);
    void (*initialize_default_args)(sim_dims *dims, void *args);
    int (*calculate_memory_size)(sim_dims *dims, void *args);
    void *(*assign_memory)(sim_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(sim_dims *dims, void *args);
} sim_solver_fcn_ptrs;



int sim_dims_calculate_size();

sim_dims *assign_sim_dims(void *raw_memory);

int sim_in_calculate_size(sim_dims *dims);

sim_in *assign_sim_in(sim_dims *dims, void *raw_memory);

int sim_out_calculate_size(sim_dims *dims);

sim_out *assign_sim_out(sim_dims *dims, void *raw_memory);

#endif  // ACADOS_SIM_SIM_COMMON_H_
