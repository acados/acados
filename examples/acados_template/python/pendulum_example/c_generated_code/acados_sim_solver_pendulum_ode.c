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

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados/ocp_nlp/ocp_nlp_cost_ls.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

// example specific

#include "pendulum_ode_model/pendulum_ode_model.h"
#include "acados_sim_solver_pendulum_ode.h"

#define NX_    4
#define NZ_    0
#define NU_    1
#define NP_    0

#define NX   NX_
#define NZ   NZ_
#define NU   NU_
#define NP   NP_

sim_config  * pendulum_ode_sim_config;
sim_in      * pendulum_ode_sim_in;
sim_out     * pendulum_ode_sim_out;
void        * pendulum_ode_sim_dims;
sim_opts    * pendulum_ode_sim_opts;
sim_solver  * pendulum_ode_sim_solver;

// ** global data **



external_function_casadi * sim_forw_vde_casadi;
external_function_casadi * sim_expl_ode_fun_casadi;





int pendulum_ode_acados_sim_create() {

	// initialize

    int ii;
    int jj;

    int nx = NX;
    int nu = NU;
    int nz = NZ;

    double Td = 2.0/ 50;

    
    // explicit ode
    
    sim_forw_vde_casadi = (external_function_casadi *) malloc(sizeof(external_function_casadi));
    sim_expl_ode_fun_casadi = (external_function_casadi *) malloc(sizeof(external_function_casadi));
    

    sim_forw_vde_casadi->casadi_fun = &pendulum_ode_expl_vde_forw;
    sim_forw_vde_casadi->casadi_n_in = &pendulum_ode_expl_vde_forw_n_in;
    sim_forw_vde_casadi->casadi_n_out = &pendulum_ode_expl_vde_forw_n_out;
    sim_forw_vde_casadi->casadi_sparsity_in = &pendulum_ode_expl_vde_forw_sparsity_in;
    sim_forw_vde_casadi->casadi_sparsity_out = &pendulum_ode_expl_vde_forw_sparsity_out;
    sim_forw_vde_casadi->casadi_work = &pendulum_ode_expl_vde_forw_work;
    external_function_casadi_create(sim_forw_vde_casadi);

    sim_expl_ode_fun_casadi->casadi_fun = &pendulum_ode_expl_ode_fun;
    sim_expl_ode_fun_casadi->casadi_n_in = &pendulum_ode_expl_ode_fun_n_in;
    sim_expl_ode_fun_casadi->casadi_n_out = &pendulum_ode_expl_ode_fun_n_out;
    sim_expl_ode_fun_casadi->casadi_sparsity_in = &pendulum_ode_expl_ode_fun_sparsity_in;
    sim_expl_ode_fun_casadi->casadi_sparsity_out = &pendulum_ode_expl_ode_fun_sparsity_out;
    sim_expl_ode_fun_casadi->casadi_work = &pendulum_ode_expl_ode_fun_work;
    external_function_casadi_create(sim_expl_ode_fun_casadi);
    

    // sim plan & config

    // choose plan
    sim_solver_plan plan;

    plan.sim_solver = ERK;

    // create correct config based on plan
    pendulum_ode_sim_config = sim_config_create(plan);

    // sim dims

    pendulum_ode_sim_dims = sim_dims_create(pendulum_ode_sim_config);
    sim_dims_set(pendulum_ode_sim_config, pendulum_ode_sim_dims, "nx", &nx);
    sim_dims_set(pendulum_ode_sim_config, pendulum_ode_sim_dims, "nu", &nu);
    sim_dims_set(pendulum_ode_sim_config, pendulum_ode_sim_dims, "nz", &nz);

    // sim opts

    pendulum_ode_sim_opts = sim_opts_create(pendulum_ode_sim_config, pendulum_ode_sim_dims);

    pendulum_ode_sim_opts->ns = 1; // number of stages in rk integrator
    pendulum_ode_sim_opts->num_steps = 1; // number of integration steps
    pendulum_ode_sim_opts->sens_adj = false;
    pendulum_ode_sim_opts->sens_forw = true;
    

    // sim in / out

    pendulum_ode_sim_in  = sim_in_create(pendulum_ode_sim_config, pendulum_ode_sim_dims);
    pendulum_ode_sim_out = sim_out_create(pendulum_ode_sim_config, pendulum_ode_sim_dims);

    pendulum_ode_sim_in->T = Td;

    // external functions
    
     
    pendulum_ode_sim_config->model_set(pendulum_ode_sim_in->model, "expl_vde_for", sim_forw_vde_casadi);
    pendulum_ode_sim_config->model_set(pendulum_ode_sim_in->model, "expl_ode_fun", sim_expl_ode_fun_casadi);
    
    

    // sim solver

    pendulum_ode_sim_solver = sim_solver_create(pendulum_ode_sim_config, pendulum_ode_sim_dims, pendulum_ode_sim_opts);
    
    // initialize state and input to zero
    // x
    for (ii = 0; ii < NX; ii++)
        pendulum_ode_sim_in->x[ii] = 0.0;
    
    // u
    for (ii = 0; ii < NU; ii++)
        pendulum_ode_sim_in->u[ii] = 0.0;
    int status = 0;

    return status;
}

int pendulum_ode_acados_sim_solve() {

    // integrate dynamics using acados sim_solver 
    int status = sim_solve(pendulum_ode_sim_solver, pendulum_ode_sim_in, pendulum_ode_sim_out);
    if (status != 0)
        printf("error in pendulum_ode_acados_sim_solve()! Exiting.\n");

    return status;
}

int pendulum_ode_acados_sim_free() {

    // free memory
    sim_solver_destroy(pendulum_ode_sim_solver);
    sim_in_destroy(pendulum_ode_sim_in);
    sim_out_destroy(pendulum_ode_sim_out);
    sim_opts_destroy(pendulum_ode_sim_opts);
    sim_dims_destroy(pendulum_ode_sim_dims);
    sim_config_destroy(pendulum_ode_sim_config);

    // free external function 
    
        
        external_function_casadi_free(sim_forw_vde_casadi);
        external_function_casadi_free(sim_expl_ode_fun_casadi);
        
    
    
    return 0;
}

sim_config  * pendulum_ode_acados_get_sim_config() {  
    return pendulum_ode_sim_config; };

sim_in      * pendulum_ode_acados_get_sim_in(){       
    return pendulum_ode_sim_in; };

sim_out     * pendulum_ode_acados_get_sim_out(){      
    return pendulum_ode_sim_out; };

void        * pendulum_ode_acados_get_sim_dims(){     
    return pendulum_ode_sim_dims; };

sim_opts    * pendulum_ode_acados_get_sim_opts(){     
    return pendulum_ode_sim_opts; };

sim_solver  * pendulum_ode_acados_get_sim_solver(){   
    return pendulum_ode_sim_solver; };
