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

#include "acados_c/sim/sim_lifted_irk_integrator.h"

#include <assert.h>
#include <string.h>

// void *sim_lifted_irk_integrator_copy_args(sim_solver_config *config, sim_dims *dims, void *raw_memory, void *source_)
// {
//     sim_lifted_irk_integrator_args *source = (sim_lifted_irk_integrator_args *) source_;
//     sim_lifted_irk_integrator_args *dest;

//     sim_lifted_irk_integrator_submodules submodules = {
//         source->forward_vde,
//         source->jacobian_ode
//     };

//     dest = sim_lifted_irk_integrator_assign_args(dims, &submodules, raw_memory);

//     dest->interval = source->interval;
//     dest->num_stages = source->num_stages;
//     dest->num_steps = source->num_steps;
//     dest->num_forw_sens = source->num_forw_sens;
//     dest->sens_forw = source->sens_forw;
//     dest->sens_adj = source->sens_adj;
//     dest->sens_hess = source->sens_hess;
//     dest->newton_iter = source->newton_iter;

//     int ns = dims->num_stages;

//     memcpy(dest->A_mat, source->A_mat, ns*ns*sizeof(double));
//     memcpy(dest->c_vec, source->c_vec, ns*sizeof(double));
//     memcpy(dest->b_vec, source->b_vec, ns*sizeof(double));

//     dest->scheme->type = source->scheme->type;
//     dest->scheme->single = source->scheme->single;
//     dest->scheme->freeze = source->scheme->freeze;
//     memcpy(dest->scheme->eig, source->scheme->eig, ns*sizeof(double));
//     memcpy(dest->scheme->transf1, source->scheme->transf1, ns*ns*sizeof(double));
//     memcpy(dest->scheme->transf2, source->scheme->transf2, ns*ns*sizeof(double));
//     memcpy(dest->scheme->transf1_T, source->scheme->transf1_T, ns*ns*sizeof(double));
//     memcpy(dest->scheme->transf2_T, source->scheme->transf2_T, ns*ns*sizeof(double));
//     dest->scheme->low_tria = source->scheme->low_tria;

//     extern external_function_dims sim_lifted_irk_forward_vde_dims;
//     extern external_function_dims sim_lifted_irk_jacobian_ode_dims;
//     external_function_copy_args(&config->extfun, &sim_lifted_irk_forward_vde_dims, dest->forward_vde_args, source->forward_vde_args);
//     external_function_copy_args(&config->extfun, &sim_lifted_irk_jacobian_ode_dims, dest->jacobian_ode_args, source->jacobian_ode_args);

//     return (void *)dest;
// }



int sim_lifted_irk_integrator_calculate_submodules_size(sim_solver_config *config, sim_dims *dims)
{
    int size = sizeof(sim_lifted_irk_integrator_submodules);

    extern external_function_dims sim_lifted_irk_forward_vde_dims;
    size += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_lifted_irk_forward_vde_dims);

    extern external_function_dims sim_lifted_irk_jacobian_ode_dims;
    size += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_lifted_irk_jacobian_ode_dims);

    return size;
}



void *sim_lifted_irk_integrator_assign_submodules(sim_solver_config *config, sim_dims *dims, void *raw_memory)
{
    sim_lifted_irk_integrator_submodules *submodules;

    char *c_ptr = (char *) raw_memory;

    submodules = (sim_lifted_irk_integrator_submodules *) c_ptr;
    c_ptr += sizeof(sim_lifted_irk_integrator_submodules);

    extern external_function_dims sim_lifted_irk_forward_vde_dims;
    submodules->forward_vde = assign_external_function_fcn_ptrs(&config->extfun, &sim_lifted_irk_forward_vde_dims, c_ptr);
    c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_lifted_irk_forward_vde_dims);

    extern external_function_dims sim_lifted_irk_jacobian_ode_dims;
    submodules->jacobian_ode = assign_external_function_fcn_ptrs(&config->extfun, &sim_lifted_irk_jacobian_ode_dims, c_ptr);
    c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_lifted_irk_jacobian_ode_dims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + sim_lifted_irk_integrator_calculate_submodules_size(config, dims) == c_ptr);

    return (void *)submodules;
}
