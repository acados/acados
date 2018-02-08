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

#include "acados_c/sim/sim_erk_integrator.h"

#include <assert.h>
#include <string.h>

#include "acados_c/external_function.h"



// void *sim_erk_integrator_copy_args(sim_solver_config *config, sim_dims *dims, void *raw_memory, void *source_)
// {
//     sim_erk_integrator_args *source = (sim_erk_integrator_args *)source_;
//     sim_erk_integrator_args *dest;

//     sim_erk_integrator_submodules submodules = {
//         source->forward_vde,
//         source->adjoint_vde,
//         source->hess_vde
//     };

//     dest = sim_erk_integrator_assign_args(dims, &submodules, raw_memory);

//     dest->interval = source->interval;
//     dest->num_stages = source->num_stages;
//     dest->num_steps = source->num_steps;
//     dest->num_forw_sens = source->num_forw_sens;
//     dest->sens_forw = source->sens_forw;
//     dest->sens_adj = source->sens_adj;
//     dest->sens_hess = source->sens_hess;

//     int ns = dims->num_stages;

//     memcpy(dest->A_mat, source->A_mat, ns*ns*sizeof(double));
//     memcpy(dest->c_vec, source->c_vec, ns*sizeof(double));
//     memcpy(dest->b_vec, source->b_vec, ns*sizeof(double));

//     extern external_function_dims sim_erk_forward_vde_dims;
//     extern external_function_dims sim_erk_adjoint_vde_dims;
//     extern external_function_dims sim_erk_hess_vde_dims;
//     external_function_copy_args(&config->extfun, &sim_erk_forward_vde_dims, dest->forward_vde_args, source->forward_vde_args);
//     external_function_copy_args(&config->extfun, &sim_erk_adjoint_vde_dims, dest->adjoint_vde_args, source->adjoint_vde_args);
//     external_function_copy_args(&config->extfun, &sim_erk_hess_vde_dims, dest->hess_vde_args, source->hess_vde_args);

//     return (void *)dest;
// }



int sim_erk_integrator_calculate_submodules_size(sim_solver_config *config, sim_dims *dims)
{
    int size = sizeof(sim_erk_integrator_submodules);

    extern external_function_dims sim_erk_forward_vde_dims;
    size += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_erk_forward_vde_dims);

    extern external_function_dims sim_erk_adjoint_vde_dims;
    size += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_erk_adjoint_vde_dims);

    extern external_function_dims sim_erk_hess_vde_dims;
    size += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_erk_hess_vde_dims);

    return size;
}



void *sim_erk_integrator_assign_submodules(sim_solver_config *config, sim_dims *dims, void *raw_memory)
{
    sim_erk_integrator_submodules *submodules;

    char *c_ptr = (char *) raw_memory;

    submodules = (sim_erk_integrator_submodules *) c_ptr;
    c_ptr += sizeof(sim_erk_integrator_submodules);

    extern external_function_dims sim_erk_forward_vde_dims;
    submodules->forward_vde = assign_external_function_fcn_ptrs(&config->extfun, &sim_erk_forward_vde_dims, c_ptr);
    c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_erk_forward_vde_dims);

    extern external_function_dims sim_erk_adjoint_vde_dims;
    submodules->adjoint_vde = assign_external_function_fcn_ptrs(&config->extfun, &sim_erk_adjoint_vde_dims, c_ptr);
    c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_erk_adjoint_vde_dims);

    extern external_function_dims sim_erk_hess_vde_dims;
    submodules->hess_vde = assign_external_function_fcn_ptrs(&config->extfun, &sim_erk_hess_vde_dims, c_ptr);
    c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_erk_hess_vde_dims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + sim_erk_integrator_calculate_submodules_size(config, dims) == c_ptr);

    return (void *)submodules;
}
