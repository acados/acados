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

#include "acados_c/sim/sim_irk_integrator.h"

#include <assert.h>
#include <string.h>

#include "acados_c/external_function.h"



int sim_irk_integrator_calculate_submodules_size(sim_solver_config *config, sim_dims *dims)
{
    int size = sizeof(sim_irk_integrator_submodules);

    extern external_function_dims sim_irk_impl_res_dims;
    size += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_irk_impl_res_dims);

    extern external_function_dims sim_irk_impl_jac_dims;
    size += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_irk_impl_jac_dims);

    return size;
}



void *sim_irk_integrator_assign_submodules(sim_solver_config *config, sim_dims *dims, void *raw_memory)
{
    sim_irk_integrator_submodules *submodules;

    char *c_ptr = (char *) raw_memory;

    submodules = (sim_irk_integrator_submodules *) c_ptr;
    c_ptr += sizeof(sim_irk_integrator_submodules);

    extern external_function_dims sim_irk_impl_res_dims;
    submodules->impl_res = assign_external_function_fcn_ptrs(&config->extfun, &sim_irk_impl_res_dims, c_ptr);
    c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_irk_impl_res_dims);

    extern external_function_dims sim_irk_impl_jac_dims;
    submodules->impl_jac = assign_external_function_fcn_ptrs(&config->extfun, &sim_irk_impl_jac_dims, c_ptr);
    c_ptr += calculate_external_function_fcn_ptrs_size(&config->extfun, &sim_irk_impl_jac_dims);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + sim_irk_integrator_calculate_submodules_size(config, dims) == c_ptr);

    return (void *)submodules;
}
