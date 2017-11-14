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

#include "acados/sim/sim_common.h"

int sim_in_calculate_size(sim_dims *dims)
{
    int size = sizeof(sim_in);
    size += dims->nx;  // x
    size += dims->nu;  // u
    size += dims->num_forw_sens * (dims->nx + dims->nu);  // S_forw
    size += dims->nx + dims->nu;  // S_adj
    size += dims->nx * dims->num_stages;  // grad_K
    return size;
}

sim_in *sim_in_assign(sim_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_in *in = (sim_in *) c_ptr;
    c_ptr += sizeof(sim_in);

    in->x = c_ptr;
    c_ptr += dims->nx;

    in->u = c_ptr;
    c_ptr += dims->nu;

    in->S_forw = c_ptr;
    c_ptr += dims->num_forw_sens * (dims->nx + dims->nu);

    in->S_adj = c_ptr;
    c_ptr += dims->nx + dims->nu;

    in->grad_K = c_ptr;
    c_ptr += dims->nx * dims->num_stages;

    assert((char *)raw_memory + sim_in_calculate_size(dims) == c_ptr);

    return in;
}
