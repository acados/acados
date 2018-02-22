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

#include <string.h>



void *sim_erk_copy_opts(sim_dims *dims, void *raw_memory, void *source_)
{
    sim_rk_opts *source = (sim_rk_opts *)source_;
    sim_rk_opts *dest;

    dest = sim_erk_opts_assign(dims, raw_memory);

    dest->interval = source->interval;
    dest->num_stages = source->num_stages;
    dest->num_steps = source->num_steps;
    dest->num_forw_sens = source->num_forw_sens;
    dest->sens_forw = source->sens_forw;
    dest->sens_adj = source->sens_adj;
    dest->sens_hess = source->sens_hess;
    dest->newton_iter = source->newton_iter;

    int ns = dims->num_stages;

    memcpy(dest->A_mat, source->A_mat, ns*ns*sizeof(double));
    memcpy(dest->c_vec, source->c_vec, ns*sizeof(double));
    memcpy(dest->b_vec, source->b_vec, ns*sizeof(double));

    dest->scheme = NULL;

    return (void *)dest;
}
