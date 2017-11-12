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

#include "acados/sim/allocate_sim.h"

#include <stdlib.h>

void allocate_sim_in(sim_in *sim_in, int_t num_stages) {
    int_t nx = sim_in->nx;
    int_t nu = sim_in->nu;
    sim_in->x = calloc(nx, sizeof(*sim_in->x));
    sim_in->u = calloc(nu, sizeof(*sim_in->u));
    sim_in->S_forw = calloc(nx * sim_in->num_forw_sens, sizeof(*sim_in->S_forw));
    sim_in->S_adj = calloc(nx + nu, sizeof(*sim_in->S_adj));
    sim_in->grad_K = calloc(nx * num_stages, sizeof(*sim_in->grad_K));

}

void allocate_sim_out(sim_in *sim_in, sim_out *sim_out) {
    int_t nx = sim_in->nx;
    int_t nu = sim_in->nu;
    sim_out->xn = calloc(nx, sizeof(*sim_out->xn));
    sim_out->S_forw = calloc(nx * (nx + nu), sizeof(*sim_out->S_forw));
    sim_out->S_adj = calloc(nx + nu, sizeof(*sim_out->S_adj));
    int_t nhess = (nx + nu + 1) * (nx + nu) / 2;
    sim_out->S_hess = calloc(nhess, sizeof(*sim_out->S_hess));
    sim_out->grad = calloc(sim_in->num_forw_sens, sizeof(*sim_out->grad));
    sim_out->info = malloc(sizeof(sim_info));
}
