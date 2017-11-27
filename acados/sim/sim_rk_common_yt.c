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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
// #include "sim/sim_collocation.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_rk_common_yt.h"


int sim_RK_opts_calculate_size(sim_dims *dims)
{
    // int size = sizeof(Newton_scheme);
    int size = sizeof(sim_RK_opts);

    int ns = dims->num_stages;
    size += ns * ns * sizeof(double);  // A_mat
    size += ns * sizeof(double);  // b_vec
    size += ns * sizeof(double);  // c_vec

    size = (size + 63) / 64 * 64;
    size += 1 * 64;

    return size;
}



void *assign_sim_RK_opts(sim_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_RK_opts *opts = (sim_RK_opts *) c_ptr;
    c_ptr += sizeof(sim_RK_opts);

    int ns = dims->num_stages;
    opts->num_stages = ns;

    // c_ptr += sizeof(Newton_scheme);

    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    opts->A_mat = (double *) c_ptr;
    c_ptr += ns*ns*sizeof(double);

    opts->b_vec = (double *) c_ptr;
    c_ptr += ns*sizeof(double);

    opts->c_vec = (double *) c_ptr;
    c_ptr += ns*sizeof(double);

    assert((char*)raw_memory + sim_RK_opts_calculate_size(dims) >= c_ptr);

    return (void *)opts;
}



void sim_rk_initialize_default_args(void *opts_)
{
    sim_RK_opts *opts = (sim_RK_opts *) opts_;
    int ns = opts->num_stages;

    assert(opts->num_stages == 4 && "only number of stages = 4 implemented!");

    memcpy(opts->A_mat,
        ((real_t[]){0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0}),
        sizeof(*opts->A_mat) * (ns * ns));
    memcpy(opts->b_vec, ((real_t[]){1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6}),
        sizeof(*opts->b_vec) * (ns));
    memcpy(opts->c_vec, ((real_t[]){0.0, 0.5, 0.5, 1.0}),
        sizeof(*opts->c_vec) * (ns));
}



void *create_sim_RK_opts(sim_dims *dims)
{
    int bytes = sim_RK_opts_calculate_size(dims);

    void *ptr = malloc(bytes);

    sim_RK_opts *opts = assign_sim_RK_opts(dims, ptr);

    sim_rk_initialize_default_args(opts);

    return (void *)opts;
}
