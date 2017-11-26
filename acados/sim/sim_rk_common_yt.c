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
#include "acados/sim/sim_rk_common_yt.h"


int sim_RK_opts_calculate_size(int ns)
{
    // int size = sizeof(Newton_scheme);
    int size = sizeof(sim_RK_opts);

    size += ns * ns * sizeof(double);  // A_mat
    size += ns * sizeof(double);  // b_vec
    size += ns * sizeof(double);  // c_vec

    size = (size + 63) / 64 * 64;
    size += 1 * 64;

    return size;
}



sim_RK_opts *assign_sim_RK_opts(int ns, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_RK_opts *opts = (sim_RK_opts *) c_ptr;
    c_ptr += sizeof(sim_RK_opts);

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

    assert((char*)raw_memory + sim_RK_opts_calculate_size(ns) >= c_ptr);

    return opts;
}



void sim_rk_initialize_default_args(sim_RK_opts *opts, int ns)
{
    assert(ns == 4 && "only number of stages = 4 implemented!");

    memcpy(opts->A_mat,
        ((real_t[]){0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0}),
        sizeof(*opts->A_mat) * (ns * ns));
    memcpy(opts->b_vec, ((real_t[]){1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6}),
        sizeof(*opts->b_vec) * (ns));
    memcpy(opts->c_vec, ((real_t[]){0.0, 0.5, 0.5, 1.0}),
        sizeof(*opts->c_vec) * (ns));
}



sim_RK_opts *create_sim_RK_opts(int ns)
{
    int bytes = sim_RK_opts_calculate_size(ns);

    void *ptr = malloc(bytes);

    sim_RK_opts *opts = assign_sim_RK_opts(ns, ptr);

    sim_rk_initialize_default_args(opts, ns);

    return opts;
}
