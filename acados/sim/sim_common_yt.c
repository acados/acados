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
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_common_yt.h"


// TODO(dimitris): write sim_dims
int sim_in_calculate_size(int nx, int nu, int NF)
{
    int size = sizeof(sim_in);

    size += nx * sizeof(double);  // x
    size += nu * sizeof(double);  // u
    size += nx * NF * sizeof(double);  // S_forw
    size += nx * sizeof(double);  // S_adj

    size = (size + 63) / 64 * 64;
    size += 1 * 64;

    return size;
}



sim_in *assign_sim_in(int nx, int nu, int NF, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_in *in = (sim_in *) c_ptr;
    c_ptr += sizeof(sim_in);

    in->nx = nx;
    in->nu = nu;
    in->NF = NF;

    // replace with mem.c functions
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    in->x = (double *) c_ptr;
    c_ptr += nx*sizeof(double);

    in->u = (double *) c_ptr;
    c_ptr += nu*sizeof(double);

    in->S_forw = (double *) c_ptr;
    c_ptr += nx * NF *sizeof(double);

    in->S_adj = (double *) c_ptr;
    c_ptr += nx *sizeof(double);

    assert((char*)raw_memory + sim_in_calculate_size(nx, nu, NF) >= c_ptr);

    return in;
}



// TODO(dimitris): move to create.c
sim_in *create_sim_in(int nx, int nu, int NF)
{
    int bytes = sim_in_calculate_size(nx, nu, NF);

    void *ptr = acados_malloc(bytes, 1);

    sim_in *in = assign_sim_in(nx, nu, NF, ptr);

    return in;
}



int sim_out_calculate_size(int nx, int nu, int NF)
{
    int size = sizeof(sim_out);

    size += sizeof(sim_info);

    size += nx * sizeof(double);  // xn
    size += nx * NF * sizeof(double);  // S_forw
    size += (nx + nu) * sizeof(double);  // S_adj
    size += ((NF + 1) * NF / 2) * sizeof(double);  // S_hess

    size = (size + 63) / 64 * 64;
    size += 1 * 64;

    return size;
}



sim_out *assign_sim_out(int nx, int nu, int NF, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    sim_out *out = (sim_out *) c_ptr;
    c_ptr += sizeof(sim_out);

    out->info = (sim_info *)c_ptr;
    c_ptr += sizeof(sim_info);

    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    out->xn = (double *) c_ptr;
    c_ptr += nx*sizeof(double);

    out->S_forw = (double *) c_ptr;
    c_ptr += nx * NF *sizeof(double);

    out->S_adj = (double *) c_ptr;
    c_ptr += (nx + nu) *sizeof(double);

    out->S_hess = (double *) c_ptr;
    c_ptr += ((NF + 1) * NF / 2) *sizeof(double);

    assert((char*)raw_memory + sim_out_calculate_size(nx, nu, NF) >= c_ptr);

    return out;
}



sim_out *create_sim_out(int nx, int nu, int NF)
{
    int bytes = sim_out_calculate_size(nx, nu, NF);

    void *ptr = malloc(bytes);

    sim_out *out = assign_sim_out(nx, nu, NF, ptr);

    return out;
}
