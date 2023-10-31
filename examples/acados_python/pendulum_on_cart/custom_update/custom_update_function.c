/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

#include <stdlib.h>
#include <stdio.h>

#include "acados_solver_pendulum.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados/utils/mem.h"

#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"


typedef struct custom_memory
{
    // struct blasfeo_dvec *ux;
    double *traj_buffer;

    void *raw_memory; // Pointer to allocated memory, to be used for freeing

} custom_memory;


static int custom_memory_calculate_size(ocp_nlp_dims *nlp_dims)
{
    int N = nlp_dims->N;
    int nx = nlp_dims->nx[0];
    int nu = nlp_dims->nu[0];
    int nxu = nx + nu;

    acados_size_t size = sizeof(custom_memory);
    size += (N+1) * (nxu) * sizeof(double); // traj_buffer

    make_int_multiple_of(8, &size);
    // printf("custom_memory_calculate_size %d\n", size);
    // printf("custom_memory_calculate_size sizeof double %d\n", sizeof(double));
    return size;
}


static custom_memory *custom_memory_assign(ocp_nlp_dims *nlp_dims, void *raw_memory)
{
    int N = nlp_dims->N;
    int nx = nlp_dims->nx[0];
    int nu = nlp_dims->nu[0];
    int nxu = nx + nu;

    char *c_ptr = (char *) raw_memory;
    custom_memory *mem = (custom_memory *) c_ptr;
    c_ptr += sizeof(custom_memory);

    align_char_to(8, &c_ptr);
    assign_and_advance_double(nxu*(N+1), &mem->traj_buffer, &c_ptr);
    // printf("custom_memory traj_buffer at %p\n", mem->traj_buffer);

    mem->raw_memory = raw_memory;

    return mem;
}



static void *example_custom_memory_create(pendulum_solver_capsule* capsule)
{
    ocp_nlp_dims *nlp_dims = pendulum_acados_get_nlp_dims(capsule);
    acados_size_t bytes = custom_memory_calculate_size(nlp_dims);

    void *ptr = acados_calloc(1, bytes);

    custom_memory *custom_mem = custom_memory_assign(nlp_dims, ptr);
    custom_mem->raw_memory = ptr;

    return custom_mem;
}


int custom_update_init_function(pendulum_solver_capsule* capsule)
{
    printf("\nin example_custom_update_init_function\n");
    capsule->custom_update_memory = example_custom_memory_create(capsule);

    return 1;
}


static void print_x_trajectory(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, custom_memory *custom_mem)
{
    ocp_nlp_config *nlp_config = solver->config;
    ocp_nlp_dims *nlp_dims = solver->dims;

    int N = nlp_dims->N;
    int nx = nlp_dims->nx[0];

    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "x", &custom_mem->traj_buffer[ii*nx]);

    printf("xtraj in custom function:\t");
    d_print_exp_mat( nx, N+1, custom_mem->traj_buffer, nx);
}


static void print_u_trajectory(ocp_nlp_solver *solver, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, custom_memory *custom_mem)
{
    ocp_nlp_config *nlp_config = solver->config;
    ocp_nlp_dims *nlp_dims = solver->dims;

    int N = nlp_dims->N;
    int nu = nlp_dims->nu[0];

    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "u", &custom_mem->traj_buffer[ii*nu]);

    printf("utraj in custom function:\t");
    d_print_exp_mat( nu, N+1, custom_mem->traj_buffer, nu);
}



int custom_update_function(pendulum_solver_capsule* capsule, double* data, int data_len)
{
    printf("\nin example_custom_update_function\n");
    custom_memory *custom_mem = (custom_memory *) capsule->custom_update_memory;
    ocp_nlp_config *nlp_config = pendulum_acados_get_nlp_config(capsule);
    ocp_nlp_dims *nlp_dims = pendulum_acados_get_nlp_dims(capsule);
    ocp_nlp_in *nlp_in = pendulum_acados_get_nlp_in(capsule);
    ocp_nlp_out *nlp_out = pendulum_acados_get_nlp_out(capsule);
    ocp_nlp_solver *nlp_solver = pendulum_acados_get_nlp_solver(capsule);
    void *nlp_opts = pendulum_acados_get_nlp_opts(capsule);

    printf("got data\t:");
    for (int i = 0; i < data_len; i++)
    {
        printf("%e\t", data[i]);
    }
    printf("\n");



    // EXAMPLE 1: print x trajectory
    // print_x_trajectory(nlp_solver, nlp_in, nlp_out, custom_mem);
    // EXAMPLE 2: print u trajectory
    print_u_trajectory(nlp_solver, nlp_in, nlp_out, custom_mem);

    return 1;
}


int custom_update_terminate_function(pendulum_solver_capsule* capsule)
{
    printf("\nin example_custom_update_terminate_function\n");

    custom_memory *mem = capsule->custom_update_memory;

    free(mem->raw_memory);
    return 1;
}
