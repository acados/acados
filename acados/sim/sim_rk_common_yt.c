// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
// #include "sim/sim_collocation.h"
#include "acados/sim/sim_rk_common_yt.h"


int_t sim_RK_opts_calculate_size(int_t ns)
{
    // int_t size = sizeof(Newton_scheme);
    int_t size = sizeof(sim_RK_opts);

    size += ns * ns * sizeof(real_t); // A_mat
    size += ns * sizeof(real_t); // b_vec
    size += ns * sizeof(real_t); // c_vec

    size = (size + 63) / 64 * 64;
    size += 1 * 64;

    return size;
}

char *assign_sim_RK_opts(int_t ns, sim_RK_opts **opts, void *ptr)
    {
    char *c_ptr = (char *) ptr;

    *opts = (sim_RK_opts *) c_ptr;
    c_ptr += sizeof(sim_RK_opts);

    (*opts)->num_stages = ns;

    // c_ptr += sizeof(Newton_scheme);

    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    (*opts)->A_mat = (real_t *) c_ptr;
    c_ptr += ns*ns*sizeof(real_t);

    (*opts)->b_vec = (real_t *) c_ptr;
    c_ptr += ns*sizeof(real_t);

    (*opts)->c_vec = (real_t *) c_ptr;
    c_ptr += ns*sizeof(real_t);

    return c_ptr;
}

sim_RK_opts *create_sim_RK_opts(int_t ns) {

    sim_RK_opts *opts;

    int_t bytes = sim_RK_opts_calculate_size(ns);

    void *ptr = malloc(bytes);

    char *ptr_end = assign_sim_RK_opts(ns, &opts, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    return opts;
}