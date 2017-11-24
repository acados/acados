// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
#include "acados/sim/sim_common_yt.h"

int_t sim_in_calculate_size(int_t nx, int_t nu, int_t NF)
{

int_t size = sizeof(sim_in);

size += nx * sizeof(real_t); // x
size += nu * sizeof(real_t); // u
size += nx * NF * sizeof(real_t); // S_forw
size += nx * sizeof(real_t); // S_adj

size = (size + 63) / 64 * 64;
size += 1 * 64;

return size;
}

char *assign_sim_in(int_t nx, int_t nu, int_t NF, sim_in **in, void *ptr)
{

char *c_ptr = (char *) ptr;

*in = (sim_in *) c_ptr;
c_ptr += sizeof(sim_in);

(*in)->nx = nx;
(*in)->nu = nu;
(*in)->NF = NF;

size_t s_ptr = (size_t)c_ptr;
s_ptr = (s_ptr + 63) / 64 * 64;
c_ptr = (char *)s_ptr;

(*in)->x = (real_t *) c_ptr;
c_ptr += nx*sizeof(real_t);

(*in)->u = (real_t *) c_ptr;
c_ptr += nu*sizeof(real_t);

(*in)->S_forw = (real_t *) c_ptr;
c_ptr += nx * NF *sizeof(real_t);

(*in)->S_adj = (real_t *) c_ptr;
c_ptr += nx *sizeof(real_t);

return c_ptr;
}

sim_in *create_sim_in(int_t nx, int_t nu, int_t NF) {

sim_in *in;

int_t bytes = sim_in_calculate_size(nx, nu, NF);

void *ptr = malloc(bytes);

char *ptr_end = assign_sim_in(nx, nu, NF, &in, ptr);
assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

return in;
}

int_t sim_out_calculate_size(int_t nx, int_t nu, int_t NF)
{

int_t size = sizeof(sim_out);

size += sizeof(sim_info);

size += nx * sizeof(real_t); // xn
size += nx * NF * sizeof(real_t); // S_forw
size += (nx + nu) * sizeof(real_t); // S_adj
size += ((NF + 1) * NF / 2) * sizeof(real_t); // S_hess

size = (size + 63) / 64 * 64;
size += 1 * 64;

return size;
}

char *assign_sim_out(int_t nx, int_t nu, int_t NF, sim_out **out, void *ptr)
{

char *c_ptr = (char *) ptr;

*out = (sim_out *) c_ptr;
c_ptr += sizeof(sim_out);

(*out)->info = (sim_info *)c_ptr;
c_ptr += sizeof(sim_info);

size_t s_ptr = (size_t)c_ptr;
s_ptr = (s_ptr + 63) / 64 * 64;
c_ptr = (char *)s_ptr;

(*out)->xn = (real_t *) c_ptr;
c_ptr += nx*sizeof(real_t);

(*out)->S_forw = (real_t *) c_ptr;
c_ptr += nx * NF *sizeof(real_t);

(*out)->S_adj = (real_t *) c_ptr;
c_ptr += (nx + nu) *sizeof(real_t);

(*out)->S_hess = (real_t *) c_ptr;
c_ptr += ((NF + 1) * NF / 2) *sizeof(real_t);

return c_ptr;
}

sim_out *create_sim_out(int_t nx, int_t nu, int_t NF) {

sim_out *out;

int_t bytes = sim_out_calculate_size(nx, nu, NF);

void *ptr = malloc(bytes);

char *ptr_end = assign_sim_out(nx, nu, NF, &out, ptr);
assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

return out;
}
