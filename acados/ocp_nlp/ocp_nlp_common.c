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

#include "acados/ocp_nlp/ocp_nlp_common.h"

// TODO(nielsvd): only perform assert in debug mode?
#include <assert.h>

#include <stdlib.h>

void size_of_memory_elements(const ocp_nlp_in *nlp_in, int_t *size_hess_l, 
                             int_t *size_grad_f, int_t *size_jac_h,
                             int_t *size_jac_g, int_t *size_h, int_t *size_g,
                             int_t *size_x, int_t *size_u, int_t *size_pi,
                             int_t *size_lam) {
    int_t N = nlp_in->N;
    int_t *nx = (int_t *)nlp_in->nx;
    int_t *nu = (int_t *)nlp_in->nu;
    int_t *ng = (int_t *)nlp_in->ng;

    *size_hess_l = 0;
    *size_grad_f = 0;
    *size_jac_h = 0;
    *size_jac_g = 0;
    *size_h = 0;
    *size_g = 0;
    *size_x = 0;
    *size_u = 0;
    *size_pi = 0;
    *size_lam = 0;

    // hess_l
    for (int_t i = 0; i <= N; i++) 
        *size_hess_l += (nx[i] + nu[i]) * (nx[i] + nu[i]);

    // grad_f
    for (int_t i = 0; i <= N; i++) 
        *size_grad_f += nx[i] + nu[i];
    
    // jac_h
    for (int_t i = 0; i < N; i++) 
        *size_jac_h += nx[i+1] * (nx[i] + nu[i]);

    // jac_g
    for (int_t i = 0; i <= N; i++) 
        *size_jac_g += ng[i] * (nx[i] + nu[i]);

    // h
    for (int_t i = 0; i < N; i++) 
        *size_h += nx[i+1];

    // g
    for (int_t i = 0; i <= N; i++) 
        *size_g += ng[i];

    // x
    for (int_t i = 0; i <= N; i++) 
        *size_x += nx[i];

    // u
    for (int_t i = 0; i <= N; i++)
        *size_u += nu[i];

    // pi
    for (int_t i = 0; i < N; i++)
        *size_pi += nx[i+1];

    // lam
    for (int_t i = 0; i <= N; i++)
        *size_lam += ng[i];

}
                            

int_t ocp_nlp_calculate_memory_size(const ocp_nlp_in *nlp_in) {

    int_t size = sizeof(ocp_nlp_memory);

    int_t size_hess_l, size_grad_f, size_jac_h,size_jac_g, size_h, size_g;
    int_t size_x, size_u, size_pi, size_lam;
    size_of_memory_elements(nlp_in, &size_hess_l, &size_grad_f, &size_jac_h,
                            &size_jac_g, &size_h, &size_g, &size_x, &size_u, 
                            &size_pi, &size_lam);

    size += size_hess_l + size_grad_f + size_jac_h + size_jac_g + size_h + size_g;
    size += size_x + size_u + size_pi + size_lam;

    return size;
}

char *ocp_nlp_assign_memory(const ocp_nlp_in *nlp_in, void **mem_, void *raw_memory) {
    ocp_nlp_memory **nlp_memory = (ocp_nlp_memory **)mem_;
    char *c_ptr = (char *)raw_memory;

    *nlp_memory = (ocp_nlp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_memory);

    int_t size_hess_l, size_grad_f, size_jac_h, size_jac_g, size_h, size_g;
    int_t size_x, size_u, size_pi, size_lam;
    size_of_memory_elements(nlp_in, &size_hess_l, &size_grad_f, &size_jac_h,
                            &size_jac_g, &size_h, &size_g, &size_x, &size_u,
                            &size_pi, &size_lam);

    (*nlp_memory)->hess_l = (const real_t **)c_ptr;
    c_ptr += size_hess_l;

    (*nlp_memory)->grad_f = (const real_t **)c_ptr;
    c_ptr += size_grad_f;

    (*nlp_memory)->jac_h = (const real_t **)c_ptr;
    c_ptr += size_jac_h;

    (*nlp_memory)->jac_g = (const real_t **)c_ptr;
    c_ptr += size_jac_g;

    (*nlp_memory)->h = (const real_t **)c_ptr;
    c_ptr += size_h;

    (*nlp_memory)->g = (const real_t **)c_ptr;
    c_ptr += size_g;

    (*nlp_memory)->x = (const real_t **)c_ptr;
    c_ptr += size_x;

    (*nlp_memory)->u = (const real_t **)c_ptr;
    c_ptr += size_u;
   
    (*nlp_memory)->pi = (const real_t **)c_ptr;
    c_ptr += size_pi;

    (*nlp_memory)->lam = (const real_t **)c_ptr;
    c_ptr += size_lam;

    return c_ptr;
}

ocp_nlp_memory *ocp_nlp_create_memory(const ocp_nlp_in *nlp_in) {
    
    ocp_nlp_memory *mem;

    int_t memory_size = ocp_nlp_calculate_memory_size(nlp_in);
    void *raw_memory_ptr = malloc(memory_size);

    char *ptr_end = ocp_nlp_assign_memory(nlp_in, (void **)&mem, raw_memory_ptr);
    assert((char *)raw_memory_ptr + memory_size >= ptr_end); (void)ptr_end;

    return mem;
}

void ocp_nlp_free_memory(int_t N, void *mem_) {

    // TODO(nielsvd): implement memory clean-up

}