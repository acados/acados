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

#include <assert.h>

#include <stdlib.h>

void size_of_memory_elements(const ocp_nlp_in *nlp_in, const int_t stage,
                             int_t *size_hess_l, int_t *size_grad_f,
                             int_t *size_jac_h, int_t *size_jac_g,
                             int_t *size_h, int_t *size_g, int_t *size_x,
                             int_t *size_u, int_t *size_pi, int_t *size_lam) {
    int_t N = nlp_in->N;
    int_t *nx = (int_t *)nlp_in->nx;
    int_t *nu = (int_t *)nlp_in->nu;
    int_t *ng = (int_t *)nlp_in->ng;
    int_t *nb = (int_t *)nlp_in->nb;

    *size_hess_l =
        ((nx[stage] + nu[stage]) * (nx[stage] + nu[stage])) * sizeof(real_t);
    *size_grad_f = (nx[stage] + nu[stage]) * sizeof(real_t);
    *size_jac_h =
        (stage < N) ? (nx[stage + 1] * (nx[stage] + nu[stage])) * sizeof(real_t)
                    : 0;
    *size_jac_g = (ng[stage] * (nx[stage] + nu[stage])) * sizeof(real_t);
    *size_h = (stage < N) ? nx[stage + 1] * sizeof(real_t) : 0;
    *size_g = ng[stage] * sizeof(real_t);
    *size_x = nx[stage] * sizeof(real_t);
    *size_u = nu[stage] * sizeof(real_t);
    *size_pi = (stage < N) ? nx[stage+1] * sizeof(real_t) : 0;
    *size_lam = (2 * nb[stage] + 2 * ng[stage]) * sizeof(real_t);
}

int_t ocp_nlp_calculate_memory_size(const ocp_nlp_in *nlp_in) {
    int_t N = nlp_in->N;

    int_t size = sizeof(ocp_nlp_memory);

    size += sizeof(real_t *) * (N + 1);  // hess_l
    size += sizeof(real_t *) * (N + 1);  // grad_f
    size += sizeof(real_t *) * (N + 1);  // jac_h
    size += sizeof(real_t *) * (N + 1);  // jac_g
    size += sizeof(real_t *) * (N + 1);  // h
    size += sizeof(real_t *) * (N + 1);  // g
    size += sizeof(real_t *) * (N + 1);  // x
    size += sizeof(real_t *) * (N + 1);  // u
    size += sizeof(real_t *) * (N + 1);  // pi
    size += sizeof(real_t *) * (N + 1);  // lam

    for (int_t i = 0; i <= nlp_in->N; i++) {
        int_t size_hess_l, size_grad_f, size_jac_h, size_jac_g, size_h, size_g;
        int_t size_x, size_u, size_pi, size_lam;
        size_of_memory_elements(nlp_in, i, &size_hess_l, &size_grad_f,
                                &size_jac_h, &size_jac_g, &size_h, &size_g,
                                &size_x, &size_u, &size_pi, &size_lam);
        size += size_hess_l + size_grad_f + size_jac_h + size_jac_g + size_h +
                size_g;
        size += size_x + size_u + size_pi + size_lam;
    }

    return size;
}

char *ocp_nlp_assign_memory(const ocp_nlp_in *nlp_in, void **mem_,
                            void *raw_memory) {
    int_t N = nlp_in->N;

    ocp_nlp_memory **nlp_memory = (ocp_nlp_memory **)mem_;
    char *c_ptr = (char *)raw_memory;

    *nlp_memory = (ocp_nlp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_memory);

    (*nlp_memory)->hess_l = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->grad_f = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->jac_h = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->jac_g = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->h = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->g = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->x = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->u = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->pi = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    (*nlp_memory)->lam = (real_t **)c_ptr;
    c_ptr += sizeof(const real_t *) * (N + 1);

    for (int_t i = 0; i <= nlp_in->N; i++) {
        int_t size_hess_l, size_grad_f, size_jac_h, size_jac_g, size_h, size_g;
        int_t size_x, size_u, size_pi, size_lam;
        size_of_memory_elements(nlp_in, i, &size_hess_l, &size_grad_f,
                                &size_jac_h, &size_jac_g, &size_h, &size_g,
                                &size_x, &size_u, &size_pi, &size_lam);

        (*nlp_memory)->hess_l[i] = (real_t *)c_ptr;
        c_ptr += size_hess_l;

        (*nlp_memory)->grad_f[i] = (real_t *)c_ptr;
        c_ptr += size_grad_f;

        (*nlp_memory)->jac_h[i] = (real_t *)c_ptr;
        c_ptr += size_jac_h;

        (*nlp_memory)->jac_g[i] = (real_t *)c_ptr;
        c_ptr += size_jac_g;

        (*nlp_memory)->h[i] = (real_t *)c_ptr;
        c_ptr += size_h;

        (*nlp_memory)->g[i] = (real_t *)c_ptr;
        c_ptr += size_g;

        (*nlp_memory)->x[i] = (real_t *)c_ptr;
        c_ptr += size_x;

        (*nlp_memory)->u[i] = (real_t *)c_ptr;
        c_ptr += size_u;

        (*nlp_memory)->pi[i] = (real_t *)c_ptr;
        c_ptr += size_pi;

        (*nlp_memory)->lam[i] = (real_t *)c_ptr;
        c_ptr += size_lam;
    }

    return c_ptr;
}

ocp_nlp_memory *ocp_nlp_create_memory(const ocp_nlp_in *nlp_in) {
    ocp_nlp_memory *mem;

    int_t memory_size = ocp_nlp_calculate_memory_size(nlp_in);
    void *raw_memory_ptr = malloc(memory_size);

    char *ptr_end =
        ocp_nlp_assign_memory(nlp_in, (void **)&mem, raw_memory_ptr);
    assert((char *)raw_memory_ptr + memory_size >= ptr_end);
    (void)ptr_end;

    return mem;
}

void ocp_nlp_destroy(void *mem_) {
    ocp_nlp_memory *mem = (ocp_nlp_memory *)mem_;

    free(mem);
}
