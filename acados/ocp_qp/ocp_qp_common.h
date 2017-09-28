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

#ifndef ACADOS_OCP_QP_OCP_QP_COMMON_H_
#define ACADOS_OCP_QP_OCP_QP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/utils/types.h"

typedef struct {
    int_t N;
    const int_t *nx;
    const int_t *nu;
    const int_t *nb;
    const int_t *nc;
    const real_t **A;
    const real_t **B;
    const real_t **b;
    const real_t **Q;
    const real_t **S;
    const real_t **R;
    const real_t **q;
    const real_t **r;
    const int_t **idxb;
    const real_t **lb;
    const real_t **ub;
    const real_t **Cx;
    const real_t **Cu;
    const real_t **lc;
    const real_t **uc;
} ocp_qp_in;

typedef struct {
    real_t **x;
    real_t **u;
    real_t **pi;
    real_t **lam;
    real_t **t;  // TODO(roversch): remove!
} ocp_qp_out;

typedef struct {
    int_t (*fun)(const ocp_qp_in *qp_in, ocp_qp_out *qp_out, void *args, void *mem, void *work);
    void (*initialize)(const ocp_qp_in *qp_in, void *args, void **mem, void **work);
    void (*destroy)(void *mem, void *work);
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;
    void *args;
    void *mem;
    void *work;
} ocp_qp_solver;

int_t ocp_qp_in_calculate_size(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                               const int_t *nc);

char *assign_ocp_qp_in(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                       const int_t *nc, ocp_qp_in **qp_in, void *ptr);

ocp_qp_in *create_ocp_qp_in(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                            const int_t *nc);

int_t ocp_qp_out_calculate_size(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                                const int_t *nc);

char *assign_ocp_qp_out(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                        const int_t *nc, ocp_qp_out **qp_out, void *ptr);

ocp_qp_out *create_ocp_qp_out(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                              const int_t *nc);

void ocp_qp_in_copy_dynamics(const real_t *A, const real_t *B, const real_t *b, ocp_qp_in *qp_in,
                             int_t stage);

void ocp_qp_in_copy_objective(const real_t *Q, const real_t *S, const real_t *R, const real_t *q,
                              const real_t *r, ocp_qp_in *qp_in, int_t stage);

ocp_qp_solver *create_ocp_qp_solver(const ocp_qp_in *qp_in, const char *name, void *options);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_QP_OCP_QP_COMMON_H_
