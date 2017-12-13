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

#include "acados/ocp_qp/ocp_qp_common.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <assert.h>

#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#ifdef ACADOS_WITH_HPMPC
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#endif
#ifdef OOQP
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#endif
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/utils/types.h"

int_t ocp_qp_in_calculate_size(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                               const int_t *nc) {

    int_t bytes = sizeof(ocp_qp_in);

    bytes += 4*(N+1)*sizeof(int_t);  // nx, nu, nb, nc
    bytes += 3*N*sizeof(real_t *);  // A, B, b
    bytes += 11*(N+1)*sizeof(real_t *);  // ...
    bytes += 1*(N+1)*sizeof(int_t *);  // idxb

    for (int_t k = 0; k < N+1; k++) {

        if (k < N) {
            bytes += nx[k+1]*nx[k]*sizeof(real_t);  // A
            bytes += nx[k+1]*nu[k]*sizeof(real_t);  // B
            bytes += nx[k+1]*sizeof(real_t);  // b
        }

        bytes += nx[k]*nx[k]*sizeof(real_t);  // Q
        bytes += nu[k]*nx[k]*sizeof(real_t);  // S
        bytes += nu[k]*nu[k]*sizeof(real_t);  // R
        bytes += nx[k]*sizeof(real_t);  // q
        bytes += nu[k]*sizeof(real_t);  // r
        bytes += nb[k]*sizeof(int_t);  // idxb
        bytes += 2*nb[k]*sizeof(real_t);  // lb, ub
        bytes += nc[k]*nx[k]*sizeof(real_t);  // Cx
        bytes += nc[k]*nu[k]*sizeof(real_t);  // Cu
        bytes += 2*nc[k]*sizeof(real_t);  // lc, uc
    }

    bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    bytes += ALIGNMENT;

    return bytes;
}


char *assign_ocp_qp_in(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                       const int_t *nc, ocp_qp_in **qp_in, void *ptr) {

    // pointer to initialize QP data to zero
    char *c_ptr_QPdata;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_in = (ocp_qp_in *) c_ptr;
    c_ptr += sizeof(ocp_qp_in);

    // copy dimensions to workspace
    (*qp_in)->N = N;

    (*qp_in)->nx = (int_t *) c_ptr;
    memcpy(c_ptr, nx, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    (*qp_in)->nu = (int_t *) c_ptr;
    memcpy(c_ptr, nu, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    (*qp_in)->nb = (int_t *) c_ptr;
    memcpy(c_ptr, nb, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    (*qp_in)->nc = (int_t *) c_ptr;
    memcpy(c_ptr, nc, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    // assign double pointers
    (*qp_in)->A = (const real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_in)->B = (const real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_in)->b = (const real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_in)->Q = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->S = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->R = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->q = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->r = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->idxb = (const int_t **) c_ptr;
    c_ptr += (N+1)*sizeof(int_t *);

    (*qp_in)->lb = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->ub = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->Cx = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->Cu = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->lc = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->uc = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    // assign pointers to ints
    for (int_t k = 0; k < N+1; k++) {
        (*qp_in)->idxb[k] = (int_t *) c_ptr;
        c_ptr += nb[k]*sizeof(int_t);
    }

    // align data
    size_t l_ptr = (size_t) c_ptr;
    l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    c_ptr = (char *) l_ptr;

    // assign pointers to doubles
    c_ptr_QPdata = c_ptr;

    for (int_t k = 0; k < N+1; k++) {
        // printf("%zu MODULO %d = %zu\n", (size_t)c_ptr, 8, (size_t)c_ptr % 8);
        assert((size_t)c_ptr % 8 == 0);

        if (k < N) {
            (*qp_in)->A[k] = (real_t *) c_ptr;
            c_ptr += nx[k+1]*nx[k]*sizeof(real_t);

            (*qp_in)->B[k] = (real_t *) c_ptr;
            c_ptr += nx[k+1]*nu[k]*sizeof(real_t);

            (*qp_in)->b[k] = (real_t *) c_ptr;
            c_ptr += nx[k+1]*sizeof(real_t);
        }

        (*qp_in)->Q[k] = (real_t *) c_ptr;
        c_ptr += nx[k]*nx[k]*sizeof(real_t);

        (*qp_in)->S[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*nx[k]*sizeof(real_t);

        (*qp_in)->R[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*nu[k]*sizeof(real_t);

        (*qp_in)->q[k] = (real_t *) c_ptr;
        c_ptr += nx[k]*sizeof(real_t);

        (*qp_in)->r[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*sizeof(real_t);

        (*qp_in)->lb[k] = (real_t *) c_ptr;
        c_ptr += nb[k]*sizeof(real_t);

        (*qp_in)->ub[k] = (real_t *) c_ptr;
        c_ptr += nb[k]*sizeof(real_t);

        (*qp_in)->Cx[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*nx[k]*sizeof(real_t);

        (*qp_in)->Cu[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*nu[k]*sizeof(real_t);

        (*qp_in)->lc[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*sizeof(real_t);

        (*qp_in)->uc[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*sizeof(real_t);
    }

    // set QP data to zero (mainly for valgrind)
    for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
        *idx = 0;

    return c_ptr;
}


ocp_qp_in *create_ocp_qp_in(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                            const int_t *nc) {

    ocp_qp_in *qp_in;

    int_t bytes = ocp_qp_in_calculate_size(N, nx, nu, nb, nc);

    // TODO(dimitris): replace with acados_malloc to replace malloc at one place if not supported
    void *ptr = malloc(bytes);

    // // set a value for debugging
    // char *c_ptr = (char *) ptr;
    // for (int_t i = 0; i < bytes; i++) c_ptr[i] = 13;

    char *ptr_end = assign_ocp_qp_in(N, nx, nu, nb, nc, &qp_in, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    // for (int_t i = 0; i < bytes; i++) printf("%d - ", c_ptr[i]);
    // exit(1);

    return qp_in;
}


int_t ocp_qp_out_calculate_size(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                                const int_t *nc) {

    int_t bytes = sizeof(ocp_qp_out);

    bytes += 3*(N+1)*sizeof(real_t *);  // u, x, lam
    bytes += N*sizeof(real_t *);  // pi

    for (int_t k = 0; k < N+1; k++) {
        bytes += (nx[k] + nu[k])*sizeof(real_t);  // u, x
        if (k < N)
            bytes += (nx[k+1])*sizeof(real_t);  // pi
        bytes += 2*(nb[k] + nc[k])*sizeof(real_t);  // lam
        }

    bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    bytes += ALIGNMENT;

    return bytes;
}


char *assign_ocp_qp_out(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                        const int_t *nc, ocp_qp_out **qp_out,
    void *ptr) {

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_out = (ocp_qp_out *) c_ptr;
    c_ptr += sizeof(ocp_qp_out);

    // assign double pointers
    (*qp_out)->x = (real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_out)->u = (real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_out)->pi = (real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_out)->lam = (real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    // align data
    size_t l_ptr = (size_t) c_ptr;
    l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    c_ptr = (char *) l_ptr;

    // NOTE(dimitris): splitted the loops below to be able to print primal/dual solution at once

    // assign pointers to QP solution
    for (int_t k = 0; k < N+1; k++) {
        assert((size_t)c_ptr % 8 == 0);

        (*qp_out)->x[k] = (real_t *) c_ptr;
        c_ptr += nx[k]*sizeof(real_t);

        (*qp_out)->u[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*sizeof(real_t);
    }

    for (int_t k = 0; k < N; k++) {
        assert((size_t)c_ptr % 8 == 0);
        (*qp_out)->pi[k] = (real_t *) c_ptr;
        c_ptr += nx[k+1]*sizeof(real_t);
    }

    for (int_t k = 0; k < N+1; k++) {
        assert((size_t)c_ptr % 8 == 0);
        (*qp_out)->lam[k] = (real_t *) c_ptr;
        c_ptr += 2*(nb[k] + nc[k])*sizeof(real_t);
    }
    return c_ptr;
}

ocp_qp_out *create_ocp_qp_out(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                              const int_t *nc) {

    ocp_qp_out *qp_out;

    int_t bytes = ocp_qp_out_calculate_size(N, nx, nu, nb, nc);
    void *ptr = malloc(bytes);
    char *ptr_end = assign_ocp_qp_out(N, nx, nu, nb, nc, &qp_out, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    return qp_out;
}


void ocp_qp_in_copy_dynamics(const real_t *A, const real_t *B, const real_t *b,
                             ocp_qp_in *qp_in, int_t stage) {

    real_t **hA = (real_t **) qp_in->A;
    real_t **hB = (real_t **) qp_in->B;
    real_t **hb = (real_t **) qp_in->b;

    memcpy(hA[stage], A, qp_in->nx[stage+1]*qp_in->nx[stage]*sizeof(real_t));
    memcpy(hB[stage], B, qp_in->nx[stage+1]*qp_in->nu[stage]*sizeof(real_t));
    memcpy(hb[stage], b, qp_in->nx[stage+1]*sizeof(real_t));
}

void ocp_qp_in_copy_objective(const real_t *Q, const real_t *S, const real_t *R, const real_t *q,
                              const real_t *r, ocp_qp_in *qp_in, int_t stage) {

        real_t **hQ = (real_t **) qp_in->Q;
        real_t **hS = (real_t **) qp_in->S;
        real_t **hR = (real_t **) qp_in->R;
        real_t **hq = (real_t **) qp_in->q;
        real_t **hr = (real_t **) qp_in->r;

        memcpy(hQ[stage], Q, qp_in->nx[stage]*qp_in->nx[stage]*sizeof(real_t));
        memcpy(hS[stage], S, qp_in->nu[stage]*qp_in->nx[stage]*sizeof(real_t));
        memcpy(hR[stage], R, qp_in->nu[stage]*qp_in->nu[stage]*sizeof(real_t));
        memcpy(hq[stage], q, qp_in->nx[stage]*sizeof(real_t));
        memcpy(hr[stage], r, qp_in->nu[stage]*sizeof(real_t));
}

ocp_qp_solver *create_ocp_qp_solver(const ocp_qp_in *qp_in, const char *solver_name,
                                    void *solver_options) {
    ocp_qp_solver *qp_solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver));

    qp_solver->qp_in = (ocp_qp_in *) qp_in;
    qp_solver->qp_out = create_ocp_qp_out(qp_in->N, qp_in->nx, qp_in->nu, qp_in->nb, qp_in->nc);
    qp_solver->args = solver_options;

    if (!strcmp(solver_name, "qpdunes")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_qpdunes_create_arguments(QPDUNES_NONLINEAR_MPC);
        qp_solver->fun = &ocp_qp_qpdunes;
        qp_solver->initialize = &ocp_qp_qpdunes_initialize;
        qp_solver->destroy = &ocp_qp_qpdunes_destroy;
#ifdef OOQP
    } else if (!strcmp(solver_name, "ooqp")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_ooqp_create_arguments();
        qp_solver->fun = &ocp_qp_ooqp;
        qp_solver->initialize = &ocp_qp_ooqp_initialize;
        qp_solver->destroy = &ocp_qp_ooqp_destroy;
#endif
    } else if (!strcmp(solver_name, "condensing_qpoases")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_condensing_qpoases_create_arguments(qp_in);
        qp_solver->fun = &ocp_qp_condensing_qpoases;
        qp_solver->initialize = &ocp_qp_condensing_qpoases_initialize;
        qp_solver->destroy = &ocp_qp_condensing_qpoases_destroy;
#ifdef ACADOS_WITH_HPMPC
    } else if (!strcmp(solver_name, "hpmpc")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_hpmpc_create_arguments(qp_in, HPMPC_DEFAULT_ARGUMENTS);
        qp_solver->fun = &ocp_qp_hpmpc;
        qp_solver->initialize = &ocp_qp_hpmpc_initialize;
        qp_solver->destroy = &ocp_qp_hpmpc_destroy;
#endif
    } else if (!strcmp(solver_name, "condensing_hpipm")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_condensing_hpipm_create_arguments(qp_in);
        qp_solver->fun = &ocp_qp_condensing_hpipm;
        qp_solver->initialize = &ocp_qp_condensing_hpipm_initialize;
        qp_solver->destroy = &ocp_qp_condensing_hpipm_destroy;
    } else if (!strcmp(solver_name, "hpipm")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_hpipm_create_arguments(qp_in);
        qp_solver->fun = &ocp_qp_hpipm;
        qp_solver->initialize = &ocp_qp_hpipm_initialize;
        qp_solver->destroy = &ocp_qp_hpipm_destroy;
    } else {
        printf("Chosen QP solver not available\n");
        exit(1);
    }
    qp_solver->initialize(qp_solver->qp_in, qp_solver->args, &qp_solver->mem, &qp_solver->work);

    return qp_solver;
}
