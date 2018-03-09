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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"


#define NN 10
#define NX 2
#define NU 1

int main() {
    int_t nMPC = 10;

    acados_timer timer;
    real_t *cputimes;
    d_zeros(&cputimes, 1, nMPC);

    // MPC trajectories
    real_t *uMPC;
    real_t *xMPC;
    d_zeros(&uMPC, NU, nMPC);
    d_zeros(&xMPC, NX, nMPC + 1);

    // LTI dynamics
    real_t A[NX*NX] = {0.66, -1.2, 0.6, 0.6};  // NOTE(dimitris): column major format
    real_t B[NX*NU] = {1, 0.5};
    real_t b[NX] = {0, 0};

    // Quadratic objective
    real_t Q[NX*NX] = {0};
    real_t S[NU*NX] = {0};
    real_t R[NU*NU] = {0};
    real_t q[NX] = {0};
    real_t r[NU] = {0};

    for (int i = 0; i < NX; i++) Q[i * NX + i] = 1.0;
    for (int i = 0; i < NU; i++) R[i * NU + i] = 1.0;

    // NOTE(dimitris): uncomment to use qpOASES for stage QPs instead of clipping
    // Q[1] = 0.0000001;
    // Q[NX] = 0.000001;
    // S[0] = 0.0000001;
    // S[1] = 0.0000001;

    // Bounds
    real_t umin[NU] = {-1};
    real_t umax[NU] = {1};

    real_t xmin[NX] = {-3, -5};
    real_t xmax[NX] = {5, 5};

    // Initial condition
    real_t x0[NX] = {-2, 5};

    // Setup vectors that define problem size
    int_t nx[NN + 1];
    int_t nu[NN + 1];
    int_t nb[NN + 1];
    int_t nc[NN + 1];

    for (int_t k = 0; k < NN+1; k++) {
        nx[k] = NX;
        if (k < NN) {
            nu[k] = NU;
        } else {
            nu[k] = 0;
        }
        nb[k] = nx[k] + nu[k];
        nc[k] = 0;
    }

    ocp_qp_in *qp_in = ocp_qp_in_create(NN, nx, nu, nb, nc);

    // Copy LTI dynamics and constraints to QP memory
    for (int_t k = 0; k < NN+1; k++) {
        ocp_qp_in_copy_objective(Q, S, R, q, r, qp_in, k);
        if (k < NN)
            ocp_qp_in_copy_dynamics(A, B, b, qp_in, k);
    }

    int_t **hidxb = (int_t **) qp_in->idxb;
    real_t **hlb = (real_t **) qp_in->lb;
    real_t **hub = (real_t **) qp_in->ub;

    // Set up bounds
    for (int_t k = 0; k < NN+1; k++) {
        for (int_t i = 0; i < nb[k]; i++) {
            hidxb[k][i] = i;
#ifdef FLIP_BOUNDS
            if (i < nu[k]) {
                hlb[k][i] = umin[i];
                hub[k][i] = umax[i];
            } else {
                hlb[k][i] = xmin[i - nu[k]];
                hub[k][i] = xmax[i - nu[k]];
            }
#else
            if (i < nx[k]) {
                hlb[k][i] = xmin[i];
                hub[k][i] = xmax[i];
            } else {
                hlb[k][i] = umin[i - nx[k]];
                hub[k][i] = umax[i - nx[k]];
            }
#endif
        }
    }

    ocp_qp_out *qp_out = ocp_qp_out_create(NN, nx, nu, nb, nc);

    ocp_qp_qpdunes_opts *args = ocp_qp_qpdunes_create_arguments(QPDUNES_LINEAR_MPC);

    int_t work_space_size = ocp_qp_qpdunes_calculate_workspace_size(qp_in, args);
    void *work = (void *)malloc(work_space_size);

    ocp_qp_qpdunes_memory *mem = ocp_qp_qpdunes_create_memory(qp_in, args);

    // store initial condition
    for (int_t i = 0; i < NX; i++) xMPC[i] = x0[i];

    // linear MPC loop
    for (int_t kk = 0; kk < nMPC; kk++) {
        // update constraint on x0
        for (int_t i = 0; i < NX; i++) {
#ifdef FLIP_BOUNDS
            hlb[0][NU + i] = x0[i];
            hub[0][NU + i] = x0[i];
#else
            hlb[0][i] = x0[i];
            hub[0][i] = x0[i];
#endif
        }

        // solve QP
        acados_tic(&timer);
        ocp_qp_qpdunes(qp_in, qp_out, args, mem, work);
        cputimes[kk] = 1000*acados_toc(&timer);
        // print_ocp_qp_in(qp_in);
        // exit(1);
        // simulate system
        for (int ii = 0; ii < NX; ii++) {
            x0[ii] = b[ii];
            for (int jj = 0; jj < NX; jj++)
                x0[ii] += A[ii + jj * NX] * qp_out->x[0][jj];
            for (int jj = 0; jj < NU; jj++)
                x0[ii] += B[ii + jj * NX] * qp_out->u[0][jj];
        }

        // store MPC trajectories
        for (int ii = 0; ii < NX; ii++) xMPC[ii + (kk + 1) * NX] = x0[ii];
        for (int ii = 0; ii < NU; ii++) uMPC[ii + kk * NU] = qp_out->u[0][ii];
    }

    // print trajectories
    printf("state trajectory:\n");
    d_print_mat(NX, nMPC, xMPC, NX);
    printf("control trajectory:\n");
    d_print_mat(NU, nMPC, uMPC, NU);
    printf("cpu times [ms]:\n");
    d_print_mat(1, nMPC, cputimes, 1);

    // d_print_mat(NX+NU, NN, qp_out.x[0][0], NX+NU);

    // free dynamically allocated memory
    ocp_qp_qpdunes_free_memory(mem);
    free(work);
    free(qp_in);
    free(qp_out);
    d_free(xMPC);
    d_free(uMPC);
    d_free(cputimes);
}
