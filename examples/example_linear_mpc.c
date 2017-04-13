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

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_i_aux.h"

#include "acados/utils/types.h"
#include "acados/utils/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"

// #include "test/test_utils/read_ocp_qp_in.h"

#define NN 10
#define NX 2
#define NU 1

int main( ) {

    int_t nMPC = 10;

    // MPC trajectories
    real_t *uMPC;
    real_t *xMPC;
    d_zeros(&uMPC, NU, nMPC);
    d_zeros(&xMPC, NX, nMPC+1);

    // LTI dynamics
    real_t A[NX*NX] = {0.66, -1.2, 0.6, 0.6};  // NOTE: column major format
    real_t B[NX*NU] = {1, 0.5};
    real_t b[NX] = {0, 0};

    // Quadratic objective
    real_t *Q;
    real_t *R;
    real_t *S;
    real_t *q;
    real_t *r;

    d_zeros(&Q, NX, NX);
    d_zeros(&R, NU, NU);
    d_zeros(&S, NU, NX);
    d_zeros(&q, NX, 1);
    d_zeros(&r, NU, 1);

    for (int ii = 0; ii < NX; ii++) Q[ii*NX+ii] = 1.0;
    for (int ii = 0; ii < NU; ii++) R[ii*NU+ii] = 1.0;

    // Bounds
    real_t umin[NU] = {-1};
    real_t umax[NU] = {1};

    real_t xmin[NX] = {-3, -5};
    real_t xmax[NX] = {5, 5};

    // Initial condition
    real_t x0[NX] = {-2, 5};

    // Concatenate bounds
    real_t zmin[NX+NU];
    real_t zmax[NX+NU];
    real_t z0min[NX+NU];
    real_t z0max[NX+NU];
    int_t idxb[NX+NU];

    for (int ii = 0; ii < NX; ii++) {
        zmin[ii] = xmin[ii];
        zmax[ii] = xmax[ii];
        idxb[ii] = ii;
    }
    for (int ii = 0; ii < NU; ii++) {
        zmin[NX+ii] = umin[ii];
        zmax[NX+ii] = umax[ii];
        z0min[NX+ii] = umin[ii];  // NOTE: First NX elements are updated in MPC loop
        z0max[NX+ii] = umax[ii];
        idxb[NX+ii] = NX+ii;
    }

    // Setup ocp_qp_in
    // TODO(dimitris): maybe function in utils for LTI systems?
    real_t *hQ[NN+1];
    real_t *hR[NN];
    real_t *hS[NN];
    real_t *hq[NN+1];
    real_t *hr[NN];
    real_t *hA[NN];
    real_t *hB[NN];
    real_t *hb[NN];
    real_t *hlb[NN+1];
    real_t *hub[NN+1];
    int_t *hidxb[NN+1];
    int_t nx[NN+1];
    int_t nu[NN+1];
    int_t nb[NN+1];
    int_t *nc;
    int_zeros(&nc, NN+1, 1);

    for (int kk = 0; kk < NN; kk++) {
        nx[kk] = NX;
        nu[kk] = NU;
        nb[kk] = NX+NU;
        hQ[kk] = Q;
        hR[kk] = R;
        hS[kk] = S;
        hq[kk] = q;
        hr[kk] = r;
        hA[kk] = A;
        hB[kk] = B;
        hb[kk] = b;
        if (kk > 0) {
            hlb[kk] = zmin;
            hub[kk] = zmax;
        } else {
            hlb[kk] = z0min;
            hub[kk] = z0max;
        }
        hidxb[kk] = idxb;
    }
    nx[NN] = NX;
    nu[NN] = 0;
    nb[NN] = NX;
    hQ[NN] = Q;
    hq[NN] = q;
    hlb[NN] = xmin;
    hub[NN] = xmax;
    hidxb[NN] = idxb;  // NOTE: the first nb[N]=NX will be read anyway

    ocp_qp_in qp_in;
    qp_in.N  = NN;
    qp_in.nx = nx;
    qp_in.nu = nu;
    qp_in.nb = nb;
    qp_in.nc = nc;
    qp_in.Q = (const real_t **)hQ;
    qp_in.R = (const real_t **)hR;
    qp_in.S = (const real_t **)hS;
    qp_in.q = (const real_t **)hq;
    qp_in.r = (const real_t **)hr;
    qp_in.A = (const real_t **)hA;
    qp_in.B = (const real_t **)hB;
    qp_in.b = (const real_t **)hb;
    qp_in.lb = (const real_t **)hlb;
    qp_in.ub = (const real_t **)hub;
    qp_in.idxb = (const int_t **)hidxb;

    // Setup ocp_qp_out, ocp_qp_args and ocp_qp_mem
    ocp_qp_out qp_out;
    allocate_ocp_qp_out(&qp_in, &qp_out);

    ocp_qp_qpdunes_args args;
    ocp_qp_qpdunes_memory mem;

    // TODO(dimitris): use LINEAR_MPC options instead and add timings
    ocp_qp_qpdunes_create_arguments(&args, QPDUNES_DEFAULT_ARGUMENTS);

    int_t work_space_size = ocp_qp_qpdunes_calculate_workspace_size(&qp_in, &args);
    void *work = (void*)malloc(work_space_size);

    ocp_qp_qpdunes_create_memory(&qp_in, &args, &mem);

    // print_ocp_qp_in(qp_in);

    // store initial condition
    for (int ii = 0; ii < NX; ii++) xMPC[ii] = x0[ii];

    for (int kk = 0; kk < nMPC; kk++) {

        // update constraint on x0
        for (int ii = 0; ii < NX; ii++) {
            z0min[ii] = x0[ii];
            z0max[ii] = x0[ii];
        }

        // solve QP
        ocp_qp_qpdunes(&qp_in, &qp_out, &args, &mem, work);

        // simulate system
        for (int ii = 0; ii < NX; ii++) {
            x0[ii] = b[ii];
            for (int jj = 0; jj < NX; jj++) x0[ii]+= A[ii+jj*NX]*qp_out.x[0][jj];
            for (int jj = 0; jj < NU; jj++) x0[ii]+= B[ii+jj*NX]*qp_out.u[0][jj];
        }

        // store MPC trajectories
        for (int ii = 0; ii < NX; ii++) xMPC[ii+(kk+1)*NX] = x0[ii];
        for (int ii = 0; ii < NU; ii++) uMPC[ii+kk*NU] = qp_out.u[0][ii];
    }

    // print trajectories
    printf("state trajectory:\n");
    d_print_mat(NX, nMPC, xMPC, NX);
    printf("control trajectory:\n");
    d_print_mat(NU, nMPC, uMPC, NU);

    // free dynamically allocated memory
    ocp_qp_qpdunes_free_memory(&mem);
    free(work);
    int_free(nc);
    d_free(r);
    d_free(q);
    d_free(S);
    d_free(R);
    d_free(Q);
    d_free(xMPC);
    d_free(uMPC);
}
