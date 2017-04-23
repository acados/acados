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

#include "hpmpc/include/aux_d.h"

#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "examples/c/Chen_model/chen_model.h"

#define NN 13
#define NX 2
#define NU 1

#ifdef DEBUG
static void print_states_controls(real_t *w, int_t N) {
    printf("node\tx\t\t\t\tu\n");
    for (int_t i = 0; i < N; i++) {
        printf("%4d\t%+e %+e\t%+e\n", i, w[i*(NX+NU)], w[i*(NX+NU)+1], w[i*(NX+NU)+2]);
    }
    printf("%4d\t%+e %+e\n", N, w[N*(NX+NU)], w[N*(NX+NU)+1]);
}
#endif  // DEBUG

static void shift_states(real_t *w, real_t *x_end, int_t N) {
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = w[(i+1)*(NX+NU)+j];
    }
    for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] = x_end[j];
}

static void shift_controls(real_t *w, real_t *u_end, int_t N) {
    for (int_t i = 0; i < N-1; i++) {
        for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = w[(i+1)*(NX+NU)+NX+j];
    }
    for (int_t j = 0; j < NU; j++) w[(N-1)*(NX+NU)+NX+j] = u_end[j];
}

// Simple SQP example for acados
int main() {
    // Problem data
    int_t   N                   = NN;
    real_t  x0[NX]              = {0.5, 0};
    real_t  w[NN*(NX+NU)+NX]    = {0};  // States and controls stacked
    real_t  Q[NX*NX]            = {0};
    real_t  R[NU*NU]            = {0};
    real_t  xref[NX]            = {0};
    real_t  uref[NX]            = {0};
    int_t   max_sqp_iters       = 1;
    int_t   max_iters           = 10000;
    real_t  x_end[NX]           = {0};
    real_t  u_end[NU]           = {0};

    for (int_t i = 0; i < NX; i++) w[i] = x0[i];
    for (int_t i = 0; i < NX; i++) Q[i*(NX+1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i*(NU+1)] = 0.05;

    // Integrator structs
    real_t T = 0.1;
    sim_in  sim_in;
    sim_out sim_out;
    sim_in.nSteps = 10;
    sim_in.step = T/sim_in.nSteps;
    sim_in.VDE_forw = &VDE_fun;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.sens_forw = true;
    sim_in.sens_adj = false;
    sim_in.sens_hess = false;
    sim_in.nsens_forw = NX+NU;

    sim_in.x = malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = malloc(sizeof(*sim_in.u) * (NU));
    sim_in.S_forw = malloc(sizeof(*sim_in.S_forw) * (NX*(NX+NU)));
    for (int_t i = 0; i < NX*(NX+NU); i++) sim_in.S_forw[i] = 0.0;
    for (int_t i = 0; i < NX; i++) sim_in.S_forw[i*(NX+1)] = 1.0;

    sim_out.xn = malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.S_forw = malloc(sizeof(*sim_out.S_forw) * (NX*(NX+NU)));

    sim_info erk_info;
    sim_out.info = &erk_info;

    sim_RK_opts rk_opts;
    sim_erk_create_opts(4, &rk_opts);
    sim_in.opts = &rk_opts;
    sim_erk_workspace erk_work;
    sim_erk_create_workspace(&sim_in, &erk_work);

    int_t nx[NN+1] = {0};
    int_t nu[NN] = {0};
    int_t nb[NN+1] = {0};
    int_t nc[NN+1] = {0};
    for (int_t i = 0; i < N; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[N] = NX;

    real_t *pA[N];
    real_t *pB[N];
    real_t *pb[N];
    real_t *pQ[N+1];
    real_t *pS[N];
    real_t *pR[N];
    real_t *pq[N+1];
    real_t *pr[N];
    real_t *px[N+1];
    real_t *pu[N];
    real_t *ppi[N];
    real_t *plam[N+1];
    real_t *px0[1];
    d_zeros(&px0[0], nx[0], 1);
    for (int_t i = 0; i < N; i++) {
        d_zeros(&pA[i], nx[i+1], nx[i]);
        d_zeros(&pB[i], nx[i+1], nu[i]);
        d_zeros(&pb[i], nx[i+1], 1);
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
        d_zeros(&ppi[i], nx[i], 1);
        d_zeros(&plam[i], nb[i]+nc[i], 1);
    }
    d_zeros(&pq[N], nx[N], 1);
    d_zeros(&px[N], nx[N], 1);
    d_zeros(&plam[N], nb[N]+nc[N], 1);

    // Allocate OCP QP variables
    ocp_qp_in qp_in;
    qp_in.N = N;
    ocp_qp_out qp_out;
    ocp_qp_condensing_qpoases_args args;
    real_t *work = NULL;
    qp_in.nx = nx;
    qp_in.nu = nu;
    qp_in.nb = nb;
    qp_in.nc = nc;
    for (int_t i = 0; i < N; i++) {
        pQ[i] = Q;
        pR[i] = R;
    }
    pQ[N] = Q;
    qp_in.Q = (const real_t **) pQ;
    qp_in.S = (const real_t **) pS;
    qp_in.R = (const real_t **) pR;
    qp_in.q = (const real_t **) pq;
    qp_in.r = (const real_t **) pr;
    qp_in.A = (const real_t **) pA;
    qp_in.B = (const real_t **) pB;
    qp_in.b = (const real_t **) pb;
    qp_in.lb = (const real_t **) px0;
    qp_out.x = px;
    qp_out.u = pu;
    qp_out.pi = ppi;
    qp_out.lam = plam;

    acado_timer timer;
    real_t total_time = 0;
    acado_tic(&timer);
    for (int_t iter = 0; iter < max_iters; iter++) {
        // printf("\n------ ITERATION %d ------\n", iter);
        for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
            for (int_t i = 0; i < N; i++) {
                // Pass state and control to integrator
                for (int_t j = 0; j < NX; j++) sim_in.x[j] = w[i*(NX+NU)+j];
                for (int_t j = 0; j < NU; j++) sim_in.u[j] = w[i*(NX+NU)+NX+j];
                sim_erk(&sim_in, &sim_out, 0, &erk_work);
                // Construct QP matrices
                for (int_t j = 0; j < NX; j++) {
                    pq[i][j] = Q[j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
                }
                for (int_t j = 0; j < NU; j++) {
                    pr[i][j] = R[j*(NX+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
                }
                for (int_t j = 0; j < NX; j++) {
                    pb[i][j] = sim_out.xn[j] - w[(i+1)*(NX+NU)+j];
                    for (int_t k = 0; k < NX; k++) pA[i][j*NX+k] = sim_out.S_forw[j*(NX)+k];
                }
                for (int_t j = 0; j < NU; j++)
                    for (int_t k = 0; k < NX; k++) pB[i][j*NX+k] = sim_out.S_forw[NX*NX + NX*j+k];
            }
            for (int_t j = 0; j < NX; j++) {
                px0[0][j] = (x0[j]-w[j]);
            }
            for (int_t j = 0; j < NX; j++) {
                pq[N][j] = Q[j*(NX+1)]*(w[N*(NX+NU)+j]-xref[j]);
            }
            int status = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, NULL, work);
            if (status) {
                printf("qpOASES returned error status %d\n", status);
                return -1;
            }
            for (int_t i = 0; i < N; i++) {
                for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] += qp_out.x[i][j];
                for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] += qp_out.u[i][j];
            }
            for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] += qp_out.x[N][j];
        }
        for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
        shift_states(w, x_end, N);
        shift_controls(w, u_end, N);
    }
    #ifdef DEBUG
    print_states_controls(&w[0], N);
    #endif  // DEBUG
    total_time = acado_toc(&timer);  // in seconds
    printf("Average of %.3f ms per iteration.\n", 1e3*total_time/max_iters);

    free(sim_in.x);
    free(sim_in.u);
    free(sim_in.S_forw);
    free(sim_out.xn);
    free(sim_out.S_forw);

    d_free(px0[0]);
    for (int_t i = 0; i < N; i++) {
        d_free(pA[i]);
        d_free(pB[i]);
        d_free(pb[i]);
        d_free(pS[i]);
        d_free(pq[i]);
        d_free(pr[i]);
        d_free(px[i]);
        d_free(pu[i]);
    }
    d_free(pq[N]);
    d_free(px[N]);

    free(rk_opts.A_mat);
    free(rk_opts.b_vec);
    free(rk_opts.c_vec);

    free(erk_work.rhs_forw_in);
    free(erk_work.K_traj);
    free(erk_work.out_forw_traj);

    return 0;
}
