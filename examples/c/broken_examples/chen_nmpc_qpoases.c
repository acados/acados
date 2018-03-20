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
#include "examples/c/chen_model/chen_model.h"

#define NN 13
#define NX 2
#define NU 1

#ifdef DEBUG
static void print_states_controls(real_t *w, int_t N) {
    printf("node\tx\t\t\t\tu\n");
    for (int_t i = 0; i < N; i++) {
        printf("%4d\t%+e %+e\t%+e\n", i, w[i * (NX + NU)], w[i * (NX + NU) + 1],
               w[i * (NX + NU) + 2]);
    }
    printf("%4d\t%+e %+e\n", N, w[N * (NX + NU)], w[N * (NX + NU) + 1]);
}
#endif  // DEBUG

static void shift_states(real_t *w, real_t *x_end, int_t N) {
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < NX; j++)
            w[i * (NX + NU) + j] = w[(i + 1) * (NX + NU) + j];
    }
    for (int_t j = 0; j < NX; j++) w[N * (NX + NU) + j] = x_end[j];
}

static void shift_controls(real_t *w, real_t *u_end, int_t N) {
    for (int_t i = 0; i < N - 1; i++) {
        for (int_t j = 0; j < NU; j++)
            w[i * (NX + NU) + NX + j] = w[(i + 1) * (NX + NU) + NX + j];
    }
    for (int_t j = 0; j < NU; j++) w[(N - 1) * (NX + NU) + NX + j] = u_end[j];
}

// Simple SQP example for acados
int main() {
    // problem size
    int_t nx[NN + 1];
    int_t nu[NN + 1];
    int_t nb[NN + 1];
    int_t nc[NN + 1];
    nx[0] = NX;
    nu[0] = NU;
    nb[0] = NX;
    nc[0] = 0;
    for (int_t i = 1; i < NN; i++) {
        nx[i] = NX;
        nu[i] = NU;
        nb[i] = 0;
        nc[i] = 0;
    }
    nu[NN] = 0;
    nx[NN] = NX;
    nb[NN] = 0;
    nc[NN] = 0;

    // Problem data
    real_t x0[NX] = {0.5, 0};
    real_t w[NN * (NX + NU) + NX] = {0};  // States and controls stacked
    real_t Q[NX * NX] = {0};
    real_t R[NU * NU] = {0};
    real_t xref[NX] = {0};
    real_t uref[NX] = {0};
    int_t max_sqp_iters = 1;
    int_t max_iters = 1000;
    real_t x_end[NX] = {0};
    real_t u_end[NU] = {0};
    real_t px0[NX];
    int_t idxb0[NX];

    for (int_t i = 0; i < NX; i++) w[i] = x0[i];
    for (int_t i = 0; i < NX; i++) Q[i * (NX + 1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i * (NU + 1)] = 0.05;
    for (int_t i = 0; i < NX; i++) idxb0[i] = nu[0] + i;

    // Integrator structs
    real_t T = 0.1;
    sim_in sim_in;
    sim_out sim_out;
    sim_in.num_steps = 10;
    sim_in.step = T / sim_in.num_steps;
    sim_in.forward_vde_wrapper = &VDE_fun;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.sens_forw = true;
    sim_in.sens_adj = false;
    sim_in.sens_hess = false;
    sim_in.num_forw_sens = NX + NU;

    sim_in.x = malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = malloc(sizeof(*sim_in.u) * (NU));
    sim_in.S_forw = malloc(sizeof(*sim_in.S_forw) * (NX * (NX + NU)));
    for (int_t i = 0; i < NX * (NX + NU); i++) sim_in.S_forw[i] = 0.0;
    for (int_t i = 0; i < NX; i++) sim_in.S_forw[i * (NX + 1)] = 1.0;

    sim_out.xn = malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.S_forw = malloc(sizeof(*sim_out.S_forw) * (NX * (NX + NU)));

    sim_info erk_info;
    sim_out.info = &erk_info;

    sim_rk_opts rk_opts;
    sim_erk_create_arguments(&rk_opts, 4);
    void *erk_work;
    int_t work_space_size = sim_erk_calculate_workspace_size(&sim_in, &rk_opts);
    erk_work = (void *) malloc(work_space_size);

    real_t *pA[NN];
    real_t *pB[NN];
    real_t *pCx[NN + 1];
    real_t *pCu[NN];
    real_t *pb[NN];
    real_t *pQ[NN + 1];
    real_t *pS[NN];
    real_t *pR[NN];
    real_t *pq[NN + 1];
    real_t *pr[NN];
    real_t *px[NN + 1];
    real_t *pu[NN + 1];
    real_t *ppi[NN];
    real_t *plam[NN + 1];
    real_t *plb[NN + 1];
    real_t *pub[NN + 1];
    real_t *plg[NN + 1];
    real_t *pug[NN + 1];
    int_t *pidxb[NN + 1];

    d_zeros(&pA[0], nx[1], nx[0]);  // TODO(dimitris): do we need this? It's the same in the loop
    d_zeros(&pB[0], nx[1], nu[0]);
    d_zeros(&pb[0], nx[1], 1);
    pQ[0] = Q;
    pR[0] = R;
    d_zeros(&pS[0], nu[0], nx[0]);
    d_zeros(&pq[0], nx[0], 1);
    d_zeros(&pr[0], nu[0], 1);
    d_zeros(&px[0], nx[0], 1);
    d_zeros(&pu[0], nu[0], 1);
    d_zeros(&ppi[0], nx[1], 1);
    d_zeros(&plam[0], 2 * nb[0] + 2 * nc[0], 1);
    plb[0] = px0;
    pub[0] = px0;
    pidxb[0] = idxb0;
    for (int_t i = 0; i < NN; i++) {
        d_zeros(&pA[i], nx[i + 1], nx[i]);
        d_zeros(&pB[i], nx[i + 1], nu[i]);
        d_zeros(&pCx[i], nc[i], nx[i]);
        d_zeros(&pCu[i], nc[i], nu[i]);
        d_zeros(&pb[i], nx[i + 1], 1);
        pQ[i] = Q;
        pR[i] = R;
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
        d_zeros(&ppi[i], nx[i + 1], 1);
        d_zeros(&plam[i], 2*nb[i] + 2*nc[i], 1);
    }
    pQ[NN] = Q;
    d_zeros(&pq[NN], nx[NN], 1);
    d_zeros(&px[NN], nx[NN], 1);
    d_zeros(&pu[NN], nu[NN], 1);
    d_zeros(&plam[NN], 2*nb[NN] + 2*nc[NN], 1);
    d_zeros(&pCx[NN], nc[NN], nx[NN]);

    // Allocate OCP QP variables
    ocp_qp_in qp_in;
    qp_in.N = NN;
    qp_in.nx = nx;
    qp_in.nu = nu;
    qp_in.nb = nb;
    qp_in.nc = nc;
    qp_in.idxb = (const int_t **)pidxb;
    qp_in.Q = (const real_t **)pQ;
    qp_in.S = (const real_t **)pS;
    qp_in.R = (const real_t **)pR;
    qp_in.q = (const real_t **)pq;
    qp_in.r = (const real_t **)pr;
    qp_in.A = (const real_t **)pA;
    qp_in.B = (const real_t **)pB;
    qp_in.Cx = (const real_t **)pCx;
    qp_in.Cu = (const real_t **)pCu;
    qp_in.b = (const real_t **)pb;
    qp_in.lb = (const real_t **)plb;
    qp_in.ub = (const real_t **)pub;
    qp_in.lc = (const real_t **)plg;
    qp_in.uc = (const real_t **)pug;

    ocp_qp_out qp_out;
    qp_out.x = px;
    qp_out.u = pu;
    qp_out.pi = ppi;
    qp_out.lam = plam;

    ocp_qp_condensing_qpoases_args *qpoases_args =
        ocp_qp_condensing_qpoases_create_arguments(&qp_in);
    qpoases_args->cputime = 1000.0;  // maximum cpu time in seconds
    qpoases_args->nwsr = 1000;  // maximum number of working set recalculations
    qpoases_args->warm_start = 0;  // wam start with dual_sol in memory

    int workspace_size = ocp_qp_condensing_qpoases_calculate_workspace_size(&qp_in, qpoases_args);
    void *workspace = malloc(workspace_size);

    int memory_size = ocp_qp_condensing_qpoases_calculate_memory_size(&qp_in, qpoases_args);
    void *memory = malloc(memory_size);
    ocp_qp_condensing_qpoases_memory *qpoases_memory;
    ocp_qp_condensing_qpoases_assign_memory(&qp_in, qpoases_args,
                                            (void **) &qpoases_memory, memory);

    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);
    for (int_t iter = 0; iter < max_iters; iter++) {
        // printf("\n------ ITERATION %d ------\n", iter);
        for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
            for (int_t i = 0; i < NN; i++) {
                // Pass state and control to integrator
                for (int_t j = 0; j < NX; j++)
                    sim_in.x[j] = w[i * (NX + NU) + j];
                for (int_t j = 0; j < NU; j++)
                    sim_in.u[j] = w[i * (NX + NU) + NX + j];
                sim_erk(&sim_in, &sim_out, &rk_opts, 0, erk_work);
                // Construct QP matrices
                for (int_t j = 0; j < NX; j++) {
                    pq[i][j] =
                        Q[j * (NX + 1)] * (w[i * (NX + NU) + j] - xref[j]);
                }
                for (int_t j = 0; j < NU; j++) {
                    pr[i][j] =
                        R[j * (NX + 1)] * (w[i * (NX + NU) + NX + j] - uref[j]);
                }
                for (int_t j = 0; j < NX; j++) {
                    pb[i][j] = sim_out.xn[j] - w[(i + 1) * (NX + NU) + j];
                    for (int_t k = 0; k < NX; k++)
                        pA[i][j * NX + k] = sim_out.S_forw[j * (NX) + k];
                }
                for (int_t j = 0; j < NU; j++)
                    for (int_t k = 0; k < NX; k++)
                        pB[i][j * NX + k] =
                            sim_out.S_forw[NX * NX + NX * j + k];
            }
            for (int_t j = 0; j < NX; j++) {
                px0[j] = (x0[j] - w[j]);
            }
            for (int_t j = 0; j < NX; j++) {
                pq[NN][j] = Q[j * (NX + 1)] * (w[NN * (NX + NU) + j] - xref[j]);
            }
            int status = ocp_qp_condensing_qpoases(
                &qp_in, &qp_out, qpoases_args, qpoases_memory, workspace);
            if (status) {
                printf("qpOASES returned error status %d\n", status);
                return -1;
            }
            for (int_t i = 0; i < NN; i++) {
                for (int_t j = 0; j < NX; j++)
                    w[i * (NX + NU) + j] += qp_out.x[i][j];
                for (int_t j = 0; j < NU; j++)
                    w[i * (NX + NU) + NX + j] += qp_out.u[i][j];
            }
            for (int_t j = 0; j < NX; j++)
                w[NN * (NX + NU) + j] += qp_out.x[NN][j];
        }
        for (int_t i = 0; i < NX; i++) x0[i] = w[NX + NU + i];
        shift_states(w, x_end, NN);
        shift_controls(w, u_end, NN);
    }
#ifdef DEBUG
    print_states_controls(&w[0], NN);
#endif  // DEBUG
    total_time = acados_toc(&timer);  // in seconds
    printf("Average of %.3f ms per iteration.\n", 1e3*total_time/max_iters);

    free(sim_in.x);
    free(sim_in.u);
    free(sim_in.S_forw);
    free(sim_out.xn);
    free(sim_out.S_forw);

    for (int_t i = 0; i < NN; i++) {
        d_free(pA[i]);
        d_free(pB[i]);
        d_free(pb[i]);
        d_free(pS[i]);
        d_free(pq[i]);
        d_free(pr[i]);
        d_free(px[i]);
        d_free(pu[i]);
        d_free(ppi[i]);
        d_free(plam[i]);
    }
    d_free(pq[NN]);
    d_free(px[NN]);
    d_free(plam[NN]);

    free(rk_opts.A_mat);
    free(rk_opts.b_vec);
    free(rk_opts.c_vec);

    free(erk_work);

    return 0;
}
