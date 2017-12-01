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
#include "acados/utils/math.h"
#include "acados/utils/types.h"
#include "examples/chen/chen_model.h"

#define NN 13
#define NX 2
#define NU 1
#define NBX 0
#define NBU 0

#ifdef DEBUG
static void print_states_controls(real_t **w, int_t N) {
    printf("node\tx\t\t\t\tu\n");
    for (int_t i = 0; i < N; i++) {
        printf("%4d\t%+e %+e\t%+e\n", i, w[i][0], w[i][1], w[i][2]);
    }
    printf("%4d\t%+e %+e\n", N, w[N][0], w[N][1]);
}
#endif  // DEBUG

// Simple SQP example for acados
int main() {
    // Problem data
    real_t x0[NX] = {0.5, 0};
    real_t Q[NX * NX] = {0};
    real_t R[NU * NU] = {0};
    real_t xref[NX] = {0};
    real_t uref[NX] = {0};
    int_t max_sqp_iters = 30;
    int_t timing_iters = 1;
    real_t x_end[NX] = {0};
    real_t u_end[NU] = {0};

    real_t *w[NN + 1];      // nlp states and controls stacked
    real_t *pi_n[NN + 1];   // nlp eq. mult
    real_t *lam_n[NN + 1];  // nlp ineq. mult

    for (int_t i = 0; i < NN; i++) {
        d_zeros(&w[i], NX + NU, 1);
        d_zeros(&pi_n[i], NX, 1);
        d_zeros(&lam_n[i], 2 * (NBX + NBU), 1);
    }

    d_zeros(&w[NN], NX, 1);
    d_zeros(&pi_n[NN], NX, 1);
    d_zeros(&lam_n[NN], NBX, 1);

    for (int_t i = 0; i < NX; i++) w[0][i] = x0[i];
    for (int_t i = 0; i < NX; i++) Q[i * (NX + 1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i * (NU + 1)] = 0.05;

    // Integrator structs
    real_t T = 0.1;
    sim_in sim_in;
    sim_out sim_out;
    sim_in.num_steps = 30;
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

    sim_erk_workspace erk_work;
    sim_rk_opts rk_opts;
    sim_erk_create_arguments(4, &rk_opts);
    sim_erk_create_workspace(&sim_in, &rk_opts, &erk_work);

    int_t nx[NN + 1] = {0};
    int_t nu[NN + 1] = {0};
    int_t nb[NN + 1] = {0};
    int_t nc[NN + 1] = {0};
    for (int_t i = 0; i < NN; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[NN] = NX;

    real_t *pA[NN];
    real_t *pB[NN];
    real_t *pb[NN];
    real_t *pQ[NN + 1];
    real_t *pS[NN];
    real_t *pR[NN];
    real_t *pq[NN + 1];
    real_t *pr[NN];
    real_t *px[NN + 1];
    real_t *pu[NN];
    real_t *px0[1];
    d_zeros(&px0[0], nx[0], 1);
    for (int_t i = 0; i < NN; i++) {
        d_zeros(&pA[i], nx[i + 1], nx[i]);
        d_zeros(&pB[i], nx[i + 1], nu[i]);
        d_zeros(&pb[i], nx[i + 1], 1);
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
    }
    d_zeros(&pq[NN], nx[NN], 1);
    d_zeros(&px[NN], nx[NN], 1);

    real_t *ppi[NN];
    real_t *plam[NN + 1];

    int ii;
    ii = 0;
    d_zeros(&ppi[ii], nx[ii + 1], 1);
    d_zeros(&plam[ii], 2 * nb[ii], 1);
    for (ii = 1; ii < NN; ii++) {
        d_zeros(&ppi[ii], nx[ii + 1], 1);
        d_zeros(&plam[ii], 2 * nb[ii], 1);
    }
    d_zeros(&plam[NN], 2 * nb[NN], 1);

    // Allocate OCP QP variables
    ocp_qp_in qp_in;
    qp_in.N = NN;
    ocp_qp_out qp_out;
    ocp_qp_condensing_qpoases_args args;
    real_t *work = NULL;
    qp_in.nx = nx;
    qp_in.nu = nu;
    qp_in.nb = nb;
    qp_in.nc = nc;
    for (int_t i = 0; i < NN; i++) {
        pQ[i] = Q;
        pR[i] = R;
    }
    pQ[NN] = Q;
    qp_in.Q = (const real_t **)pQ;
    qp_in.S = (const real_t **)pS;
    qp_in.R = (const real_t **)pR;
    qp_in.q = (const real_t **)pq;
    qp_in.r = (const real_t **)pr;
    qp_in.A = (const real_t **)pA;
    qp_in.B = (const real_t **)pB;
    qp_in.b = (const real_t **)pb;
    qp_in.lb = (const real_t **)px0;
    qp_in.ub = (const real_t **)
        px0;  // Andrea: we need to set the lower bound too here, right?
    qp_out.x = px;
    qp_out.u = pu;
    qp_out.pi = ppi;
    qp_out.lam = plam;

    initialise_qpoases(&qp_in);

    acados_timer timer;
    real_t total_time = 0;

    // Define residuals
    real_t *res_stat[NN + 1];
    real_t *res_eq[NN];
    real_t *res_ineq[NN + 1];
    real_t *res_compl[NN + 1];

    // Allocate memory for the residuals
    for (ii = 0; ii < NN; ii++) {
        d_zeros(&res_stat[ii], nx[ii] + nu[ii], 1);
        d_zeros(&res_eq[ii], nx[ii + 1], 1);
        d_zeros(&res_ineq[ii], 2 * nb[ii], 1);
        d_zeros(&res_compl[ii], 2 * nb[ii], 1);
    }
    ii = NN;
    d_zeros(&res_stat[ii], nx[ii] + nu[ii], 1);
    d_zeros(&res_ineq[ii], 2 * nb[ii], 1);
    d_zeros(&res_compl[ii], 2 * nb[ii], 1);

    acados_tic(&timer);
    for ( int_t iter = 0; iter < timing_iters; iter++ ) {
        for ( int_t i = 0; i < NX; i++ ) w[0][i] = x0[i];

        for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
            printf("\n------ ITERATION %d ------\n", sqp_iter);
            for (int_t i = 0; i < NN; i++) {
                // Pass state and control to integrator
                for (int_t j = 0; j < NX; j++) sim_in.x[j] = w[i][j];
                for (int_t j = 0; j < NU; j++) sim_in.u[j] = w[i][j + NX];
                sim_erk(&sim_in, &sim_out, &rk_opts, &erk_work);

                // Construct QP matrices
                for (int_t j = 0; j < NX; j++) {
                    pq[i][j] = Q[j * (NX + 1)] * (w[i][j] - xref[j]);
                }
                for (int_t j = 0; j < NU; j++) {
                    pr[i][j] = R[j * (NX + 1)] * (w[i][j + NX] - uref[j]);
                }
                for (int_t j = 0; j < NX; j++) {
                    pb[i][j] = sim_out.xn[j] - w[i + 1][j];
                    for (int_t k = 0; k < NX; k++)
                        pA[i][j * NX + k] = sim_out.S_forw[j * (NX) + k];
                }
                for (int_t j = 0; j < NU; j++)
                    for (int_t k = 0; k < NX; k++)
                        pB[i][j * NX + k] =
                            sim_out.S_forw[NX * NX + NX * j + k];
            }

            for (int_t j = 0; j < NX; j++) {
                px0[0][j] = (x0[j] - w[0][j]);
            }
            for (int_t j = 0; j < NX; j++) {
                pq[NN][j] = Q[j * (NX + 1)] * (w[NN][j] - xref[j]);
            }

            // Compute residuals
            int_t i;
            for (i = 0; i < NN - 1; i++) {
                for (int_t j = 0; j < nx[i] + nu[i]; j++)
                    res_stat[i][j] = -pi_n[i][j];
                dgemv_n_3l(NX, NX, pQ[i], NX, w[i], res_stat[i]);
                dgemv_t_3l(NX, NX, pA[i], NX, pi_n[i + 1], res_stat[i]);
                double temp_stat[NX];
                for (int_t j = 0; j < nx[i]; j++) temp_stat[j] = 0.0;
                daxpy_3l(NX, -1.0, lam_n[i], temp_stat);
                dgemv_n_3l(NX, NX, pQ[i], NX, w[i], res_stat[i]);

                // Print residuals
                printf("Stationarity residuals = %f\n",
                       twonormv(NX, res_stat[i]));
            }
            i = NN;
            for (int_t j = 0; j < nx[i] + nu[i]; j++)
                res_stat[i][j] = -pi_n[i][j];
            dgemv_n_3l(NX, NX, pQ[i], NX, w[i], res_stat[i]);
            double temp_stat[NX];
            for (int_t j = 0; j < nx[i]; j++) temp_stat[j] = 0.0;
            daxpy_3l(NX, -1.0, lam_n[i], temp_stat);
            dgemv_n_3l(NX, NX, pQ[i], NX, w[i], res_stat[i]);

            // Print residuals
            printf("Stationarity residuals = %f\n", twonormv(NX, res_stat[i]));

            int status =
                ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, work);
            if (status) {
                printf("qpOASES returned error status %d\n", status);
                return -1;
            }
            i = 0;

            for (int_t j = 0; j < NU; j++) w[i][j + NX] += qp_out.u[i][j];
            for (int_t j = 0; j < NX; j++) pi_n[i][j] += qp_out.pi[i][j];
            for (int_t j = 0; j < 2 * (NBU); j++)
                lam_n[i][j] += qp_out.lam[i][j];
            // printf("checkpoint=%i\n",status);

            for (i = 1; i < NN; i++) {
                for (int_t j = 0; j < NX; j++) w[i][j] += qp_out.x[i][j];
                for (int_t j = 0; j < NU; j++) w[i][j + NX] += qp_out.u[i][j];
                printf("x_step=%f\n", qp_out.x[i][0]);

                for (int_t j = 0; j < NX; j++) pi_n[i][j] += qp_out.pi[i][j];
                for (int_t j = 0; j < 2 * (NBX + NBU); j++)
                    lam_n[i][j] += qp_out.lam[i][j];
            }
            i = NN;
            for (int_t j = 0; j < 2 * (NBU); j++)
                lam_n[i][j] += qp_out.lam[i][j];
            for (int_t j = 0; j < NX; j++) w[i][j] += qp_out.x[i][j];
            printf("x_step=%f\n", qp_out.x[i][0]);
        }
    }
#ifdef DEBUG
    print_states_controls(&w[0], NN);
#endif  // DEBUG
    total_time = acados_toc(&timer);
    printf("Average of %.3f ms per iteration.\n", 1e3*total_time/timing_iters);

    free(sim_in.x);
    free(sim_in.u);
    free(sim_in.S_forw);
    free(sim_out.xn);
    free(sim_out.S_forw);

    return 0;
}
