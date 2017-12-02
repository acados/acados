
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

// system headers
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// flush denormals to zero
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) ||  \
    defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || \
    defined(TARGET_AMD_SSE3)
#include <xmmintrin.h>  // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

#include "hpmpc/include/aux_d.h"

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/utils/timing.h"
#include "acados/utils/math.h"
#include "examples/c/chen_model/chen_model.h"

// define IP solver arguments && number of repetitions
#define NREP 1000
#define MAXITER 10
#define TOL 1e-8
#define MINSTEP 1e-8

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
    // Problem data
    real_t x0[NX] = {0.05, 0};
    real_t w[NN * (NX + NU) + NX] = {0};  // States and controls stacked
    real_t Q[NX * NX] = {0};
    real_t R[NU * NU] = {0};
    real_t xref[NX] = {0};
    real_t uref[NX] = {0};
    int_t max_sqp_iters = 1;
    int_t max_iters = 10000;
    real_t x_end[NX] = {0};
    real_t u_end[NU] = {0};

    for (int_t i = 0; i < NX; i++) Q[i * (NX + 1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i * (NU + 1)] = 0.05;

    // Integrator structs
    real_t T = 0.01;
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

    sim_erk_workspace erk_work;
    sim_rk_opts rk_opts;
    sim_erk_create_arguments(&rk_opts, 4);
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
    int *hidxb[NN + 1];
    real_t *pC[NN + 1];
    real_t *pD[NN];
    real_t *plg[NN + 1];
    real_t *pug[NN + 1];

    /************************************************
     * box constraints
     ************************************************/
    int ii;

    nb[0] = 0;
    for (ii = 1; ii < NN; ii++) nb[ii] = 0;
    nb[NN] = 0;

    int *idxb0;
    int_zeros(&idxb0, nb[0], 1);
    idxb0[0] = 0;
    idxb0[1] = 1;

    int *idxb1;
    int_zeros(&idxb1, nb[1], 1);

    int *idxbN;
    int_zeros(&idxbN, nb[NN], 1);

    for (ii = 1; ii < NN; ii++) {
        hidxb[ii] = idxb1;
    }

    hidxb[NN] = idxbN;

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
    // hidxb[NN] = idxbN;

    nx[0] = 0;
    //    d_print_mat(nx, 1, b, nx);
    //    d_print_mat(nx, 1, b0, nx);

    // then A0 is a matrix of size 0x0
    double *A0;
    d_zeros(&A0, 0, 0);

    int ng = 0;   // 4;  // number of general constraints
    int ngN = 0;  // 4;  // number of general constraints at the last stage

    int ngg[NN + 1];
    for (ii = 0; ii < NN; ii++) ngg[ii] = ng;
    ngg[NN] = ngN;

    /************************************************
     * general constraints
     ************************************************/

    double *C0;
    d_zeros(&C0, 0, NX);
    double *D0;
    d_zeros(&D0, 0, NU);
    double *lg0;
    d_zeros(&lg0, 0, 1);
    double *ug0;
    d_zeros(&ug0, 0, 1);

    double *C;
    d_zeros(&C, 0, NX);
    double *D;
    d_zeros(&D, 0, NU);
    double *lg;
    d_zeros(&lg, 0, 1);
    double *ug;
    d_zeros(&ug, 0, 1);

    double *CN;
    d_zeros(&CN, ngN, NX);
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

    pC[0] = C0;
    pD[0] = D0;
    plg[0] = lg0;
    pug[0] = ug0;
    real_t *ppi[NN];
    real_t *plam[NN + 1];

    ii = 0;
    d_zeros(&ppi[ii], nx[ii + 1], 1);
    d_zeros(&plam[ii], 2 * nb[ii] + 2 * nb[ii], 1);
    for (ii = 1; ii < NN; ii++) {
        pC[ii] = C;
        pD[ii] = D;
        plg[ii] = lg;
        pug[ii] = ug;
        d_zeros(&ppi[ii], nx[ii + 1], 1);
        d_zeros(&plam[ii], 2 * nb[ii] + 2 * nb[ii], 1);
    }
    pC[NN] = CN;
    plg[NN] = lgN;
    pug[NN] = ugN;

    /************************************************
     * solver arguments
     ************************************************/

    // solver arguments
    ocp_qp_hpmpc_args hpmpc_args;
    hpmpc_args.tol = TOL;
    hpmpc_args.max_iter = MAXITER;
    //  hpmpc_args.min_step = MINSTEP;
    hpmpc_args.mu0 = 0.0;
    //  hpmpc_args.sigma_min = 1e-3;
    hpmpc_args.warm_start = 0;
    hpmpc_args.N2 = NN;

    /************************************************
     * work space
     ************************************************/

    int work_space_size = ocp_qp_hpmpc_workspace_size_bytes(NN, nx, nu, nb, ngg,
                                                            hidxb, &hpmpc_args);
    printf("\nwork space size: %d bytes\n", work_space_size);

    // void *workspace = malloc(work_space_size);
    void *workspace = malloc(500000);

    // double workspace[500000];
    // Allocate OCP QP variables
    ocp_qp_in qp_in;
    qp_in.N = NN;
    ocp_qp_out qp_out;
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
    // qp_in.lb = (const real_t **) px0;
    // qp_in.ub = (const real_t **) px0;
    // qp_in.idxb = (const int_t **) hidxb;
    qp_in.Cx = (const real_t **)pC;
    qp_in.Cu = (const real_t **)pD;
    qp_in.lc = (const real_t **)plg;
    qp_in.uc = (const real_t **)pug;

    qp_out.x = px;
    qp_out.u = pu;
    qp_out.pi = ppi;
    qp_out.lam = plam;

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
                sim_erk(&sim_in, &sim_out, &rk_opts, 0, &erk_work);

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

            dgemv_n_3l(NX, NX, pA[0], NX, x0, pb[0]);

            for (int_t j = 0; j < NX; j++) {
                pq[NN][j] = Q[j] * (w[NN * (NX + NU) + j] - xref[j]);
            }
            int status = ocp_qp_hpmpc(&qp_in, &qp_out, &hpmpc_args, workspace);
            // int status = 0;
            // printf("hpmpc_status=%i\n",status);
            if (status == 1) printf("status = ACADOS_MAXITER\n");

            if (status == 2) printf("status = ACADOS_MINSTEP\n");

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
    total_time = acados_toc(&timer);
    printf("Average of %.3f ms per iteration.\n", 1e3*total_time/max_iters);
    // free(workspace);
    return 0;
}
