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

#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#define PARALLEL 0

// Simple fixed-step Gauss-Newton based SQP routine
int_t ocp_nlp_gn_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out,
    void *nlp_mem_, void *nlp_work_) {

    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost*) nlp_in->cost;
    sim_solver *sim = nlp_in->sim;
    ocp_nlp_mem *nlp_mem = (ocp_nlp_mem*) nlp_mem_;
    ocp_nlp_work *work = (ocp_nlp_work*) nlp_work_;

    int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;

    real_t *w = work->w;
    real_t **y_ref = cost->y_ref;

    real_t **qp_Q = work->Q;
    real_t **qp_S = work->S;
    real_t **qp_R = work->R;
    // TODO(rien): only for least squares cost with state and control reference atm
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            for (int_t k = 0; k < nx[i]; k++) {
                qp_Q[i][j*nx[i]+k] = cost->W[i][j*(nx[i]+nu[i])+k];
            }
            for (int_t k = 0; k < nu[i]; k++) {
                qp_S[i][j*nu[i]+k] = cost->W[i][j*(nx[i]+nu[i])+nx[i]+k];
            }
        }
        for (int_t j = 0; j < nu[i]; j++) {
            for (int_t k = 0; k < nu[i]; k++) {
                qp_R[i][j*nu[i]+k] = cost->W[i][(nx[i]+j)*(nx[i]+nu[i])+nx[i]+k];
            }
        }
    }
    for (int_t j = 0; j < nx[N]; j++) {
        for (int_t k = 0; k < nx[N]; k++) {
            qp_Q[N][j*nx[N]+k] = cost->W[N][j*nx[N]+k];
        }
    }
    real_t **qp_q = work->q;
    real_t **qp_r = work->r;
    real_t **qp_lb = work->lb;
    real_t **qp_ub = work->ub;

    // Initialization of trajectories:
    int_t w_idx = 0;
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx+j] = nlp_mem->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            w[w_idx+nx[i]+j] = nlp_mem->u[i][j];
        }
        w_idx += nx[i]+nu[i];
    }
    for (int_t j = 0; j < nx[N]; j++) {
        w[w_idx+j] = nlp_mem->x[N][j];
    }

    acado_timer timer;
    real_t timings = 0;
    real_t timings_sim = 0;
    real_t timings_la = 0;
    real_t timings_ad = 0;

    real_t feas, stepX, stepU;
    int_t status;

    acado_tic(&timer);
    for (int_t sqp_iter = 0; sqp_iter < nlp_in->maxIter; sqp_iter++) {
        feas = stepX = stepU = -1e10;
#if PARALLEL
#pragma omp parallel for
#endif
        w_idx = 0;
        for (int_t i = 0; i < N; i++) {
            // Pass state and control to integrator
            for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = w[w_idx+j];
            for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = w[w_idx+nx[i]+j];
            sim[i].fun(sim[i].in, sim[i].out, sim[i].mem, sim[i].work);

            // TODO(rien): transition functions for changing dimensions not yet implemented!
            for (int_t j = 0; j < nx[i]; j++) {
                work->b[i][j] = sim[i].out->xn[j] - w[w_idx+nx[i]+nu[i]+j];
                if (fabs(work->b[i][j]) > feas)
                    feas = fabs(work->b[i][j]);
                for (int_t k = 0; k < nx[i]; k++)
                    work->A[i][j*nx[i]+k] = sim[i].out->S_forw[j*nx[i]+k];  // COLUMN MAJOR
            }
            for (int_t j = 0; j < nu[i]; j++)
                for (int_t k = 0; k < nx[i]; k++)
                    work->B[i][j*nx[i]+k] = sim[i].out->S_forw[(nx[i]+j)*nx[i]+k];  // COLUMN MAJOR

            timings_sim += sim[i].out->info->CPUtime;
            timings_la += sim[i].out->info->LAtime;
            timings_ad += sim[i].out->info->ADtime;
            w_idx += nx[i]+nu[i];
        }
        w_idx = 0;
        for (int_t i = 0; i < N+1; i++) {
            // Update bounds:
            for (int_t j = 0; j < nlp_in->nb[i]; j++) {
                qp_lb[i][j] = nlp_in->lb[i][j] - w[w_idx+nlp_in->idxb[i][j]];
                qp_ub[i][j] = nlp_in->ub[i][j] - w[w_idx+nlp_in->idxb[i][j]];
            }
//            print_matrix_name((char*)"stdout", (char*)"qp_lb: ", work->solver->qp_in->lb[i],
//            1, work->solver->qp_in->nb[i]);
//            print_matrix_name((char*)"stdout", (char*)"qp_ub: ", work->solver->qp_in->ub[i],
//            1, work->solver->qp_in->nb[i]);
//
//            print_matrix_name((char*)"stdout", (char*)"nlp_lb: ", nlp_in->lb[i],
//            1, nlp_in->nb[i]);
//            print_matrix_name((char*)"stdout", (char*)"nlp_ub: ", nlp_in->ub[i],
//            1, nlp_in->nb[i]);

            // Update gradients
            // TODO(rien): only for diagonal Q, R matrices atm
            // TODO(rien): only for least squares cost with state and control reference atm
            if (i < N) {
                for (int_t j = 0; j < nx[i]; j++) {
                    qp_q[i][j] = cost->W[i][j*(nx[i]+nu[i]+1)]*(w[w_idx+j]-y_ref[i][j]);
                }
                for (int_t j = 0; j < nu[i]; j++) {
                    qp_r[i][j] = cost->W[i][(nx[i]+j)*(nx[i]+nu[i]+1)]
                                                    *(w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j]);
                }
                w_idx += nx[i]+nu[i];
            } else {
                for (int_t j = 0; j < nx[i]; j++) {
                    qp_q[N][j] = cost->W[N][j*(nx[i]+1)]*(w[w_idx+j]-y_ref[N][j]);
                }
                w_idx += nx[i];
            }
        }
//        printf("nb[N]: %d \n", work->solver->qp_in->nb[N]);
//        print_matrix_name((char*)"stdout", (char*)"qp_lb[N]", qp_lb[N], 1, nx[N]);
//        print_matrix_name((char*)"stdout", (char*)"qp_ub[N]", qp_ub[N], 1, nx[N]);

        status = work->solver->fun(work->solver->qp_in, work->solver->qp_out,
                work->solver->mem, NULL, work->solver->work);
        if (status) {
            printf("qpOASES returned error status %d\n", status);
            return -1;
        }
        w_idx = 0;
        for (int_t i = 0; i < N; i++) {
            for (int_t j = 0; j < nx[i]; j++) {
                w[w_idx+j] += work->solver->qp_out->x[i][j];
                if (fabs(work->solver->qp_out->x[i][j]) > stepX)
                    stepX = fabs(work->solver->qp_out->x[i][j]);
            }
            for (int_t j = 0; j < nu[i]; j++) {
                w[w_idx+nx[i]+j] += work->solver->qp_out->u[i][j];
                if (fabs(work->solver->qp_out->u[i][j]) > stepU)
                    stepU = fabs(work->solver->qp_out->u[i][j]);
            }
            w_idx += nx[i]+nu[i];
//            print_matrix_name((char*)"stdout", (char*)"solver->qp_out->x[i]: ",
//              work->solver->qp_out->x[i], 1, nx[i]);
//            print_matrix_name((char*)"stdout", (char*)"solver->qp_out->u[i]: ",
//              work->solver->qp_out->u[i], 1, nu[i]);
        }
        for (int_t j = 0; j < nx[N]; j++) {
            w[w_idx+j] += work->solver->qp_out->x[N][j];
            if (fabs(work->solver->qp_out->x[N][j]) > stepX)
                stepX = fabs(work->solver->qp_out->x[N][j]);
        }
//        print_matrix_name((char*)"stdout", (char*)"solver->qp_out->x[N]: ",
    //        work->solver->qp_out->x[N], 1, nx[N]);
//        w_idx += nx[N];
//        print_matrix_name((char*)"stdout", (char*)"w_cur: ", w, 1, w_idx);

        fprintf(stdout, "--- ITERATION %d, Infeasibility: %+.3e , step X: %+.3e, "
                        "step U: %+.3e \n", sqp_iter, feas, stepX, stepU);
    }
    timings += acado_toc(&timer);

    printf("\nAverage of %.3f ms in the integrator,\n",
            1e3*timings_sim/(nlp_in->maxIter));
    printf("  of which %.3f ms spent in CasADi and\n",
            1e3*timings_ad/(nlp_in->maxIter));
    printf("  of which %.3f ms spent in BLASFEO.\n",
            1e3*timings_la/(nlp_in->maxIter));
    printf("--Total of %.3f ms per SQP iteration.--\n",
            1e3*timings/(nlp_in->maxIter));

    // Store trajectories:
    w_idx = 0;
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            nlp_mem->x[i][j] = w[w_idx+j];
            nlp_out->x[i][j] = w[w_idx+j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            nlp_mem->u[i][j] = w[w_idx+nx[i]+j];
            nlp_out->u[i][j] = w[w_idx+nx[i]+j];
        }
        w_idx += nx[i]+nu[i];
    }
    for (int_t j = 0; j < nx[N]; j++) {
        nlp_mem->x[N][j] = w[w_idx+j];
        nlp_out->x[N][j] = w[w_idx+j];
    }

    return 0;
}

void ocp_nlp_sqp_create_workspace(const ocp_nlp_in *in, ocp_nlp_work *work) {
    int_t num_vars = 0;
    work->A = (real_t **) malloc(sizeof(*work->A) * (in->N));
    work->B = (real_t **) malloc(sizeof(*work->B) * (in->N));
    work->b = (real_t **) malloc(sizeof(*work->b) * (in->N));

    work->Q = (real_t **) malloc(sizeof(*work->Q) * (in->N+1));
    work->q = (real_t **) malloc(sizeof(*work->q) * (in->N+1));
    work->S = (real_t **) malloc(sizeof(*work->S) * (in->N));
    work->R = (real_t **) malloc(sizeof(*work->R) * (in->N));
    work->r = (real_t **) malloc(sizeof(*work->r) * (in->N));

    work->lb = (real_t **) malloc(sizeof(*work->lb) * (in->N+1));
    work->ub = (real_t **) malloc(sizeof(*work->ub) * (in->N+1));
    work->lc = (real_t **) malloc(sizeof(*work->lc) * (in->N+1));
    work->uc = (real_t **) malloc(sizeof(*work->uc) * (in->N+1));

    for (int_t i = 0; i < in->N; i++) {
        work->A[i] = (real_t *) malloc(sizeof(*work->A[i]) * (in->nx[i+1]*in->nx[i]));
        work->B[i] = (real_t *) malloc(sizeof(*work->B[i]) * (in->nx[i+1]*in->nu[i]));
        work->b[i] = (real_t *) malloc(sizeof(*work->b[i]) * (in->nx[i+1]));

        work->Q[i] = (real_t *) malloc(sizeof(*work->Q[i]) * (in->nx[i]*in->nx[i]));
        work->q[i] = (real_t *) malloc(sizeof(*work->q[i]) * (in->nx[i]));
        work->S[i] = (real_t *) malloc(sizeof(*work->S[i]) * (in->nx[i]*in->nu[i]));
        work->R[i] = (real_t *) malloc(sizeof(*work->R[i]) * (in->nu[i]*in->nu[i]));
        work->r[i] = (real_t *) malloc(sizeof(*work->r[i]) * (in->nu[i]));

        if (in->nb[i]) {
            work->lb[i] = (real_t *) malloc(sizeof(*work->lb[i]) * (in->nb[i]));
            work->ub[i] = (real_t *) malloc(sizeof(*work->ub[i]) * (in->nb[i]));
        }
        if (in->nc[i]) {
            work->lc[i] = (real_t *) malloc(sizeof(*work->lc[i]) * (in->nc[i]));
            work->uc[i] = (real_t *) malloc(sizeof(*work->uc[i]) * (in->nc[i]));
        }

        num_vars += in->nx[i] + in->nu[i];
    }
    num_vars += in->nx[in->N];
    work->Q[in->N] = (real_t *) malloc(sizeof(*work->Q[in->N]) * (in->nx[in->N]*in->nx[in->N]));
    work->q[in->N] = (real_t *) malloc(sizeof(*work->q[in->N]) * (in->nx[in->N]));
    if (in->nb[in->N]) {
        work->lb[in->N] = (real_t *) malloc(sizeof(*work->lb[in->N]) * (in->nb[in->N]));
        work->ub[in->N] = (real_t *) malloc(sizeof(*work->ub[in->N]) * (in->nb[in->N]));
    }
    if (in->nc[in->N]) {
        work->lc[in->N] = (real_t *) malloc(sizeof(*work->lc[in->N]) * (in->nc[in->N]));
        work->uc[in->N] = (real_t *) malloc(sizeof(*work->uc[in->N]) * (in->nc[in->N]));
    }

    work->w = (real_t *) malloc(sizeof(*work->w) * num_vars);

    // Set OCP QP variables
    work->solver->qp_in->N = in->N;
    work->solver->qp_in->nx = in->nx;
    work->solver->qp_in->nu = in->nu;
    work->solver->qp_in->nb = in->nb;
    work->solver->qp_in->nc = in->nc;

    work->solver->qp_in->lb = (const real_t **) work->lb;
    work->solver->qp_in->ub = (const real_t **) work->ub;
    work->solver->qp_in->idxb = in->idxb;

    work->solver->qp_in->Q = (const real_t **) work->Q;
    work->solver->qp_in->S = (const real_t **) work->S;
    work->solver->qp_in->R = (const real_t **) work->R;
    work->solver->qp_in->q = (const real_t **) work->q;
    work->solver->qp_in->r = (const real_t **) work->r;

    work->solver->qp_in->A = (const real_t **) work->A;
    work->solver->qp_in->B = (const real_t **) work->B;
    work->solver->qp_in->b = (const real_t **) work->b;

    work->solver->qp_out->x = (real_t **) malloc(sizeof(*work->solver->qp_out->x) * (in->N+1));
    work->solver->qp_out->u = (real_t **) malloc(sizeof(*work->solver->qp_out->u) * (in->N));
    for (int_t i = 0; i < in->N; i++) {
        work->solver->qp_out->x[i] = \
            (real_t *) malloc(sizeof(*work->solver->qp_out->x[i]) * (in->nx[i]));
        work->solver->qp_out->u[i] = \
            (real_t *) malloc(sizeof(*work->solver->qp_out->u[i]) * (in->nu[i]));
    }
    work->solver->qp_out->x[in->N] = (real_t *)
            malloc(sizeof(*work->solver->qp_out->x[in->N]) * (in->nx[in->N]));
}
