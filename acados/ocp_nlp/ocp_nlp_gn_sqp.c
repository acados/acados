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
#include <string.h>

#include "acados/ocp_qp/ocp_qp_common.h"
#ifdef OOQP
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#endif
#include "acados/ocp_qp/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#define PARALLEL 0

int_t ocp_nlp_gn_sqp_calculate_workspace_size(const ocp_nlp_in *in, void *args_) {
    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args*) args_;

    int_t size;

    size = sizeof(ocp_nlp_gn_sqp_work);
    size += ocp_nlp_calculate_workspace_size(in, args->common);
    return size;
}


static void ocp_nlp_gn_sqp_cast_workspace(ocp_nlp_gn_sqp_work *work,
    ocp_nlp_gn_sqp_memory *mem) {

    char *ptr = (char *)work;

    ptr += sizeof(ocp_nlp_gn_sqp_work);
    work->common = (ocp_nlp_work*) ptr;
    ocp_nlp_cast_workspace(work->common, mem->common);
}


// Simple fixed-step Gauss-Newton based SQP routine
int_t ocp_nlp_gn_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, void *nlp_args_,
    void *nlp_mem_, void *nlp_work_) {

    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost*) nlp_in->cost;
    sim_solver *sim = nlp_in->sim;
    ocp_nlp_gn_sqp_args *gn_sqp_args = (ocp_nlp_gn_sqp_args *) nlp_args_;
    ocp_nlp_gn_sqp_memory *gn_sqp_mem = (ocp_nlp_gn_sqp_memory *) nlp_mem_;
    ocp_nlp_gn_sqp_work *work = (ocp_nlp_gn_sqp_work*) nlp_work_;

    ocp_nlp_gn_sqp_cast_workspace(work, gn_sqp_mem);

    int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;

    real_t *w = work->common->w;
    real_t **y_ref = cost->y_ref;

    int_t **qp_idxb = (int_t **) gn_sqp_mem->qp_solver->qp_in->idxb;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nlp_in->nb[i]; j++) {
            qp_idxb[i][j] = nlp_in->idxb[i][j];
        }
    }
    real_t **qp_A = (real_t **) gn_sqp_mem->qp_solver->qp_in->A;
    real_t **qp_B = (real_t **) gn_sqp_mem->qp_solver->qp_in->B;
    real_t **qp_b = (real_t **) gn_sqp_mem->qp_solver->qp_in->b;
    real_t **qp_Q = (real_t **) gn_sqp_mem->qp_solver->qp_in->Q;
    real_t **qp_S = (real_t **) gn_sqp_mem->qp_solver->qp_in->S;
    real_t **qp_R = (real_t **) gn_sqp_mem->qp_solver->qp_in->R;
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
    real_t **qp_q = (real_t **) gn_sqp_mem->qp_solver->qp_in->q;
    real_t **qp_r = (real_t **) gn_sqp_mem->qp_solver->qp_in->r;
    real_t **qp_lb = (real_t **) gn_sqp_mem->qp_solver->qp_in->lb;
    real_t **qp_ub = (real_t **) gn_sqp_mem->qp_solver->qp_in->ub;

    // Initialization of trajectories:
    int_t w_idx = 0;
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx+j] = gn_sqp_mem->common->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            w[w_idx+nx[i]+j] = gn_sqp_mem->common->u[i][j];
        }
        w_idx += nx[i]+nu[i];
    }
    for (int_t j = 0; j < nx[N]; j++) {
        w[w_idx+j] = gn_sqp_mem->common->x[N][j];
    }

#ifdef MEASURE_TIMINGS
    acados_timer timer;
    real_t timings = 0;
    real_t timings_sim = 0;
    real_t timings_la = 0;
    real_t timings_ad = 0;
#endif

    real_t feas, stepX, stepU;
    int_t status;

#ifdef MEASURE_TIMINGS
    acados_tic(&timer);
#endif

    for (int_t sqp_iter = 0; sqp_iter < gn_sqp_args->common->maxIter; sqp_iter++) {
        feas = stepX = stepU = -1e10;
#if PARALLEL
#pragma omp parallel for
#endif
        w_idx = 0;
        for (int_t i = 0; i < N; i++) {
            sim_RK_opts *opts = (sim_RK_opts*) sim[i].args;
            // Pass state and control to integrator
            for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = w[w_idx+j];
            for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = w[w_idx+nx[i]+j];
            sim[i].in->sens_adj = (opts->scheme.type != exact && sqp_iter > 0);
            if (nlp_in->freezeSens) {  // freeze inexact sensitivities after first SQP iteration !!
                opts->scheme.freeze = sqp_iter > 0;
            }
            sim[i].fun(sim[i].in, sim[i].out, sim[i].args, sim[i].mem, sim[i].work);

            // TODO(rien): transition functions for changing dimensions not yet implemented!
            for (int_t j = 0; j < nx[i]; j++) {
                qp_b[i][j] = sim[i].out->xn[j] - w[w_idx+nx[i]+nu[i]+j];
                if (fabs(qp_b[i][j]) > feas)
                    feas = fabs(qp_b[i][j]);
                for (int_t k = 0; k < nx[i]; k++)
                    qp_A[i][j*nx[i]+k] = sim[i].out->S_forw[j*nx[i]+k];  // COLUMN MAJOR
            }
            for (int_t j = 0; j < nu[i]; j++)
                for (int_t k = 0; k < nx[i]; k++)
                    qp_B[i][j*nx[i]+k] = sim[i].out->S_forw[(nx[i]+j)*nx[i]+k];  // COLUMN MAJOR

            // printf("w\n");
            // print_matrix("stdout", &w[w_idx], nx[i]+nu[i], 1);
            // printf("A[%d]\n",i);
            // d_print_mat(nx[i], nx[i], qp_A[i], nx[i]);
            // printf("B[%d]\n",i);
            // d_print_mat(nx[i], nu[i], qp_B[i], nx[i]);

#ifdef MEASURE_TIMINGS
            timings_sim += sim[i].out->info->CPUtime;
            timings_la += sim[i].out->info->LAtime;
            timings_ad += sim[i].out->info->ADtime;
#endif

            w_idx += nx[i]+nu[i];
        }
        w_idx = 0;
        for (int_t i = 0; i < N+1; i++) {
            // Update bounds:
            for (int_t j = 0; j < nlp_in->nb[i]; j++) {
                qp_lb[i][j] = nlp_in->lb[i][j] - w[w_idx+nlp_in->idxb[i][j]];
                qp_ub[i][j] = nlp_in->ub[i][j] - w[w_idx+nlp_in->idxb[i][j]];
            }
        //    print_matrix_name((char*)"stdout", (char*)"qp_lb: ",
        //    gn_sqp_mem->qp_solver->qp_in->lb[i], 1, gn_sqp_mem->qp_solver->qp_in->nb[i]);
        //    print_matrix_name((char*)"stdout", (char*)"qp_ub: ",
        //    gn_sqp_mem->qp_solver->qp_in->ub[i], 1, gn_sqp_mem->qp_solver->qp_in->nb[i]);

        //    print_matrix_name((char*)"stdout", (char*)"nlp_lb: ", nlp_in->lb[i],
        //    1, nlp_in->nb[i]);
        //    print_matrix_name((char*)"stdout", (char*)"nlp_ub: ", nlp_in->ub[i],
        //    1, nlp_in->nb[i]);

            // Update gradients
            // TODO(rien): only for diagonal Q, R matrices atm
            // TODO(rien): only for least squares cost with state and control reference atm
            if (i < N) {
                sim_RK_opts *opts = (sim_RK_opts*) sim[i].args;
                for (int_t j = 0; j < nx[i]; j++) {
                    qp_q[i][j] = cost->W[i][j*(nx[i]+nu[i]+1)]*(w[w_idx+j]-y_ref[i][j]);
                    // adjoint-based gradient correction:
                    if (opts->scheme.type != exact) qp_q[i][j] += sim[i].out->grad[j];
                }
                for (int_t j = 0; j < nu[i]; j++) {
                    qp_r[i][j] = cost->W[i][(nx[i]+j)*(nx[i]+nu[i]+1)]
                                                    *(w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j]);
                    // adjoint-based gradient correction:
                    if (opts->scheme.type != exact) qp_r[i][j] += sim[i].out->grad[nx[i]+j];
                }
                w_idx += nx[i]+nu[i];
            } else {
                for (int_t j = 0; j < nx[i]; j++) {
                    qp_q[N][j] = cost->W[N][j*(nx[i]+1)]*(w[w_idx+j]-y_ref[N][j]);
                }
                w_idx += nx[i];
            }
        }
//        printf("nb[N]: %d \n", gn_sqp_mem->qp_solver->qp_in->nb[N]);
//        print_matrix_name((char*)"stdout", (char*)"qp_lb[N]", qp_lb[N], 1, nx[N]);
//        print_matrix_name((char*)"stdout", (char*)"qp_ub[N]", qp_ub[N], 1, nx[N]);

    //    gn_sqp_mem->qp_solver->initialize(gn_sqp_mem->qp_solver->qp_in,
    //     gn_sqp_mem->qp_solver->args, gn_sqp_mem->qp_solver->mem, gn_sqp_mem->qp_solver->work);
#ifdef DEBUG
        print_ocp_qp(gn_sqp_mem->qp_solver->qp_in);
#endif
        status = gn_sqp_mem->qp_solver->fun(gn_sqp_mem->qp_solver->qp_in,
            gn_sqp_mem->qp_solver->qp_out, gn_sqp_mem->qp_solver->args, gn_sqp_mem->qp_solver->mem,
            gn_sqp_mem->qp_solver->work);
        if (status) {
            printf("QP solver returned error status %d\n", status);
            return -1;
        }
        w_idx = 0;
        for (int_t i = 0; i < N; i++) {
            for (int_t j = 0; j < nx[i]; j++) {
                sim[i].in->S_adj[j] = -gn_sqp_mem->qp_solver->qp_out->pi[i][j];
            }
            for (int_t j = 0; j < nx[i]; j++) {
                w[w_idx+j] += gn_sqp_mem->qp_solver->qp_out->x[i][j];
                if (fabs(gn_sqp_mem->qp_solver->qp_out->x[i][j]) > stepX)
                    stepX = fabs(gn_sqp_mem->qp_solver->qp_out->x[i][j]);
            }
            for (int_t j = 0; j < nu[i]; j++) {
                w[w_idx+nx[i]+j] += gn_sqp_mem->qp_solver->qp_out->u[i][j];
                if (fabs(gn_sqp_mem->qp_solver->qp_out->u[i][j]) > stepU)
                    stepU = fabs(gn_sqp_mem->qp_solver->qp_out->u[i][j]);
            }
            w_idx += nx[i]+nu[i];
           print_matrix_name((char*)"stdout", (char*)"solver->qp_out->x[i]: ",
             gn_sqp_mem->qp_solver->qp_out->x[i], 1, nx[i]);
           print_matrix_name((char*)"stdout", (char*)"solver->qp_out->u[i]: ",
             gn_sqp_mem->qp_solver->qp_out->u[i], 1, nu[i]);
        }
        for (int_t j = 0; j < nx[N]; j++) {
            w[w_idx+j] += gn_sqp_mem->qp_solver->qp_out->x[N][j];
            if (fabs(gn_sqp_mem->qp_solver->qp_out->x[N][j]) > stepX)
                stepX = fabs(gn_sqp_mem->qp_solver->qp_out->x[N][j]);
        }
       print_matrix_name((char*)"stdout", (char*)"solver->qp_out->x[N]: ",
           gn_sqp_mem->qp_solver->qp_out->x[N], 1, nx[N]);

        fprintf(stdout, "--- ITERATION %d, Infeasibility: %+.3e , step X: %+.3e, "
                        "step U: %+.3e \n", sqp_iter, feas, stepX, stepU);
    }

#ifdef MEASURE_TIMINGS
    timings += acados_toc(&timer);
    printf("\nAverage of %.3f ms in the integrator,\n",
            1e3*timings_sim/(gn_sqp_args->common->maxIter));
    printf("  of which %.3f ms spent in CasADi and\n",
            1e3*timings_ad/(gn_sqp_args->common->maxIter));
    printf("  of which %.3f ms spent in BLASFEO.\n",
            1e3*timings_la/(gn_sqp_args->common->maxIter));
    printf("--Total of %.3f ms per SQP iteration.--\n",
            1e3*timings/(gn_sqp_args->common->maxIter));
#endif

    // Store trajectories:
    w_idx = 0;
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            gn_sqp_mem->common->x[i][j] = w[w_idx+j];
            nlp_out->x[i][j] = w[w_idx+j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            gn_sqp_mem->common->u[i][j] = w[w_idx+nx[i]+j];
            nlp_out->u[i][j] = w[w_idx+nx[i]+j];
        }
        w_idx += nx[i]+nu[i];
    }
    for (int_t j = 0; j < nx[N]; j++) {
        gn_sqp_mem->common->x[N][j] = w[w_idx+j];
        nlp_out->x[N][j] = w[w_idx+j];
    }

    return 0;
}


void ocp_nlp_gn_sqp_create_memory(const ocp_nlp_in *in, void *args_, void *memory_) {
    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args *) args_;
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *) memory_;

    mem->qp_solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver));
    ocp_qp_in *qp_in = (ocp_qp_in *) malloc(sizeof(ocp_qp_in));
    allocate_ocp_qp_in(in->N, in->nx, in->nu, in->nb, in->nc, qp_in);
    ocp_qp_out *qp_out = (ocp_qp_out *) malloc(sizeof(ocp_qp_out));
    allocate_ocp_qp_out(qp_in, qp_out);
    void *qp_args = NULL, *qp_mem = NULL, *qp_work = NULL;
    if (!strcmp(args->qp_solver_name, "qpdunes")) {
        mem->qp_solver->fun = &ocp_qp_qpdunes;
        mem->qp_solver->initialize = &ocp_qp_qpdunes_initialize;
        mem->qp_solver->destroy = &ocp_qp_qpdunes_destroy;
        qp_args = (void *) malloc(sizeof(ocp_qp_qpdunes_args));
        qp_mem = (void *) malloc(sizeof(ocp_qp_qpdunes_memory));
    #ifdef OOQP
    } else if (!strcmp(args->qp_solver_name, "ooqp")) {
        mem->qp_solver->fun = &ocp_qp_ooqp;
        mem->qp_solver->initialize = &ocp_qp_ooqp_initialize;
        mem->qp_solver->destroy = &ocp_qp_ooqp_destroy;
        qp_args = (void *) malloc(sizeof(ocp_qp_ooqp_args));
        qp_mem = (void *) malloc(sizeof(ocp_qp_ooqp_memory));
    #endif
    } else if (!strcmp(args->qp_solver_name, "condensing_qpoases")) {
        mem->qp_solver->fun = &ocp_qp_condensing_qpoases;
        mem->qp_solver->initialize = &ocp_qp_condensing_qpoases_initialize;
        mem->qp_solver->destroy = &ocp_qp_condensing_qpoases_destroy;
        qp_args = (void *) malloc(sizeof(ocp_qp_condensing_qpoases_args));
    } else if (!strcmp(args->qp_solver_name, "hpmpc")) {
        mem->qp_solver->fun = &ocp_qp_hpmpc;
        mem->qp_solver->initialize = &ocp_qp_hpmpc_initialize;
        mem->qp_solver->destroy = &ocp_qp_hpmpc_destroy;
        qp_args = (void *) malloc(sizeof(ocp_qp_hpmpc_args));
    } else {
        printf("CHOSEN QP SOLVER FOR SQP METHOD NOT AVAILABLE!\n");
        exit(1);
    }
    mem->qp_solver->initialize(qp_in, qp_args, qp_mem, &qp_work);
    mem->qp_solver->qp_in = qp_in;
    mem->qp_solver->qp_out = qp_out;
    mem->qp_solver->args = qp_args;
    mem->qp_solver->mem = qp_mem;
    mem->qp_solver->work = qp_work;

    ocp_nlp_create_memory(in, mem->common);
}


void ocp_nlp_gn_sqp_free_memory(void *mem_) {
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *) mem_;

    int_t N = mem->qp_solver->qp_in->N;
    ocp_nlp_free_memory(N, mem->common);

    mem->qp_solver->destroy(mem->qp_solver->mem, mem->qp_solver->work);
    free_ocp_qp_in(mem->qp_solver->qp_in);
    free_ocp_qp_out(mem->qp_solver->qp_out);
    free(mem->qp_solver->qp_in);
    free(mem->qp_solver->qp_out);
    free(mem->qp_solver->args);
    free(mem->qp_solver->mem);
    free(mem->qp_solver);
    // TODO(dimitris): where do we free the integrators?
}
