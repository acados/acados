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
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_qp/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

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
    work->common = (ocp_nlp_work *)ptr;
    ocp_nlp_cast_workspace(work->common, mem->common);
}


static void initialize_objective(
    const ocp_nlp_in *nlp_in,
    ocp_nlp_gn_sqp_memory *gn_sqp_mem,
    ocp_nlp_gn_sqp_work *work) {

    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost*) nlp_in->cost;

    real_t **qp_Q = (real_t **) gn_sqp_mem->qp_solver->qp_in->Q;
    real_t **qp_S = (real_t **) gn_sqp_mem->qp_solver->qp_in->S;
    real_t **qp_R = (real_t **) gn_sqp_mem->qp_solver->qp_in->R;
    // TODO(rien): only for least squares cost with state and control reference atm
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            for (int_t k = 0; k < nx[i]; k++) {
                qp_Q[i][j * nx[i] + k] = cost->W[i][j * (nx[i] + nu[i]) + k];
            }
            for (int_t k = 0; k < nu[i]; k++) {
                qp_S[i][j * nu[i] + k] =
                    cost->W[i][j * (nx[i] + nu[i]) + nx[i] + k];
            }
        }
        for (int_t j = 0; j < nu[i]; j++) {
            for (int_t k = 0; k < nu[i]; k++) {
                qp_R[i][j * nu[i] + k] =
                    cost->W[i][(nx[i] + j) * (nx[i] + nu[i]) + nx[i] + k];
            }
        }
    }
}


static void initialize_trajectories(
    const ocp_nlp_in *nlp_in,
    ocp_nlp_gn_sqp_memory *gn_sqp_mem,
    ocp_nlp_gn_sqp_work *work) {

    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
    real_t *w = work->common->w;

    int_t w_idx = 0;
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx + j] = gn_sqp_mem->common->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            w[w_idx + nx[i] + j] = gn_sqp_mem->common->u[i][j];
        }
        w_idx += nx[i] + nu[i];
    }
    for (int_t j = 0; j < nx[N]; j++) {
        w[w_idx + j] = gn_sqp_mem->common->x[N][j];
    }
}


static void multiple_shooting(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *mem, real_t *w) {

    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;
    sim_solver *sim = nlp->sim;
    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *) nlp->cost;
    real_t **y_ref = cost->y_ref;

    real_t **qp_A = (real_t **) mem->qp_solver->qp_in->A;
    real_t **qp_B = (real_t **) mem->qp_solver->qp_in->B;
    real_t **qp_b = (real_t **) mem->qp_solver->qp_in->b;
    real_t **qp_q = (real_t **) mem->qp_solver->qp_in->q;
    real_t **qp_r = (real_t **) mem->qp_solver->qp_in->r;
    real_t **qp_lb = (real_t **) mem->qp_solver->qp_in->lb;
    real_t **qp_ub = (real_t **) mem->qp_solver->qp_in->ub;

    int_t w_idx = 0;

    for (int_t i = 0; i < N; i++) {
        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = w[w_idx+j];
        for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = w[w_idx+nx[i]+j];
        sim[i].fun(sim[i].in, sim[i].out, sim[i].args, sim[i].mem, sim[i].work);

        // TODO(rien): transition functions for changing dimensions not yet implemented!
        for (int_t j = 0; j < nx[i]; j++) {
            qp_b[i][j] = sim[i].out->xn[j] - w[w_idx+nx[i]+nu[i]+j];
            for (int_t k = 0; k < nx[i]; k++)
                qp_A[i][j*nx[i]+k] = sim[i].out->S_forw[j*nx[i]+k];
        }
        for (int_t j = 0; j < nu[i]; j++)
            for (int_t k = 0; k < nx[i]; k++)
                qp_B[i][j*nx[i]+k] = sim[i].out->S_forw[(nx[i]+j)*nx[i]+k];

        // Update bounds:
        for (int_t j = 0; j < nlp->nb[i]; j++) {
            if (nlp->idxb[i][j] < nu[i]) {
                qp_lb[i][j] = nlp->lb[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
                qp_ub[i][j] = nlp->ub[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
            } else {
                qp_lb[i][j] = nlp->lb[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
                qp_ub[i][j] = nlp->ub[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
            }
        }

        // Update gradients
        // TODO(rien): only for diagonal Q, R matrices atm
        // TODO(rien): only for least squares cost with state and control reference atm
        sim_RK_opts *opts = (sim_RK_opts*) sim[i].args;
        for (int_t j = 0; j < nx[i]; j++) {
            qp_q[i][j] = cost->W[i][j*(nx[i]+nu[i]+1)]*(w[w_idx+j]-y_ref[i][j]);
            // adjoint-based gradient correction:
            if (opts->scheme.type != exact) qp_q[i][j] += sim[i].out->grad[j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            qp_r[i][j] = cost->W[i][(nx[i]+j)*(nx[i]+nu[i]+1)]*(w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j]);
            // adjoint-based gradient correction:
            if (opts->scheme.type != exact) qp_r[i][j] += sim[i].out->grad[nx[i]+j];
        }
        w_idx += nx[i]+nu[i];
    }

    for (int_t j = 0; j < nlp->nb[N]; j++) {
        if (nlp->idxb[N][j] < nu[N]) {
            qp_lb[N][j] = nlp->lb[N][j] - w[w_idx + nx[N] + nlp->idxb[N][j]];
            qp_ub[N][j] = nlp->ub[N][j] - w[w_idx + nx[N] + nlp->idxb[N][j]];
        } else {
            qp_lb[N][j] = nlp->lb[N][j] - w[w_idx - nu[N] + nlp->idxb[N][j]];
            qp_ub[N][j] = nlp->ub[N][j] - w[w_idx - nu[N] + nlp->idxb[N][j]];
        }
    }

    for (int_t j = 0; j < nx[N]; j++)
        qp_q[N][j] = cost->W[N][j*(nx[N]+nu[N]+1)]*(w[w_idx+j]-y_ref[N][j]);
}


static void update_variables(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *mem, real_t *w) {
    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;
    sim_solver *sim = nlp->sim;

    int_t w_idx = 0;
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            sim[i].in->S_adj[j] = -mem->qp_solver->qp_out->pi[i][j];
            w[w_idx+j] += mem->qp_solver->qp_out->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++)
            w[w_idx+nx[i]+j] += mem->qp_solver->qp_out->u[i][j];
        w_idx += nx[i]+nu[i];
    }
    for (int_t j = 0; j < nx[N]; j++)
        w[w_idx+j] += mem->qp_solver->qp_out->x[N][j];
}


static void store_trajectories(const ocp_nlp_in *nlp, ocp_nlp_memory *memory, ocp_nlp_out *out,
    real_t *w) {

    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            memory->x[i][j] = w[w_idx+j];
            out->x[i][j] = w[w_idx+j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            memory->u[i][j] = w[w_idx+nx[i]+j];
            out->u[i][j] = w[w_idx+nx[i]+j];
        }
        w_idx += nx[i] + nu[i];
    }
}


// Simple fixed-step Gauss-Newton based SQP routine
int_t ocp_nlp_gn_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, void *nlp_args_,
    void *nlp_mem_, void *nlp_work_) {

    ocp_nlp_gn_sqp_memory *gn_sqp_mem = (ocp_nlp_gn_sqp_memory *) nlp_mem_;
    ocp_nlp_gn_sqp_work *work = (ocp_nlp_gn_sqp_work*) nlp_work_;
    ocp_nlp_gn_sqp_cast_workspace(work, gn_sqp_mem);

    initialize_objective(nlp_in, gn_sqp_mem, work);
    initialize_trajectories(nlp_in, gn_sqp_mem, work);

    int_t **qp_idxb = (int_t **) gn_sqp_mem->qp_solver->qp_in->idxb;
    for (int_t i = 0; i <= nlp_in->N; i++) {
        for (int_t j = 0; j < nlp_in->nb[i]; j++) {
            qp_idxb[i][j] = nlp_in->idxb[i][j];
        }
    }

    int_t max_sqp_iterations = ((ocp_nlp_gn_sqp_args *) nlp_args_)->common->maxIter;

    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iterations; sqp_iter++) {

        multiple_shooting(nlp_in, gn_sqp_mem, work->common->w);

        int_t qp_status = gn_sqp_mem->qp_solver->fun(
            gn_sqp_mem->qp_solver->qp_in,
            gn_sqp_mem->qp_solver->qp_out,
            gn_sqp_mem->qp_solver->args,
            gn_sqp_mem->qp_solver->mem,
            gn_sqp_mem->qp_solver->work);

        if (qp_status != 0) {
            printf("QP solver returned error status %d\n", qp_status);
            return -1;
        }

        update_variables(nlp_in, gn_sqp_mem, work->common->w);

        for (int_t i = 0; i < nlp_in->N; i++) {
            sim_RK_opts *opts = (sim_RK_opts*) nlp_in->sim[i].args;
            nlp_in->sim[i].in->sens_adj = (opts->scheme.type != exact);
            if (nlp_in->freezeSens) {  // freeze inexact sensitivities after first SQP iteration !!
                opts->scheme.freeze = true;
            }
        }
    }

    total_time += acados_toc(&timer);
    store_trajectories(nlp_in, gn_sqp_mem->common, nlp_out, work->common->w);
    return 0;
}

void ocp_nlp_gn_sqp_create_memory(const ocp_nlp_in *in, void *args_,
                                  void *memory_) {
    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args *)args_;
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)memory_;

    mem->qp_solver = (ocp_qp_solver *)malloc(sizeof(ocp_qp_solver));
    ocp_qp_in *qp_in = (ocp_qp_in *)malloc(sizeof(ocp_qp_in));
    allocate_ocp_qp_in(in->N, in->nx, in->nu, in->nb, in->nc, qp_in);
    ocp_qp_out *qp_out = (ocp_qp_out *)malloc(sizeof(ocp_qp_out));
    allocate_ocp_qp_out(qp_in, qp_out);
    void *qp_args = NULL, *qp_mem = NULL, *qp_work = NULL;
    if (!strcmp(args->qp_solver_name, "qpdunes")) {
        mem->qp_solver->fun = &ocp_qp_qpdunes;
        mem->qp_solver->initialize = &ocp_qp_qpdunes_initialize;
        mem->qp_solver->destroy = &ocp_qp_qpdunes_destroy;
        qp_args = (void *)malloc(sizeof(ocp_qp_qpdunes_args));
        qp_mem = (void *)malloc(sizeof(ocp_qp_qpdunes_memory));
#ifdef OOQP
    } else if (!strcmp(args->qp_solver_name, "ooqp")) {
        mem->qp_solver->fun = &ocp_qp_ooqp;
        mem->qp_solver->initialize = &ocp_qp_ooqp_initialize;
        mem->qp_solver->destroy = &ocp_qp_ooqp_destroy;
        qp_args = (void *)malloc(sizeof(ocp_qp_ooqp_args));
        qp_mem = (void *)malloc(sizeof(ocp_qp_ooqp_memory));
#endif
    } else if (!strcmp(args->qp_solver_name, "condensing_qpoases")) {
        mem->qp_solver->fun = &ocp_qp_condensing_qpoases;
        mem->qp_solver->initialize = &ocp_qp_condensing_qpoases_initialize;
        mem->qp_solver->destroy = &ocp_qp_condensing_qpoases_destroy;
        qp_args = (void *)malloc(sizeof(ocp_qp_condensing_qpoases_args));
    } else if (!strcmp(args->qp_solver_name, "hpmpc")) {
        mem->qp_solver->fun = &ocp_qp_hpmpc;
        mem->qp_solver->initialize = &ocp_qp_hpmpc_initialize;
        mem->qp_solver->destroy = &ocp_qp_hpmpc_destroy;
        qp_args = (void *)malloc(sizeof(ocp_qp_hpmpc_args));
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
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)mem_;

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
