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

#include "acados/ocp_nlp/ocp_nlp_sqp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

void prepare_qp(const ocp_nlp_in *nlp_in, ocp_nlp_sqp_args *sqp_args, ocp_nlp_sqp_memory *sqp_mem) {
    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
    const int_t *nb = nlp_in->nb;
    const int_t **idxb = nlp_in->idxb;
    const int_t *ng = nlp_in->ng;

    real_t **qp_A = (real_t **)sqp_args->qp_solver->qp_in->A;
    real_t **qp_B = (real_t **)sqp_args->qp_solver->qp_in->B;
    real_t **qp_b = (real_t **)sqp_args->qp_solver->qp_in->b;
    real_t **qp_q = (real_t **)sqp_args->qp_solver->qp_in->q;
    real_t **qp_r = (real_t **)sqp_args->qp_solver->qp_in->r;
    real_t **qp_R = (real_t **)sqp_args->qp_solver->qp_in->R;
    real_t **qp_Q = (real_t **)sqp_args->qp_solver->qp_in->Q;
    real_t **qp_S = (real_t **)sqp_args->qp_solver->qp_in->S;
    real_t **qp_lb = (real_t **)sqp_args->qp_solver->qp_in->lb;
    real_t **qp_ub = (real_t **)sqp_args->qp_solver->qp_in->ub;
    real_t **qp_Cx = (real_t **)sqp_args->qp_solver->qp_in->Cx;
    real_t **qp_Cu = (real_t **)sqp_args->qp_solver->qp_in->Cu;
    real_t **qp_lc = (real_t **)sqp_args->qp_solver->qp_in->lc;
    real_t **qp_uc = (real_t **)sqp_args->qp_solver->qp_in->uc;

    real_t **hess_l = (real_t **)sqp_mem->common->hess_l;
    real_t **grad_f = (real_t **)sqp_mem->common->grad_f;
    real_t **jac_h = (real_t **)sqp_mem->common->jac_h;
    real_t **jac_g = (real_t **)sqp_mem->common->jac_g;
    real_t **h = (real_t **)sqp_mem->common->h;
    real_t **g = (real_t **)sqp_mem->common->g;

    real_t **nlp_x = (real_t **)sqp_mem->common->x;
    real_t **nlp_u = (real_t **)sqp_mem->common->u;
    real_t **nlp_lg = (real_t **)nlp_in->lg;
    real_t **nlp_ug = (real_t **)nlp_in->ug;

    // Objective
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            for (int_t k = 0; k < nx[i]; k++) {
                qp_Q[i][j * nx[i] + k] = hess_l[i][j * (nx[i] + nu[i]) + k];
            }
            for (int_t k = 0; k < nu[i]; k++) {
                qp_S[i][j * nu[i] + k] = hess_l[i][j * (nx[i] + nu[i]) + nx[i] + k];
            }
        }
        for (int_t j = 0; j < nu[i]; j++) {
            for (int_t k = 0; k < nu[i]; k++) {
                qp_R[i][j * nu[i] + k] = hess_l[i][(nx[i] + j) * (nx[i] + nu[i]) + nx[i] + k];
            }
        }
        for (int_t j = 0; j < nx[i]; j++) {
            qp_q[i][j] = grad_f[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            qp_r[i][j] = grad_f[i][nx[i] + j];
        }
    }

    // State-continuity constraints, and state/control bounds
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            qp_b[i][j] = h[i][j] - nlp_x[i + 1][j];
            for (int_t k = 0; k < nx[i]; k++) {
                qp_A[i][k * nx[i] + j] = jac_h[i][k * nx[i] + j];
            }
            for (int_t k = 0; k < nu[i]; k++) {
                qp_B[i][k * nx[i] + j] = jac_h[i][(nx[i] + k) * nx[i] + j];
            }
        }
        for (int_t j = 0; j < nb[i]; j++) {
#ifdef FLIP_BOUNDS
            if (idxb[i][j] < nu[i]) {
                qp_lb[i][j] = nlp_in->lb[i][j] - nlp_u[i][idxb[i][j]];
                qp_ub[i][j] = nlp_in->ub[i][j] - nlp_u[i][idxb[i][j]];
            } else {
                qp_lb[i][j] = nlp_in->lb[i][j] - nlp_x[i][idxb[i][j] - nu[i]];
                qp_ub[i][j] = nlp_in->ub[i][j] - nlp_x[i][idxb[i][j] - nu[i]];
            }
#else
            if (idxb[i][j] < nx[i]) {
                qp_lb[i][j] = nlp_in->lb[i][j] - nlp_x[i][idxb[i][j]];
                qp_ub[i][j] = nlp_in->ub[i][j] - nlp_x[i][idxb[i][j]];
            } else {
                qp_lb[i][j] = nlp_in->lb[i][j] - nlp_u[i][idxb[i][j] - nx[i]];
                qp_ub[i][j] = nlp_in->ub[i][j] - nlp_u[i][idxb[i][j] - nx[i]];
            }
#endif
        }
    }

    // Path constraints
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < ng[i]; j++) {
            qp_lc[i][j] = nlp_lg[i][j] - g[i][j];
            qp_uc[i][j] = nlp_ug[i][j] - g[i][j];
            for (int_t k = 0; k < nx[i]; k++) qp_Cx[i][k * ng[i] + j] = jac_g[i][k * ng[i] + j];
            for (int_t k = 0; k < nu[i]; k++)
                qp_Cu[i][k * ng[i] + j] = jac_g[i][(nx[i] + k) * ng[i] + j];
        }
    }
}

void update_variables(const ocp_nlp_in *nlp_in, ocp_nlp_sqp_args *sqp_args,
                      ocp_nlp_sqp_memory *sqp_mem) {
    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
    const int_t *nb = nlp_in->nb;
    const int_t *ng = nlp_in->ng;

    real_t **nlp_x = (real_t **)sqp_mem->common->x;
    real_t **nlp_u = (real_t **)sqp_mem->common->u;
    real_t **nlp_pi = (real_t **)sqp_mem->common->pi;
    real_t **nlp_lam = (real_t **)sqp_mem->common->lam;

    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            nlp_x[i][j] += sqp_args->qp_solver->qp_out->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            nlp_u[i][j] += sqp_args->qp_solver->qp_out->u[i][j];
        }
        for (int_t j = 0; j < 2 * nb[i] + 2 * ng[i]; j++) {
            nlp_lam[i][j] = sqp_args->qp_solver->qp_out->lam[i][j];
        }
    }

    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i + 1]; j++) {
            nlp_pi[i][j] = sqp_args->qp_solver->qp_out->pi[i][j];
        }
    }
}

void store_variables(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_sqp_memory *sqp_mem) {
    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
    const int_t *nb = nlp_in->nb;
    const int_t *ng = nlp_in->ng;

    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            nlp_out->x[i][j] = sqp_mem->common->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            nlp_out->u[i][j] = sqp_mem->common->u[i][j];
        }
        for (int_t j = 0; j < 2 * nb[i] + 2 * ng[i]; j++) {
            nlp_out->lam[i][j] = sqp_mem->common->lam[i][j];
        }
    }

    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < nx[i + 1]; j++) {
            nlp_out->pi[i][j] = sqp_mem->common->pi[i][j];
        }
    }
}

ocp_nlp_sqp_args *ocp_nlp_sqp_create_arguments() {
    ocp_nlp_sqp_args *args = (ocp_nlp_sqp_args *)malloc(sizeof(ocp_nlp_sqp_args));
    args->maxIter = 10;

    return args;
}

int_t ocp_nlp_sqp_calculate_memory_size(const ocp_nlp_in *nlp_in, void *args_) {
    int_t size = sizeof(ocp_nlp_sqp_memory);

    size += ocp_nlp_calculate_memory_size(nlp_in);
    size += sizeof(ocp_nlp_sm_in);
    size += sizeof(ocp_nlp_sm_out);

    return size;
}

char *ocp_nlp_sqp_assign_memory(const ocp_nlp_in *nlp_in, void *args_, void **mem_,
                                void *raw_memory) {
    ocp_nlp_sqp_memory **sqp_memory = (ocp_nlp_sqp_memory **)mem_;
    char *c_ptr = (char *)raw_memory;

    *sqp_memory = (ocp_nlp_sqp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_memory);

    c_ptr = ocp_nlp_assign_memory(nlp_in, (void **)(&(*sqp_memory)->common), (void *)c_ptr);

    (*sqp_memory)->sm_in = (ocp_nlp_sm_in *)c_ptr;
    c_ptr += sizeof(ocp_nlp_sm_in);

    (*sqp_memory)->sm_out = (ocp_nlp_sm_out *)c_ptr;
    c_ptr += sizeof(ocp_nlp_sm_out);

    return c_ptr;
}

ocp_nlp_sqp_memory *ocp_nlp_sqp_create_memory(const ocp_nlp_in *nlp_in, void *args_) {
    ocp_nlp_sqp_memory *mem;

    int_t memory_size = ocp_nlp_sqp_calculate_memory_size(nlp_in, args_);
    void *raw_memory_ptr = malloc(memory_size);

    char *ptr_end = ocp_nlp_sqp_assign_memory(nlp_in, args_, (void **)&mem, raw_memory_ptr);
    assert((char *)raw_memory_ptr + memory_size >= ptr_end);
    (void)ptr_end;

    return mem;
}

int_t ocp_nlp_sqp_calculate_workspace_size(const ocp_nlp_in *nlp_in, void *args_) {
    int_t size = sizeof(ocp_nlp_sqp_workspace);

    return size;
}

char *ocp_nlp_sqp_assign_workspace(const ocp_nlp_in *nlp_in, void *args_, void **work_,
                                   void *raw_memory) {
    ocp_nlp_sqp_workspace **sqp_workspace = (ocp_nlp_sqp_workspace **)work_;
    char *c_ptr = (char *)raw_memory;

    *sqp_workspace = (ocp_nlp_sqp_workspace *)c_ptr;
    c_ptr += sizeof(ocp_nlp_sqp_workspace);

    return c_ptr;
}

ocp_nlp_sqp_workspace *ocp_nlp_sqp_create_workspace(const ocp_nlp_in *nlp_in, void *args_) {
    ocp_nlp_sqp_workspace *work;

    int_t workspace_size = ocp_nlp_sqp_calculate_workspace_size(nlp_in, args_);
    void *raw_memory_ptr = malloc(workspace_size);

    char *ptr_end = ocp_nlp_sqp_assign_workspace(nlp_in, args_, (void **)&work, raw_memory_ptr);
    assert((char *)raw_memory_ptr + workspace_size >= ptr_end);
    (void)ptr_end;

    return work;
}

int_t ocp_nlp_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, void *args_, void *memory_,
                  void *workspace_) {
    int_t return_status = 0;

    ocp_nlp_sqp_args *sqp_args = (ocp_nlp_sqp_args *)args_;
    ocp_nlp_sqp_memory *sqp_mem = (ocp_nlp_sqp_memory *)memory_;

    // SQP iterations
    int_t max_sqp_iterations = sqp_args->maxIter;

    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iterations; sqp_iter++) {
        // Compute/update quadratic approximation
        sqp_args->sensitivity_method->fun(
            sqp_mem->sm_in, sqp_mem->sm_out, sqp_args->sensitivity_method->args,
            sqp_args->sensitivity_method->mem, sqp_args->sensitivity_method->work);

        // Prepare QP
        prepare_qp(nlp_in, sqp_args, sqp_mem);

        // Solve QP
        int_t qp_status = sqp_args->qp_solver->fun(
            sqp_args->qp_solver->qp_in, sqp_args->qp_solver->qp_out, sqp_args->qp_solver->args,
            sqp_args->qp_solver->mem, sqp_args->qp_solver->work);
        if (qp_status) return_status = qp_status;

        // Update optimization variables (globalization)
        update_variables(nlp_in, sqp_args, sqp_mem);

        // TODO(nielsvd): debug, remove... Norm of step-size
        // real_t norm_step = 0;
        // for (int_t i = 0; i <= nlp_in->N; i++) {
        //     for (int_t j = 0; j < nlp_in->nx[i]; j++) {
        //         norm_step += (sqp_args->qp_solver->qp_out->x[i][j]) *
        //                      (sqp_args->qp_solver->qp_out->x[i][j]);
        //     }
        //     for (int_t j = 0; j < nlp_in->nu[i]; j++) {
        //         norm_step += (sqp_args->qp_solver->qp_out->u[i][j]) *
        //                      (sqp_args->qp_solver->qp_out->u[i][j]);
        //     }
        // }
        // norm_step = sqrt(norm_step);
        // printf("Norm_step = %.10e\n", norm_step);
    }

    // Post-process solution
    store_variables(nlp_in, nlp_out, sqp_mem);

    return return_status;
}

void ocp_nlp_sqp_initialize(const ocp_nlp_in *nlp_in, void *args_, void **mem_, void **work_) {
    ocp_nlp_sqp_args *args = (ocp_nlp_sqp_args *)args_;
    ocp_nlp_sqp_memory **mem = (ocp_nlp_sqp_memory **)mem_;
    ocp_nlp_sqp_workspace **work = (ocp_nlp_sqp_workspace **)work_;

    *mem = ocp_nlp_sqp_create_memory(nlp_in, args);
    *work = ocp_nlp_sqp_create_workspace(nlp_in, args);

    ocp_nlp_sm_in *sm_in = (*mem)->sm_in;
    ocp_nlp_sm_out *sm_out = (*mem)->sm_out;

    ocp_nlp_sm *nlp_sm = args->sensitivity_method;
    ocp_qp_solver *qp_solver = args->qp_solver;

    // Sensitivity method input
    sm_in->N = nlp_in->N;
    sm_in->nx = nlp_in->nx;
    sm_in->nu = nlp_in->nu;
    sm_in->nb = nlp_in->nb;
    sm_in->ng = nlp_in->ng;
    sm_in->cost = nlp_in->cost;
    sm_in->sim = (sim_solver **)nlp_in->sim;
    sm_in->path_constraints = (ocp_nlp_function **)nlp_in->path_constraints;
    sm_in->x = (const real_t **)(*mem)->common->x;
    sm_in->u = (const real_t **)(*mem)->common->u;
    sm_in->pi = (const real_t **)(*mem)->common->pi;
    sm_in->lam = (const real_t **)(*mem)->common->lam;

    // Sensitivity method output
    sm_out->hess_l = (const real_t **)(*mem)->common->hess_l;
    sm_out->grad_f = (const real_t **)(*mem)->common->grad_f;
    sm_out->jac_h = (const real_t **)(*mem)->common->jac_h;
    sm_out->jac_g = (const real_t **)(*mem)->common->jac_g;
    sm_out->h = (const real_t **)(*mem)->common->h;
    sm_out->g = (const real_t **)(*mem)->common->g;

    // // Initialize sensitivity method and QP solver
    nlp_sm->initialize(sm_in, nlp_sm->args, &nlp_sm->mem, &nlp_sm->work);
    qp_solver->initialize(qp_solver->qp_in, qp_solver->args, &qp_solver->mem, &qp_solver->work);
}

void ocp_nlp_sqp_destroy(void *mem_, void *work_) {
    ocp_nlp_sqp_memory *mem = (ocp_nlp_sqp_memory *)mem_;
    ocp_nlp_sqp_workspace *work = (ocp_nlp_sqp_workspace *)work_;

    free(work);
    free(mem);
}
