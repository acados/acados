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

#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"

#include <stdlib.h>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
/* Ignore compiler warnings from qpOASES */
#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wtautological-pointer-compare"
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wunused-function"
    #include "qpOASES_e/QProblemB.h"
    #include "qpOASES_e/QProblem.h"
    #pragma clang diagnostic pop
#elif defined(__GNUC__)
    #if __GNUC__ >= 6
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
        #pragma GCC diagnostic pop
    #else
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/QProblemB.h"
        #include "qpOASES_e/QProblem.h"
    #endif
#else
    #include "qpOASES_e/QProblemB.h"
    #include "qpOASES_e/QProblem.h"
#endif

#include "acados/ocp_qp/condensing.h"
#include "acados/utils/print.h"

QProblemB QPB;
QProblem QP;
real_t *A_row_major;
real_t *primal_solution;
real_t *dual_solution;
condensing_in in;
condensing_out out;
condensing_workspace work;

// #ifdef DEBUG
// static void print_ocp_qp(ocp_qp_in *in) {
//     int_t N = in->N;
//     for (int_t i = 0; i < N; i++) {
//         print_int_array("nx.txt", in->nx, N+1);
//         print_int_array("nu.txt", in->nu, N);
//         print_int_array("nb.txt", in->nb, N+1);
//         print_int_array("nc.txt", in->nc, N+1);
//     }
// }

// static void print_condensed_QP(const int_t ncv, const int_t nc,
//     condensing_out *out) {

//     print_matrix("hessian.txt", out->H, ncv, ncv);
//     print_array("gradient.txt", out->h, ncv);
//     print_matrix("A.txt", out->A, nc, ncv);
//     print_array("lbA.txt", out->lbA, nc);
//     print_array("ubA.txt", out->ubA, nc);
//     print_array("lb.txt", out->lb, ncv);
//     print_array("ub.txt", out->ub, ncv);
// }
// #endif

static int_t get_num_condensed_vars(ocp_qp_in *in) {
    int_t num_condensed_vars = 0;
    // TODO(robin): this only holds for MPC, not MHE
    num_condensed_vars += 0*(in->nx[1]);
    for (int_t i = 0; i < in->N; i++)
        num_condensed_vars += in->nu[i];
    return num_condensed_vars;
}

static void calculate_num_state_bounds(ocp_qp_in *in) {
    int_zeros(&work.nstate_bounds, in->N+1, 1);
    int_t num_state_bounds;
    for (int_t i = 1; i <= in->N; i++) {
        num_state_bounds = 0;
        for (int_t j = 0; j < in->nb[i]; j++) {
            if (in->idxb[i][j] < in->nx[i])
                num_state_bounds++;
        }
        work.nstate_bounds[i] = num_state_bounds;
    }
}

static int_t get_num_constraints(ocp_qp_in *in) {
    int_t num_constraints = in->nc[0];
    for (int_t i = 1; i <= in->N; i++) {
        num_constraints += in->nc[i] + work.nstate_bounds[i];
    }
    return num_constraints;
}

static void fill_in_condensing_structs(ocp_qp_in *qp_in) {
    int_t N = qp_in->N;
    const int_t *nc = qp_in->nc;

    // condensing input
    in.qp_input = qp_in;

    // condensing output
    int_t nconvars = get_num_condensed_vars(qp_in);
    int_t nconstraints = get_num_constraints(qp_in);
    d_zeros(&out.H, nconvars, nconvars);
    d_zeros(&out.h, nconvars, 1);
    d_zeros(&out.lb, nconvars, 1);
    d_zeros(&out.ub, nconvars, 1);
    d_zeros(&out.A, nconstraints, nconvars);
    d_zeros(&out.lbA, nconstraints, 1);
    d_zeros(&out.ubA, nconstraints, 1);

    for (int_t i = 0; i < nconvars; i++) {
        out.lb[i] = -QPOASES_INFTY;
        out.ub[i] = +QPOASES_INFTY;
    }

    for (int_t i = 0; i < nconstraints; i++) {
        out.lbA[i] = -QPOASES_INFTY;
        out.ubA[i] = +QPOASES_INFTY;
    }

    // condensing workspace
    work.nconvars = nconvars;
    work.nconstraints = nconstraints;
    work.G = calloc(N, sizeof(*work.G));
    work.g = calloc(N, sizeof(*work.g));
    work.D = calloc(N+1, sizeof(*work.D));
    for (int_t i = 0; i < N; i++) {
        work.G[i] = calloc(i+1, sizeof(*(work.G[i])));
        work.D[i] = calloc(i+1, sizeof(*(work.D[i])));
        d_zeros(&work.g[i], qp_in->nx[i], 1);
        for (int_t j = 0; j <= i; j++) {
            d_zeros(&work.G[i][j], qp_in->nx[i], qp_in->nu[j]);
            d_zeros(&work.D[i][j], nc[i], qp_in->nu[j]);
        }
    }
    work.D[N] = calloc(N, sizeof(*(work.D[N])));
    for (int_t i = 0; i < N; i++) {
        d_zeros(&work.D[N][i], nc[N], qp_in->nu[i]);
    }
    d_zeros(&work.W1_x, qp_in->nx[0], qp_in->nx[0]);
    d_zeros(&work.W2_x, qp_in->nx[0], qp_in->nx[0]);
    d_zeros(&work.W1_u, qp_in->nx[0], qp_in->nu[0]);
    d_zeros(&work.W2_u, qp_in->nx[0], qp_in->nu[0]);
    d_zeros(&work.w1, qp_in->nx[0], 1);
    d_zeros(&work.w2, qp_in->nx[0], 1);
    d_zeros(&work.Sx0, qp_in->nu[0], 1);
}

static int_t solve_condensed_QPB(const int_t ncv, QProblemB *QP,
    real_t* primal_solution_, real_t* dual_solution_) {

    int_t nwsr = 1000;
    real_t cpu_time = 100.0;

    QProblemBCON(QP, ncv, HST_POSDEF);
    QProblemB_setPrintLevel(QP, PL_MEDIUM);
    QProblemB_printProperties(QP);

    int_t return_flag = QProblemB_initW(QP, out.H, out.h, out.lb,
            out.ub, &nwsr, &cpu_time, NULL,
            dual_solution_, NULL, NULL);
    QProblemB_getPrimalSolution(QP, primal_solution_);
    QProblemB_getDualSolution(QP, dual_solution_);
    return return_flag;
}

static int_t solve_condensed_QP(const int_t ncv, const int_t ncon, QProblem *QP,
    real_t* primal_solution_, real_t* dual_solution_) {

    int_t nwsr = 1000;
    real_t cpu_time = 100.0;

    QProblemCON(QP, ncv, ncon, HST_POSDEF);
    QProblem_setPrintLevel(QP, PL_MEDIUM);
    QProblem_printProperties(QP);

    int_t return_flag = QProblem_initW(QP, out.H, out.h, A_row_major, out.lb,
            out.ub, out.lbA, out.ubA, &nwsr, &cpu_time, NULL,
            dual_solution_, NULL, NULL, NULL);
    QProblem_getPrimalSolution(QP, primal_solution_);
    QProblem_getDualSolution(QP, dual_solution_);
    return return_flag;
}

static void recover_state_trajectory(ocp_qp_in *qp_in, real_t **x, real_t **u, real_t **pi,
    real_t *primal_solution, real_t *dual_solution, const real_t *x0,
    int_t nconstraints, int_t nconvars) {

    for (int_t i = 0; i < qp_in->nx[0]; i++) {
        x[0][i] = x0[i];
    }
    for (int_t i = 0; i < qp_in->N; i++) {
        for (int_t j = 0; j < qp_in->nx[0]; j++) {
            x[i+1][j] = work.g[i][j];
            for (int_t k = 0; k <= i; k++) {
                for (int_t l = 0; l < qp_in->nu[0]; l++) {
                    x[i+1][j] += work.G[i][k][j+qp_in->nx[0]*l]*primal_solution[qp_in->nu[0]*k+l];
                }
            }
        }
        for (int_t j = 0; j < qp_in->nu[0]; j++) u[i][j] = primal_solution[i*qp_in->nu[0]+j];
        // TODO(robin): this only holds for MPC, not MHE
        // for (int_t j = 0; j < NU; j++) u[i][j] = primal_solution[NX+i*NU+j];
    }

    // Recover multipliers in pi
    // pi_k = lambda_k + Q_k x_k
    int_t cindex = nconvars+nconstraints;
    for (int_t i = qp_in->N; i > 0; i--) {
        int_t idxb = 0;
        cindex -= (work.nstate_bounds[i]+qp_in->nc[i]);
        // pi_k = lambda_k
        for (int_t j = 0; j < qp_in->nx[0]; j++) {
            if (qp_in->nb[i] > idxb && qp_in->idxb[i][idxb] == j) {  // bound on x[j]
                pi[(i-1)][j] = -dual_solution[cindex+idxb];
                idxb++;
            } else {
                pi[(i-1)][j] = 0.0;
            }
        }
        if (qp_in->nc[i] > 0) {
            // pi_k = pi_k + Cx_k+1^T lambda_k+1
            for (int_t j = 0; j < qp_in->nx[0]; j++) {
                for (int_t k = 0; k < qp_in->nc[i]; k++) {
                    pi[(i-1)][j] -= qp_in->Cx[i][j*qp_in->nc[i]+k]
                                                 *dual_solution[cindex+work.nstate_bounds[i]+k];
                }
            }
        }

        // pi_k = pi_k + Q_k x_k + q_k
        for (int_t j = 0; j < qp_in->nx[0]; j++) {
            pi[(i-1)][j] += qp_in->q[i][j];
            for (int_t k = 0; k < qp_in->nx[0]; k++) {
                pi[(i-1)][j] += qp_in->Q[i][k*qp_in->nx[0]+j]*x[i][k];
            }
        }

        if (i < qp_in->N) {
            // pi_k = pi_k + S_k u_k
            for (int_t j = 0; j < qp_in->nx[0]; j++) {
                for (int_t k = 0; k < qp_in->nu[0]; k++) {
                    pi[(i-1)][j] += qp_in->S[i][j*qp_in->nu[0]+k]*u[i][k];
                }
            }

            // pi_k = pi_k + A_k^T pi_k+1
            for (int_t j = 0; j < qp_in->nx[0]; j++) {
                for (int_t k = 0; k < qp_in->nx[0]; k++) {
                    pi[(i-1)][j] += qp_in->A[i][j*qp_in->nx[0]+k]*pi[i][k];
                }
            }
        }
    }
}

static void convert_to_row_major(const real_t *input, real_t *output, const int_t nrows,
    const int_t ncols) {

    for (int_t i = 0; i < nrows; i++) {
        for (int_t j = 0; j < ncols; j++) {
            output[i*ncols+j] = input[j*nrows+i];
        }
    }
}

int_t ocp_qp_condensing_qpoases(ocp_qp_in *qp_in, ocp_qp_out *qp_out,
    void *args_, void *mem_, void *workspace_) {

    int_t ncv = get_num_condensed_vars(qp_in);
    calculate_num_state_bounds(qp_in);
    int_t nconstraints = get_num_constraints(qp_in);
    d_zeros(&primal_solution, ncv, 1);
    d_zeros(&dual_solution, ncv+nconstraints, 1);

    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args*) args_;
    double *workspace = (double*) workspace_;
    fill_in_condensing_structs(qp_in);
    // #ifdef DEBUG
    // print_ocp_qp(qp_in);
    // #endif
    condensing_N2_fixed_initial_state(&in, &out, &work);

    // Process arguments
    args->dummy = 1.0;
    workspace++;
    workspace = 0;
    mem_ = 0;
    (void) mem_;

    d_zeros(&A_row_major, work.nconstraints, work.nconvars);
    convert_to_row_major(out.A, A_row_major, work.nconstraints, work.nconvars);

    // #ifdef DEBUG
    // print_condensed_QP(work.nconvars, work.nconstraints, &out);
    // #endif

    int_t return_flag;
    if (work.nconstraints) {
        return_flag = solve_condensed_QP(work.nconvars, work.nconstraints, &QP,
                primal_solution, dual_solution);
    } else {
        return_flag = solve_condensed_QPB(work.nconvars, &QPB, primal_solution, dual_solution);
    }
    recover_state_trajectory(qp_in, qp_out->x, qp_out->u, qp_out->pi,
            primal_solution, dual_solution, qp_in->lb[0], work.nconstraints, work.nconvars);

    d_free(out.H);
    d_free(out.h);
    d_free(out.lb);
    d_free(out.ub);
    d_free(out.A);
    d_free(out.lbA);
    d_free(out.ubA);

    for (int_t i = 0; i < qp_in->N; i++) {
        for (int_t j = 0; j <= i; j++) {
            d_free(work.G[i][j]);
            d_free(work.D[i][j]);
        }
        free(work.G[i]);
        free(work.D[i]);
        d_free(work.g[i]);
    }
    for (int_t i = 0; i < qp_in->N; i++) {
        d_free(work.D[qp_in->N][i]);
    }
    free(work.D[qp_in->N]);
    free(work.G);
    free(work.g);
    free(work.D);
    d_free(work.W1_x);
    d_free(work.W2_x);
    d_free(work.W1_u);
    d_free(work.W2_u);
    d_free(work.w1);
    d_free(work.w2);
    d_free(work.Sx0);
    int_free(work.nstate_bounds);
    d_free(A_row_major);
    d_free(primal_solution);
    d_free(dual_solution);

    return return_flag;
}

int_t ocp_qp_condensing_qpoases_workspace_size(ocp_qp_in *in,
    ocp_qp_condensing_qpoases_args *args) {

    int_t ws_size = 0;
    args->dummy = 0.0;

    int_t N = in->N;
    int_t ncv = get_num_condensed_vars(in);
    int_t ncs = get_num_constraints(in);

    ws_size += 0*(ncv*ncv + ncv + ncv + ncv + ncs*ncv + ncs + ncs);  // condensing output
    ws_size += sizeof(*work.G)*N + sizeof(*(work.G[0]))*(N+1)*N/2
                + in->nx[0]*in->nu[0]*sizeof(*(work.G[0][0]))*(N+1)*N/2;

    return ws_size;
}

void ocp_qp_condensing_qpoases_initialize(ocp_qp_in *qp_in, void *args_, void *mem_, void **work) {
    ocp_qp_condensing_qpoases_args *args = (ocp_qp_condensing_qpoases_args*) args_;

    // TODO(dimitris): replace dummy commands once interface completed
    args->dummy = 42.0;
    if (qp_in->nx[0] > 0)
        (void) mem_;
    (void) work;
}

void ocp_qp_condensing_qpoases_destroy(void *mem_, void *work) {
    // TODO(dimitris): replace dummy commands once interface completed
    (void) mem_;
    (void) work;
}
