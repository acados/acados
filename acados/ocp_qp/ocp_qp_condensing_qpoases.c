/*
 *    This file is part of ACADOS.
 *
 *    ACADOS is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADOS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADOS; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <stdlib.h>
#include "hpmpc/include/aux_d.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/condensing.h"
#include "acados/utils/print.h"
#include "blasfeo/include/blasfeo_i_aux.h"

/* Ignore compiler warnings from qpOASES */
#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wtypedef-redefinition"
    #pragma clang diagnostic ignored "-Wtautological-pointer-compare"
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wunused-function"
    #include "qpOASES_e/Constants.h"
    #include "qpOASES_e/QProblemB.h"
    #pragma clang diagnostic pop
#elif defined(__GNUC__)
    #if __GNUC__ == 4
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/Constants.h"
        #include "qpOASES_e/QProblemB.h"
        #pragma GCC diagnostic error "-Wunused-parameter"
        #pragma GCC diagnostic error "-Wunused-function"
    #elif __GNUC__ == 6
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
        #pragma GCC diagnostic ignored "-Wunused-parameter"
        #pragma GCC diagnostic ignored "-Wunused-function"
        #include "qpOASES_e/Constants.h"
        #include "qpOASES_e/QProblemB.h"
        #pragma GCC diagnostic pop
    #endif
#endif

QProblemB QP;
real_t *A_row_major;
real_t *primal_solution;
real_t *dual_solution;
condensing_in in;
condensing_out out;
condensing_workspace work;

#ifdef DEBUG
static void print_condensed_QP(const int_t ncv, const int_t nc,
    condensing_out *out) {

    print_matrix("../experimental/robin/hessian.txt", out->H, ncv, ncv);
    print_array("../experimental/robin/gradient.txt", out->h, ncv);
    print_matrix("../experimental/robin/A.txt", out->A, nc, ncv);
    print_array("../experimental/robin/lbA.txt", out->lbA, nc);
    print_array("../experimental/robin/ubA.txt", out->ubA, nc);
    print_array("../experimental/robin/lb.txt", out->lb, ncv);
    print_array("../experimental/robin/ub.txt", out->ub, ncv);
}
#endif

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
    work.G = malloc(sizeof(*work.G) * N);
    work.g = malloc(sizeof(*work.g) * N);
    work.D = malloc(sizeof(*work.D) * (N+1));
    for (int_t i = 0; i < N; i++) {
        work.G[i] = malloc(sizeof(*(work.G[i])) * (i+1));
        work.D[i] = malloc(sizeof(*(work.D[i])) * (i+1));
        d_zeros(&work.g[i], qp_in->nx[i], 1);
        for (int_t j = 0; j <= i; j++) {
            d_zeros(&work.G[i][j], qp_in->nx[i], qp_in->nu[j]);
            d_zeros(&work.D[i][j], nc[i], qp_in->nu[j]);
        }
    }
    work.D[N] = malloc(sizeof(*(work.D[N])) * N);
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

static int_t solve_condensed_QP(const int_t ncv, QProblemB *QP,
    real_t* primal_solution, real_t* dual_solution) {

    int_t nwsr = 1000;
    real_t cpu_time = 100.0;

    QProblemBCON(QP, ncv, HST_POSDEF);
    QProblemB_setPrintLevel(QP, PL_NONE);
    QProblemB_printProperties(QP);

    int_t return_flag = QProblemB_initW(QP, out.H, out.h, out.lb,
                        out.ub, &nwsr, &cpu_time, NULL,
                        dual_solution, NULL, NULL);
    QProblemB_getPrimalSolution(QP, primal_solution);
    QProblemB_getDualSolution(QP, dual_solution);
    return return_flag;
}

static void recover_state_trajectory(ocp_qp_in *qp_in, real_t **x, real_t **u,
    real_t *primal_solution, const real_t *x0) {

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
    ocp_qp_condensing_qpoases_args *args, double *workspace) {

    fill_in_condensing_structs(qp_in);
    condensing_N2_fixed_initial_state(&in, &out, &work);

    // Process arguments
    args->dummy = 1.0;
    workspace++;
    workspace = 0;

    d_zeros(&A_row_major, work.nconstraints, work.nconvars);
    convert_to_row_major(out.A, A_row_major, work.nconstraints, work.nconvars);

    #ifdef DEBUG
    print_condensed_QP(work.nconvars, work.nconstraints, &out);
    #endif

    int_t return_flag = solve_condensed_QP(work.nconvars, &QP, primal_solution, dual_solution);
    recover_state_trajectory(qp_in, qp_out->x, qp_out->u, primal_solution, qp_in->lb[0]);

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

    d_free(A_row_major);
    // d_free(primal_solution);
    // d_free(dual_solution);

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

void initialise_qpoases(ocp_qp_in *in) {
    int_t ncv = get_num_condensed_vars(in);
    calculate_num_state_bounds(in);
    int_t nconstraints = get_num_constraints(in);
    d_zeros(&primal_solution, ncv, 1);
    d_zeros(&dual_solution, ncv+nconstraints, 1);
}
