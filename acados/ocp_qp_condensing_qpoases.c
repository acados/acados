#include <stdlib.h>
#include "hpmpc/include/aux_d.h"
#include "acados/ocp_qp_condensing_qpoases.h"
#include "acados/condensing.h"
#include "acados/print.h"

/* Ignore compiler warnings from qpOASES */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtypedef-redefinition"
#pragma clang diagnostic ignored "-Wtautological-pointer-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "qpOASES_e/Constants.h"
#include "qpOASES_e/QProblem.h"
#pragma clang diagnostic pop

QProblem QP;
real_t *A_row_major;
real_t *primal_solution;
real_t *dual_solution;
condensing_input in;
condensing_output out;
condensing_workspace work;

#ifdef DEBUG
static void print_condensed_QP(const int_t ncv, const int_t nc,
    condensing_output *out) {

    print_matrix("../experimental/robin/hessian.txt", out->H, ncv, ncv);
    print_array("../experimental/robin/gradient.txt", out->h, ncv);
    print_matrix("../experimental/robin/A.txt", out->A, nc, ncv);
    print_array("../experimental/robin/lbA.txt", out->lbA, nc);
    print_array("../experimental/robin/ubA.txt", out->ubA, nc);
    print_array("../experimental/robin/lb.txt", out->lb, ncv);
    print_array("../experimental/robin/ub.txt", out->ub, ncv);
}
#endif

static int_t get_num_condensed_vars(ocp_qp_input *in) {
    int_t num_condensed_vars = 0;
    // TODO(robin): this only holds for MPC, not MHE
    num_condensed_vars += 0*(in->nx[1]);
    for (int_t i = 0; i < in->N; i++)
        num_condensed_vars += in->nu[i];
    return num_condensed_vars;
}

static void calculate_num_state_bounds(ocp_qp_input *in) {
    i_zeros(&work.nstate_bounds, in->N+1, 1);
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

static int_t get_num_constraints(ocp_qp_input *in) {
    int_t num_constraints = in->nc[0];
    for (int_t i = 1; i <= in->N; i++) {
        num_constraints += in->nc[i] + work.nstate_bounds[i];
    }
    return num_constraints;
}

static void fill_in_condensing_structs(ocp_qp_input *qp_in) {
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
        d_zeros(&work.g[i], NX, 1);
        for (int_t j = 0; j <= i; j++) {
            d_zeros(&work.G[i][j], NX, NU);
            d_zeros(&work.D[i][j], nc[i], NU);
        }
    }
    work.D[N] = malloc(sizeof(*(work.D[N])) * N);
    for (int_t i = 0; i < N; i++) {
        d_zeros(&work.D[N][i], nc[N], NU);
    }
    d_zeros(&work.W1_x, NX, NX);
    d_zeros(&work.W2_x, NX, NX);
    d_zeros(&work.W1_u, NX, NU);
    d_zeros(&work.W2_u, NX, NU);
    d_zeros(&work.w1, NX, 1);
    d_zeros(&work.w2, NX, 1);
}

static int_t solve_condensed_QP(QProblem QP, real_t* primal_solution, real_t* dual_solution) {
    int_t nwsr = 1000;
    real_t cpu_time = 100.0;

    int_t return_flag = QProblem_initW(&QP, out.H, out.h, A_row_major, out.lb,
                        out.ub, out.lbA, out.ubA, &nwsr, &cpu_time, NULL,
                        dual_solution, NULL, NULL, NULL);
    QProblem_getPrimalSolution(&QP, primal_solution);
    QProblem_getDualSolution(&QP, dual_solution);
    return return_flag;
}

static void recover_state_trajectory(int_t N, real_t **x, real_t **u,
    real_t *primal_solution, const real_t *x0) {

    for (int_t i = 0; i < NX; i++) {
        x[0][i] = x0[i];
    }
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < NX; j++) {
            x[i+1][j] = work.g[i][j];
            for (int_t k = 0; k <= i; k++) {
                for (int_t l = 0; l < NU; l++) {
                    x[i+1][j] += work.G[i][k][j+NX*l]*primal_solution[NU*k+l];
                }
            }
        }
        for (int_t j = 0; j < NU; j++) u[i][j] = primal_solution[i*NU+j];
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

int_t ocp_qp_condensing_qpoases(ocp_qp_input *qp_in, ocp_qp_output *qp_out,
    ocp_qp_condensing_qpoases_args *args, double *workspace) {

    fill_in_condensing_structs(qp_in);
    condensingN2_fixed_initial_state(&in, &out, &work);

    // Process arguments
    args->dummy = 1.0;
    workspace = 0;

    d_zeros(&A_row_major, work.nconstraints, work.nconvars);
    convert_to_row_major(out.A, A_row_major, work.nconstraints, work.nconvars);

    #ifdef DEBUG
    print_condensed_QP(work.nconvars, work.nconstraints, &out);
    #endif

    int_t return_flag = solve_condensed_QP(QP, primal_solution, dual_solution);
    recover_state_trajectory(qp_in->N, qp_out->x, qp_out->u, primal_solution, qp_in->lb[0]);

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

int_t ocp_qp_condensing_qpoases_workspace_size(ocp_qp_input *in,
    ocp_qp_condensing_qpoases_args *args) {

    int_t ws_size = 0;
    args->dummy = 0.0;

    int_t N = in->N;
    int_t ncv = get_num_condensed_vars(in);
    int_t ncs = get_num_constraints(in);

    ws_size += 0*(ncv*ncv + ncv + ncv + ncv + ncs*ncv + ncs + ncs);  // condensing output
    ws_size += sizeof(*work.G)*N + sizeof(*(work.G[0]))*(N+1)*N/2
                + NX*NU*sizeof(*(work.G[0][0]))*(N+1)*N/2;

    return ws_size;
}

void initialise_qpoases(ocp_qp_input *in) {
    int_t ncv = get_num_condensed_vars(in);
    calculate_num_state_bounds(in);
    int_t nconstraints = get_num_constraints(in);
    d_zeros(&primal_solution, ncv, 1);
    d_zeros(&dual_solution, ncv+nconstraints, 1);
    QProblemCON(&QP, ncv, nconstraints, HST_POSDEF);
    QProblem_setPrintLevel(&QP, PL_NONE);
    QProblem_printProperties(&QP);
}
