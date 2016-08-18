#include <stdlib.h>
#include "ocp_qp_condensing_qpoases.h"
#include "condensing.h"
#include "hpmpc/include/aux_d.h"

/* qpOASES specifics */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtypedef-redefinition"
#pragma clang diagnostic ignored "-Wtautological-pointer-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblem.h"
#pragma clang diagnostic pop

QProblem QP;
real_t *A_row_major;
real_t *primal_solution;
real_t *dual_solution;
condensing_input in;
condensing_output out;
condensing_workspace work;

// static int_t get_num_opt_vars(int_t N, int_t *nx, int_t *nu) {
//     int_t num_opt_vars = 0;
//     for (int_t i = 0; i < N; i++)
//         num_opt_vars += nx[i] + nu[i];
//     num_opt_vars += nx[N];
//     return num_opt_vars;
// }

static int_t get_num_condensed_vars(ocp_qp_input *in) {
    int_t num_condensed_vars = 0;
    // TODO(robin): this only holds for MPC, not MHE
    num_condensed_vars += 0*(in->nx[1]);
    for (int_t i = 0; i < in->N; i++)
        num_condensed_vars += in->nu[i];
    return num_condensed_vars;
}

static int_t get_num_constraints(ocp_qp_input *in) {
    int_t num_constraints = 0;
    for (int_t i = 0; i < in->N; i++) {
        // TODO(robin): count actual simple bounds on states
        num_constraints += in->nc[i] + in->nx[i];
    }
    num_constraints += in->nc[in->N];
    return num_constraints;
}

// static void write_array_to_file(FILE *outputFile, real_t *array, int_t size) {
//     for (int_t i = 0; i < size; i++) fprintf(outputFile, "%g ", array[i]);
//     fprintf(outputFile, "\n");
// }
//
// static void write_condensed_QP_to_file(const int_t ncv, const int_t nc,
//     condensing_output *out) {
//
//     FILE *outFile = fopen("../experimental/robin/QP_data.txt", "w");
//     if (outFile == NULL) {
//         fprintf(stderr, "%s\n", "OPEN FILE FAILED!");
//     }
//     write_array_to_file(outFile, out->H, ncv*ncv);
//     write_array_to_file(outFile, out->h, ncv);
//     write_array_to_file(outFile, out->A, nc*ncv);
//     write_array_to_file(outFile, out->lb, ncv);
//     write_array_to_file(outFile, out->ub, ncv);
//     write_array_to_file(outFile, out->lbA, nc);
//     write_array_to_file(outFile, out->ubA, nc);
//     fclose(outFile);
// }

static void fill_in_condensing_structs(ocp_qp_input *qp_in) {
    // Input
    in.qp_input = qp_in;
    int_t N = qp_in->N;
    const int_t *nc = qp_in->nc;

    // Output
    int_t nconvars = get_num_condensed_vars(qp_in);
    int_t nconstraints = get_num_constraints(qp_in);
    d_zeros(&out.H, nconvars, nconvars);
    d_zeros(&out.h, nconvars, 1);
    d_zeros(&out.lb, nconvars, 1);
    d_zeros(&out.ub, nconvars, 1);
    d_zeros(&out.A, nconstraints, nconvars);
    d_zeros(&out.lbA, nconstraints, 1);
    d_zeros(&out.ubA, nconstraints, 1);

    // Workspace
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
    real_t cput = 100.0;

    int_t return_flag = QProblem_initW(&QP, out.H, out.h, A_row_major, out.lb,
                        out.ub, out.lbA, out.ubA,
                        &nwsr, &cput, NULL, dual_solution, NULL, NULL, NULL);
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
        // for (int_t j = 0; j < NU; j++) u[i][j] = primal_solution[NX+i*NU+j];
    }
}

int_t ocp_qp_condensing_qpoases(ocp_qp_input *qp_in, ocp_qp_output *qp_out,
    ocp_qp_condensing_qpoases_args *args, double *workspace) {

    fill_in_condensing_structs(qp_in);
    condensingN2_fixed_initial_state(&in, &out, &work);

    // Process arguments
    args->dummy = 1.0;
    workspace = 0;

    // Convert A to row major
    int_t num_condensed_vars = get_num_condensed_vars(qp_in);
    int_t num_constraints = get_num_constraints(qp_in);
    d_zeros(&A_row_major, num_constraints, num_condensed_vars);
    for (int_t i = 0; i < num_constraints; i++) {
        for (int_t j = 0; j < num_condensed_vars; j++) {
            A_row_major[i*num_condensed_vars+j] = out.A[j*num_constraints+i];
        }
    }
    // write_condensed_QP_to_file(num_condensed_vars, num_constraints, &out);
    d_zeros(&primal_solution, num_condensed_vars, 1);
    d_zeros(&dual_solution, num_condensed_vars+num_constraints, 1);
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
    d_free(primal_solution);
    d_free(dual_solution);

    return return_flag;
}

void initialise_qpoases(ocp_qp_input *in) {
    int_t ncv = get_num_condensed_vars(in);
    int_t nconstraints = get_num_constraints(in);
    QProblemCON(&QP, ncv, nconstraints, HST_POSDEF);
    QProblem_setPrintLevel(&QP, PL_NONE);
    QProblem_printProperties(&QP);
}
