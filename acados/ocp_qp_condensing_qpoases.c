#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "ocp_qp_condensing_qpoases.h"
#include "condensing.h"

/* qpOASES specifics */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtypedef-redefinition"
#pragma clang diagnostic ignored "-Wtautological-pointer-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-function"
#include "qpOASES_e/QProblem.h"
#pragma clang diagnostic pop
QProblem    QP;
real_t      _A[(NNN*(NX+NA)+NA)*NVC] = {0};
real_t      cput;
int_t       nwsr;
real_t      primal_solution[NVC]                     = {0};  // QP primal solution vector
real_t      dual_solution[(NNN+1)*NX+NNN*(NX+NU)+NX]    = {0};  // QP dual solution vector
/* condensing specifics */
data_struct data = {{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, \
                    {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};

int_t get_num_opt_vars(int_t NN, int_t *nx, int_t *nu) {
    int_t num_opt_vars = 0;
    for (int_t i = 0; i < NN; i++)
        num_opt_vars += nx[i] + nu[i];
    num_opt_vars += nx[NN];
    return num_opt_vars;
}

int_t get_num_condensed_vars(int_t NN, int_t *nx, int_t *nu) {
    int_t num_condensed_vars = 0;
    #if FIXED_INITIAL_STATE == 0
    num_condensed_vars += nx[1];
    #endif
    for (int_t i = 0; i < NN; i++)
        num_condensed_vars += nu[i];
    return num_condensed_vars;
}

void write_array_to_file(FILE *outputFile, real_t *array, int_t size) {
    for (int_t i = 0; i < size; i++) fprintf(outputFile, "%g ", array[i]);
    fprintf(outputFile, "\n");
}

void write_QP_data_to_file() {
    FILE *outFile = fopen("../experimental/robin/QP_data.txt", "w");
    if (outFile == NULL) {
        fprintf(stderr, "%s\n", "OPEN FILE FAILED!");
    }
    write_array_to_file(outFile, data.Hc, NVC*NVC);
    write_array_to_file(outFile, data.gc, NVC);
    write_array_to_file(outFile, data.Ac, (NNN*(NX+NA)+NA)*NVC);
    write_array_to_file(outFile, data.lbU, NVC);
    write_array_to_file(outFile, data.ubU, NVC);
    write_array_to_file(outFile, data.lbA, NNN*(NX+NA)+NA);
    write_array_to_file(outFile, data.ubA, NNN*(NX+NA)+NA);
    write_array_to_file(outFile, data.G, NNN*(NX)*NVC);
    write_array_to_file(outFile, data.g, NNN*NX);
    write_array_to_file(outFile, data.D, (NNN+1)*NA*NVC);
}

static void fill_in_objective_variables(int_t NN, int_t* nx, int_t* nu,
    real_t** Q, real_t** S, real_t** R, real_t** q, real_t** r) {
    for (int_t i = 0; i < NN+1; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            for (int_t k = 0; k < nx[i]; k++) data.Q[i*NX*NX+j*NX+k] = Q[i][j*NX+k];
        }
    }
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            for (int_t k = 0; k < nu[i]; k++) data.S[i*NU*NX+j*NU+k] = S[i][j*NU+k];
        }
    }
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < nu[i]; j++) {
            for (int_t k = 0; k < nu[i]; k++) data.R[i*NU*NU+j*NU+k] = R[i][j*NU+k];
        }
    }
    int_t start_of_current_block = 0;
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < nx[i]; j++)
            data.f[start_of_current_block+j] = q[i][j];
        for (int_t j = 0; j < nu[i]; j++)
            data.f[start_of_current_block+NX+j] = r[i][j];
        start_of_current_block += (NX + NU);
    }
    for (int_t j = 0; j < nx[NNN]; j++) data.f[NNN*(NX+NU)+j] = q[NNN][j];
}

static void fill_in_dynamics(int_t NN, int_t* nx, int_t* nu,
    real_t** A, real_t** B, real_t** b) {
    int_t start_of_current_block = 0;
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            for (int_t k = 0; k < nx[i+1]; k++) {
                data.A[start_of_current_block+j*nx[i+1]+k] = A[i][j*nx[i+1]+k];
            }
        }
        start_of_current_block += nx[i+1]*nx[i+1];
    }

    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < nx[i+1]; j++) data.b[i*NX+j] = b[i][j];
    }
    start_of_current_block = 0;
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < nu[i]; j++) {
            for (int_t k = 0; k < nx[i+1]; k++) {
                data.B[start_of_current_block+j*nx[i+1]+k] = B[i][j*nx[i+1]+k];
            }
        }
        start_of_current_block += nu[i]*nx[i+1];
    }
}

static void fill_in_bounds(int_t NN, int_t* nx, int_t* nu, int_t* nb,
    int_t** idxb, real_t** lb, real_t** ub) {
    for (int_t j = 0; j < nu[0]; j++) {
        data.lb[NX+idxb[0][j]] = lb[0][idxb[0][j]];
        data.ub[NX+idxb[0][j]] = ub[0][idxb[0][j]];
    }
    for (int_t i = 1; i < NN; i++) {
        for (int_t j = 0; j < nu[i]; j++) {
            // the interface assumes [u,x] ordering
            data.lb[NX+i*(nu[i]+nx[i])+idxb[i][j]] = lb[i][idxb[i][j]];
            data.ub[NX+i*(nu[i]+nx[i])+idxb[i][j]] = ub[i][idxb[i][j]];
        }
        for (int_t j = nu[i]; j < nb[i]; j++) {
            data.lb[-NU+i*(nu[i]+nx[i])+idxb[i][j]] = lb[i][idxb[i][j]];
            data.ub[-NU+i*(nu[i]+nx[i])+idxb[i][j]] = ub[i][idxb[i][j]];
        }
    }
    for (int_t j = 0; j < nb[NN]; j++) {
        data.lb[NN*(NX+NU)+idxb[NN][j]] = lb[NN][idxb[NN][j]];
        data.ub[NN*(NX+NU)+idxb[NN][j]] = ub[NN][idxb[NN][j]];
    }
}

static void fill_in_polytopic_constraints() {
    // TODO(robin):
}

static void fill_data_for_condensing(int_t NN, int_t *nx, int_t *nu, int_t *nb, int_t *ng, \
                                double **A, double **B, double **b, \
                                double **Q, double **S, double **R, double **q, double **r, \
                                int_t **idxb, double **lb, double **ub) {
    // Condensing implicitly assumes zeros initialisation
    memset(&data, 0, sizeof(data_struct));
    fill_in_objective_variables(NN, nx, nu, Q, S, R, q, r);
    fill_in_dynamics(NN, nx, nu, A, B, b);
    fill_in_bounds(NN, nx, nu, nb, idxb, lb, ub);
    fill_in_polytopic_constraints();
    for (int_t k = 0; k < NX; k++) {
        data.Dx[NNN*NA*NX+k*(NA+1)] = 1.0;
    }
}

static int_t solve_QP(QProblem QP, real_t* primal_solution, real_t* dual_solution) {
    nwsr = 1000;
    cput = 100.0;

    int_t return_flag = QProblem_initW(&QP, &(data.Hc[0]), &(data.gc[0]), &(_A[0]), &(data.lbU[0]),
                        &(data.ubU[0]), &(data.lbA[0]), &(data.ubA[0]),
                        &nwsr, &cput, NULL, dual_solution, NULL, NULL, NULL);
    QProblem_getPrimalSolution(&QP, primal_solution);
    QProblem_getDualSolution(&QP, dual_solution);
    return return_flag;
}

static void recover_state_trajectory(int_t NN,
    real_t** x, real_t** u, real_t* primal_solution) {
    for (int_t i = 0; i < NN; i++) {
        for (int_t j = 0; j < NX; j++) {
            x[i+1][j] = 0.0;
            for (int_t k = 0; k < NVC; k++) {
                x[i+1][j] = x[i+1][j] + data.G[i*NX+j+k*NN*NX]*primal_solution[k];
            }
            x[i+1][j] = x[i+1][j] + data.g[i*NX+j];
        }
        #if FIXED_INITIAL_STATE == 1
        for (int_t j = 0; j < NU; j++) u[i][j] = primal_solution[i*NU+j];
        #else
        for (int_t j = 0; j < NU; j++) u[i][j] = primal_solution[NX+i*NU+j];
        #endif
    }
}

int_t ocp_qp_condensing_qpoases(int_t NN, int_t *nx, int_t *nu, int_t *nb, int_t *ng,
    double **A, double **B, double **b,
    double **Q, double **S, double **R, double **q, double **r,
    int_t **idxb, double **lb, double **ub,
    double **C, double **D, double **ld, double **ud,
    double **x, double **u,
    struct ocp_qp_condensing_qpoases_args *args, double *work) {

    fill_data_for_condensing(NN, nx, nu, nb, ng, A, B, b,
    Q, S, R, q, r, idxb, lb, ub);
    condensingN2_fixed_initial_state();

    // Symmetrize H
    int_t num_condensed_vars = get_num_condensed_vars(NN, nx, nu);
    for (int_t i = 1; i < num_condensed_vars; i++) {
        for (int_t j = 0; j < i; j++) {
            data.Hc[i*num_condensed_vars+j] = data.Hc[j*num_condensed_vars+i];
        }
    }
    // Convert C to row major in A
    for (int_t i = 0; i < NN*(NX+NA)+NA; i++) {
        for (int_t j = 0; j < num_condensed_vars; j++) {
            _A[i*num_condensed_vars+j] = data.Ac[j*(NN*(NX+NA)+NA)+i];
        }
    }
    write_QP_data_to_file();
    int_t return_flag = solve_QP(QP, &(primal_solution[0]), &(dual_solution[0]));
    recover_state_trajectory(NN, x, u, &(primal_solution[0]));

    return return_flag;
}

void initialise_qpoases() {
    QProblemCON(&QP, NVC, NCONSTRAINTS, HST_POSDEF);
    QProblem_setPrintLevel(&QP, PL_NONE);
    QProblem_printProperties(&QP);
}
#pragma clang diagnostic pop
