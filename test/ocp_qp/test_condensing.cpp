#include <string>

#include "catch/include/catch.hpp"
#include "test/test_utils/read_matrix.hpp"
#include "acados/ocp_qp/condensing.c"
#include "test_condensing_helper.cpp"

#define COMPARISON_TOLERANCE 1.0e-15

using Eigen::MatrixXd;
using Eigen::VectorXd;

static MatrixXd readMatrixFromFile(std::string filename, int_t rows, int_t cols) {
    MatrixXd M = readMatrix(filename);
    REQUIRE(M.rows() == rows);
    REQUIRE(M.cols() == cols);
    return M;
}

static VectorXd readVectorFromFile(std::string filename, int_t length) {
    VectorXd v = readMatrix(filename);
    REQUIRE(v.rows() == length);
    return v;
}

static void readInputDimensionsFromFile(int_t *N, int_t *nx, int_t *nu) {
    *N = (int_t) readMatrix("N.dat")(0, 0);
    REQUIRE(*N > 0);
    *nx = (int_t) readMatrix("nx.dat")(0, 0);
    REQUIRE(*nx > 0);
    *nu = (int_t) readMatrix("nu.dat")(0, 0);
    REQUIRE(*nu > 0);
}

static void readInputDataFromFile(int_t nx, int_t nu, MatrixXd *A, MatrixXd *B,
    VectorXd *b, VectorXd *x0, MatrixXd *Q, MatrixXd *S, MatrixXd *R, VectorXd *q, VectorXd *r) {
    *A = readMatrixFromFile("A.dat", nx, nx);
    *B = readMatrixFromFile("B.dat", nx, nu);
    *b = readVectorFromFile("bv.dat", nx);
    *x0 = readVectorFromFile("x0.dat", nx);
    *Q = readMatrixFromFile("Q.dat", nx, nx);
    *S = readMatrixFromFile("S.dat", nu, nx);
    *R = readMatrixFromFile("R.dat", nu, nu);
    *q = readVectorFromFile("qv.dat", nx);
    *r = readVectorFromFile("rv.dat", nu);
}

TEST_CASE("Unconstrained LTV system", "[condensing]") {
    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;

    int_t N, nx, nu;
    readInputDimensionsFromFile(&N, &nx, &nu);
    MatrixXd A, B, Q, S, R;
    VectorXd b, x0, q, r;
    readInputDataFromFile(nx, nu, &A, &B, &b, &x0, &Q, &S, &R, &q, &r);

    int_t nx_vector[N+1], nu_vector[N], nb_vector[N+1], nc_vector[N+1];
    real_t *A_vector[N], *B_vector[N], *b_vector[N], *lb_vector[1], *ub_vector[1];
    real_t *Q_vector[N+1], *S_vector[N], *R_vector[N], *q_vector[N+1], *r_vector[N];

    for (int_t i = 0; i < N; i++) {
        nx_vector[i] = nx;
        nu_vector[i] = nu;
        nb_vector[i] = 0;
        nc_vector[i] = 0;
        A_vector[i] = A.data();
        B_vector[i] = B.data();
        b_vector[i] = b.data();
        Q_vector[i] = Q.data();
        S_vector[i] = S.data();
        R_vector[i] = R.data();
        q_vector[i] = q.data();
        r_vector[i] = r.data();
    }
    // Initial state
    nb_vector[0] = nx;
    lb_vector[0] = x0.data();
    ub_vector[0] = x0.data();
    // Final state
    nx_vector[N] = nx;
    nb_vector[N] = 0;
    nc_vector[N] = 0;
    Q_vector[N] = Q.data();
    q_vector[N] = q.data();

    // Finalize the data
    qp.N = N;
    qp.nx = (const int_t *) nx_vector;
    qp.nu = (const int_t *) nu_vector;
    qp.nb = (const int_t *) nb_vector;
    qp.nc = (const int_t *) nc_vector;
    qp.A = (const real_t **) A_vector;
    qp.B = (const real_t **) B_vector;
    qp.b = (const real_t **) b_vector;
    qp.Q = (const real_t **) Q_vector;
    qp.S = (const real_t **) S_vector;
    qp.R = (const real_t **) R_vector;
    qp.q = (const real_t **) q_vector;
    qp.r = (const real_t **) r_vector;
    qp.lb = (const real_t **) lb_vector;
    qp.ub = (const real_t **) ub_vector;

    fill_in_condensing_structs(&qp, &input, &output, &work);

    SECTION("Transition vector") {
        calculate_transition_vector(&qp, &work, x0.data());
        VectorXd true_g = readMatrix("transition_vector.dat");
        VectorXd acados_g = VectorXd(true_g).setZero();
        for (int_t i = 0; i < N; i++) {
            acados_g.block(i*nx, 0, nx, 1) = Eigen::Map<VectorXd>(work.g[i], nx);
        }
        REQUIRE(acados_g.isApprox(true_g, COMPARISON_TOLERANCE));
    }

    SECTION("Transition matrix") {
        calculate_transition_matrix(&qp, &work);
        MatrixXd true_G = readMatrixFromFile("transition_matrix.dat", N*nx, N*nu);
        MatrixXd acados_G = MatrixXd(true_G).setZero();
        for (int_t i = 0; i < N; i++) {
            for (int_t j = 0; j <= i; j++) {
                acados_G.block(i*nx, j*nu, nx, nu) = Eigen::Map<MatrixXd>(work.G[i][j], nx, nu);
            }
        }
        REQUIRE(acados_G.isApprox(true_G, COMPARISON_TOLERANCE));
    }

    SECTION("Condensed gradient") {
        calculate_transition_vector(&qp, &work, x0.data());
        calculate_transition_matrix(&qp, &work);
        calculate_gradient(&qp, &output, &work, 0, x0.data());
        VectorXd true_h = readVectorFromFile("condensed_gradient.dat", N*nu);
        VectorXd acados_h = VectorXd(true_h).setZero();
        for (int_t i = 0; i < N; i++) {
            acados_h.block(i*nu, 0, nu, 1) = Eigen::Map<VectorXd>(&output.h[i*nu], nu);
        }
        REQUIRE(acados_h.isApprox(true_h, COMPARISON_TOLERANCE));
    }
}
