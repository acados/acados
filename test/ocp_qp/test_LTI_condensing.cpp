#include <string>

#include "catch/include/catch.hpp"
#include "test/test_utils/read_matrix.hpp"
#include "acados/ocp_qp/condensing.c"
#include "test_condensing_helper.cpp"

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

TEST_CASE("Unconstrained LTI system", "[condensing]") {
    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;
    VectorXd x0;

    fillWithUnconstrainedData(&qp, &input, &output, &work, &x0);
    int_t N = qp.N;
    int_t nx = qp.nx[0];
    int_t nu = qp.nu[0];
    real_t *x0d = x0.data();

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
        VectorXd acados_h = Eigen::Map<VectorXd>(&output.h[0], N*nu);
        REQUIRE(acados_h.isApprox(true_h, COMPARISON_TOLERANCE));
    }

    SECTION("Condensed Hessian") {
        calculate_transition_vector(&qp, &work, x0.data());
        calculate_transition_matrix(&qp, &work);
        calculate_hessian(&qp, &output, &work, 0);
        MatrixXd true_H = readMatrixFromFile("condensed_hessian.dat", N*nu, N*nu);
        MatrixXd acados_H = Eigen::Map<MatrixXd>(&output.H[0], N*nu, N*nu);
        REQUIRE(acados_H.isApprox(true_H, COMPARISON_TOLERANCE));
    }
}

TEST_CASE("Constrained LTI system, simple bounds", "[condensing]") {
    REQUIRE(true);
}
