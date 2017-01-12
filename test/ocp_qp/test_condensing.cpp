#include <string>
#include <vector>
#include "catch/include/catch.hpp"
#include "test/test_utils/read_matrix.h"
#include "test/ocp_qp/condensing_test_helper.h"
#include "acados/ocp_qp/condensing_helper_functions.c"

using Eigen::MatrixXd;
using Eigen::VectorXd;

extern real_t COMPARISON_TOLERANCE;
// For each scenario, there must be a corresponding subfolder with test data
std::vector<std::string> scenarios = {"LTI", "LTV"};

void readInputDimensionsFromFile(int_t *N, int_t *nx, int_t *nu, std::string folder) {
    *N = (int_t) readMatrix(folder + "/N.dat")(0, 0);
    REQUIRE(*N > 0);
    *nx = (int_t) readMatrix(folder + "/nx.dat")(0, 0);
    REQUIRE(*nx > 0);
    *nu = (int_t) readMatrix(folder + "/nu.dat")(0, 0);
    REQUIRE(*nu > 0);
}

TEST_CASE("Unconstrained system", "[condensing]") {
    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;
    VectorXd x0;

    for (std::string scenario : scenarios) {
        fillWithUnconstrainedData(&qp, &x0, scenario);
        int_t N = qp.N;
        int_t nx = qp.nx[0];
        int_t nu = qp.nu[0];

        allocateCondensingData(&qp, &input, &output, &work);

        SECTION(scenario + " " + "Transition vector") {
            calculate_transition_vector(&qp, &work, x0.data());
            VectorXd true_g = readMatrix(scenario + "/transition_vector.dat");
            VectorXd acados_g = VectorXd(true_g).setZero();
            for (int_t i = 0; i < N; i++) {
                acados_g.block(i*nx, 0, nx, 1) = Eigen::Map<VectorXd>(work.g[i], nx);
            }
            REQUIRE(acados_g.isApprox(true_g, COMPARISON_TOLERANCE));
        }

        SECTION(scenario + " " + "Transition matrix") {
            calculate_transition_matrix(&qp, &work);
            MatrixXd true_G = readMatrixFromFile(scenario + "/transition_matrix.dat", N*nx, N*nu);
            MatrixXd acados_G = MatrixXd(true_G).setZero();
            for (int_t i = 0; i < N; i++) {
                for (int_t j = 0; j <= i; j++) {
                    acados_G.block(i*nx, j*nu, nx, nu) = Eigen::Map<MatrixXd>(work.G[i][j], nx, nu);
                }
            }
            REQUIRE(acados_G.isApprox(true_G, COMPARISON_TOLERANCE));
        }

        SECTION(scenario + " " + "Condensed gradient") {
            calculate_transition_vector(&qp, &work, x0.data());
            calculate_transition_matrix(&qp, &work);
            calculate_gradient(&qp, &output, &work, 0, x0.data());
            VectorXd true_h = readVectorFromFile(scenario + "/condensed_gradient.dat", N*nu);
            VectorXd acados_h = Eigen::Map<VectorXd>(&output.h[0], N*nu);
            REQUIRE(acados_h.isApprox(true_h, COMPARISON_TOLERANCE));
        }

        SECTION(scenario + " " + "Condensed Hessian") {
            calculate_transition_vector(&qp, &work, x0.data());
            calculate_transition_matrix(&qp, &work);
            calculate_hessian(&qp, &output, &work, 0);
            MatrixXd true_H = readMatrixFromFile(scenario + "/condensed_hessian.dat", N*nu, N*nu);
            MatrixXd acados_H = Eigen::Map<MatrixXd>(&output.H[0], N*nu, N*nu);
            REQUIRE(acados_H.isApprox(true_H, COMPARISON_TOLERANCE));
        }
    }
}

TEST_CASE("Constrained LTV system, simple bounds", "[condensing]") {
    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;
    VectorXd x0;

    for (std::string scenario : scenarios) {
        fillWithUnconstrainedData(&qp, &x0, scenario);
        int_t N = qp.N;
        int_t nx = qp.nx[0];
        int_t nu = qp.nu[0];
        fillWithBoundsData(&qp, N, nx, nu, scenario);

        allocateCondensingData(&qp, &input, &output, &work);

        calculate_transition_vector(&qp, &work, x0.data());
        calculate_transition_matrix(&qp, &work);

        SECTION(scenario + " Simple bounds", "[condensing]") {
            calculate_simple_bounds(&qp, &output);
            VectorXd true_lb = readVectorFromFile(scenario + "/u_lower_bound.dat", N*nu);
            VectorXd acados_lb = Eigen::Map<VectorXd>(&output.lb[0], N*nu);
            REQUIRE(acados_lb.isApprox(true_lb, COMPARISON_TOLERANCE));
            VectorXd true_ub = readVectorFromFile(scenario + "/u_upper_bound.dat", N*nu);
            VectorXd acados_ub = Eigen::Map<VectorXd>(&output.ub[0], N*nu);
            REQUIRE(acados_ub.isApprox(true_ub, COMPARISON_TOLERANCE));
        }

        SECTION(scenario + " Constraint bounds", "[condensing]") {
            int_t nA = get_num_constraints(&qp, &work);
            calculate_constraint_bounds(&qp, &output, &work, x0.data());
            VectorXd true_lbA = readVectorFromFile(scenario + "/condensed_lower_bound.dat", nA);
            VectorXd acados_lbA = Eigen::Map<VectorXd>(&output.lbA[0], nA);
            REQUIRE(acados_lbA.isApprox(true_lbA, COMPARISON_TOLERANCE));
            VectorXd true_ubA = readVectorFromFile(scenario + "/condensed_upper_bound.dat", nA);
            VectorXd acados_ubA = Eigen::Map<VectorXd>(&output.ubA[0], nA);
            REQUIRE(acados_ubA.isApprox(true_ubA, COMPARISON_TOLERANCE));
        }

        SECTION(scenario + " Constraint matrix", "[condensing]") {
            int_t nA = get_num_constraints(&qp, &work);
            calculate_constraint_matrix(&qp, &output, &work);
            MatrixXd true_A = readMatrixFromFile(scenario + "/transition_matrix.dat", nA, N*nu);
            MatrixXd acados_A = Eigen::Map<MatrixXd>(&output.A[0], nA, N*nu);
            REQUIRE(acados_A.isApprox(true_A, COMPARISON_TOLERANCE));
        }
    }
}

TEST_CASE("Constrained LTV system, general polytopic constraints", "[condensing]") {
    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;
    VectorXd x0;

    for (std::string scenario : scenarios) {
        fillWithUnconstrainedData(&qp, &x0, scenario);
        int_t N = qp.N;
        int_t nx = qp.nx[0];
        int_t nu = qp.nu[0];
        fillWithBoundsData(&qp, N, nx, nu, scenario);
        fillWithGeneralConstraintsData(&qp, N, nx, nu, scenario);

        allocateCondensingData(&qp, &input, &output, &work);

        calculate_transition_vector(&qp, &work, x0.data());
        calculate_transition_matrix(&qp, &work);

        SECTION(scenario + " General constraint bounds", "[condensing]") {
            int_t nA = get_num_constraints(&qp, &work);
            calculate_constraint_bounds(&qp, &output, &work, x0.data());
            VectorXd true_lbA = readVectorFromFile(scenario +
                "/condensed_general_constraint_lb.dat", nA);
            VectorXd acados_lbA = Eigen::Map<VectorXd>(&output.lbA[0], nA);
            REQUIRE(acados_lbA.isApprox(true_lbA, COMPARISON_TOLERANCE));
            VectorXd true_ubA = readVectorFromFile(scenario +
                "/condensed_general_constraint_ub.dat", nA);
            VectorXd acados_ubA = Eigen::Map<VectorXd>(&output.ubA[0], nA);
            REQUIRE(acados_ubA.isApprox(true_ubA, COMPARISON_TOLERANCE));
        }

        SECTION(scenario + " General constraint matrix", "[condensing]") {
            int_t nA = get_num_constraints(&qp, &work);
            calculate_constraint_matrix(&qp, &output, &work);
            MatrixXd true_A = readMatrixFromFile(scenario +
                "/condensed_general_constraint_matrix.dat", nA, N*nu);
            MatrixXd acados_A = Eigen::Map<MatrixXd>(&output.A[0], nA, N*nu);
            REQUIRE(acados_A.isApprox(true_A, COMPARISON_TOLERANCE));
        }
    }
}
