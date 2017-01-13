#include <iostream>
#include <string>
#include <vector>
#include "catch/include/catch.hpp"
#include "test/test_utils/read_matrix.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "test/ocp_qp/condensing_test_helper.h"
// TODO(dimitris): Add CPP code in blasfeo header for this to work!
#include "blasfeo/include/blasfeo_d_aux.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

extern real_t COMPARISON_TOLERANCE;

std::vector<std::string> scenarios = {"LTI", "LTV"};

void readInputDimensionsFromFile(int_t *N, int_t *nx, int_t *nu, std::string folder) {
    *N = (int_t) readMatrix(folder + "/N.dat")(0, 0);
    REQUIRE(*N > 0);
    *nx = (int_t) readMatrix(folder + "/nx.dat")(0, 0);
    REQUIRE(*nx > 0);
    *nu = (int_t) readMatrix(folder + "/nu.dat")(0, 0);
    REQUIRE(*nu > 0);
}

// TODO(dimitris): Fix problems with bounded and constrained cases (and implement them)
// TODO(dimitris): Clean up octave code
TEST_CASE("Solve random OCP_QP", "[qpOASES]") {

    ocp_qp_in qp_in;
    ocp_qp_out qp_out;

    VectorXd x0;

    int MODE = 2; // TEMP. MODE = 1: bounds

    for (std::string scenario : scenarios) {

        // TODO(dimitris): remove line below once test is run via cmake
        scenario = "../build/test/" + scenario;

        fillWithUnconstrainedData(&qp_in, &x0, scenario);
        int_t N  = qp_in.N;
        int_t nx = qp_in.nx[0];
        int_t nu = qp_in.nu[0];

        // if (MODE > 0) fillWithBoundsData(&qp_in, N, nx, nu, scenario);
        if (MODE > 1) fillWithGeneralConstraintsData(&qp_in, N, nx, nu, scenario);

        double *hx[N + 1];
        double *hu[N + 1];

        for (int ii = 0; ii < N; ii++) {
            d_zeros(&hx[ii], nx, 1);
            d_zeros(&hu[ii], nu, 1);
        }
        d_zeros(&hx[N], nx, 1);

        qp_out.x = hx;
        qp_out.u = hu;

        ocp_qp_condensing_qpoases_args args;
        args.dummy = 42.0;

        int return_value;
        initialise_qpoases(&qp_in);

        // TODO(dimitris): Ask why NULL in example (don't we warm-start?)
        return_value = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, NULL);

        double concatenated_W[(N+1)*nx+N*nu];
        for (int ii = 0; ii < N; ii++) {
            for (int jj = 0; jj < nx; jj++) concatenated_W[ii*(nx+nu)+jj] = *(hx[ii]+jj);
            for (int jj = 0; jj < nu; jj++) concatenated_W[ii*(nx+nu)+nx+jj] = *(hu[ii]+jj);
        }
        for (int jj = 0; jj < nx; jj++) concatenated_W[N*(nx+nu)+jj] = *(hx[N]+jj);

        MatrixXd acados_W = Eigen::Map<VectorXd>(concatenated_W, (N+1)*nx + N*nu);
        MatrixXd true_W = readMatrixFromFile(scenario +
            "/w_star_ocp_unconstrained.dat", (N+1)*nx + N*nu, 1);

        // std::cout << "ACADOS output:\n" << acados_W << std::endl;
        // printf("-------------------\n");
        // std::cout << "OCTAVE output:\n" << true_W << std::endl;
        // printf("-------------------\n");
        // printf("return value = %d\n", return_value);
        // printf("-------------------\n");

        if (MODE == 0) {
            REQUIRE(return_value == 0);
            // TODO(dimitris): update qpOASES or quadprog settings to reach COMPARISON_TOLERANCE
            REQUIRE(acados_W.isApprox(true_W, 1e-10));
        } else if (MODE == 1) {
            REQUIRE(return_value == 0);
        }
    }
}
