#include <iostream>
#include <string>
#include <vector>
#include "catch/include/catch.hpp"
#include "OOQP/include/cQpGenSparse.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "test/test_utils/read_matrix.h"
#include "test/test_utils/read_ocp_qp_in.h"
#include "test/ocp_qp/condensing_test_helper.h"
#include "acados/utils/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_ooqp.h"

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

int_t TEST_OOQP = 1;
real_t TOL_OOQP = 1e-6;
int_t TEST_QPOASES = 1;
real_t TOL_QPOASES = 1e-10;

static vector<std::string> scenarios = {"LTI", "LTV"};
vector<std::string> constraints = {"UNCONSTRAINED", "ONLY_BOUNDS", "ONLY_AFFINE", "CONSTRAINED"};

// TODO(dimitris): Clean up octave code
TEST_CASE("Solve random OCP_QP", "[QP solvers]") {
    ocp_qp_in qp_in;
    ocp_qp_out qp_out;

    int_t SET_BOUNDS = 0;
    int_t SET_INEQUALITIES = 0;
    int_t SET_x0 = 1;
    int_t QUIET = 1;

    int return_value;
    VectorXd acados_W, true_W;

    for (std::string constraint : constraints) {
        SECTION(constraint) {
            if (constraint == "CONSTRAINED" || constraint == "ONLY_BOUNDS") SET_BOUNDS = 1;
            if (constraint == "CONSTRAINED" || constraint == "ONLY_AFFINE") SET_INEQUALITIES = 1;

            for (std::string scenario : scenarios) {
                SECTION(scenario) {
                    read_ocp_qp_in(&qp_in, (char*) scenario.c_str(),
                    SET_BOUNDS, SET_INEQUALITIES, SET_x0, QUIET);
                    allocate_ocp_qp_out(&qp_in, &qp_out);

                    // TODO(dimitris): extend to variable dimensions
                    int_t N = qp_in.N;
                    int_t nx = qp_in.nx[0];
                    int_t nu = qp_in.nu[0];

                    // load optimal solution from quadprog
                    if (constraint == "UNCONSTRAINED") {
                        true_W = readMatrixFromFile(scenario +
                            "/w_star_ocp_unconstrained.dat", (N+1)*nx + N*nu, 1);
                    } else if (constraint == "ONLY_BOUNDS") {
                        true_W = readMatrixFromFile(scenario +
                            "/w_star_ocp_bounds.dat", (N+1)*nx + N*nu, 1);
                    } else if (constraint == "ONLY_AFFINE") {
                        true_W = readMatrixFromFile(scenario +
                            "/w_star_ocp_no_bounds.dat", (N+1)*nx + N*nu, 1);
                    } else if (constraint == "CONSTRAINED") {
                        true_W = readMatrixFromFile(scenario +
                            "/w_star_ocp_constrained.dat", (N+1)*nx + N*nu, 1);
                    }
                    if (TEST_QPOASES) {
                        SECTION("qpOASES") {
                            std::cout <<"---> TESTING qpOASES with QP: "<< scenario <<
                                ", " << constraint << std::endl;

                            ocp_qp_condensing_qpoases_args args;
                            args.dummy = 42.0;

                            initialise_qpoases(&qp_in);
                            // TODO(dimitris): also test that qp_in has not changed
                            return_value = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, NULL);
                            acados_W = Eigen::Map<VectorXd>(qp_out.x[0], (N+1)*nx + N*nu);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_QPOASES));
                            std::cout <<"---> PASSED " << std::endl;
                        }
                    }
                    if (TEST_OOQP) {
                        SECTION("OOQP") {
                            std::cout <<"---> TESTING OOQP with QP: "<< scenario <<
                                ", " << constraint << std::endl;

                            ocp_qp_ooqp_args args;
                            ocp_qp_ooqp_memory mem;
                            ocp_qp_ooqp_workspace work;

                            args.printLevel = 0;

                            int_t mem_return = ocp_qp_ooqp_create_memory(&qp_in, &args, &mem);
                            REQUIRE(mem_return == 0);
                            int_t work_return = ocp_qp_ooqp_create_workspace(&qp_in, &args, &work);
                            REQUIRE(work_return == 0);

                            return_value = ocp_qp_ooqp(&qp_in, &qp_out, &args, &mem, &work);
                            acados_W = Eigen::Map<VectorXd>(qp_out.x[0], (N+1)*nx + N*nu);
                            // TODO(dimitris): WHY FREEING MEMORY CRASHES CATCH?!?!
                            // (IF TEST OCP_QP IS BEFORE OCP_NLP)
                            ocp_qp_ooqp_free_workspace(&work);
                            ocp_qp_ooqp_free_memory(&mem);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_OOQP));
                            std::cout <<"---> PASSED " << std::endl;
                        }
                    }
                    // std::cout << "ACADOS output:\n" << acados_W << std::endl;
                    // printf("-------------------\n");
                    // std::cout << "OCTAVE output:\n" << true_W << std::endl;
                    // printf("-------------------\n");
                    // printf("return value = %d\n", return_value);
                    // printf("-------------------\n");
                    free_ocp_qp_in(&qp_in);
                    free_ocp_qp_out(&qp_out);
                }  // END_SECTION_SCENARIOS
            }  // END_FOR_SCENARIOS
        }  // END_SECTION_CONSTRAINTS
    }  // END_FOR_CONSTRAINTS
}  // END_TEST_CASE
