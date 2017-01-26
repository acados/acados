#include <iostream>
#include <string>
#include <vector>
#include "catch/include/catch.hpp"
#include "test/test_utils/read_matrix.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "test/ocp_qp/condensing_test_helper.h"
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#include "OOQP/include/cQpGenSparse.h"
#include "blasfeo/include/blasfeo_d_aux.h"

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

// TODO(dimitris): update qpOASES/quadprog/OOQP accuracies
// to reach COMPARISON_TOLERANCE in the tests
extern real_t COMPARISON_TOLERANCE;

int_t TEST_OOQP = 1;
int_t TEST_QPOASES = 0;

// TODO(dimitris): RE-WRITE THESE TESTS AS THEY ARE NOT STABLE!

static vector<std::string> scenarios = {"LTI", "LTV"};
vector<std::string> constraints = {"CONSTRAINED"};


void concatenateSolution(int_t N, int_t nx, int_t nu, const ocp_qp_out *out, VectorXd *acados_W) {
    int ii, jj;
    double *concatenated_W;

    d_zeros(&concatenated_W, (N+1)*nx+N*nu, 1);

    for (ii = 0; ii < N; ii++) {
        for (jj = 0; jj < nx; jj++) concatenated_W[ii*(nx+nu)+jj] = out->x[ii][jj];
        for (int jj = 0; jj < nu; jj++) concatenated_W[ii*(nx+nu)+nx+jj] = out->u[ii][jj];
    }
    for (jj = 0; jj < nx; jj++) concatenated_W[N*(nx+nu)+jj] = *(out->x[N]+jj);
    *acados_W = Eigen::Map<VectorXd>(concatenated_W, (N+1)*nx + N*nu);
}

// TODO(dimitris): Clean up octave code with Robin
TEST_CASE("Solve random OCP_QP", "[QP solvers]") {
    ocp_qp_in qp_in;
    ocp_qp_out qp_out;

    VectorXd x0;

    for (std::string constraint : constraints) {
        SECTION(constraint) {
            for (std::string scenario : scenarios) {
                SECTION(scenario) {
                    // fill-in qp_in struct with data
                    fillWithUnconstrainedData(&qp_in, &x0, scenario);
                    int_t N  = qp_in.N;
                    int_t nx = qp_in.nx[0];
                    int_t nu = qp_in.nu[0];

                    // store x0 before overwritten from fillWithBoundsData
                    // TODO(dimitris): Me or Robin should fix this to avoid the workaround..
                    real_t *x0val = (real_t *) malloc(sizeof(real_t)*qp_in.nx[0]);
                    for (int ii = 0; ii < qp_in.nx[0]; ii++) {
                        x0val[ii] = qp_in.lb[0][ii];
                    }

                    if (constraint == "ONLY_BOUNDS" || constraint == "CONSTRAINED") {
                        fillWithBoundsData(&qp_in, N, nx, nu, scenario);
                    }
                    if (constraint == "ONLY_AFFINE" || constraint == "CONSTRAINED") {
                        fillWithGeneralConstraintsData(&qp_in, N, nx, nu, scenario);
                    }

                    // fix overwritten x0...
                    real_t *lb0 = (real_t *) qp_in.lb[0];
                    real_t *ub0 = (real_t *) qp_in.ub[0];
                    for (int ii = 0; ii < qp_in.nx[0]; ii++) {
                        lb0[ii] = x0val[ii];
                        ub0[ii] = x0val[ii];
                    }

                    // load optimal solution from quadprog
                    MatrixXd true_W;
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

                    // initialize qp_out struct
                    double *hx[N + 1];
                    double *hu[N + 1];

                    for (int ii = 0; ii < N; ii++) {
                        d_zeros(&hx[ii], nx, 1);
                        d_zeros(&hu[ii], nu, 1);
                    }
                    d_zeros(&hx[N], nx, 1);

                    qp_out.x = hx;
                    qp_out.u = hu;

                    int return_value;
                    VectorXd acados_W;

                    if (TEST_QPOASES) {
                        SECTION("qpOASES") {
                            std::cout <<"---> TESTING qpOASES with QP: "<< scenario <<
                                ", " << constraint << std::endl;

                            ocp_qp_condensing_qpoases_args args;
                            args.dummy = 42.0;

                            initialise_qpoases(&qp_in);

                            return_value = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, NULL);
                            concatenateSolution(N, nx, nu, &qp_out, &acados_W);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, 1e-10));
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
                            // std::cout << "ACADOS output:\n" << acados_W << std::endl;
                            // printf("-------------------\n");
                            // std::cout << "OCTAVE output:\n" << true_W << std::endl;
                            // printf("-------------------\n");
                            // printf("return value = %d\n", return_value);
                            // printf("-------------------\n");
                            concatenateSolution(N, nx, nu, &qp_out, &acados_W);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, 1e-6));
                            std::cout <<"---> PASSED " << std::endl;
                            // TODO(dimitris): also test that qp_in has not changed!!
                            // TODO(dimitris): re-write tests such that freeing the memory below
                            // does not interfere with other tests
                            // ocp_qp_ooqp_free_workspace(&work);
                            // ocp_qp_ooqp_free_memory(&mem);
                        }
                    }
                    for (int ii = 0; ii < N; ii++) {
                        d_free(hx[ii]);
                        d_free(hu[ii]);
                    }
                    d_free(hx[N]);
                }  // END_SECTION_SCENARIOS
            }  // END_FOR_SCENARIOS
        }  // END_SECTION_CONSTRAINTS
    }  // END_FOR_CONSTRAINTS
}  // END_TEST_CASE
