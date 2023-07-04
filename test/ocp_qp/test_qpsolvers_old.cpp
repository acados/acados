/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include <iostream>
#include <string>
#include <vector>

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "catch/include/catch.hpp"

#ifdef OOQP
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#endif

#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "test/test_utils/read_matrix.h"
#include "test/test_utils/read_ocp_qp_in.h"

#include "acados/utils/print.h"

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

// TODO(dimitris): enable tests of condensing solvers after updating hpipm submodule

int_t TEST_OOQP = 1;
real_t TOL_OOQP = 1e-6;
int_t TEST_QPOASES = 0;
real_t TOL_QPOASES = 1e-10;
int_t TEST_QPDUNES = 1;
real_t TOL_QPDUNES = 1e-10;
int_t TEST_HPMPC = 0;
real_t TOL_HPMPC = 1e-5;
int_t TEST_CON_HPIPM = 0;
real_t TOL_CON_HPIPM = 1e-5;
int_t TEST_HPIPM = 0;
real_t TOL_HPIPM = 1e-5;

static vector<std::string> scenarios = {"ocp_qp/LTI", "ocp_qp/LTV"};
// TODO(dimitris): add back "ONLY_AFFINE" after fixing problem
vector<std::string> constraints = {"UNCONSTRAINED", "ONLY_BOUNDS", "CONSTRAINED"};

// TODO(dimitris): Clean up octave code
TEST_CASE("Solve random OCP_QP", "[QP solvers]")
{
    ocp_qp_in *qp_in;
    ocp_qp_out *qp_out;

    int_t SET_BOUNDS = 0;
    int_t SET_INEQUALITIES = 0;
    int_t SET_x0 = 1;
    int_t QUIET = 1;

    int return_value;
    VectorXd acados_W, acados_PI, true_W, true_PI;

    for (std::string constraint : constraints)
    {
        SECTION(constraint)
        {
            if (constraint == "CONSTRAINED" || constraint == "ONLY_BOUNDS") SET_BOUNDS = 1;
            if (constraint == "CONSTRAINED" || constraint == "ONLY_AFFINE") SET_INEQUALITIES = 1;

            for (std::string scenario : scenarios)
            {
                SECTION(scenario)
                {
                    qp_in = read_ocp_qp_in((char *) scenario.c_str(), SET_BOUNDS, SET_INEQUALITIES,
                                           SET_x0, QUIET);
                    qp_out = ocp_qp_out_create(qp_in->N, (int *) qp_in->nx, (int *) qp_in->nu,
                                               (int *) qp_in->nb, (int *) qp_in->nc);

                    // TODO(dimitris): extend to variable dimensions
                    int_t N = qp_in->N;
                    int_t nx = qp_in->nx[0];
                    int_t nu = qp_in->nu[0];

                    // load optimal solution from quadprog
                    if (constraint == "UNCONSTRAINED")
                    {
                        true_W = readMatrixFromFile(scenario + "/w_star_ocp_unconstrained.dat",
                                                    (N + 1) * nx + N * nu, 1);
                        //                        true_PI = readMatrixFromFile(scenario +
                        //                            "/pi_star_ocp_unconstrained.dat", N*nx, 1);
                    }
                    else if (constraint == "ONLY_BOUNDS")
                    {
                        true_W = readMatrixFromFile(scenario + "/w_star_ocp_bounds.dat",
                                                    (N + 1) * nx + N * nu, 1);
                        //                        true_PI = readMatrixFromFile(scenario +
                        //                            "/pi_star_ocp_bounds.dat", N*nx, 1);
                    }
                    else if (constraint == "ONLY_AFFINE")
                    {
                        true_W = readMatrixFromFile(scenario + "/w_star_ocp_no_bounds.dat",
                                                    (N + 1) * nx + N * nu, 1);
                        //                        true_PI = readMatrixFromFile(scenario +
                        //                            "/pi_star_ocp_no_bounds.dat", N*nx, 1);
                    }
                    else if (constraint == "CONSTRAINED")
                    {
                        true_W = readMatrixFromFile(scenario + "/w_star_ocp_constrained.dat",
                                                    (N + 1) * nx + N * nu, 1);
                        true_PI = readMatrixFromFile(scenario + "/pi_star_ocp_constrained.dat",
                                                     N * nx, 1);
                    }
                    if (TEST_QPOASES)
                    {
                        SECTION("qpOASES")
                        {
                            std::cout << "---> TESTING qpOASES with QP: " << scenario << ", "
                                      << constraint << std::endl;

                            ocp_qp_solver *solver =
                                create_ocp_qp_solver(qp_in, "condensing_qpoases", NULL);

                            // TODO(dimitris): also test that qp_in has not changed
                            return_value = solver->fun(solver->qp_in, solver->qp_out, solver->args,
                                                       solver->mem, solver->work);

                            acados_W =
                                Eigen::Map<VectorXd>(solver->qp_out->x[0], (N + 1) * nx + N * nu);
                            acados_PI = Eigen::Map<VectorXd>(solver->qp_out->pi[0], N * nx);

                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_QPOASES));
                            // TODO(dimitris): check multipliers in other solvers too
                            if (constraint == "CONSTRAINED")
                            {
                                // for (int j = 0; j < N*nx; j++) {
                                //     printf(" %5.2e \t %5.2e\n", acados_PI(j), true_PI(j));
                                // }
                                // TODO(dimitris): re-enable this once HPIPM is updated in acados
                                // REQUIRE(acados_PI.isApprox(true_PI, TOL_QPOASES));
                            }
                            std::cout << "---> PASSED " << std::endl;
                        }
                    }
                    if (TEST_QPDUNES)
                    {
                        SECTION("qpDUNES")
                        {
                            std::cout << "---> TESTING qpDUNES with QP: " << scenario << ", "
                                      << constraint << std::endl;

                            ocp_qp_solver *solver = create_ocp_qp_solver(qp_in, "qpdunes", NULL);

                            return_value = solver->fun(solver->qp_in, solver->qp_out, solver->args,
                                                       solver->mem, solver->work);

                            acados_W =
                                Eigen::Map<VectorXd>(solver->qp_out->x[0], (N + 1) * nx + N * nu);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_OOQP));
                            std::cout << "---> PASSED " << std::endl;
                        }
                    }
#ifdef OOQP
                    if (TEST_OOQP)
                    {
                        SECTION("OOQP")
                        {
                            std::cout << "---> TESTING OOQP with QP: " << scenario << ", "
                                      << constraint << std::endl;

                            ocp_qp_solver *solver = create_ocp_qp_solver(qp_in, "ooqp", NULL);

                            return_value = solver->fun(solver->qp_in, solver->qp_out, solver->args,
                                                       solver->mem, solver->work);

                            acados_W =
                                Eigen::Map<VectorXd>(solver->qp_out->x[0], (N + 1) * nx + N * nu);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_OOQP));
                            std::cout << "---> PASSED " << std::endl;
                        }
                    }
#endif
                    if (TEST_HPMPC)
                    {
                        SECTION("HPMPC")
                        {
                            std::cout << "---> TESTING HPMPC with QP: " << scenario << ", "
                                      << constraint << std::endl;

                            ocp_qp_solver *solver = create_ocp_qp_solver(qp_in, "hpmpc", NULL);

                            // TODO(dimitris): also test that qp_in has not changed
                            return_value = solver->fun(solver->qp_in, solver->qp_out, solver->args,
                                                       solver->mem, solver->work);

                            acados_W =
                                Eigen::Map<VectorXd>(solver->qp_out->x[0], (N + 1) * nx + N * nu);

                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_HPMPC));
                            std::cout << "---> PASSED " << std::endl;
                        }
                    }
                    if (TEST_CON_HPIPM)
                    {
                        SECTION("CONDENSING_HPIPM")
                        {
                            std::cout << "---> TESTING condensing + HPIPM with QP: " << scenario
                                      << ", " << constraint << std::endl;

                            // ocp_qp_condensing_hpipm_args args;
                            // args.mu_max = 1e-8;
                            // args.iter_max = 30;
                            // args.alpha_min = 1e-8;
                            // args.mu0 = 1;

                            ocp_qp_solver *solver =
                                create_ocp_qp_solver(qp_in, "condensing_hpipm", NULL);

                            return_value = solver->fun(solver->qp_in, solver->qp_out, solver->args,
                                                       solver->mem, solver->work);

                            acados_W =
                                Eigen::Map<VectorXd>(solver->qp_out->x[0], (N + 1) * nx + N * nu);

                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_CON_HPIPM));
                            if (constraint == "CONSTRAINED")
                            {
                                // for (int j = 0; j < N*nx; j++) {
                                //     printf(" %5.2e \t %5.2e\n", acados_PI(j), true_PI(j));
                                // }
                                // REQUIRE(acados_PI.isApprox(true_PI, TOL_CON_HPIPM));
                            }
                            std::cout << "---> PASSED " << std::endl;
                        }
                    }
                    if (TEST_HPIPM)
                    {
                        SECTION("HPIPM")
                        {
                            std::cout << "---> TESTING HPIPM with QP: " << scenario << ", "
                                      << constraint << std::endl;

                            ocp_qp_solver *solver = create_ocp_qp_solver(qp_in, "hpipm", NULL);

                            // TODO(dimitris): also test that qp_in has not changed
                            return_value = solver->fun(solver->qp_in, solver->qp_out, solver->args,
                                                       solver->mem, solver->work);

                            acados_W =
                                Eigen::Map<VectorXd>(solver->qp_out->x[0], (N + 1) * nx + N * nu);
                            acados_PI = Eigen::Map<VectorXd>(solver->qp_out->pi[0], N * nx);

                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_HPIPM));
                            if (constraint == "CONSTRAINED")
                            {
                                // for (int j = 0; j < N*nx; j++) {
                                //     printf(" %5.2e \t %5.2e\n", acados_PI(j), true_PI(j));
                                // }
                                REQUIRE(acados_PI.isApprox(true_PI, TOL_HPIPM));
                            }
                            std::cout << "---> PASSED " << std::endl;
                        }
                    }
                    // std::cout << "ACADOS output:\n" << acados_W << std::endl;
                    // printf("-------------------\n");
                    // std::cout << "OCTAVE output:\n" << true_W << std::endl;
                    // printf("-------------------\n");
                    // printf("return value = %d\n", return_value);
                    // printf("-------------------\n");
                    free(qp_in);
                    free(qp_out);
                }  // END_SECTION_SCENARIOS
            }      // END_FOR_SCENARIOS
        }          // END_SECTION_CONSTRAINTS
    }              // END_FOR_CONSTRAINTS
}  // END_TEST_CASE
