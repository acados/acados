/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
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

#include "acados/ocp_qp/allocate_ocp_qp.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "test/test_utils/read_matrix.h"
#include "test/test_utils/read_ocp_qp_in.h"

#include "acados/utils/print.h"

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

int_t TEST_OOQP = 1;
real_t TOL_OOQP = 1e-6;
int_t TEST_QPOASES = 1;
real_t TOL_QPOASES = 1e-10;
int_t TEST_QPDUNES = 1;
real_t TOL_QPDUNES = 1e-10;
int_t TEST_HPMPC = 1;
real_t TOL_HPMPC = 1e-5;

static vector<std::string> scenarios = {"ocp_qp/LTI", "ocp_qp/LTV"};
// TODO(dimitris): add back "ONLY_AFFINE" after fixing problem
vector<std::string> constraints = {"UNCONSTRAINED", "ONLY_BOUNDS", "CONSTRAINED"};

// TODO(dimitris): Clean up octave code
TEST_CASE("Solve random OCP_QP", "[QP solvers]") {
    ocp_qp_in qp_in;
    ocp_qp_out qp_out;

    int_t SET_BOUNDS = 0;
    int_t SET_INEQUALITIES = 0;
    int_t SET_x0 = 1;
    int_t QUIET = 1;

    int return_value;
    VectorXd acados_W, acados_PI, true_W, true_PI;

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
//                        true_PI = readMatrixFromFile(scenario +
//                            "/pi_star_ocp_unconstrained.dat", N*nx, 1);
                    } else if (constraint == "ONLY_BOUNDS") {
                        true_W = readMatrixFromFile(scenario +
                            "/w_star_ocp_bounds.dat", (N+1)*nx + N*nu, 1);
//                        true_PI = readMatrixFromFile(scenario +
//                            "/pi_star_ocp_bounds.dat", N*nx, 1);
                    } else if (constraint == "ONLY_AFFINE") {
                        true_W = readMatrixFromFile(scenario +
                            "/w_star_ocp_no_bounds.dat", (N+1)*nx + N*nu, 1);
//                        true_PI = readMatrixFromFile(scenario +
//                            "/pi_star_ocp_no_bounds.dat", N*nx, 1);
                    } else if (constraint == "CONSTRAINED") {
                        true_W = readMatrixFromFile(scenario +
                            "/w_star_ocp_constrained.dat", (N+1)*nx + N*nu, 1);
                        true_PI = readMatrixFromFile(scenario +
                            "/pi_star_ocp_constrained.dat", N*nx, 1);
                    }
                    if (TEST_QPOASES) {
                        SECTION("qpOASES") {
                            std::cout <<"---> TESTING qpOASES with QP: "<< scenario <<
                                ", " << constraint << std::endl;

                            ocp_qp_condensing_qpoases_args args;
                            args.cputime = 100.0;  // maximum cpu time in seconds
                            args.nwsr = 1000;  // maximum number of working set recalculations
                            args.warm_start = 0;  // wam start with dual_sol in memory

                            int workspace_size =
                                ocp_qp_condensing_qpoases_calculate_workspace_size(&qp_in, &args);
                            void *work = malloc(workspace_size);

                            int memory_size =
                                ocp_qp_condensing_qpoases_calculate_memory_size(&qp_in, &args);
                            void *mem = malloc(memory_size);

                            ocp_qp_condensing_qpoases_memory memory;
                            ocp_qp_condensing_qpoases_create_memory(&qp_in, &args, &memory, mem);

                            // TODO(dimitris): also test that qp_in has not changed
                            return_value = \
                                ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, &memory, work);
                            acados_W = Eigen::Map<VectorXd>(qp_out.x[0], (N+1)*nx + N*nu);
                            acados_PI = Eigen::Map<VectorXd>(qp_out.pi[0], N*nx);
                            free(work);
                            free(mem);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_QPOASES));
                            // TODO(dimitris): check multipliers in other solvers too
                            if (constraint == "CONSTRAINED")
                                // for (int j = 0; j < N*nx; j++) {
                                //     printf(" %5.2e \t %5.2e\n", acados_PI(j), true_PI(j));
                                // }
                                // REQUIRE(acados_PI.isApprox(true_PI, TOL_QPOASES));
                            std::cout <<"---> PASSED " << std::endl;
                        }
                    }
                    if (TEST_QPDUNES) {
                        SECTION("qpDUNES") {
                            std::cout <<"---> TESTING qpDUNES with QP: "<< scenario <<
                                ", " << constraint << std::endl;

                            ocp_qp_qpdunes_args args;
                            ocp_qp_qpdunes_memory mem;
                            void *work;

                            ocp_qp_qpdunes_create_arguments(&args, QPDUNES_DEFAULT_ARGUMENTS);

                            int_t workspace_size =
                                ocp_qp_qpdunes_calculate_workspace_size(&qp_in, &args);
                            work = (void*)malloc(workspace_size);

                            int_t mem_return = ocp_qp_qpdunes_create_memory(&qp_in, &args, &mem);
                            REQUIRE(mem_return == 0);

                            return_value = ocp_qp_qpdunes(&qp_in, &qp_out, &args, &mem, work);
                            acados_W = Eigen::Map<VectorXd>(qp_out.x[0], (N+1)*nx + N*nu);
                            free(work);
                            ocp_qp_qpdunes_free_memory(&mem);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_OOQP));
                            std::cout <<"---> PASSED " << std::endl;
                        }
                    }
                    #ifdef OOQP
                    if (TEST_OOQP) {
                        SECTION("OOQP") {
                            std::cout <<"---> TESTING OOQP with QP: "<< scenario <<
                                ", " << constraint << std::endl;

                            ocp_qp_ooqp_args args;
                            ocp_qp_ooqp_memory mem;
                            void *work;

                            args.printLevel = 0;

                            int_t workspace_size =
                                ocp_qp_ooqp_calculate_workspace_size(&qp_in, &args);
                            work = (void*)malloc(workspace_size);

                            int_t mem_return = ocp_qp_ooqp_create_memory(&qp_in, &args, &mem);
                            REQUIRE(mem_return == 0);

                            return_value = ocp_qp_ooqp(&qp_in, &qp_out, &args, &mem, work);
                            acados_W = Eigen::Map<VectorXd>(qp_out.x[0], (N+1)*nx + N*nu);
                            free(work);
                            ocp_qp_ooqp_free_memory(&mem);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_OOQP));
                            std::cout <<"---> PASSED " << std::endl;
                        }
                    }
                    #endif
                    if (TEST_HPMPC) {
                        SECTION("HPMPC") {
                            std::cout <<"---> TESTING HPMPC with QP: "<< scenario <<
                                ", " << constraint << std::endl;

                            ocp_qp_hpmpc_args args;

                            args.N = N;

                            ocp_qp_hpmpc_create_arguments(&args, HPMPC_DEFAULT_ARGUMENTS);

                            double inf_norm_res[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
                            args.inf_norm_res = &inf_norm_res[0];

                            void *workspace = 0;
                            void *mem = 0;

                            int_t workspace_size = 0;
                            workspace_size = ocp_qp_hpmpc_calculate_workspace_size(&qp_in, &args);
                            // printf("workspace_size = %i", workspace_size);
                            v_zeros_align(&workspace, workspace_size);
                            ocp_qp_hpmpc_create_memory(&qp_in, &args, &mem);

                            return_value = ocp_qp_hpmpc(&qp_in, &qp_out, &args, mem, workspace);

                            acados_W = Eigen::Map<VectorXd>(qp_out.x[0], (N+1)*nx + N*nu);
                            free(workspace);
                            REQUIRE(return_value == 0);
                            REQUIRE(acados_W.isApprox(true_W, TOL_HPMPC));
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
