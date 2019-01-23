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

#include "catch/include/catch.hpp"
#include "test/test_utils/eigen.h"

#include "acados_c/ocp_qp_interface.h"
#include "acados_c/options_interface.h"

extern "C" {
ocp_qp_dims *create_ocp_qp_dims_mass_spring(int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(void *config, ocp_qp_dims *dims);
}

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

ocp_qp_solver_t hashit(std::string const &inString)
{
    if (inString == "SPARSE_HPIPM") return PARTIAL_CONDENSING_HPIPM;
    if (inString == "SPARSE_HPMPC") return PARTIAL_CONDENSING_HPMPC;
    if (inString == "SPARSE_QPDUNES") return PARTIAL_CONDENSING_QPDUNES;

    if (inString == "DENSE_HPIPM") return FULL_CONDENSING_HPIPM;
    if (inString == "DENSE_QPOASES") return FULL_CONDENSING_QPOASES;
#ifdef ACADOS_WITH_QORE
    if (inString == "DENSE_QORE") return FULL_CONDENSING_QORE;
#endif
#ifdef ACADOS_WITH_OOQP
    if (inString == "DENSE_OOQP") return FULL_CONDENSING_OOQP;
    if (inString == "SPARSE_OOQP") return PARTIAL_CONDENSING_OOQP;
#endif
#ifdef ACADOS_WITH_OSQP
    if (inString == "SPARSE_OSQP") return PARTIAL_CONDENSING_OSQP;
#endif

    return (ocp_qp_solver_t) -1;
}

double solver_tolerance(std::string const &inString)
{
    if (inString == "SPARSE_HPIPM") return 1e-8;
    if (inString == "SPARSE_HPMPC") return 1e-5;
    if (inString == "SPARSE_QPDUNES") return 1e-8;

    if (inString == "DENSE_HPIPM") return 1e-8;
    if (inString == "DENSE_QPOASES") return 1e-10;
    if (inString == "DENSE_QORE") return 1e-10;
    if (inString == "SPARSE_OOQP") return 1e-5;
    if (inString == "DENSE_OOQP") return 1e-5;
    if (inString == "SPARSE_OSQP") return 1e-8;

    return -1;
}

void set_N2(std::string const &inString, void *opts, int N2, int N)
{
    bool option_found = false;

    if (inString == "SPARSE_HPIPM")
    {
        option_found = set_option_int(opts, "sparse_hpipm.N2", N2);
        REQUIRE(option_found == true);
    }

    if (inString == "SPARSE_HPMPC")
    {
        option_found = set_option_int(opts, "hpmpc.N2", N2);
        REQUIRE(option_found == true);
    }

    if (inString == "SPARSE_QPDUNES")
    {
        option_found = set_option_int(opts, "qpdunes.N2", N2);
        REQUIRE(option_found == true);
        if (N2 == N)
        {
            option_found = set_option_int(opts, "qpdunes.clipping", 1);
            REQUIRE(option_found == true);
        }
        else
        {
            option_found = set_option_int(opts, "qpdunes.clipping", 0);
            REQUIRE(option_found == true);
        }
    }

    if (inString == "SPARSE_OOQP")
    {
        option_found = set_option_int(opts, "sparse_ooqp.N2", N2);
        REQUIRE(option_found == true);
    }

    if (inString == "SPARSE_OSQP")
    {
        option_found = set_option_int(opts, "sparse_osqp.N2", N2);
        REQUIRE(option_found == true);
    }
}

TEST_CASE("mass spring example", "[QP solvers]")
{
    vector<std::string> solvers = {"SPARSE_HPIPM",
                                   "SPARSE_HPMPC",
                                   "SPARSE_QPDUNES",
                                   "DENSE_HPIPM",
                                   "DENSE_QPOASES"
#ifdef ACADOS_WITH_OOQP
                                   // ,
                                   // "DENSE_OOQP",
                                   // "SPARSE_OOQP"
#endif
#ifdef ACADOS_WITH_OSQP
                                    ,
                                   "SPARSE_OSQP"
#endif
#ifdef ACADOS_WITH_QORE
                                   ,
                                   "DENSE_QORE"};
#else
    };
#endif

    /************************************************
     * set up dimensions
     ************************************************/

    int nx_ = 8;  // number of states (it has to be even for the mass-spring system test problem)

    int nu_ = 3;  // number of inputs (controllers) (  1 <= nu_ <= nx_/2 )

    int N = 15;  // horizon length

    int nb_ = 11;  // number of box constrained inputs and states

    int ng_ = 0;  // number of general constraints

    int ngN = 0;  // number of general constraints at last stage

    double N2_values[] = {15, 5, 3};  // horizon of partially condensed QP for sparse solvers

    ocp_qp_dims *qp_dims = create_ocp_qp_dims_mass_spring(N, nx_, nu_, nb_, ng_, ngN);

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring(NULL, qp_dims);

    ocp_qp_out *qp_out = ocp_qp_out_create(NULL, qp_dims);

    ocp_qp_solver_plan plan;

    ocp_qp_xcond_solver_config *config;

    ocp_qp_solver *qp_solver;

    void *opts;

    int acados_return;

    double res[4];

    double max_res;

    double tol;

    int N2_length;  // 3 for sparse solvers, 1 for dense solvers

    std::size_t sparse_solver;

    for (std::string solver : solvers)
    {
        SECTION(solver)
        {
            plan.qp_solver = hashit(solver);  // convert string to enum

            sparse_solver = !solver.find("SPARSE");

            tol = solver_tolerance(solver);

            config = ocp_qp_config_create(plan);

            opts = ocp_qp_opts_create(config, qp_dims);

            if (sparse_solver)
            {
                N2_length = 3;
                if (plan.qp_solver == PARTIAL_CONDENSING_HPMPC)
                    N2_length = 1;  // TODO(dimitris): fix this
            }
            else
            {
                N2_length = 1;
            }

            for (int ii = 0; ii < N2_length; ii++)
            {
                SECTION("N2 = " + std::to_string((int)N2_values[ii]))
                {
                    set_N2(solver, opts, N2_values[ii], N);

                    qp_solver = ocp_qp_create(config, qp_dims, opts);

                    acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

                    // TODO(dimitris): fix this hack for qpDUNES
                    // (it terminates one iteration before optimal solution,
                    // fixed with warm-start and calling solve twice)
                    if (plan.qp_solver == PARTIAL_CONDENSING_QPDUNES)
                    {
                        acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);
                    }

                    REQUIRE(acados_return == 0);

                    ocp_qp_inf_norm_residuals(qp_dims, qp_in, qp_out, res);

                    max_res = 0.0;
                    for (int ii = 0; ii < 4; ii++)
                    {
                        max_res = (res[ii] > max_res) ? res[ii] : max_res;
                    }

                    std::cout << "\n---> residuals of " << solver << " (N2 = "
                                                        << N2_values[ii] << ")\n";
                    printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);
                    REQUIRE(max_res <= tol);

                    free(qp_solver);
                }
            }

            free(opts);
            free(config);

        }  // END_FOR_N2

    }  // END_FOR_SOLVERS

    free(qp_out);
    free(qp_in);
    free(qp_dims);

}  // END_TEST_CASE
