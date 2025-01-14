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

#include "catch/include/catch.hpp"
//#include "test/test_utils/eigen.h"

#include "acados_c/ocp_qp_interface.h"

extern "C" {
ocp_qp_xcond_solver_dims *create_ocp_qp_dims_mass_spring(ocp_qp_xcond_solver_config *config, int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(ocp_qp_dims *dims);
}

using std::vector;



ocp_qp_solver_t hashit(std::string const &inString)
{
    if (inString == "SPARSE_HPIPM") return PARTIAL_CONDENSING_HPIPM;
    if (inString == "DENSE_HPIPM") return FULL_CONDENSING_HPIPM;
#ifdef ACADOS_WITH_HPMPC
    if (inString == "SPARSE_HPMPC") return PARTIAL_CONDENSING_HPMPC;
#endif
#ifdef ACADOS_WITH_QPOASES
    if (inString == "DENSE_QPOASES") return FULL_CONDENSING_QPOASES;
#endif
#ifdef ACADOS_WITH_DAQP
    if (inString == "DENSE_DAQP") return FULL_CONDENSING_DAQP;
#endif
#ifdef ACADOS_WITH_QPDUNES
    if (inString == "SPARSE_QPDUNES") return PARTIAL_CONDENSING_QPDUNES;
#endif
#ifdef ACADOS_WITH_OOQP
    if (inString == "DENSE_OOQP") return FULL_CONDENSING_OOQP;
    if (inString == "SPARSE_OOQP") return PARTIAL_CONDENSING_OOQP;
#endif
#ifdef ACADOS_WITH_OSQP
    if (inString == "SPARSE_OSQP") return PARTIAL_CONDENSING_OSQP;
#endif
#ifdef ACADOS_WITH_QORE
    if (inString == "DENSE_QORE") return FULL_CONDENSING_QORE;
#endif
    return (ocp_qp_solver_t) -1;
}



double solver_tolerance(std::string const &inString)
{
    if (inString == "SPARSE_HPIPM") return 1e-8;
    if (inString == "SPARSE_HPMPC") return 1e-5;
    // if (inString == "SPARSE_QPDUNES") return 1e-8;
    if (inString == "DENSE_HPIPM") return 1e-8;
    if (inString == "DENSE_QPOASES") return 1e-10;
    if (inString == "DENSE_DAQP") return 1e-10;
    if (inString == "DENSE_QORE") return 1e-10;
    if (inString == "SPARSE_OOQP") return 1e-5;
    if (inString == "DENSE_OOQP") return 1e-5;
    if (inString == "SPARSE_OSQP") return 1e-8;

    return -1;
}



void set_N2(std::string const &inString, ocp_qp_xcond_solver_config *config, void *opts, int N2, int N)
{
    bool option_found = false;

    if ( inString=="SPARSE_HPIPM" | inString=="SPARSE_HPMPC" | inString == "SPARSE_OOQP" | inString == "SPARSE_OSQP" )
    {
		config->opts_set(config, opts, "cond_N", &N2);
    }

    if (inString == "SPARSE_QPDUNES")
    {
		config->opts_set(config, opts, "cond_N", &N2);
    }

}



TEST_CASE("mass spring example", "[QP solvers]")
{
    vector<std::string> solvers = {"DENSE_HPIPM",
                                   "SPARSE_HPIPM"
#ifdef ACADOS_WITH_HPMPC
                                   ,
                                   "SPARSE_HPMPC"
#endif
#ifdef ACADOS_WITH_QPOASES
                                   ,
                                   "DENSE_QPOASES"
#endif
#ifdef ACADOS_WITH_DAQP
                                   ,
                                   "DENSE_DAQP"
#endif
// #ifdef ACADOS_WITH_QPDUNES
//                                    ,"SPARSE_QPDUNES"
// #endif
#ifdef ACADOS_WITH_OOQP
                                   ,
                                   "DENSE_OOQP",
                                   "SPARSE_OOQP"
#endif
#ifdef ACADOS_WITH_OSQP
                                   ,
                                   "SPARSE_OSQP"
#endif
#ifdef ACADOS_WITH_QORE
                                   ,
                                   "DENSE_QORE"
#endif
    };

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

    ocp_qp_xcond_solver_dims *qp_dims;

    ocp_qp_in *qp_in;

    ocp_qp_out *qp_out;

    ocp_qp_solver_plan_t plan;

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

            config = ocp_qp_xcond_solver_config_create(plan);

            qp_dims = create_ocp_qp_dims_mass_spring(config, N, nx_, nu_, nb_, ng_, ngN);

            qp_in = create_ocp_qp_in_mass_spring(qp_dims->orig_dims);

            qp_out = ocp_qp_out_create(qp_dims->orig_dims);

            opts = ocp_qp_xcond_solver_opts_create(config, qp_dims);

            if (sparse_solver)
            {
                N2_length = 3;
#ifdef ACADOS_WITH_HPMPC
                if (plan.qp_solver == PARTIAL_CONDENSING_HPMPC)
                    N2_length = 1;  // TODO(dimitris): fix this
#endif
            }
            else
            {
                N2_length = 1;
            }

            for (int ii = 0; ii < N2_length; ii++)
            {
                SECTION("N2 = " + std::to_string((int) N2_values[ii]))
                {
                    set_N2(solver, config, opts, N2_values[ii], N);

                    qp_solver = ocp_qp_create(config, qp_dims, opts);

                    acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

                    REQUIRE(acados_return == 0);

                    ocp_qp_inf_norm_residuals(qp_dims->orig_dims, qp_in, qp_out, res);

                    max_res = 0.0;
                    for (int ii = 0; ii < 4; ii++)
                    {
                        max_res = (res[ii] > max_res) ? res[ii] : max_res;
                    }

                    std::cout << "\n---> residuals of " << solver << " (N2 = " << N2_values[ii]
                              << ")\n";
                    printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);
                    REQUIRE(max_res <= tol);

                    free(qp_solver);
                }
            }

            free(qp_out);
            free(qp_in);
            free(qp_dims);

            free(opts);
            free(config);

        }  // END_FOR_N2

    }  // END_FOR_SOLVERS

}  // END_TEST_CASE
