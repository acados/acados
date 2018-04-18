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
#include "test/test_utils/eigen.h"

// #include "blasfeo/include/blasfeo_target.h"
// #include "blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "catch/include/catch.hpp"

// #include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
// #include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
// #include "acados/ocp_qp/ocp_qp_hpipm.h"
// #include "acados/ocp_qp/ocp_qp_hpmpc.h"
// #include "acados/ocp_qp/ocp_qp_qpdunes.h"
// #include "test/test_utils/read_matrix.h"
// #include "test/test_utils/read_ocp_qp_in.h"
#include "acados/utils/print.h"

#include "acados_c/ocp_qp_interface.h"
#include "acados_c/options.h"

extern "C"
{
ocp_qp_dims *create_ocp_qp_dims_mass_spring(int N, int nx_, int nu_, int nb_, int ng_, int ngN);
ocp_qp_in *create_ocp_qp_in_mass_spring(void *config, int N, int nx_, int nu_, int nb_, int ng_, int ngN);
}

using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

double TOL_HPIPM = 1e-8;
double TOL_HPMPC = 1e-8;
double TOL_QPDUNES = 1e-10;
double TOL_QPOASES = 1e-10;
double TOL_QORE = 1e-10;



ocp_qp_solver_t hashit(std::string const& inString)
{
    if (inString == "SPARSE_HPIPM") return PARTIAL_CONDENSING_HPIPM;
    if (inString == "SPARSE_HPMPC") return PARTIAL_CONDENSING_HPMPC;
    if (inString == "SPARSE_QPDUNES") return PARTIAL_CONDENSING_QPDUNES;

    if (inString == "DENSE_HPIPM") return FULL_CONDENSING_HPIPM;
    if (inString == "DENSE_QPOASES") return FULL_CONDENSING_QPOASES;
    if (inString == "DENSE_QORE") return FULL_CONDENSING_QORE;

    return (ocp_qp_solver_t) -1;
}


double solver_tolerance(std::string const& inString)
{
    if (inString == "SPARSE_HPIPM") return 1e-8;
    if (inString == "SPARSE_HPMPC") return 1e-5;
    if (inString == "SPARSE_QPDUNES") return 1e-10;

    if (inString == "DENSE_HPIPM") return 1e-8;
    if (inString == "DENSE_QPOASES") return 1e-10;
    if (inString == "DENSE_QORE") return 1e-10;

    return -1;
}


TEST_CASE("mass spring example", "[QP solvers]")
{
    vector<std::string> solvers = {"SPARSE_HPIPM", "SPARSE_HPMPC", "DENSE_HPIPM", "DENSE_QPOASES", "DENSE_QORE"};
    // TODO(dimitris): FIX QPDUNES CASE!!!
    // vector<std::string> solvers = {"SPARSE_HPIPM", "SPARSE_HPMPC", "SPARSE_QPDUNES", "DENSE_HPIPM", "DENSE_QPOASES", "DENSE_QORE"};

    /************************************************
     * set up dimensions
     ************************************************/

    int nx_ = 8;   // number of states (it has to be even for the mass-spring system test problem)

    int nu_ = 3;   // number of inputs (controllers) (  1 <= nu_ <= nx_/2 )

    int N = 15;    // horizon length

    int nb_ = 11;  // number of box constrained inputs and states

    int ng_ = 0;   // number of general constraints

    int ngN = 0;   // number of general constraints at last stage


	ocp_qp_dims *qp_dims = create_ocp_qp_dims_mass_spring(N, nx_, nu_, nb_, ng_, ngN);

    ocp_qp_in *qp_in = create_ocp_qp_in_mass_spring(NULL, N, nx_, nu_, nb_, ng_, ngN);

    ocp_qp_out *qp_out = ocp_qp_out_create(NULL, qp_dims);

    ocp_qp_solver_plan plan;

    ocp_qp_xcond_solver_config *config;

    ocp_qp_solver *qp_solver;

    void *opts;

    int acados_return;

    double res[4];

    double max_res;

    double tol;

    for (std::string solver : solvers)
    {

        SECTION(solver)
        {
            plan.qp_solver = hashit(solver);  // convert string to enum

            tol = solver_tolerance(solver);

            config = ocp_qp_config_create(plan);

            opts = ocp_qp_opts_create(config, qp_dims);

            // TODO(dimitris): the t slacks calculated here (after solving with qpDUNES) are different than in the example!!

            // TODO(dimitris): fix low accuracy in qpDUNES

            // TODO(dimitris): test qpdunes both with clipping and qpOASES

            qp_solver = ocp_qp_create(config, qp_dims, opts);

            acados_return = ocp_qp_solve(qp_solver, qp_in, qp_out);

            REQUIRE(acados_return == 0);

            ocp_qp_inf_norm_residuals(qp_dims, qp_in, qp_out, res);

            // if (plan.qp_solver == PARTIAL_CONDENSING_QPDUNES)
            //     print_ocp_qp_out(qp_out);

            max_res = 0.0;
            for (int ii = 0; ii < 4; ii++)
            {
                max_res = (res[ii] > max_res) ? res[ii] : max_res;
            }

            std::cout << "\n---> residuals of " << solver << "\n";
            printf("\ninf norm res: %e, %e, %e, %e\n", res[0], res[1], res[2], res[3]);
            REQUIRE(max_res <= tol);

            free(qp_solver);
            free(opts);
            free(config);
        }
    }  // END_FOR_SOLVERS

    free(qp_out);
    free(qp_in);
    free(qp_dims);

}  // END_TEST_CASE
