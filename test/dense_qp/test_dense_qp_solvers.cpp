/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "catch/include/catch.hpp"

#include "acados_c/dense_qp_interface.h"

using std::vector;


struct dense_qp_data
{
    int nv, ne, nb, ng, ns;
    vector<double> H, g, A, b;
    vector<int> idxb;
    vector<double> d_lb, d_ub, C, d_lg, d_ug;
    vector<int> idxs;
    vector<double> Zl, Zu, zl, zu, d_ls, d_us;
};


struct dense_qp_result
{
    int status;
    vector<double> v, pi;
    double res[4];
};


static dense_qp_result solve_dense_qp(dense_qp_solver_t solver_type, dense_qp_data &qp)
{
    dense_qp_solver_plan plan;
    plan.qp_solver = solver_type;
    qp_solver_config *config = dense_qp_config_create(&plan);

    dense_qp_dims *dims = dense_qp_dims_create();
    d_dense_qp_dim_set_all(qp.nv, qp.ne, qp.nb, qp.ng, qp.ns, dims);

    dense_qp_in *qp_in = dense_qp_in_create(config, dims);
    d_dense_qp_set_all(qp.H.data(), qp.g.data(), qp.A.data(), qp.b.data(), qp.idxb.data(),
                       qp.d_lb.data(), qp.d_ub.data(), qp.C.data(), qp.d_lg.data(), qp.d_ug.data(),
                       qp.Zl.data(), qp.Zu.data(), qp.zl.data(), qp.zu.data(), qp.idxs.data(), NULL,
                       qp.d_ls.data(), qp.d_us.data(), qp_in);

    void *opts = dense_qp_opts_create(config, dims);
    dense_qp_out *qp_out = dense_qp_out_create(config, dims);
    dense_qp_solver *solver = dense_qp_create(config, dims, opts);

    dense_qp_result result;
    result.status = dense_qp_solve(solver, qp_in, qp_out);
    result.v.resize(qp.nv);
    result.pi.resize(qp.ne);
    d_dense_qp_sol_get_v(qp_out, result.v.data());
    d_dense_qp_sol_get_pi(qp_out, result.pi.data());
    dense_qp_inf_norm_residuals(dims, qp_in, qp_out, result.res);

    free(solver);
    free(qp_out);
    free(opts);
    free(qp_in);
    free(dims);
    free(config);

    return result;
}


// min 1/2 |v|^2  s.t.  v0 + v1 = 1  ->  v = (0.5, 0.5), pi = 0.5
static dense_qp_data equality_only_qp()
{
    dense_qp_data qp;
    qp.nv = 2;
    qp.ne = 1;
    qp.nb = 0;
    qp.ng = 0;
    qp.ns = 0;
    qp.H = {1.0, 0.0, 0.0, 1.0};
    qp.g = {0.0, 0.0};
    qp.A = {1.0, 1.0};
    qp.b = {1.0};
    return qp;
}


// min 1/2 |v|^2 - v2  s.t.  v0 + v1 = 1, 0.2 <= v0 - v1 <= 10, -10 <= v2 <= 0.5  ->  v = (0.6, 0.4, 0.5), pi = 0.5
static dense_qp_data equality_bound_general_qp()
{
    dense_qp_data qp;
    qp.nv = 3;
    qp.ne = 1;
    qp.nb = 1;
    qp.ng = 1;
    qp.ns = 0;
    qp.H = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    qp.g = {0.0, 0.0, -1.0};
    qp.A = {1.0, 1.0, 0.0};
    qp.b = {1.0};
    qp.idxb = {2};
    qp.d_lb = {-10.0};
    qp.d_ub = {0.5};
    qp.C = {1.0, -1.0, 0.0};
    qp.d_lg = {0.2};
    qp.d_ug = {10.0};
    return qp;
}


// soft general constraint with 1/2 * 100 * s^2  ->  s = 0.1 / 100.5, v0 - v1 = 0.2 - s, pi = 0.5
static dense_qp_data equality_soft_general_qp()
{
    dense_qp_data qp = equality_bound_general_qp();
    qp.ns = 1;
    qp.idxs = {qp.nb + 0};
    qp.Zl = {100.0};
    qp.Zu = {100.0};
    qp.zl = {0.0};
    qp.zu = {0.0};
    qp.d_ls = {0.0};
    qp.d_us = {0.0};
    return qp;
}


// without the equality constraint  ->  v0 = -v1 = (0.2 - s) / 2, v2 = 0.5
static dense_qp_data soft_general_qp()
{
    dense_qp_data qp = equality_soft_general_qp();
    qp.ne = 0;
    qp.A = {};
    qp.b = {};
    return qp;
}


TEST_CASE("dense QP solvers with equality and soft constraints", "[dense QP solvers]")
{
    vector<std::string> solver_names = {"HPIPM"};
    vector<dense_qp_solver_t> solver_types = {DENSE_QP_HPIPM};
#ifdef ACADOS_WITH_QPOASES
    solver_names.push_back("QPOASES");
    solver_types.push_back(DENSE_QP_QPOASES);
#endif
#ifdef ACADOS_WITH_DAQP
    solver_names.push_back("DAQP");
    solver_types.push_back(DENSE_QP_DAQP);
#endif

    const double tol = 1e-6;
    const double slack = 0.1 / 100.5;

    for (size_t ii = 0; ii < solver_types.size(); ii++)
    {
        SECTION(solver_names[ii] + ", equality constraints only")
        {
            dense_qp_data qp = equality_only_qp();
            dense_qp_result result = solve_dense_qp(solver_types[ii], qp);

            REQUIRE(result.status == 0);
            REQUIRE(std::fabs(result.v[0] - 0.5) <= tol);
            REQUIRE(std::fabs(result.v[1] - 0.5) <= tol);
            REQUIRE(std::fabs(result.pi[0] - 0.5) <= tol);
            for (int jj = 0; jj < 4; jj++)
                REQUIRE(result.res[jj] <= tol);
        }

        SECTION(solver_names[ii] + ", equality, bound and general constraints")
        {
            dense_qp_data qp = equality_bound_general_qp();
            dense_qp_result result = solve_dense_qp(solver_types[ii], qp);

            REQUIRE(result.status == 0);
            REQUIRE(std::fabs(result.v[0] - 0.6) <= tol);
            REQUIRE(std::fabs(result.v[1] - 0.4) <= tol);
            REQUIRE(std::fabs(result.v[2] - 0.5) <= tol);
            REQUIRE(std::fabs(result.pi[0] - 0.5) <= tol);
            for (int jj = 0; jj < 4; jj++)
                REQUIRE(result.res[jj] <= tol);
        }

        SECTION(solver_names[ii] + ", soft general constraint and a bound")
        {
            dense_qp_data qp = soft_general_qp();
            dense_qp_result result = solve_dense_qp(solver_types[ii], qp);

            REQUIRE(result.status == 0);
            REQUIRE(std::fabs(result.v[0] - (0.2 - slack) / 2.0) <= tol);
            REQUIRE(std::fabs(result.v[1] + (0.2 - slack) / 2.0) <= tol);
            REQUIRE(std::fabs(result.v[2] - 0.5) <= tol);
            for (int jj = 0; jj < 4; jj++)
                REQUIRE(result.res[jj] <= tol);
        }

        SECTION(solver_names[ii] + ", equality and soft general constraint")
        {
            dense_qp_data qp = equality_soft_general_qp();
            dense_qp_result result = solve_dense_qp(solver_types[ii], qp);

            REQUIRE(result.status == 0);
            REQUIRE(std::fabs(result.v[0] - (1.0 + 0.2 - slack) / 2.0) <= tol);
            REQUIRE(std::fabs(result.v[1] - (1.0 - 0.2 + slack) / 2.0) <= tol);
            REQUIRE(std::fabs(result.v[2] - 0.5) <= tol);
            REQUIRE(std::fabs(result.pi[0] - 0.5) <= tol);
            for (int jj = 0; jj < 4; jj++)
                REQUIRE(result.res[jj] <= tol);
        }
    }
}  // END_TEST_CASE
