#include "catch/include/catch.hpp"
#include "test/sim/pendulum/casadi/casadi_pendulum.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_erk_integrator.h"
#include "test/test_utils/eigen.h"

extern real_t COMPARISON_TOLERANCE;
using Eigen::VectorXd;

TEST_CASE("ERK simulation with adjoint sensitivities", "[simulation]") {
    sim_in  sim_in;
    sim_out sim_out;
    sim_info info;

    int_t NX = 4;
    int_t NU = 1;
    real_t T = 0.1;

    sim_RK_opts rk_opts;
    sim_erk_workspace erk_work;

    real_t adj[NX+NU];
    real_t seed[NX+NU];

    sim_in.nSteps = 2;
    sim_in.step = T/sim_in.nSteps;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.sens_forw = true;
    sim_in.sens_adj = true;
    sim_in.sens_hess = false;
    sim_in.nsens_forw = NX+NU;

    sim_in.VDE_forw = &VDE_forw_pendulum;
    sim_in.VDE_adj = &VDE_adj_pendulum;
    sim_in.jac_fun = &jac_fun_pendulum;

    sim_in.x = (real_t*) malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = (real_t*) malloc(sizeof(*sim_in.u) * (NU));
    sim_in.S_forw = (real_t*) malloc(sizeof(*sim_in.S_forw) * (NX*(NX+NU)));
    sim_in.S_adj = (real_t*) malloc(sizeof(*sim_in.S_adj) * (NX+NU));
    for (int_t i = 0; i < NX*(NX+NU); i++) sim_in.S_forw[i] = 0.0;
    for (int_t i = 0; i < NX; i++) sim_in.S_forw[i*(NX+1)] = 1.0;

    // adjoint seed:
    for (int_t i = 0; i < NX; i++) seed[i] = 1.0;
    for (int_t i = 0; i < NU; i++) seed[NX+i] = 0.0;
    for (int_t i = 0; i < NX+NU; i++) sim_in.S_adj[i] = seed[i];

    sim_out.xn = (real_t*) malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.S_forw = (real_t*) malloc(sizeof(*sim_out.S_forw) * (NX*(NX+NU)));
    sim_out.info = &info;

    sim_erk_create_opts(4, &rk_opts);
    sim_erk_create_workspace(&sim_in, &rk_opts, &erk_work);

    SECTION("Adjoint sensitivities") {
        for (int_t i = 0; i < NX; i++) sim_in.x[i] = 0.0;
        sim_in.u[0] = 0.1;

        sim_erk(&sim_in, &sim_out, &rk_opts, &erk_work);

        for (int_t i = 0; i < NX+NU; i++) adj[i] = 0.0;
        for (int_t j = 0; j < NX+NU; j++) {
            for (int_t i = 0; i < NX; i++) {
                adj[j] += seed[i]*sim_in.S_forw[j*NX+i];
            }
        }

        print_matrix_name((char*)"stdout", (char*)"adj", adj, 1, NX+NU);
        print_matrix_name((char*)"stdout", (char*)"adj_test", sim_in.S_adj, 1, NX+NU);

        VectorXd true_adj = Eigen::Map<VectorXd>(&adj[0], NX+NU);
        VectorXd test_adj = Eigen::Map<VectorXd>(&sim_in.S_adj[0], NX+NU);
        REQUIRE(true_adj.isApprox(true_adj, COMPARISON_TOLERANCE));
    }
}
