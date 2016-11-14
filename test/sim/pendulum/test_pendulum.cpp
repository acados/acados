#include "catch/include/catch.hpp"
#include "test/sim/pendulum/casadi/casadi_pendulum.h"
#include "acados/utils/types.h"
#include "acados/utils/timing.h"
#include "acados/utils/print.h"
#include "acados/sim/sim_erk_integrator.h"

extern real_t COMPARISON_TOLERANCE;

TEST_CASE("ERK simulation without sensitivities", "[simulation]") {
    sim_in  sim_in;
    sim_out sim_out;
    sim_info info;

    sim_RK_opts rk_opts;
    sim_erk_workspace erk_work;
}
