#include "catch/include/catch.hpp"

#include "acados/ocp_qp/condensing.c"

TEST_CASE("Transition vector", "[condensing]") {

    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;
    real_t *x0;

    calculate_transition_vector(&qp, &work, x0);


    REQUIRE(true);
}
