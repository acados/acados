#include "catch/include/catch.hpp"

#include "test/test_utils/read_matrix.hpp"
#include "acados/ocp_qp/condensing.c"

TEST_CASE("Transition vector", "[condensing]") {
    ocp_qp_in qp;
    condensing_in input;
    condensing_out output;
    condensing_workspace work;
    real_t *x0;

    Eigen::MatrixXd A = readMatrix("A.dat");
    std::cout << A << std::endl;

    // calculate_transition_vector(&qp, &work, x0);


    REQUIRE(true);
}
