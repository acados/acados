file(GLOB UNIT_TESTS_SRC_TEST_UTILS
    ${PROJECT_SOURCE_DIR}/test/test_utils/read_matrix.cpp
    ${PROJECT_SOURCE_DIR}/test/test_utils/zeros.cpp
)

file(GLOB UNIT_TESTS_SRC_OCP_QP
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/test_condensing.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/condensing_test_helper.cpp
)

file(GLOB UNIT_TESTS_SRC_SIM
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/test_pendulum.cpp
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/casadi/casadi_pendulum.c
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/casadi/jac_pendulum.c
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/casadi/vde_forw_pendulum.c
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/casadi/vde_adj_pendulum.c
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/casadi/vde_hess_pendulum.c
)

set(UNIT_TESTS_SRC ${PROJECT_SOURCE_DIR}/test/all_tests.cpp)
