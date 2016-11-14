file(GLOB UNIT_TESTS_SRC_TEST_UTILS
    ${PROJECT_SOURCE_DIR}/test/test_utils/read_matrix.cpp
    ${PROJECT_SOURCE_DIR}/test/test_utils/zeros.cpp
)

file(GLOB UNIT_TESTS_SRC_OCP_QP
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/test_LTI_condensing.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/test_LTV_condensing.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/condensing_test_helper.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/LTI_condensing_test_helper.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/LTV_condensing_test_helper.cpp
)

file(GLOB UNIT_TESTS_SRC_SIM
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/casadi/*.c
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/test_pendulum.cpp
)

set(UNIT_TESTS_SRC ${PROJECT_SOURCE_DIR}/test/all_tests.cpp)
