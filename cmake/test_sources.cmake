set(UNIT_TESTS_SRC_TEST_UTILS ${PROJECT_SOURCE_DIR}/test/test_utils/read_matrix.cpp
    ${PROJECT_SOURCE_DIR}/test/test_utils/zeros.cpp
) # TODO (dimitris): eliminate this at some point

set(UNIT_TESTS_SRC_SIM ${PROJECT_SOURCE_DIR}/test/sim/pendulum/test_pendulum.cpp
    ${PROJECT_SOURCE_DIR}/test/sim/pendulum/casadi/casadi_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/jac_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_forw_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_adj_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_hess_pendulum.c
)

set(UNIT_TESTS_SRC_OCP_NLP ${PROJECT_SOURCE_DIR}/test/ocp_nlp/chain/test_chain.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_nlp/chain/Chain_model.c
    ${PROJECT_SOURCE_DIR}/build/test/jac_chain_nm2.c
    ${PROJECT_SOURCE_DIR}/build/test/jac_chain_nm3.c
    ${PROJECT_SOURCE_DIR}/build/test/jac_chain_nm4.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_chain_nm2.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_chain_nm3.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_chain_nm4.c
)

set(UNIT_TESTS_SRC_CASADI ${PROJECT_SOURCE_DIR}/build/test/jac_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_forw_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_adj_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_hess_pendulum.c
    ${PROJECT_SOURCE_DIR}/build/test/jac_chain_nm2.c
    ${PROJECT_SOURCE_DIR}/build/test/jac_chain_nm3.c
    ${PROJECT_SOURCE_DIR}/build/test/jac_chain_nm4.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_chain_nm2.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_chain_nm3.c
    ${PROJECT_SOURCE_DIR}/build/test/vde_chain_nm4.c
)

set(UNIT_TESTS_SRC_OCP_QP ${PROJECT_SOURCE_DIR}/test/ocp_qp/test_condensing.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/condensing_test_helper.cpp
    ${PROJECT_SOURCE_DIR}/test/ocp_qp/test_qpsolvers.cpp
    ${PROJECT_SOURCE_DIR}/test/test_utils/read_ocp_qp_in.c
)

set(UNIT_TESTS_SRC ${PROJECT_SOURCE_DIR}/test/all_tests.cpp)
