file(GLOB ACADOS_SRC_SIM
    ${PROJECT_SOURCE_DIR}/acados/sim/sim_erk_integrator.c
    ${PROJECT_SOURCE_DIR}/acados/sim/sim_lifted_irk_integrator.c
)

file(GLOB ACADOS_SRC_OCP_QP
    ${PROJECT_SOURCE_DIR}/acados/ocp_qp/condensing.c
    ${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_condensing_qpoases.c
    ${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_hpmpc.c
)

file(GLOB ACADOS_SRC_OCP_NLP
    ${PROJECT_SOURCE_DIR}/acados/ocp_nlp/ocp_nlp_gn_sqp.c
)

file(GLOB ACADOS_SRC_UTILS
    ${PROJECT_SOURCE_DIR}/acados/utils/print.c
    ${PROJECT_SOURCE_DIR}/acados/utils/timing.c
    ${PROJECT_SOURCE_DIR}/acados/utils/tools.c
)

set(TEST_HPMPC_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_hpmpc.c)
set(TEST_CONDENSING_QPOASES_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_condensing_qpoases.c)
set(TEST_CHEN_SRC ${PROJECT_SOURCE_DIR}/examples/Chen/Chen_model.c)
set(TEST_NMPC_SRC ${PROJECT_SOURCE_DIR}/examples/test_nmpc.c)
set(TEST_CHAIN_SRC ${PROJECT_SOURCE_DIR}/examples/casadi_chain/Chain_model.c)

file(GLOB TEST_VDE_CHAIN_SRC
    ${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde*.c
)

file(GLOB TEST_JAC_CHAIN_SRC
    ${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac*.c
)

set(TEST_CHAIN_OCP_SRC ${PROJECT_SOURCE_DIR}/examples/test_chain.c)
