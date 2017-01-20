file(GLOB ACADOS_SRC_SIM
    ${PROJECT_SOURCE_DIR}/acados/sim/sim_erk_integrator.c
    ${PROJECT_SOURCE_DIR}/acados/sim/sim_lifted_irk_integrator.c
)

file(GLOB ACADOS_SRC_OCP_QP
    ${PROJECT_SOURCE_DIR}/acados/ocp_qp/condensing.c
    ${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_condensing_qpoases.c
    ${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_hpmpc.c
)

file(GLOB ACADOS_SRC_UTILS
    ${PROJECT_SOURCE_DIR}/acados/utils/print.c
    ${PROJECT_SOURCE_DIR}/acados/utils/timing.c
    ${PROJECT_SOURCE_DIR}/acados/utils/tools.c
)

set(TEST_HPMPC_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_hpmpc.c)
set(TEST_HPMPC_LIBSTR_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_hpmpc_libstr.c)
set(TEST_HPMPC_PT_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_hpmpc_libstr_partial_tightening.c)
set(TEST_OCP_NLP_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_nlp.c)
set(TEST_OCP_NLP_HPNMPC_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_nlp_hpnmpc.c)
set(TEST_CONDENSING_QPOASES_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_condensing_qpoases.c)
set(TEST_CHEN_SRC ${PROJECT_SOURCE_DIR}/examples/Chen/Chen_model.c)
set(TEST_NMPC_SRC ${PROJECT_SOURCE_DIR}/examples/test_nmpc.c)
set(TEST_NMPC_HPMPC_SRC ${PROJECT_SOURCE_DIR}/examples/test_nmpc_hpmpc.c)
set(TEST_CHAIN_SRC ${PROJECT_SOURCE_DIR}/examples/casadi_chain/Chain_model.c)

file(GLOB TEST_VDE_CHAIN_SRC
    ${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde*.c
)

file(GLOB TEST_JAC_CHAIN_SRC
    ${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac*.c
)

set(TEST_CHAIN_OCP_SRC ${PROJECT_SOURCE_DIR}/examples/test_chain.c)




set(TEST_PENDULUM_SRC ${PROJECT_SOURCE_DIR}/examples/casadi_pendulum/pendulum_model.c)

file(GLOB TEST_VDE_PENDULUM_SRC
    ${PROJECT_SOURCE_DIR}/examples/casadi_pendulum/vde_forw_pendulum.c
)

file(GLOB TEST_JAC_PENDULUM_SRC
    ${PROJECT_SOURCE_DIR}/examples/casadi_pendulum/jac_pendulum.c
)

set(TEST_PENDULUM_PT_OCP_SRC ${PROJECT_SOURCE_DIR}/examples/test_nmpc_pendulum_hpmpc_libstr_pt.c)
set(TEST_PENDULUM_OCP_SRC ${PROJECT_SOURCE_DIR}/examples/test_nmpc_pendulum_hpmpc_libstr.c)
