# Build list with all source files to go into the acados library
file(GLOB_RECURSE ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/*.c)
# Exclude helper files
list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/ocp_qp/condensing_helper_functions.c)

# Sources for examples
set(HPMPC_EXAMPLE_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_hpmpc.c)
set(CONDENSING_QPOASES_EXAMPLE_SRC ${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_condensing_qpoases.c)
set(CHEN_MODEL_SRC ${PROJECT_SOURCE_DIR}/examples/Chen/Chen_model.c)
set(NMPC_EXAMPLE_SRC ${PROJECT_SOURCE_DIR}/examples/test_nmpc.c)
file(GLOB CHAIN_EXAMPLE_SRC ${PROJECT_SOURCE_DIR}/examples/casadi_chain/Chain_model.c
    ${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde*.c
    ${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac*.c
    ${PROJECT_SOURCE_DIR}/examples/test_chain.c
)
