# Build list with all source files to go into the acados library
file(GLOB_RECURSE ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/*.c)
# Exclude helper files
list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/ocp_qp/condensing_helper_functions.c)

if (NOT EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_ooqp.c)
endif()