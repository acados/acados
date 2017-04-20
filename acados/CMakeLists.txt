# Build list with all source files to go into the acados library
file(GLOB_RECURSE ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/*.c)
# Exclude helper files
list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/ocp_qp/condensing_helper_functions.c)

if (NOT EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP)
    list(REMOVE_ITEM ACADOS_SRC ${PROJECT_SOURCE_DIR}/acados/ocp_qp/ocp_qp_ooqp.c)
endif()

# Define acados library
add_library(acados STATIC ${ACADOS_SRC})
target_include_directories(acados
    PUBLIC
        "${ACADOS_SRC_DIR}/.."
        "${EXTERNAL_SRC_DIR}"
        "${EXTERNAL_SRC_DIR}/qpDUNES-dev/include/"
        "${EXTERNAL_SRC_DIR}/qpOASES/include/"
)
target_link_libraries(acados
    PRIVATE
        qpoases
        hpmpc
        blasfeo
        qpdunes
        m
)
if(EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP AND NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
    target_compile_definitions(acados PRIVATE OOQP)
    target_link_libraries(acados
        PUBLIC 
            ooqpbase
            ooqpgensparse
            ooqpsparse
            ooqpgondzio
            ma27
            blas
            ${FORTRAN_LIBRARY}
            m
    )
endif()

target_compile_definitions(acados PRIVATE USE_ACADOS_TYPES) # Needed for qpDUNES
