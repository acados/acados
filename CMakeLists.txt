cmake_minimum_required ( VERSION 2.6 )
project ( ACADOS )

set ( CMAKE_C_FLAGS "-std=c99 -O2 -Wall -Wextra -Werror" )

# Define include folders
include_directories ( "${PROJECT_SOURCE_DIR}" )

# Define .c source files
set ( TEST_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_solver.c" )

# Define libraries
set ( ACADOS_LIB "${PROJECT_SOURCE_DIR}/acados/libacados.a" )
set ( HPMPC_LIB "${PROJECT_SOURCE_DIR}/hpmpc/libhpmpc.a" )

add_executable ( test_acados ${TEST_SRC} )
target_link_libraries ( test_acados ${ACADOS_LIB} ${HPMPC_LIB} )
