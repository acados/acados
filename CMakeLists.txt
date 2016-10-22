cmake_minimum_required ( VERSION 2.6 )
project ( ACADOS )

set ( CMAKE_VERBOSE_MAKEFILE ON )

set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 -Wall -Wextra -Werror" )

if ( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Debug )
endif ( NOT CMAKE_BUILD_TYPE )

if ( CMAKE_BUILD_TYPE MATCHES Debug )
    add_definitions ( -DDEBUG )
endif ( CMAKE_BUILD_TYPE MATCHES Debug )

# Define include folders
include_directories ( "${PROJECT_SOURCE_DIR}" )
include_directories ( "${PROJECT_SOURCE_DIR}/external/" )
include_directories ( "${PROJECT_SOURCE_DIR}/external/qpOASES/include/" )

# Define .c source files
file (GLOB ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/*.c")
file (GLOB ACADOS_SRC_SIM "${PROJECT_SOURCE_DIR}/acados/sim/*.c")
file (GLOB ACADOS_SRC_QP "${PROJECT_SOURCE_DIR}/acados/ocp_qp/*.c")
file (GLOB ACADOS_SRC_UTILS "${PROJECT_SOURCE_DIR}/acados/utils/*.c")
set ( TEST_HPMPC_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_hpmpc" )
set ( TEST_CONDENSING_QPOASES_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_condensing_qpoases" )
set ( TEST_CHEN_SRC "${PROJECT_SOURCE_DIR}/test/Chen/Chen_model" )
set ( TEST_NMPC_SRC "${PROJECT_SOURCE_DIR}/test/test_nmpc" )
set ( TEST_CHAIN_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/Chain_model" )
set ( TEST_VDE_CHAIN2_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm2" )
set ( TEST_VDE_CHAIN3_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm3" )
set ( TEST_VDE_CHAIN4_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm4" )
set ( TEST_VDE_CHAIN5_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm5" )
set ( TEST_VDE_CHAIN6_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm6" )
set ( TEST_VDE_CHAIN7_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm7" )
set ( TEST_VDE_CHAIN8_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm8" )
set ( TEST_VDE_CHAIN9_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/vde_chain_nm9" )
set ( TEST_JAC_CHAIN2_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm2" )
set ( TEST_JAC_CHAIN3_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm3" )
set ( TEST_JAC_CHAIN4_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm4" )
set ( TEST_JAC_CHAIN5_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm5" )
set ( TEST_JAC_CHAIN6_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm6" )
set ( TEST_JAC_CHAIN7_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm7" )
set ( TEST_JAC_CHAIN8_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm8" )
set ( TEST_JAC_CHAIN9_SRC "${PROJECT_SOURCE_DIR}/test/casadi_chain/jac_chain_nm9" )
set ( TEST_CHAIN_OCP_SRC "${PROJECT_SOURCE_DIR}/test/test_chain" )

# Define external libraries
# BLASFEO
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/blasfeo make static_library )
add_library ( blasfeo STATIC IMPORTED )
set_property ( TARGET blasfeo PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/blasfeo/libblasfeo.a")

# HPMPC
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/hpmpc make static_library )
add_library ( hpmpc STATIC IMPORTED )
set_property ( TARGET hpmpc PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/hpmpc/libhpmpc.a")

# qpOASES
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/qpOASES make )
add_library ( qpoases STATIC IMPORTED )
set_property ( TARGET qpoases PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/qpOASES/bin/libqpOASES_e.a")

# CasADi
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi mkdir build )
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi cd build )
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build cmake .. )
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build make )

# Define acados library
add_library ( acados STATIC ${ACADOS_SRC} ${ACADOS_SRC_SIM} ${ACADOS_SRC_QP} ${ACADOS_SRC_UTILS} )

# Define test
add_executable ( test_hpmpc ${TEST_HPMPC_SRC} )
target_link_libraries ( test_hpmpc acados hpmpc blasfeo m )
add_executable ( test_condensing_qpoases ${TEST_CONDENSING_QPOASES_SRC} )
target_link_libraries ( test_condensing_qpoases acados qpoases hpmpc blasfeo m )
add_executable ( test_nmpc ${TEST_NMPC_SRC} ${TEST_CHEN_SRC} )
target_link_libraries ( test_nmpc acados qpoases blasfeo m )
add_executable ( test_chain ${TEST_CHAIN_OCP_SRC} ${TEST_CHAIN_SRC} ${TEST_VDE_CHAIN2_SRC} ${TEST_VDE_CHAIN3_SRC} ${TEST_VDE_CHAIN4_SRC} ${TEST_VDE_CHAIN5_SRC} ${TEST_VDE_CHAIN6_SRC} ${TEST_VDE_CHAIN7_SRC} ${TEST_VDE_CHAIN8_SRC} ${TEST_VDE_CHAIN9_SRC} ${TEST_JAC_CHAIN2_SRC} ${TEST_JAC_CHAIN3_SRC} ${TEST_JAC_CHAIN4_SRC} ${TEST_JAC_CHAIN5_SRC} ${TEST_JAC_CHAIN6_SRC} ${TEST_JAC_CHAIN7_SRC} ${TEST_JAC_CHAIN8_SRC} ${TEST_JAC_CHAIN9_SRC} )
target_link_libraries ( test_chain acados qpoases blasfeo m )

# Lint acados
find_package ( PythonInterp )
if ( PYTHONINTERP_FOUND )
    set ( ACADOS_STYLE_FILTER
    	"-legal/copyright,"
    	"-readability/casting" )
    set ( FIND_COMMAND find acados test external/casadi/casadi/interfaces/acados -type f -name "*.c" ! -path "test/casadi_chain/*" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" )
    set ( FIND_ACADOS_FILES ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_COMMAND} )
    execute_process ( COMMAND ${FIND_ACADOS_FILES} RESULT_VARIABLE RESULT OUTPUT_VARIABLE ACADOS_FILES_TO_LINT )
    string ( REPLACE "\n" " " ACADOS_FILES_TO_LINT ${ACADOS_FILES_TO_LINT} )
    separate_arguments ( ACADOS_FILES_TO_LINT )
    set ( LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --filter="${ACADOS_STYLE_FILTER}" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 ${ACADOS_FILES_TO_LINT} )
    add_custom_target ( lint ${LINT_COMMAND} )
endif()
