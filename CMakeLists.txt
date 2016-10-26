cmake_minimum_required ( VERSION 2.6 )
project ( ACADOS )

set ( CMAKE_INSTALL_PREFIX "$ENV{HOME}/local")
set ( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set ( CMAKE_VERBOSE_MAKEFILE ON )
set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror" )

if ( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Release )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops")
endif ( NOT CMAKE_BUILD_TYPE )

if ( CMAKE_BUILD_TYPE MATCHES Debug )
    add_definitions ( -DDEBUG )
endif ( CMAKE_BUILD_TYPE MATCHES Debug )

if ( CMAKE_BUILD_TYPE MATCHES Profile )
    set( CMAKE_BUILD_TYPE Release )
	set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fprofile-use" )
endif ( CMAKE_BUILD_TYPE MATCHES Profile )

# Define include folders
include_directories ( "${PROJECT_SOURCE_DIR}" )
include_directories ( "${PROJECT_SOURCE_DIR}/external/" )
include_directories ( "${PROJECT_SOURCE_DIR}/external/qpOASES/include/" )

# Define .c source files
file (GLOB ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/*.c")
file (GLOB ACADOS_SRC_SIM "${PROJECT_SOURCE_DIR}/acados/sim/*.c")
file (GLOB ACADOS_SRC_QP "${PROJECT_SOURCE_DIR}/acados/ocp_qp/*.c")
file (GLOB ACADOS_SRC_UTILS "${PROJECT_SOURCE_DIR}/acados/utils/*.c")
set ( TEST_HPMPC_SRC "${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_hpmpc" )
set ( TEST_CONDENSING_QPOASES_SRC "${PROJECT_SOURCE_DIR}/examples/test_ocp_qp_condensing_qpoases" )
set ( TEST_CHEN_SRC "${PROJECT_SOURCE_DIR}/examples/Chen/Chen_model" )
set ( TEST_NMPC_SRC "${PROJECT_SOURCE_DIR}/examples/test_nmpc" )
set ( TEST_CHAIN_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/Chain_model" )
set ( TEST_VDE_CHAIN2_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm2" )
set ( TEST_VDE_CHAIN3_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm3" )
set ( TEST_VDE_CHAIN4_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm4" )
set ( TEST_VDE_CHAIN5_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm5" )
set ( TEST_VDE_CHAIN6_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm6" )
set ( TEST_VDE_CHAIN7_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm7" )
set ( TEST_VDE_CHAIN8_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm8" )
set ( TEST_VDE_CHAIN9_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/vde_chain_nm9" )
set ( TEST_JAC_CHAIN2_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm2" )
set ( TEST_JAC_CHAIN3_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm3" )
set ( TEST_JAC_CHAIN4_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm4" )
set ( TEST_JAC_CHAIN5_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm5" )
set ( TEST_JAC_CHAIN6_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm6" )
set ( TEST_JAC_CHAIN7_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm7" )
set ( TEST_JAC_CHAIN8_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm8" )
set ( TEST_JAC_CHAIN9_SRC "${PROJECT_SOURCE_DIR}/examples/casadi_chain/jac_chain_nm9" )
set ( TEST_CHAIN_OCP_SRC "${PROJECT_SOURCE_DIR}/examples/test_chain" )
set ( UNIT_TESTS_SRC "${PROJECT_SOURCE_DIR}/test/all_tests" )

# Define external libraries
# BLASFEO
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/blasfeo
    make static_library RESULT_VARIABLE blasfeo_built )
if ( NOT "${blasfeo_built}" STREQUAL "0")
    message( FATAL_ERROR "build of blasfeo failed with status ${blasfeo_built}!" )
endif ( )
add_library ( blasfeo STATIC IMPORTED )
set_property ( TARGET blasfeo PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/blasfeo/libblasfeo.a")

# HPMPC
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/hpmpc
    make static_library RESULT_VARIABLE hpmpc_built )
if ( NOT "${hpmpc_built}" STREQUAL "0")
    message( FATAL_ERROR "build of hpmpc failed with status ${hpmpc_built}!" )
endif ( )
add_library ( hpmpc STATIC IMPORTED )
set_property ( TARGET hpmpc PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/hpmpc/libhpmpc.a")

# qpOASES
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/qpOASES
    make RESULT_VARIABLE qpoases_built )
if ( NOT "${qpoases_built}" STREQUAL "0")
    message( FATAL_ERROR "build of qpoases failed with status ${qpoases_built}!" )
endif ( )
add_library ( qpoases STATIC IMPORTED )
set_property ( TARGET qpoases PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/qpOASES/bin/libqpOASES_e.a")

# CasADi
# execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi mkdir build )
# execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi cd build )
# execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build cmake .. )
# execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build
    # make RESULT_VARIABLE casadi_built )
# if ( NOT "${casadi_built}" STREQUAL "0")
    # message( FATAL_ERROR "build of casadi failed with status ${casadi_built}!" )
# endif ( )

# eigen
find_package ( Eigen3 )
if ( NOT EIGEN3_FOUND )
	execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/eigen mkdir build)
	execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/eigen/build cmake ..)
	execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/eigen/build make install)
endif ( )

# catch
add_library ( catch STATIC IMPORTED )
set_property ( TARGET catch PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/catch" )

# Define acados library
add_library ( acados STATIC ${ACADOS_SRC} ${ACADOS_SRC_SIM} ${ACADOS_SRC_QP} ${ACADOS_SRC_UTILS} )

# Define tests
add_executable ( test_hpmpc ${TEST_HPMPC_SRC} )
target_link_libraries ( test_hpmpc acados hpmpc blasfeo m )

add_executable ( test_condensing_qpoases ${TEST_CONDENSING_QPOASES_SRC} )
target_link_libraries ( test_condensing_qpoases acados qpoases hpmpc blasfeo m )

add_executable ( test_nmpc ${TEST_NMPC_SRC} ${TEST_CHEN_SRC} )
target_link_libraries ( test_nmpc acados qpoases blasfeo m )

add_executable ( test_chain ${TEST_CHAIN_OCP_SRC} ${TEST_CHAIN_SRC} ${TEST_VDE_CHAIN2_SRC} ${TEST_VDE_CHAIN3_SRC} ${TEST_VDE_CHAIN4_SRC} ${TEST_VDE_CHAIN5_SRC} ${TEST_VDE_CHAIN6_SRC} ${TEST_VDE_CHAIN7_SRC} ${TEST_VDE_CHAIN8_SRC} ${TEST_VDE_CHAIN9_SRC} ${TEST_JAC_CHAIN2_SRC} ${TEST_JAC_CHAIN3_SRC} ${TEST_JAC_CHAIN4_SRC} ${TEST_JAC_CHAIN5_SRC} ${TEST_JAC_CHAIN6_SRC} ${TEST_JAC_CHAIN7_SRC} ${TEST_JAC_CHAIN8_SRC} ${TEST_JAC_CHAIN9_SRC} )
target_link_libraries ( test_chain acados qpoases blasfeo m )

add_executable ( unit_tests ${UNIT_TESTS_SRC} )
target_link_libraries ( unit_tests acados qpoases blasfeo m )

# Lint acados
find_package ( PythonInterp )
if ( PYTHONINTERP_FOUND )
    set ( ACADOS_STYLE_FILTER
    	"-legal/copyright,"
    	"-readability/casting" )
    set ( FIND_COMMAND find acados test examples external/casadi/casadi/interfaces/acados -type f -name "*.c" ! -path "examples/casadi_chain/*" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" )
    set ( FIND_ACADOS_FILES ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_COMMAND} )
    execute_process ( COMMAND ${FIND_ACADOS_FILES} RESULT_VARIABLE RESULT OUTPUT_VARIABLE ACADOS_FILES_TO_LINT )
    string ( REPLACE "\n" " " ACADOS_FILES_TO_LINT ${ACADOS_FILES_TO_LINT} )
    separate_arguments ( ACADOS_FILES_TO_LINT )
    set ( LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --filter="${ACADOS_STYLE_FILTER}" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 ${ACADOS_FILES_TO_LINT} )
    add_custom_target ( lint ${LINT_COMMAND} )
endif ( )
