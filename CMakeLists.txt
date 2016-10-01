cmake_minimum_required ( VERSION 2.6 )
project ( ACADOS )

set ( CMAKE_VERBOSE_MAKEFILE ON )

set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror" )

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
file (GLOB ACADOS_SRC_SIM "${PROJECT_SOURCE_DIR}/acados/sim/*.c")
file (GLOB ACADOS_SRC_QP "${PROJECT_SOURCE_DIR}/acados/ocp_qp/*.c")
file (GLOB ACADOS_SRC_UTILS "${PROJECT_SOURCE_DIR}/acados/utils/*.c")
set ( TEST_HPMPC_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_hpmpc" )
set ( TEST_CONDENSING_QPOASES_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_condensing_qpoases" )
set ( TEST_NMPC_SRC "${PROJECT_SOURCE_DIR}/test/test_nmpc" )

# Define external libraries
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/hpmpc make )
add_library ( hpmpc STATIC IMPORTED )
set_property ( TARGET hpmpc PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/hpmpc/libhpmpc.a")
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/qpOASES make )
add_library ( qpoases STATIC IMPORTED )
set_property ( TARGET qpoases PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/external/qpOASES/bin/libqpOASES_e.a")
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build make )

# Define acados library
add_library ( acados STATIC ${ACADOS_SRC_SIM} ${ACADOS_SRC_QP} ${ACADOS_SRC_UTILS} )

# Define test
add_executable ( test_hpmpc ${TEST_HPMPC_SRC} )
target_link_libraries ( test_hpmpc acados hpmpc m )
add_executable ( test_condensing_qpoases ${TEST_CONDENSING_QPOASES_SRC} )
target_link_libraries ( test_condensing_qpoases acados qpoases hpmpc m )
add_executable ( test_nmpc ${TEST_NMPC_SRC} )
target_link_libraries ( test_nmpc acados qpoases hpmpc m )

# Lint acados
find_package ( PythonInterp )
if ( PYTHONINTERP_FOUND )
    set ( ACADOS_STYLE_FILTER
    	"-legal/copyright,"
    	"-readability/casting" )
    set ( FIND_COMMAND find acados test external/casadi/casadi/interfaces/acados -type f -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" )
    set ( FIND_ACADOS_FILES ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_COMMAND} )
    execute_process ( COMMAND ${FIND_ACADOS_FILES} RESULT_VARIABLE RESULT OUTPUT_VARIABLE ACADOS_FILES_TO_LINT )
    string ( REPLACE "\n" " " ACADOS_FILES_TO_LINT ${ACADOS_FILES_TO_LINT} )
    separate_arguments ( ACADOS_FILES_TO_LINT )
    set ( LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --filter="${ACADOS_STYLE_FILTER}" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 ${ACADOS_FILES_TO_LINT} )
    add_custom_target ( lint ${LINT_COMMAND} )
endif()
