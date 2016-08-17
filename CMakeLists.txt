cmake_minimum_required ( VERSION 2.6 )
project ( ACADOS )

set ( CMAKE_VERBOSE_MAKEFILE ON )
set ( CMAKE_C_FLAGS "-std=c99 -g -Wall -Wextra -Werror" )

# Define include folders
include_directories ( "${PROJECT_SOURCE_DIR}" )
include_directories ( "${PROJECT_SOURCE_DIR}/qpOASES/include/" )

# Define .c source files
file (GLOB ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/*.c")
set ( TEST_HPMPC_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_hpmpc" )
set ( TEST_CONDENSING_QPOASES_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_condensing_qpoases" )

# Make libraries
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/hpmpc make )
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/qpOASES make )

# Define libraries
add_library ( acados STATIC ${ACADOS_SRC} )
add_library ( hpmpc STATIC IMPORTED )
set_property ( TARGET hpmpc PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/hpmpc/libhpmpc.a")
add_library ( qpoases STATIC IMPORTED )
set_property ( TARGET qpoases PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/qpOASES/bin/libqpOASES_e.a")

# Define test
add_executable ( test_hpmpc ${TEST_HPMPC_SRC} )
add_executable ( test_condensing_qpoases ${TEST_CONDENSING_QPOASES_SRC} )
target_link_libraries ( test_hpmpc acados hpmpc m )
target_link_libraries ( test_condensing_qpoases acados qpoases hpmpc m )

# Lint acados
find_package ( PythonInterp )
if ( PYTHONINTERP_FOUND )
    set ( ACADOS_STYLE_FILTER
    	"-build/include,"
    	"-legal/copyright,"
    	"-readability/casting" )
    set ( FIND_COMMAND find acados -type f -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" )
    set ( FIND_ACADOS_FILES ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_COMMAND} )
    execute_process ( COMMAND ${FIND_ACADOS_FILES} RESULT_VARIABLE RESULT OUTPUT_VARIABLE ACADOS_FILES_TO_LINT )
    string ( REPLACE "\n" " " ACADOS_FILES_TO_LINT ${ACADOS_FILES_TO_LINT} )
    separate_arguments ( ACADOS_FILES_TO_LINT )
    set ( LINT_COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --filter="${ACADOS_STYLE_FILTER}" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 ${ACADOS_FILES_TO_LINT} )
    add_custom_target ( lint ${LINT_COMMAND} )
endif()
