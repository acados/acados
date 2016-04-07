cmake_minimum_required ( VERSION 2.6 )
project ( ACADOS )

set ( CMAKE_VERBOSE_MAKEFILE ON )
set ( CMAKE_C_FLAGS "-std=c99 -O2 -Wall -Wextra -Werror" )

# Define include folders
include_directories ( "${PROJECT_SOURCE_DIR}" )

# Define .c source files
file (GLOB ACADOS_SRC "${PROJECT_SOURCE_DIR}/acados/*.c")
set ( TEST_SRC "${PROJECT_SOURCE_DIR}/test/test_ocp_qp_solver.c" )

# Make hpmpc
execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/hpmpc make )

# Define libraries
add_library ( acados STATIC ${ACADOS_SRC} )
add_library ( hpmpc STATIC IMPORTED )
set_property ( TARGET hpmpc PROPERTY IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/hpmpc/libhpmpc.a")

# Add test
add_executable ( test_acados ${TEST_SRC} )
target_link_libraries ( test_acados acados hpmpc m )

# Set up cpplint
find_package ( PythonInterp )
if ( PYTHONINTERP_FOUND )
    set ( ACADOS_STYLE_FILTER
    	"-build/header_guard,"
    	"-build/include,"
    	"-legal/copyright,"
    	"-readability/casting" )
    set ( FIND_ACADOS find acados -type f -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp" )
    set ( FIND_LINT_FILES ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_ACADOS} )
    execute_process ( COMMAND ${FIND_LINT_FILES} RESULT_VARIABLE RESULT OUTPUT_VARIABLE ACADOS_LINT_SRC )
    string ( REPLACE "\n" " " ACADOS_LINT_SRC ${ACADOS_LINT_SRC} )
    separate_arguments ( ACADOS_LINT_SRC )
    set ( RUN_LINT ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE} ./cpplint.py --filter="${ACADOS_STYLE_FILTER}" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 ${ACADOS_LINT_SRC} )
    add_custom_target ( lint ${RUN_LINT} )
endif()
