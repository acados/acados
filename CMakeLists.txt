cmake_minimum_required ( VERSION 2.6 )
project ( ACADOS )

set ( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")
set ( CMAKE_C_COMPILER "gcc-6" )

set ( CMAKE_VERBOSE_MAKEFILE ON )
set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -v" )
set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -Wextra -Werror -Wno-unused-function" )
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wno-unused-function" )

if ( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Release )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops")
endif ( )

message ( "\n\n --- The CMake build type is ${CMAKE_BUILD_TYPE} ---\n\n" )

if ( CMAKE_BUILD_TYPE MATCHES Debug )
    add_definitions ( -DDEBUG )
endif ( )

if ( CMAKE_BUILD_TYPE MATCHES Profile )
	set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -fprofile-use" )
endif ( )

if ( CMAKE_BUILD_TYPE MATCHES Test )
	set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g" )
	set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
endif ( )

# Define include folders
include_directories ( "${PROJECT_SOURCE_DIR}" )
include_directories ( "${PROJECT_SOURCE_DIR}/external/" )
include_directories ( "${PROJECT_SOURCE_DIR}/external/qpOASES/include/" )
if ( CMAKE_BUILD_TYPE MATCHES Test )
    if ( NOT EXISTS "${PROJECT_SOURCE_DIR}/external/eigen/Eigen" )
        message ( FATAL_ERROR "The build type is ${CMAKE_BUILD_TYPE} (Test) but Eigen can not be found" )
    endif ( )
    include_directories ( "${PROJECT_SOURCE_DIR}/external/eigen/" )
endif ( )

# Define source files
include ( acados_sources )
if ( CMAKE_BUILD_TYPE MATCHES Test )
    file ( GLOB UNIT_TESTS_SRC_TEST_UTILS ${PROJECT_SOURCE_DIR}/test/test_utils/*.cpp )
    file ( GLOB UNIT_TESTS_SRC_OCP_QP_HELPER ${PROJECT_SOURCE_DIR}/test/ocp_qp/*.cpp )
    set ( UNIT_TESTS_SRC "${PROJECT_SOURCE_DIR}/test/all_tests" )
endif ( )

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
if ( BUILD_WITH_CASADI )
	execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi mkdir build )
	execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi cd build )
	execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build cmake .. )
	execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/external/casadi/build
	    make RESULT_VARIABLE casadi_built )
	if ( NOT "${casadi_built}" STREQUAL "0")
	    message( FATAL_ERROR "build of casadi failed with status ${casadi_built}!" )
	endif ( )
endif ( )

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

if ( CMAKE_BUILD_TYPE MATCHES Test )
    # Find test data generation files
	set ( FIND_COMMAND find test -type f -name "generate_*.m")
	set ( FIND_GENERATION_FILES ${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR} ${FIND_COMMAND}  )
	execute_process ( COMMAND ${FIND_GENERATION_FILES} RESULT_VARIABLE RESULT OUTPUT_VARIABLE TEST_DATA_GENERATION_FILES )
	string ( REPLACE "\n" " " TEST_DATA_GENERATION_FILES ${TEST_DATA_GENERATION_FILES} )
	separate_arguments ( TEST_DATA_GENERATION_FILES )
	# Generate the test data
	foreach ( GENERATION_FILE IN ITEMS ${TEST_DATA_GENERATION_FILES} )
		set ( GENERATION_FILE "${PROJECT_SOURCE_DIR}/${GENERATION_FILE}" )
		message ( "Current file: ${GENERATION_FILE}" )
		if ( NOT EXISTS ${PROJECT_BINARY_DIR}/test )
			execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR} mkdir test )
		endif ( )
		set ( GENERATE_TEST_DATA octave --no-gui --path "${PROJECT_SOURCE_DIR}/test/ocp_qp/octave" ${GENERATION_FILE} )
		execute_process ( COMMAND ${CMAKE_COMMAND} -E chdir ${PROJECT_BINARY_DIR}/test ${GENERATE_TEST_DATA} )
	endforeach ( )
    add_executable ( unit_tests ${UNIT_TESTS_SRC} ${UNIT_TESTS_SRC_OCP_QP_HELPER} ${UNIT_TESTS_SRC_TEST_UTILS} )
    target_link_libraries ( unit_tests acados qpoases blasfeo hpmpc m )
	set_target_properties ( unit_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/test )
endif ( )

# Linter
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
