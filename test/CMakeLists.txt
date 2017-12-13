#
#    This file is part of acados.
#
#    acados is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 3 of the License, or (at your option) any later version.
#
#    acados is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with acados; if not, write to the Free Software Foundation,
#    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#

# Check if external libraries are present; these are needed for the unit tests
find_package(Eigen3)
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/casadi-octave-v3.2.3")
    message(WARNING "You chose to run unit tests but casadi-octave-v3.2.3 was not found in: ${PROJECT_SOURCE_DIR}/external/casadi-octave-v3.2.3")
endif()
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(ocp_nlp)
add_subdirectory(ocp_qp)
add_subdirectory(sim)

set(TEST_OCP_NLP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp/chain/test_chain.cpp
)

set(TEST_OCP_QP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/test_qpsolvers.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/../test_utils/read_ocp_qp_in.c
)

set(TEST_SIM_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/sim/pendulum/test_pendulum.cpp
)

set(TEST_UTILS_SRC
    ${PROJECT_SOURCE_DIR}/test/test_utils/read_matrix.cpp
    ${PROJECT_SOURCE_DIR}/test/test_utils/zeros.cpp
) # TODO (dimitris): eliminate this at some point

# Unit test executable
add_executable(unit_tests
    ${PROJECT_SOURCE_DIR}/test/all_tests.cpp
    $<TARGET_OBJECTS:ocp_nlp_gen>
    ${TEST_OCP_NLP_SRC}
    $<TARGET_OBJECTS:ocp_qp_gen>
    ${TEST_OCP_QP_SRC}
    $<TARGET_OBJECTS:sim_gen>
    ${TEST_SIM_SRC}
    ${TEST_UTILS_SRC}
)

target_include_directories(unit_tests PRIVATE "${EXTERNAL_SRC_DIR}/eigen")
target_link_libraries(unit_tests acados)
if(ACADOS_WITH_OOQP)
    target_compile_definitions(unit_tests PRIVATE OOQP)
endif()

# Add as test in ctest
add_test(NAME unit_tests COMMAND "${CMAKE_COMMAND}" -E chdir ${CMAKE_BINARY_DIR}/test ./unit_tests -a)

file(COPY "${PROJECT_SOURCE_DIR}/acados/sim/simplified/" DESTINATION "${PROJECT_BINARY_DIR}/test/simplified/")
