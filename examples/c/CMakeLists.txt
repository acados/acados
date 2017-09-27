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

configure_file(chain_model/chain_model.h.in ${CMAKE_CURRENT_SOURCE_DIR}/chain_model/chain_model.h @ONLY)

# Define sources for models used in the examples
file(GLOB CHAIN_MODEL_SRC
    chain_model/chain_model.c
    chain_model/vde*.c
    chain_model/jac*.c)

set(CHEN_MODEL_SRC Chen_model/chen_model.c)

file(GLOB PENDULUM_MODEL_SRC
    pendulum_model/pendulum_model.c
    pendulum_model/vde_forw_pendulum.c
    pendulum_model/jac_pendulum.c)

file(GLOB QUADCOPTER_MODEL_SRC
    quadcopter_model/quadcopter_model.c
    quadcopter_model/vde_forw_quadcotper.c
    quadcopter_model/ls_res_quadcopter.c
    quadcopter_model/ls_res_end_quadcopter.c)

# Define examples
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
   add_executable(mass_spring_hpmpc_example mass_spring_hpmpc.c)
   target_link_libraries(mass_spring_hpmpc_example acados)
   add_test(mass_spring_hpmpc_example mass_spring_hpmpc_example)

    add_executable(mass_spring_hpipm_example mass_spring_hpipm.c)
    target_link_libraries(mass_spring_hpipm_example acados)
    add_test(mass_spring_hpipm_example mass_spring_hpipm_example)

    add_executable(mass_spring_condensing_hpipm_example mass_spring_condensing_hpipm.c)
    target_link_libraries(mass_spring_condensing_hpipm_example acados)
    add_test(mass_spring_condensing_hpipm_example mass_spring_condensing_hpipm_example)

    # Remove from test because of weird GCC bug (Segfault)
    add_executable(mass_spring_qpoases_example mass_spring_qpoases.c)
    target_link_libraries(mass_spring_qpoases_example acados)
    # add_test(mass_spring_qpoases_example mass_spring_qpoases_example)

    add_executable(mass_spring_qpdunes_example mass_spring_qpdunes.c)
    target_link_libraries(mass_spring_qpdunes_example acados)
    target_compile_definitions(mass_spring_qpdunes_example PRIVATE USE_ACADOS_TYPES)
    add_test(mass_spring_qpdunes_example mass_spring_qpdunes_example)

    if (EXISTS ${PROJECT_SOURCE_DIR}/external/OOQP AND NOT CMAKE_BUILD_TYPE MATCHES NoExternalLibs)
        add_executable(mass_spring_ooqp_example mass_spring_ooqp.c)
        target_link_libraries(mass_spring_ooqp_example acados)
        add_test(mass_spring_ooqp_example mass_spring_ooqp_example)
    endif()

    add_executable(nonlinear_chain_ocp_nlp_example
    nonlinear_chain_ocp_nlp.c
    chain_model/chain_model.c
    chain_model/jac_chain_nm2.c
    chain_model/jac_chain_nm3.c
    chain_model/jac_chain_nm4.c
    chain_model/jac_chain_nm5.c
    chain_model/jac_chain_nm6.c
    chain_model/jac_chain_nm7.c
    chain_model/jac_chain_nm8.c
    chain_model/jac_chain_nm9.c

    chain_model/vde_chain_nm2.c
    chain_model/vde_chain_nm3.c
    chain_model/vde_chain_nm4.c
    chain_model/vde_chain_nm5.c
    chain_model/vde_chain_nm6.c
    chain_model/vde_chain_nm7.c
    chain_model/vde_chain_nm8.c
    chain_model/vde_chain_nm9.c
    )
    target_link_libraries(nonlinear_chain_ocp_nlp_example acados)
    add_test(nonlinear_chain_ocp_nlp_example nonlinear_chain_ocp_nlp_example)

# add_executable(mass_spring_partial_tightening_example
#    mass_spring_partial_tightening.c
# )
# target_link_libraries(mass_spring_partial_tightening_example acados)
# add_test(
#     NAME mass_spring_partial_tightening_example
#     COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./mass_spring_partial_tightening_example
# )

    add_executable(chen_nmpc_qpoases_example
        chen_nmpc_qpoases.c ${CHEN_MODEL_SRC})
    target_link_libraries(chen_nmpc_qpoases_example acados)
    add_test(chen_nmpc_qpoases_example chen_nmpc_qpoases_example)

    add_executable(nonlinear_chain_sqp_example
        nonlinear_chain_sqp.c ${CHAIN_MODEL_SRC})
    target_include_directories(nonlinear_chain_sqp_example PRIVATE "${CMAKE_BINARY_DIR}")
    target_link_libraries(nonlinear_chain_sqp_example acados)
    add_test(nonlinear_chain_sqp_example nonlinear_chain_sqp_example)
endif()

# add_executable(pendulum_hpmpc_libstr_example pendulum_hpmpc_libstr.c ${PENDULUM_MODEL_SRC})
# target_link_libraries(pendulum_hpmpc_libstr_example acados)
# add_test(NAME pendulum_hpmpc_libstr_example COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./pendulum_hpmpc_libstr_example)

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_executable(pendulum_partial_tightening_example
      pendulum_partial_tightening.c ${PENDULUM_MODEL_SRC})
    target_link_libraries(pendulum_partial_tightening_example acados)
    add_test(pendulum_partial_tightening_example pendulum_partial_tightening_example)
endif()

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_executable(quadcopter_ocp_nlp
      quadcopter_ocp_nlp.c ${QUADCOPTER_MODEL_SRC})
    target_link_libraries(quadcopter_ocp_nlp acados)
    add_test(quadcopter_ocp_nlp quadcopter_ocp_nlp)
endif()

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_executable(linear_mpc_example linear_mpc.c)
    target_link_libraries(linear_mpc_example acados)
    target_compile_definitions(linear_mpc_example PRIVATE USE_ACADOS_TYPES)
    add_test(linear_mpc_example linear_mpc_example)
endif()
