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
    chain_model/jac_chain_nm2.c
    chain_model/jac_chain_nm3.c
    chain_model/jac_chain_nm4.c
    chain_model/vde_chain_nm2.c
    chain_model/vde_chain_nm3.c
    chain_model/vde_chain_nm4.c
    chain_model/ls_cost_nm2.c
    chain_model/ls_cost_nm3.c
    chain_model/ls_cost_nm4.c
    chain_model/ls_costN_nm2.c
    chain_model/ls_costN_nm3.c
    chain_model/ls_costN_nm4.c
    chain_model/pathcon_nm2.c
    chain_model/pathcon_nm3.c
    chain_model/pathcon_nm4.c
    chain_model/pathconN_nm2.c
    chain_model/pathconN_nm3.c
    chain_model/pathconN_nm4.c
    chain_model/vde_hess_chain_nm2.c
    chain_model/vde_hess_chain_nm3.c
    chain_model/vde_hess_chain_nm4.c
    implicit_chain_model/impl_jacFun_u_chain_nm2.c
    implicit_chain_model/impl_jacFun_u_chain_nm3.c
    implicit_chain_model/impl_jacFun_u_chain_nm4.c
    implicit_chain_model/impl_jacFun_x_chain_nm2.c
    implicit_chain_model/impl_jacFun_x_chain_nm3.c
    implicit_chain_model/impl_jacFun_x_chain_nm4.c
    implicit_chain_model/impl_jacFun_xdot_chain_nm2.c
    implicit_chain_model/impl_jacFun_xdot_chain_nm3.c
    implicit_chain_model/impl_jacFun_xdot_chain_nm4.c
    implicit_chain_model/impl_odeFun_chain_nm2.c
    implicit_chain_model/impl_odeFun_chain_nm3.c
    implicit_chain_model/impl_odeFun_chain_nm4.c
)

set(CHEN_MODEL_SRC chen_model/chen_model.c)

file(GLOB PENDULUM_MODEL_SRC
    pendulum_model/pendulum_model.c
    pendulum_model/vde_forw_pendulum.c
    pendulum_model/jac_pendulum.c
)

file(GLOB CRANE_MODEL_SRC
    crane_model/impl_jac_u.c
    crane_model/impl_jac_x.c
    crane_model/impl_jac_xdot.c
    crane_model/impl_ode.c
    crane_model/jac_model.c
    crane_model/ode_model.c
    crane_model/vde_adj_model.c
    crane_model/vde_forw_model.c
    crane_model/vde_hess_model.c
)

# Define examples
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

    # -------------------- mass_spring

    add_executable(mass_spring_example mass_spring_example.c mass_spring_model/mass_spring_qp.c)
    target_link_libraries(mass_spring_example acados)
    add_test(mass_spring_example mass_spring_example)

    if (1) #(NOT ACADOS_WITH_C_INTERFACE)
        add_executable(
            mass_spring_example_no_interface
            mass_spring_example_no_interface.c
            mass_spring_model/mass_spring_qp.c
        )
        target_link_libraries(
            mass_spring_example_no_interface
            acados
        )
        add_test(
            mass_spring_example_no_interface
            mass_spring_example_no_interface
        )
    endif()

    add_executable(mass_spring_fcond_split mass_spring_fcond_split.c mass_spring_model/mass_spring_qp.c)
    target_link_libraries(mass_spring_fcond_split acados)
    add_test(mass_spring_fcond_split mass_spring_example)

    add_executable(mass_spring_pcond_split mass_spring_pcond_split.c mass_spring_model/mass_spring_qp.c)
    target_link_libraries(mass_spring_pcond_split acados)
    add_test(mass_spring_pcond_split mass_spring_example)

    add_executable(mass_spring_offline_fcond_qpoases_split mass_spring_offline_fcond_qpoases_split.c mass_spring_model/mass_spring_qp.c)
    target_link_libraries(mass_spring_offline_fcond_qpoases_split acados)
    add_test(mass_spring_offline_fcond_qpoases_split mass_spring_example)

    # add_executable(mass_spring_partial_tightening_example
    #   mass_spring_partial_tightening.c
    # )
    # target_link_libraries(mass_spring_partial_tightening_example acados)
    # add_test(
    #   NAME mass_spring_partial_tightening_example
    #   COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH}
    #   ./mass_spring_partial_tightening_example
    # )

    # -------------------- nonlinear_chain

    if(ACADOS_WITH_C_INTERFACE)
        # add_executable(nonlinear_chain_ocp_nlp_example nonlinear_chain_ocp_nlp.c ${CHAIN_MODEL_SRC})
        # target_link_libraries(nonlinear_chain_ocp_nlp_example acados)
        # add_test(nonlinear_chain_ocp_nlp_example nonlinear_chain_ocp_nlp_example)
    else()
        add_executable(
            nonlinear_chain_ocp_nlp_no_interface_example
            nonlinear_chain_ocp_nlp_no_interface.c
            ${CHAIN_MODEL_SRC}
        )
        target_link_libraries(
            nonlinear_chain_ocp_nlp_no_interface_example
            acados
        )
        add_test(
            nonlinear_chain_ocp_nlp_no_interface_example
            nonlinear_chain_ocp_nlp_no_interface_example
        )
    endif()

    # -------------------- sim_crane

    add_executable(sim_crane_example sim_crane_example.c ${CRANE_MODEL_SRC})
    target_link_libraries(sim_crane_example acados)
    add_test(sim_crane_example sim_crane_no_interface_example)

    if(ACADOS_WITH_C_INTERFACE)
    else()
        add_executable(
            sim_crane_no_interface_example
            sim_crane_no_interface.c
            ${CRANE_MODEL_SRC}
        )
        target_link_libraries(
            sim_crane_no_interface_example
            acados
        )
        add_test(
            sim_crane_no_interface_example
            sim_crane_no_interface_example
        )
    endif()

    # -------------------- ocp_qp

    add_executable(ocp_qp ocp_qp.c)
    target_link_libraries(ocp_qp acados)
    add_test(ocp_qp ocp_qp)

    # -------------------- dense_qp

    add_executable(dense_qp dense_qp.c)
    target_link_libraries(dense_qp acados)
    add_test(dense_qp dense_qp)

    # -------------------- chen

    # add_executable(chen_nmpc_qpoases_example
    #   chen_nmpc_qpoases.c ${CHEN_MODEL_SRC})
    #   target_link_libraries(chen_nmpc_qpoases_example acados)
    #   add_test(chen_nmpc_qpoases_example chen_nmpc_qpoases_example)
    #   install(TARGETS chen_nmpc_qpoases_example RUNTIME DESTINATION bin)
    #   set_target_properties(chen_nmpc_qpoases_example PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib/")

    # -------------------- pendulum

    # add_executable(pendulum_partial_tightening_example
    #   pendulum_partial_tightening.c ${PENDULUM_MODEL_SRC})
    # target_link_libraries(pendulum_partial_tightening_example acados)
    # add_test(pendulum_partial_tightening_example pendulum_partial_tightening_example)

    # # add_executable(linear_mpc_example linear_mpc.c)
    # # target_link_libraries(linear_mpc_example acados)
    # # target_compile_definitions(linear_mpc_example PRIVATE USE_ACADOS_TYPES)
    # # add_test(linear_mpc_example linear_mpc_example)

    # add_executable(discrete_model
    #   discrete_model.c
    #   discrete_model/discrete_model.c
    #   discrete_model/discrete_model_cost.c
    #   discrete_model/discrete_model_costN.c
    # )
    # target_link_libraries(discrete_model acados)
    # add_test(discrete_model discrete_model)

    # add_executable(pendulum_hpmpc_libstr_example pendulum_hpmpc_libstr.c ${PENDULUM_MODEL_SRC})
    # target_link_libraries(pendulum_hpmpc_libstr_example acados)
    # add_test(NAME pendulum_hpmpc_libstr_example
    #   COMMAND "${CMAKE_COMMAND}" -E chdir ${EXECUTABLE_OUTPUT_PATH} ./pendulum_hpmpc_libstr_example)

endif()
