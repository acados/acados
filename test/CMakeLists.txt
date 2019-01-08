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

include_directories(${EIGEN3_INCLUDE_DIR})

# add_subdirectory(ocp_nlp)
add_subdirectory(ocp_qp)
add_subdirectory(sim)

set(TEST_OCP_NLP_SRC
    ######################
    ### CHAIN EXAMPLE ####
    ######################

    ### explicit ###
    # 2 masses
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_chain_nm2.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_adj_chain_nm2.c
    #chain_model/vde_hess_chain_nm2.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_cost_nm2.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_costN_nm2.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/casadi_erk4_chain_nm2.c
    # 3 masses
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_chain_nm3.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_adj_chain_nm3.c
    #chain_model/vde_hess_chain_nm3.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_cost_nm3.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_costN_nm3.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/casadi_erk4_chain_nm3.c
    # 4 masses
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_chain_nm4.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_adj_chain_nm4.c
    #chain_model/vde_hess_chain_nm4.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_cost_nm4.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_costN_nm4.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/casadi_erk4_chain_nm4.c
    # 5 masses
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_chain_nm5.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_adj_chain_nm5.c
    #chain_model/vde_hess_chain_nm5.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_cost_nm5.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_costN_nm5.c
    #${CMAKE_SOURCE_DIR}/examples/c/chain_model/casadi_erk4_chain_nm5.c
    # 6 masses
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_chain_nm6.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/vde_adj_chain_nm6.c
    #chain_model/vde_hess_chain_nm6.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_cost_nm6.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/ls_costN_nm6.c
    #${CMAKE_SOURCE_DIR}/examples/c/chain_model/casadi_erk4_chain_nm6.c

    # external cost
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/chain_nm_2_external_cost.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/chain_nm_3_external_cost.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/chain_nm_4_external_cost.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/chain_nm_5_external_cost.c
    ${CMAKE_SOURCE_DIR}/examples/c/chain_model/chain_nm_6_external_cost.c

    ### implicit ###
    # 2 masses
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_chain_nm2.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm2.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm2.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm2.c
    #${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_u_chain_nm2.c
    # 3 masses
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_chain_nm3.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm3.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm3.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm3.c
    #${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_u_chain_nm3.c
    # 4 masses
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_chain_nm4.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm4.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm4.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm4.c
    #${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_u_chain_nm4.c
    # 5 masses
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_chain_nm5.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm5.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm5.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm5.c
    #${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_u_chain_nm5.c
    # 6 masses
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_chain_nm6.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm6.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm6.c
    ${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm6.c
    #${CMAKE_SOURCE_DIR}/examples/c/implicit_chain_model/impl_ode_jac_x_u_chain_nm6.c

    ############################
    ### WIND TURBINE EXAMPLE ###
    ############################

    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_expl_ode_fun.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_expl_vde_for.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun_jac_x_xdot.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_jac_x_xdot_u.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun_jac_x_xdot_u.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_phi_fun.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_phi_fun_jac_y.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_phi_jac_y_uhat.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_f_lo_fun_jac_x1k1uz.c
    ${CMAKE_SOURCE_DIR}/examples/c/wt_model_nx6/nx6p2/wt_nx6p2_get_matrices_fun.c

    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp/test_chain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_nlp/test_wind_turbine.cpp
)

set(TEST_OCP_QP_SRC
    ${PROJECT_SOURCE_DIR}/examples/c/mass_spring_model/mass_spring_qp.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/test_qpsolvers.cpp
    # ${CMAKE_CURRENT_SOURCE_DIR}/ocp_qp/../test_utils/read_ocp_qp_in.c
)


set(TEST_SIM_ODE_SRC
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/expl_ode_fun.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/expl_vde_for.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/expl_vde_adj.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/impl_ode_fun.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/impl_ode_fun_jac_x_xdot.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/impl_ode_jac_x_xdot_u.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/impl_ode_fun_jac_x_xdot_u.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/phi_fun.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/phi_fun_jac_y.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/phi_jac_y_uhat.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/f_lo_fun_jac_x1k1uz.c
    ${PROJECT_SOURCE_DIR}/examples/c/wt_model_nx3/get_matrices_fun.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sim/sim_test_ode.cpp
)

set(TEST_SIM_DAE_SRC
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_impl_ode_fun.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_impl_ode_fun_jac_x_xdot.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_impl_ode_jac_x_xdot_u.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_impl_ode_fun_jac_x_xdot_u.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_phi_fun.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_phi_fun_jac_y.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_phi_jac_y_uhat.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_f_lo_fun_jac_x1k1uz.c
    ${PROJECT_SOURCE_DIR}/examples/c/crane_dae_model/crane_dae_get_matrices_fun.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sim/sim_test_dae.cpp
)

set(TEST_SIM_HESS_SRC
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_expl_ode_fun.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_expl_vde_forw.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_expl_vde_adj.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_expl_ode_hess.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_impl_ode_fun.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_impl_ode_fun_jac_x_xdot.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_impl_ode_jac_x_xdot_u.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_impl_ode_fun_jac_x_xdot_u.c
    ${PROJECT_SOURCE_DIR}/examples/c/pendulum_model/pendulum_ode_impl_ode_hess.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sim/sim_test_hessian.cpp
)


# Unit test executable
add_executable(unit_tests
    ${PROJECT_SOURCE_DIR}/test/all_tests.cpp
    # $<TARGET_OBJECTS:ocp_nlp_gen>
    # $<TARGET_OBJECTS:ocp_qp_gen>
    ${TEST_SIM_HESS_SRC}
    # ${TEST_SIM_DAE_SRC}
    ${TEST_SIM_ODE_SRC}
    ${TEST_OCP_QP_SRC}
    ${TEST_OCP_NLP_SRC}
    # $<TARGET_OBJECTS:sim_gen>
    # ${TEST_UTILS_SRC}
)

target_include_directories(unit_tests PRIVATE "${EXTERNAL_SRC_DIR}/eigen")
target_link_libraries(unit_tests acados)

# if(ACADOS_WITH_OOQP)
#     target_compile_definitions(unit_tests PRIVATE OOQP)
# endif()

# Add as test in ctest
add_test(NAME unit_tests COMMAND "${CMAKE_COMMAND}" -E chdir ${CMAKE_BINARY_DIR}/test ./unit_tests -a)

file(COPY "${PROJECT_SOURCE_DIR}/acados/sim/simplified/" DESTINATION "${PROJECT_BINARY_DIR}/test/simplified/")
