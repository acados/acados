#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#


if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")

configure_file(chain_model/chain_model.h.in ${CMAKE_CURRENT_SOURCE_DIR}/chain_model/chain_model.h @ONLY)

# Define sources for models used in the examples
file(GLOB CHAIN_MODEL_SRC
    ### explicit ###
    # 2 masses
    chain_model/vde_chain_nm2.c
    chain_model/vde_adj_chain_nm2.c
    #chain_model/vde_hess_chain_nm2.c
    chain_model/ls_cost_nm2.c
    chain_model/ls_costN_nm2.c
    chain_model/casadi_erk4_chain_nm2.c
    # 3 masses
    chain_model/vde_chain_nm3.c
    chain_model/vde_adj_chain_nm3.c
    #chain_model/vde_hess_chain_nm3.c
    chain_model/ls_cost_nm3.c
    chain_model/ls_costN_nm3.c
    chain_model/casadi_erk4_chain_nm3.c
    # 4 masses
    chain_model/vde_chain_nm4.c
    chain_model/vde_adj_chain_nm4.c
    #chain_model/vde_hess_chain_nm4.c
    chain_model/ls_cost_nm4.c
    chain_model/ls_costN_nm4.c
    chain_model/casadi_erk4_chain_nm4.c
    # 5 masses
    chain_model/vde_chain_nm5.c
    chain_model/vde_adj_chain_nm5.c
    #chain_model/vde_hess_chain_nm5.c
    chain_model/ls_cost_nm5.c
    chain_model/ls_costN_nm5.c
    # chain_model/casadi_erk4_chain_nm5.c
    # NOTE: chain_model/casadi_erk4_chain_nm5.c (and 6) commented, because compilation takes quite a while
    # 6 masses
    # chain_model/vde_chain_nm6.c
    # chain_model/vde_adj_chain_nm6.c
    # chain_model/vde_hess_chain_nm6.c
    # chain_model/ls_cost_nm6.c
    # chain_model/ls_costN_nm6.c
    # chain_model/casadi_erk4_chain_nm6.c

    # external cost
    chain_model/chain_nm_2_external_cost.c
    chain_model/chain_nm_3_external_cost.c
    chain_model/chain_nm_4_external_cost.c
    chain_model/chain_nm_5_external_cost.c
    chain_model/chain_nm_6_external_cost.c

    ### implicit ###
    # 2 masses
    implicit_chain_model/impl_ode_fun_chain_nm2.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm2.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm2.c
    implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm2.c
    implicit_chain_model/impl_ode_jac_x_u_chain_nm2.c
    # 3 masses
    implicit_chain_model/impl_ode_fun_chain_nm3.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm3.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm3.c
    implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm3.c
    implicit_chain_model/impl_ode_jac_x_u_chain_nm3.c
    # 4 masses
    implicit_chain_model/impl_ode_fun_chain_nm4.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm4.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm4.c
    implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm4.c
    implicit_chain_model/impl_ode_jac_x_u_chain_nm4.c
    # 5 masses
    implicit_chain_model/impl_ode_fun_chain_nm5.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm5.c
    implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm5.c
    implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm5.c
    implicit_chain_model/impl_ode_jac_x_u_chain_nm5.c
    # 6 masses
    # implicit_chain_model/impl_ode_fun_chain_nm6.c
    # implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm6.c
    # implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm6.c
    # implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm6.c
    # implicit_chain_model/impl_ode_jac_x_u_chain_nm6.c

)

file(GLOB CRANE_MODEL_SRC
    crane_model/impl_ode_fun.c
    crane_model/impl_ode_fun_jac_x_xdot.c
    crane_model/impl_ode_jac_x_xdot_u.c
    crane_model/jac_model.c
    crane_model/ode_model.c
    crane_model/vde_adj_model.c
    crane_model/vde_forw_model.c
    crane_model/vde_hess_model.c
)

file(GLOB WT_MODEL_NX3_SRC
    wt_model_nx3/expl_ode_fun.c
    wt_model_nx3/expl_vde_for.c
    wt_model_nx3/expl_vde_adj.c
    wt_model_nx3/impl_ode_fun.c
    wt_model_nx3/impl_ode_fun_jac_x_xdot.c
    wt_model_nx3/impl_ode_jac_x_xdot_u.c
    wt_model_nx3/impl_ode_fun_jac_x_xdot_u.c
    wt_model_nx3/phi_fun.c
    wt_model_nx3/phi_fun_jac_y.c
    wt_model_nx3/phi_jac_y_uhat.c
    wt_model_nx3/f_lo_fun_jac_x1k1uz.c
    wt_model_nx3/get_matrices_fun.c
)

file(GLOB WT_MODEL_NX6_SRC
    wt_model_nx6/expl_ode_fun.c
    wt_model_nx6/expl_vde_for.c
    wt_model_nx6/expl_vde_adj.c
    wt_model_nx6/impl_ode_fun.c
    wt_model_nx6/impl_ode_fun_jac_x_xdot.c
    wt_model_nx6/impl_ode_jac_x_xdot_u.c
    wt_model_nx6/impl_ode_fun_jac_x_xdot_u.c
    # gnsf functions
    wt_model_nx6/phi_fun.c
    wt_model_nx6/phi_fun_jac_y.c
    wt_model_nx6/phi_jac_y_uhat.c
    wt_model_nx6/f_lo_fun_jac_x1k1uz.c
    wt_model_nx6/get_matrices_fun.c

)

file(GLOB WT_MODEL_NX6P2_SRC
    wt_model_nx6/nx6p2/wt_nx6p2_expl_ode_fun.c
    wt_model_nx6/nx6p2/wt_nx6p2_expl_vde_for.c
    wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun.c
    wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun_jac_x_xdot.c
    wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_jac_x_xdot_u.c
    wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun_jac_x_xdot_u.c
    wt_model_nx6/nx6p2/wt_nx6p2_phi_fun.c
    wt_model_nx6/nx6p2/wt_nx6p2_phi_fun_jac_y.c
    wt_model_nx6/nx6p2/wt_nx6p2_phi_jac_y_uhat.c
    wt_model_nx6/nx6p2/wt_nx6p2_f_lo_fun_jac_x1k1uz.c
    wt_model_nx6/nx6p2/wt_nx6p2_get_matrices_fun.c
)

file(GLOB INV_PENDULUM_SRC
    pendulum_dae_model/pendulum_dae_dyn_impl_ode_fun.c
    pendulum_dae_model/pendulum_dae_dyn_impl_ode_fun_jac_x_xdot.c
    pendulum_dae_model/pendulum_dae_dyn_impl_ode_jac_x_xdot_u.c
    pendulum_dae_model/pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u.c
    pendulum_dae_model/pendulum_dae_dyn_impl_ode_hess.c
    pendulum_dae_model/pendulum_dae_dyn_gnsf_phi_fun.c
    pendulum_dae_model/pendulum_dae_dyn_gnsf_phi_fun_jac_y.c
    pendulum_dae_model/pendulum_dae_dyn_gnsf_phi_jac_y_uhat.c
    pendulum_dae_model/pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz.c
    pendulum_dae_model/pendulum_dae_dyn_gnsf_get_matrices_fun.c
)

file(GLOB GNSF_CRANE_SRC
    crane_nx9_model/crane_nx9_phi_fun.c
    crane_nx9_model/crane_nx9_phi_fun_jac_y.c
    crane_nx9_model/crane_nx9_phi_jac_y_uhat.c
    crane_nx9_model/crane_nx9_f_lo_fun_jac_x1k1uz.c
    crane_nx9_model/crane_nx9_get_matrices_fun.c
)

file(GLOB ENGINE_SRC
    engine_model/engine_impl_dae_fun.c
    engine_model/engine_impl_dae_fun_jac_x_xdot_z.c
    engine_model/engine_impl_dae_jac_x_xdot_u_z.c
    engine_model/engine_ls_cost.c
    engine_model/engine_ls_cost_N.c
)

# Define examples

# -------------------- two-staged turbocharged engine
add_executable(engine_example engine_example.c ${ENGINE_SRC})
target_link_libraries(engine_example acados)

# -------------------- nonlinear_chain
add_executable(nonlinear_chain_ocp_nlp_example nonlinear_chain_ocp_nlp.c ${CHAIN_MODEL_SRC})
target_link_libraries(nonlinear_chain_ocp_nlp_example acados)
add_test(nonlinear_chain_ocp_nlp_example nonlinear_chain_ocp_nlp_example)

# -------------------- wind turbine nmpc
add_executable(wind_turbine_nmpc_example wind_turbine_nmpc.c ${WT_MODEL_NX6P2_SRC})
target_link_libraries(wind_turbine_nmpc_example acados)
add_test(wind_turbine_nmpc_example wind_turbine_nmpc_example)

# -------------------- sim_pendulum_dae
add_executable(sim_pendulum_dae sim_pendulum_dae.c ${INV_PENDULUM_SRC})
target_link_libraries(sim_pendulum_dae acados)
add_test(sim_pendulum_dae sim_pendulum_dae)

# -------------------- sim_wt_nx6
add_executable(sim_wt_model_nx6 sim_wt_model_nx6.c ${WT_MODEL_NX6_SRC})
target_link_libraries(sim_wt_model_nx6 acados)
add_test(sim_wt_model_nx6 sim_wt_model_nx6)

# -------------------- sim_gnsf_crane
add_executable(sim_gnsf_crane sim_gnsf_crane.c ${GNSF_CRANE_SRC})
target_link_libraries(sim_gnsf_crane acados)
add_test(sim_gnsf_crane sim_gnsf_crane)

# -------------------- simple dae_example
add_executable(simple_dae_example simple_dae_example.c
    simple_dae_model/simple_dae_impl_ode_fun.c
    simple_dae_model/simple_dae_impl_ode_jac_x_xdot_u_z.c
    simple_dae_model/simple_dae_impl_ode_fun_jac_x_xdot_z.c
    simple_dae_model/simple_dae_constr_h_fun_jac_ut_xt.c
    simple_dae_model/simple_dae_constr_h_fun_jac_ut_xt_hess.c)
    
target_link_libraries(simple_dae_example acados)
add_test(simple_dae_example simple_dae_example)

# -------------------- mass_spring
add_executable(mass_spring_example no_interface_examples/mass_spring_example.c no_interface_examples/mass_spring_model/mass_spring_qp.c)
target_link_libraries(mass_spring_example acados)
add_test(mass_spring_example mass_spring_example)

#add_executable(mass_spring_fcond_split no_interface_examples/mass_spring_fcond_split.c no_interface_examples/mass_spring_model/mass_spring_qp.c)
#target_link_libraries(mass_spring_fcond_split acados)
#add_test(mass_spring_fcond_split mass_spring_example)

#add_executable(mass_spring_pcond_split no_interface_examples/mass_spring_pcond_split.c no_interface_examples/mass_spring_model/mass_spring_qp.c)
#target_link_libraries(mass_spring_pcond_split acados)
#add_test(mass_spring_pcond_split mass_spring_example)

#add_executable(mass_spring_offline_fcond_qpoases_split no_interface_examples/mass_spring_offline_fcond_qpoases_split.c no_interface_examples/mass_spring_model/mass_spring_qp.c)
#target_link_libraries(mass_spring_offline_fcond_qpoases_split acados)
#add_test(mass_spring_offline_fcond_qpoases_split mass_spring_example)

add_executable(mass_spring_nmpc_example no_interface_examples/mass_spring_nmpc_example.c no_interface_examples/mass_spring_model/mass_spring_qp.c)
target_link_libraries(mass_spring_nmpc_example acados)
add_test(mass_spring_nmpc_example mass_spring_example)

add_executable(sim_crane_example sim_crane_example.c ${CRANE_MODEL_SRC})
target_link_libraries(sim_crane_example acados)
add_test(sim_crane_example sim_crane_example)

# -------------------- sim_wt
add_executable(sim_wt_model_nx3 sim_wt_model_nx3.c ${WT_MODEL_NX3_SRC})
target_link_libraries(sim_wt_model_nx3 acados)
add_test(sim_wt_model_nx3 sim_wt_model_nx3)

add_executable(ocp_qp ocp_qp.c)
target_link_libraries(ocp_qp acados)
add_test(ocp_qp ocp_qp)

# -------------------- dense_qp
add_executable(dense_qp dense_qp.c)
target_link_libraries(dense_qp acados)
add_test(dense_qp dense_qp)

# regularization
add_executable(regularization regularization.c)
target_link_libraries(regularization acados)
add_test(regularization regularization)


endif()
