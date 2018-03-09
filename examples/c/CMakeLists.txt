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
file(GLOB PENDULUM_MODEL_SRC pendulum_model/pendulum_model.c
    pendulum_model/vde_forw_pendulum.c
    pendulum_model/jac_pendulum.c)

# Define examples

if(NOT CMAKE_C_COMPILER_ID MATCHES "MSVC")
    add_executable(mass_spring_example mass_spring_example.c)
    target_link_libraries(mass_spring_example acados)
    add_test(mass_spring_example mass_spring_example)

    add_executable(nonlinear_chain_ocp_nlp_example nonlinear_chain_ocp_nlp.c ${CHAIN_MODEL_SRC})
    target_link_libraries(nonlinear_chain_ocp_nlp_example acados)
    add_test(nonlinear_chain_ocp_nlp_example nonlinear_chain_ocp_nlp_example)
endif()

add_executable(ocp_qp ocp_qp.c)
target_link_libraries(ocp_qp acados)
add_test(ocp_qp ocp_qp)
