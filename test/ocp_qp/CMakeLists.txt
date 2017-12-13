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

set(GENERATION_DIR "${CMAKE_SOURCE_DIR}/test/ocp_qp/octave")

add_custom_command(OUTPUT ocp_qp.generated
    COMMAND "octave-cli" "--no-gui" "--path" "${GENERATION_DIR}:${EXTERNAL_SRC_DIR}/casadi-octave-v3.2.3"
        "${GENERATION_DIR}/generate_condensing_data.m"
    COMMAND ${CMAKE_COMMAND} -E touch ocp_qp.generated
    COMMENT "Generating data for OCP QP test"
)

add_library(ocp_qp_gen OBJECT ocp_qp.generated)