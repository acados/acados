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

set(TEST_SIM_GENERATED_FILES
    "${PROJECT_SOURCE_DIR}/build/test/sim/jac_pendulum.c"
    "${PROJECT_SOURCE_DIR}/build/test/sim/vde_forw_pendulum.c"
    "${PROJECT_SOURCE_DIR}/build/test/sim/vde_adj_pendulum.c"
    "${PROJECT_SOURCE_DIR}/build/test/sim/vde_hess_pendulum.c"
)

set(GENERATION_DIR "${CMAKE_SOURCE_DIR}/test/sim/pendulum/casadi")
set(SED_COMMAND "'s/real_t\ a0/mem\ =\ 0;\ mem\ +=\ 0;\ w\ =\ 0;\ w\ +=\ 0;\ iw\ =\ 0;\ iw\ +=\ 0;  real_t\ a0/g'")

add_custom_command(OUTPUT ${TEST_SIM_GENERATED_FILES}
    COMMAND "octave-cli" "--no-gui" "--path" "${GENERATION_DIR}:${EXTERNAL_SRC_DIR}/casadi-octave-v3.2.3"
        "${GENERATION_DIR}/generate_pendulum.m"
    COMMAND "sed" "-i.bak" "-e" "${SED_COMMAND}" ${TEST_SIM_GENERATED_FILES}
    COMMENT "Generating data for numerical integration test"
)

add_library(sim_gen OBJECT ${TEST_SIM_GENERATED_FILES})
