#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


# set(GENERATION_DIR "${CMAKE_SOURCE_DIR}/test/ocp_qp/octave")

# add_custom_command(OUTPUT ocp_qp.generated
#     COMMAND "octave-cli" "--no-gui" "--path" "${GENERATION_DIR}:${EXTERNAL_SRC_DIR}/casadi-octave-v3.2.3"
#         "${GENERATION_DIR}/generate_condensing_data.m"
#     COMMAND ${CMAKE_COMMAND} -E touch ocp_qp.generated
#     COMMENT "Generating data for OCP QP test"
# )

# add_library(ocp_qp_gen OBJECT ocp_qp.generated)