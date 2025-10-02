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


if(ACADOS_OCTAVE)
    # new minimal examples
    add_test(NAME octave_new_ocp
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/getting_started
        octave --no-gui --no-window-system ./minimal_example_ocp.m)
    add_test(NAME octave_new_sim
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/getting_started
        octave --no-gui --no-window-system ./minimal_example_sim.m)
    # dim check test
    add_test(NAME octave_test_dim_check
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_dim_check.m)

    # SIM
    add_test(NAME octave_test_sim_forw
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_sim_forw.m)

    add_test(NAME octave_test_sim_adj
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_sim_adj.m)

    add_test(NAME octave_test_sim_hess
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_sim_hess.m)

    add_test(NAME octave_test_sim_dae
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_sim_dae.m)

    # OCP
    add_test(NAME octave_test_ocp_pendulum
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_ocp_pendulum.m)
    add_test(NAME octave_test_ocp_pendulum_code_reuse
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./test_code_reuse.m)
    add_test(NAME octave_test_ocp_wtnx6
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_ocp_wtnx6.m)
    add_test(NAME octave_test_ocp_mass_spring
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_test_ocp_mass_spring.m)
    add_test(NAME octave_test_ocp_simple_dae
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/simple_dae_model
        octave --no-gui --no-window-system ./example_ocp.m)

    # AS-RTI closed loop example
    add_test(NAME octave_closed_loop_as_rti_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/pendulum_on_cart_model
        octave --no-gui --no-window-system ./closed_loop_as_rti.m)

    add_test(NAME octave_test_mhe_lorentz
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/lorentz
        octave --no-gui --no-window-system ./example_mhe.m)

    # OSQP
    if(ACADOS_WITH_OSQP)
        add_test(NAME octave_test_OSQP
            COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
            octave --no-gui --no-window-system ./test_ocp_OSQP.m)
    endif()

    # qpDUNES
    if(ACADOS_WITH_QPDUNES)
    add_test(NAME octave_test_qpDUNES
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./test_ocp_qpdunes.m)
    endif()
    #
    add_test(NAME octave_test_target_selector
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./test_target_selector.m)
    add_test(NAME octave_test_generic_ext_cost_disc_dyn
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/generic_dyn_disc
        octave --no-gui --no-window-system ./run_generic_disc_dyn_example.m)

    # MOCP
    add_test(NAME octave_mocp_transition_test
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/mocp_transition_example
        octave --no-gui --no-window-system ./main_multiphase_ocp.m)

    # Force serial execution for tests in the same directory
    # Directory test
    set_tests_properties(octave_new_ocp PROPERTIES DEPENDS octave_new_sim)

    set_tests_properties(octave_test_target_selector PROPERTIES DEPENDS octave_test_generic_ext_cost_disc_dyn)
    set_tests_properties(octave_test_generic_ext_cost_disc_dyn PROPERTIES DEPENDS octave_test_dim_check)
    set_tests_properties(octave_test_dim_check PROPERTIES DEPENDS octave_test_sim_forw)
    set_tests_properties(octave_test_sim_forw PROPERTIES DEPENDS octave_test_sim_adj)
    set_tests_properties(octave_test_sim_adj PROPERTIES DEPENDS octave_test_sim_hess)
    set_tests_properties(octave_test_sim_hess PROPERTIES DEPENDS octave_test_sim_dae)
    set_tests_properties(octave_test_sim_dae PROPERTIES DEPENDS octave_test_ocp_pendulum)
    set_tests_properties(octave_test_ocp_pendulum PROPERTIES DEPENDS octave_test_ocp_wtnx6)
    set_tests_properties(octave_test_ocp_wtnx6 PROPERTIES DEPENDS octave_test_ocp_mass_spring)
    set_tests_properties(octave_test_ocp_mass_spring PROPERTIES DEPENDS octave_test_ocp_simple_dae)
    set_tests_properties(octave_test_ocp_simple_dae PROPERTIES DEPENDS octave_test_target_selector)
    set_tests_properties(octave_test_target_selector PROPERTIES DEPENDS octave_test_ocp_pendulum_code_reuse)
    if(ACADOS_WITH_OSQP)
        set_tests_properties(octave_test_ocp_pendulum_code_reuse PROPERTIES DEPENDS octave_test_OSQP)
    endif()
    if(ACADOS_WITH_QPDUNES)
        set_tests_properties(octave_test_ocp_pendulum_code_reuse PROPERTIES DEPENDS "octave_test_qpDUNES;octave_test_OSQP")
        set_tests_properties(octave_test_OSQP PROPERTIES DEPENDS octave_test_qpDUNES)
    endif()
endif()
