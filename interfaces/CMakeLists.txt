#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


if(ACADOS_OCTAVE)
    # new minimal examples
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

    # Force serial execution for tests in the same directory
    # Directory test
    set_tests_properties(octave_test_generic_ext_cost_disc_dyn PROPERTIES DEPENDS octave_test_target_selector)
    set_tests_properties(octave_test_dim_check PROPERTIES DEPENDS octave_test_generic_ext_cost_disc_dyn)
    set_tests_properties(octave_test_sim_forw PROPERTIES DEPENDS octave_test_dim_check)
    set_tests_properties(octave_test_sim_adj PROPERTIES DEPENDS octave_test_sim_forw)
    set_tests_properties(octave_test_sim_hess PROPERTIES DEPENDS octave_test_sim_adj)
    set_tests_properties(octave_test_sim_dae PROPERTIES DEPENDS octave_test_sim_hess)
    set_tests_properties(octave_test_ocp_pendulum PROPERTIES DEPENDS octave_test_sim_dae)
    set_tests_properties(octave_test_ocp_wtnx6 PROPERTIES DEPENDS octave_test_ocp_pendulum)
    set_tests_properties(octave_test_ocp_mass_spring PROPERTIES DEPENDS octave_test_ocp_wtnx6)
    set_tests_properties(octave_test_ocp_simple_dae PROPERTIES DEPENDS octave_test_ocp_mass_spring)

    if(ACADOS_WITH_OSQP)
        set_tests_properties(octave_test_OSQP PROPERTIES DEPENDS octave_test_target_selector)
        set_tests_properties(octave_test_target_selector PROPERTIES DEPENDS octave_test_OSQP)
    endif()

    if(ACADOS_WITH_QPDUNES)
        set_tests_properties(octave_test_qpDUNES PROPERTIES DEPENDS octave_test_target_selector)
        set_tests_properties(octave_test_OSQP PROPERTIES DEPENDS octave_test_qpDUNES)
    endif()

endif()
