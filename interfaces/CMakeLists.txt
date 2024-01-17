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


# MATLAB / OCTAVE
if(ACADOS_MATLAB OR ACADOS_OCTAVE)
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
        add_test(NAME octave_test_ocp_pendulum_dae
                COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
                octave --no-gui --no-window-system ./run_test_ocp_pendulum_dae.m)
        add_test(NAME octave_test_ocp_mass_spring
                COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
                octave --no-gui --no-window-system ./run_test_ocp_mass_spring.m)
        add_test(NAME octave_test_ocp_simple_dae
                COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
                octave --no-gui --no-window-system ./run_test_ocp_simple_dae.m)

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
endif()



### PYTHON ###
if(ACADOS_PYTHON)

# Minimal examples
add_test(NAME python_pendulum_closed_loop_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started
        python minimal_example_closed_loop.py)
add_test(NAME python_pendulum_sim_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started
        python minimal_example_sim.py)
add_test(NAME python_pendulum_ocp_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started
        python minimal_example_ocp.py)

add_test(NAME python_test_generic_impl_dyn
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/generic_impl_dyn
        python minimal_example_ocp_generic_impl_dyn.py)
add_test(NAME python_test_reset
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python reset_test.py)
add_test(NAME python_test_reset_timing
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/timing_example
        python reset_timing.py)
add_test(NAME python_test_cython_vs_ctypes
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python test_cython_ctypes.py)
add_test(NAME python_test_ocp
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python main_test.py)

add_test(NAME python_test_cost_integration_euler
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python test_cost_integration_euler.py)

add_test(NAME python_test_cost_integration_value
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python test_cost_integration_value.py)

add_test(NAME python_pendulum_ocp_example_reuse_code
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python minimal_example_ocp_reuse_code.py)

add_test(NAME python_solution_sensitivities_and_exact_hess
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/solution_sensitivities
        python test_solution_sens_and_exact_hess.py)

add_test(NAME python_time_varying_irk
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/time_varying
        python test_time_varying_irk.py)

add_test(NAME test_polynomial_controls_and_penalties
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/time_varying
        python test_polynomial_controls_and_penalties.py)

add_test(NAME py_mocp_qp_test
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/time_varying
        python test_mocp_qp.py)


# CMake test
add_test(NAME python_pendulum_ocp_example_cmake
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python minimal_example_ocp_cmake.py)
add_test(NAME python_nonuniform_discretization_ocp_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python nonuniform_discretization_example.py)
add_test(NAME python_rti_loop_ocp_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python example_sqp_rti_loop.py)
# Python Simulink
add_test(NAME python_render_simulink_wrapper
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python simulink_example.py)
add_test(NAME python_constraints_expression_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python ocp_example_h_init_contraints.py)

add_test(NAME python_chain_ocp
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/chain_mass
        python main.py)

add_test(NAME python_chain_sim
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/chain_mass
        python minimal_example_sim.py)

# Marathos test problem with different globalization options
add_test(NAME python_marathos_test_problem_globalization
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/non_ocp_nlp
        python marathos_test_problem.py)

# Simple OCP with Marathos effect
add_test(NAME python_OCP_marathos_test_problem_globalization
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/linear_mass_model
        python linear_mass_test_problem.py)

# Armijo test problem
add_test(NAME python_armijo_test
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python armijo_test.py)

# Multiphase nonlinear constraint test problem
add_test(NAME python_multiphase_nonlinear_constraints
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/multiphase_nonlinear_constraints
        python main.py)

# qpDUNES example
if(ACADOS_WITH_QPDUNES)
        add_test(NAME python_qpDUNES_test
                COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
                python example_sqp_qpDUNES.py)
endif()

# OSQP test
if(ACADOS_WITH_OSQP)
        add_test(NAME python_OSQP_test
                COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
                python test_osqp.py)
endif()

# MHE examples
add_test(NAME python_pendulum_mhe_example_minimal
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/mhe
        python minimal_example_mhe.py)
add_test(NAME python_pendulum_mhe_example_noisy_param
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/mhe
        python minimal_example_mhe_with_noisy_param.py)
add_test(NAME python_pendulum_mhe_example_param
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/mhe
        python minimal_example_mhe_with_param.py)
add_test(NAME python_pendulum_mhe_ocp_closed_loop_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/mhe
        python closed_loop_mhe_ocp.py)

add_test(NAME python_pendulum_ocp_cost_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python ocp_example_cost_formulations.py)

add_test(NAME python_custom_update_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/custom_update
        python example_custom_rti_loop.py)

add_test(NAME python_fast_zoro_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/zoRO_example
        python pendulum_on_cart/minimal_example_zoro.py)

add_test(NAME python_zoro_diff_drive_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/zoRO_example
        python diff_drive/main.py)

# CMake test
add_test(NAME python_pendulum_sim_example_cmake
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/sim
        python minimal_example_sim_cmake.py)
add_test(NAME python_pendulum_ext_sim_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/sim
        python extensive_example_sim.py)

add_test(NAME cython_pendulum_closed_loop_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart
        python cython_example_closed_loop.py)
add_test(NAME py_time_optimal_cython_ctypes
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/crane
        python time_optimal_example.py)
add_test(NAME pendulum_optimal_value_gradient
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python example_optimal_value_derivative.py)

# example_ocp_dynamics_formulations all versions
add_test(NAME python_pendulum_ocp_IRK
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python example_ocp_dynamics_formulations.py --INTEGRATOR_TYPE=IRK)
add_test(NAME python_pendulum_ocp_ERK
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python example_ocp_dynamics_formulations.py --INTEGRATOR_TYPE=ERK)
# add_test(NAME python_pendulum_ocp_DISCRETE
#         COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
#         python example_ocp_dynamics_formulations.py --INTEGRATOR_TYPE=DISCRETE)
add_test(NAME python_pendulum_ocp_GNSF
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python example_ocp_dynamics_formulations.py --INTEGRATOR_TYPE=GNSF)
# CMake and solver=DISCRETE test
add_test(NAME python_example_ocp_dynamics_formulations_cmake
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pendulum_on_cart/ocp
        python example_ocp_dynamics_formulations.py --INTEGRATOR_TYPE DISCRETE --BUILD_SYSTEM cmake)

add_test(NAME python_pendulum_soft_constraints_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python soft_constraint_test.py)
add_test(NAME python_pendulum_parametric_nonlinear_constraint_h_test
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python test_parametric_nonlinear_constraint_h.py)
add_test(NAME python_test_sim_dae
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python test_sim_dae.py)

add_test(NAME python_pmsm_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pmsm_example
        python main.py)
add_test(NAME python_race_cars
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/race_cars
        python main.py)
add_test(NAME python_rsm_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/rsm_example
        python main.py)
endif()
