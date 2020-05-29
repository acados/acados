#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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


if(ACADOS_OCTAVE_TEMPLATE)
add_test(NAME octave_test_templates
         COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
         octave --no-gui --no-window-system ./run_tests_templates.m)
endif()


if(ACADOS_MATLAB OR ACADOS_OCTAVE)
add_test(NAME octave_test_sim
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
        octave --no-gui --no-window-system ./run_tests_sim.m)

add_test(NAME octave_test_ocp
         COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_matlab_octave/test
         octave --no-gui --no-window-system ./run_tests_ocp.m)
endif()


if(ACADOS_PYTHON)
add_test(NAME python_test_ocp
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python main_test.py)
add_test(NAME python_pendulum_ocp_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started/ocp
        python minimal_example_ocp.py)
add_test(NAME python_nonuniform_discretization_ocp_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started/ocp
        python nonuniform_discretization_example.py)
add_test(NAME python_pendulum_mhe_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started/mhe
        python minimal_example_mhe.py)
add_test(NAME python_pendulum_ocp_cost_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started/ocp
        python ocp_example_cost_formulations.py)
add_test(NAME python_pendulum_sim_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started/sim
        python minimal_example_sim.py)
add_test(NAME python_pendulum_ext_sim_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started/sim
        python extensive_example_sim.py)
add_test(NAME python_pendulum_closed_loop_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started
        python minimal_example_closed_loop.py)
add_test(NAME python_pendulum_GNSF_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/getting_started/ocp
        python example_gnsf_ocp.py)
add_test(NAME python_pendulum_soft_constraints_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python soft_constraint_test.py)
add_test(NAME python_pendulum_parametric_nonlinear_constraint_h_test
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/tests
        python test_parametric_nonlinear_constraint_h.py)
add_test(NAME python_pmsm_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/pmsm_example
        python generate_c_code.py)
add_test(NAME python_race_cars
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/race_cars
        python main.py)
add_test(NAME python_rsm_example
        COMMAND "${CMAKE_COMMAND}" -E chdir ${PROJECT_SOURCE_DIR}/examples/acados_python/rsm_example
        python generate_c_code.py)
endif()
