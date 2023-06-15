# -*- coding: future_fstrings -*-
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

import sys
sys.path.insert(0, '../common')

import json
from acados_template import AcadosOcp, AcadosOcpSolver, ocp_get_default_cmake_builder
from pendulum_model import export_pendulum_ode_model, export_pendulum_ode_model_with_discrete_rk4
import numpy as np
import scipy.linalg
from utils import plot_pendulum
import argparse


INTEGRATOR_TYPES = ['ERK', 'IRK', 'GNSF', 'DISCRETE']
BUILD_SYSTEMS = ['cmake', 'make']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='test Python interface on pendulum example.')
    parser.add_argument('--INTEGRATOR_TYPE', dest='INTEGRATOR_TYPE',
                        help=f'INTEGRATOR_TYPE: supports {INTEGRATOR_TYPES}')
    parser.add_argument('--BUILD_SYSTEM', dest='BUILD_SYSTEM',
                        default='make',
                        help=f'BUILD_SYSTEM: supports {BUILD_SYSTEMS}')

    args = parser.parse_args()

    integrator_type = args.INTEGRATOR_TYPE
    if integrator_type not in INTEGRATOR_TYPES:
        msg = f'Invalid unit test value {integrator_type} for parameter INTEGRATOR_TYPE. Possible values are' \
              f' {INTEGRATOR_TYPES}, got {integrator_type}.'
        raise Exception(msg)

    build_system = args.BUILD_SYSTEM
    if build_system not in BUILD_SYSTEMS:
        msg = f'Invalid unit test value {build_system} for parameter INTEGRATOR_TYPE. Possible values are' \
              f' {BUILD_SYSTEMS}, got {build_system}.'
        raise Exception(msg)

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    Tf = 1.0
    N = 20

    # set model
    if integrator_type == 'DISCRETE':
        model = export_pendulum_ode_model_with_discrete_rk4(Tf/N)
    else:
        model = export_pendulum_ode_model()

    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx

    # set dimensions
    ocp.dims.N = N

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
    ocp.constraints.idxbu = np.array([0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'  # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = integrator_type
    ocp.solver_options.print_level = 1

    if ocp.solver_options.integrator_type == 'GNSF':
        from acados_template import acados_dae_model_json_dump
        import os
        acados_dae_model_json_dump(model)
        # Set up Octave to be able to run the following:
        ## if using a virtual python env, the following lines can be added to the env/bin/activate script:
        # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/external/casadi-octave
        # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/
        # export OCTAVE_PATH=$OCTAVE_PATH:$ACADOS_INSTALL_DIR/interfaces/acados_matlab_octave/acados_template_mex/
        # echo
        # echo "OCTAVE_PATH=$OCTAVE_PATH"

        # status = os.system(
        #     "octave --eval \"convert_dae2gnsf({})\"".format("\'"+model.name+"_acados_dae.json\'")
        # )
        # if status == 0:
        #     print("\nsuccessfully detected GNSF structure in Octave\n")
        # else:
        #     Exception("Failed to detect GNSF structure in Octave")

        # load gnsf from json
        with open(model.name + '_gnsf_functions.json', 'r') as f:
            import json
            gnsf_dict = json.load(f)
        ocp.gnsf_model = gnsf_dict

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'  # SQP_RTI

    if build_system == 'cmake':
        print("\nusing the CMake build system")
        cmake_builder = ocp_get_default_cmake_builder()
        ocp_solver = AcadosOcpSolver(ocp, json_file='acados_ocp.json', cmake_builder=cmake_builder)
    elif build_system == 'make':
        print("\nusing the make build system")
        ocp_solver = AcadosOcpSolver(ocp, json_file='acados_ocp.json')

    simX = np.ndarray((N+1, nx))
    simU = np.ndarray((N, nu))

    status = ocp_solver.solve()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX[i, :] = ocp_solver.get(i, "x")
        simU[i, :] = ocp_solver.get(i, "u")
    simX[N, :] = ocp_solver.get(N, "x")

    plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX)
