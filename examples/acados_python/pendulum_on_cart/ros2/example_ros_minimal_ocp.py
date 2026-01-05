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

import numpy as np
import scipy.linalg

from acados_template import AcadosOcp, AcadosOcpSolver
from acados_template.ros2 import AcadosOcpRosOptions

import sys
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
common_path = os.path.join(script_dir, '..', 'common')
sys.path.insert(0, os.path.abspath(common_path))
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum

def create_minimal_ocp(export_dir: str, N: int = 20, Tf: float = 1.0, Fmax: float = 80):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # set dimensions
    ocp.solver_options.N_horizon = N

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints

    # bound on u
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    # initial state
    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'

    # set prediction horizon
    ocp.solver_options.tf = Tf

    # Ros stuff
    ocp.ros_opts = AcadosOcpRosOptions()
    ocp.ros_opts.package_name = "pendulum_on_cart_ocp"
    ocp.ros_opts.generated_code_dir = export_dir
    ocp.ros_opts.publish_control_sequence = True

    ocp.code_export_directory = str(os.path.join(export_dir, "c_generated_code"))
    return ocp



def main():
    Fmax = 80
    Tf = 1.0
    N = 20

    export_dir = os.path.join(script_dir, 'generated_ocp')
    ocp = create_minimal_ocp(export_dir, N, Tf, Fmax)
    ocp_solver = AcadosOcpSolver(ocp, json_file = str(os.path.join(export_dir, 'acados_ocp.json')))

    simX = np.zeros((N+1, ocp.dims.nx))
    simU = np.zeros((N, ocp.dims.nu))

    # call SQP_RTI solver in the loop:
    tol = 1e-6

    for i in range(20):
        status = ocp_solver.solve()
        ocp_solver.print_statistics()
        residuals = ocp_solver.get_residuals(recompute=True)
        print("residuals after ", i, "SQP_RTI iterations:\n", residuals)
        if max(residuals) < tol:
            break

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N,:] = ocp_solver.get(N, "x")

    expected_u_file = os.path.join(export_dir, 'expected_control_sequence.npy')
    np.save(expected_u_file, simU)

    ocp_solver.print_statistics()

    # plot
    plot_pendulum(np.linspace(0, Tf, N+1), Fmax, simU, simX, latexify=False)


if __name__ == "__main__":
    main()
