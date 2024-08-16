# This test is an extension of the 'minimal_example_ocp_reuse_code.py' example.
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

sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum

PLOT = False

def main(interface_type='ctypes'):

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    N_horizon = 20  # number of shooting nodes
    Tf = 1.0

    # set dimensions
    ocp.solver_options.N_horizon = N_horizon

    # set cost
    Q = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2 * np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx, :nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4, 0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP'

    # set prediction horizon
    ocp.solver_options.tf = Tf

    print(80*'-')
    print('generate code and compile...')

    if interface_type == 'cython':
        AcadosOcpSolver.generate(ocp, json_file='acados_ocp.json')
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        ocp_solver = AcadosOcpSolver.create_cython_solver('acados_ocp.json')
    elif interface_type == 'ctypes':
        ocp_solver = AcadosOcpSolver(ocp, json_file='acados_ocp.json')
    elif interface_type == 'cython_prebuilt':
        from c_generated_code.acados_ocp_solver_pyx import AcadosOcpSolverCython
        ocp_solver = AcadosOcpSolverCython(ocp.model.name, ocp.solver_options.nlp_solver_type, ocp.solver_options.N_horizon)


    # test setting HPIPM options
    ocp_solver.options_set('qp_tol_ineq', 1e-8)
    ocp_solver.options_set('qp_tau_min', 1e-10)
    ocp_solver.options_set('qp_mu0', 1e0)

    # solve the problem defined here (original from code export), analog to 'minimal_example_ocp.py'
    nvariant = 0
    simX0 = np.zeros((N_horizon + 1, nx))
    simU0 = np.zeros((N_horizon, nu))

    print(80*'-')
    print(f'solve original code with N = {N_horizon} and Tf = {Tf} s:')
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N_horizon):
        simX0[i, :] = ocp_solver.get(i, "x")
        simU0[i, :] = ocp_solver.get(i, "u")
    simX0[N_horizon, :] = ocp_solver.get(N_horizon, "x")

    ocp_solver.store_iterate(filename=f'final_iterate_{interface_type}_variant{nvariant}.json', overwrite=True)

    if PLOT:# plot but don't halt
        plot_pendulum(np.linspace(0, Tf, N_horizon + 1), Fmax, simU0, simX0, latexify=False, plt_show=False, X_true_label=f'original: N={N_horizon}, Tf={Tf}')


if __name__ == "__main__":
    for interface_type in ['ctypes', 'cython', 'cython_prebuilt']:
        main(interface_type=interface_type)

    import json
    # compare iterates
    for nvariant in [0]:
        iterate_filename = f'final_iterate_ctypes_variant{nvariant}.json'
        with open(iterate_filename, 'r') as f:
            iterate_ctypes = json.load(f)

        for interface_type in ['cython', 'cython_prebuilt']:
            iterate_filename = f'final_iterate_{interface_type}_variant{nvariant}.json'
            with open(iterate_filename, 'r') as f:
                iterate = json.load(f)

            assert iterate.keys() == iterate_ctypes.keys()
            assert(all ([iterate[k] == iterate_ctypes[k] for k in iterate]))
