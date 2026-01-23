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
import casadi as ca

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosMultiphaseOcp, AcadosSim, AcadosSimSolver
from create_mocp import create_mocp, export_double_integrator_model


def create_ocp():
    ocp = AcadosOcp()

    dim_q = 1
    nx = 2 * dim_q
    nu = dim_q
    ny = nx + nu
    # set model
    model = export_double_integrator_model(dim_q=dim_q, dt=0.1)
    ocp.model = model

    # set dimensions
    ocp.solver_options.N_horizon = 3
    ocp.solver_options.tf = 1

    # set cost
    Q = np.diag([1.0, 1.0])
    R = np.diag([0.1])
    ocp.cost.W = ca.diagcat(Q, R).full()
    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx, :] = np.eye(nx)
    ocp.cost.Vu = np.zeros((ny, nu))
    ocp.cost.Vu[nx:, :] = np.eye(nu)
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.cost.yref = np.zeros((model.x.size()[0] + model.u.size()[0],))

    # set constraints
    ocp.constraints.lbu = np.array([-1.0])
    ocp.constraints.ubu = np.array([1.0])
    ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = np.ones((nx,))

    return ocp

def main(problem_class):
    json_file = 'mocp.json'
    creation_modes = ['standard', 'precompiled']
    for i, creation_mode in enumerate(creation_modes):
        # ocp = create_ocp()
        if problem_class == 'mocp':
            ocp = create_mocp()
        elif problem_class == 'ocp':
            ocp = create_ocp()
        else:
            raise Exception(f'Unknown problem class {problem_class}')
        ocp.code_gen_opts.json_file = json_file

        if creation_mode == 'standard':
            ocp_solver = AcadosOcpSolver(ocp)
        elif creation_mode == 'precompiled':
            ocp_solver = AcadosOcpSolver(ocp, build=False, generate=False)

        status = ocp_solver.solve()
        ocp_solver.print_statistics()

        if i == 0:
            ref_sol = ocp_solver.get_iterate()
        else:
            test_sol = ocp_solver.get_iterate()
            assert ref_sol.allclose(test_sol, atol=1e-6), 'Solutions do not match!'

        if creation_mode == 'precompiled':
            assert not ocp_solver.generated, 'Expected reused code, but code was generated!'

        assert status == 0, f'acados returned status {status}'
        del ocp_solver
    

def main_sim():
    json_file = 'acados_sim.json'
    creation_modes = ['standard', 'precompiled']
    for i, creation_mode in enumerate(creation_modes):

        sim = AcadosSim()
        dim_q = 2
        sim.model = export_double_integrator_model(dim_q=dim_q, dt=0.1)
        sim.solver_options.T = 0.1
        sim.solver_options.integrator_type = 'ERK'
        sim.code_gen_opts.json_file = json_file

        if creation_mode == 'standard':
            sim_solver = AcadosSimSolver(sim)
        elif creation_mode == 'precompiled':
            sim_solver = AcadosSimSolver(sim, build=False, generate=False)

        sim_solver.solve()

        if i == 0:
            ref_sol = sim_solver.get('x')
        else:
            test_sol = sim_solver.get('x')
            assert np.allclose(ref_sol, test_sol, atol=1e-6), 'Solutions do not match!'

        if creation_mode == 'precompiled':
            assert not sim_solver.generated, 'Expected reused code, but code was generated!'

        del sim_solver

if __name__ == '__main__':
    main_sim()
    main('ocp')
    main('mocp')
