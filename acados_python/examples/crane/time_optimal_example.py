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

import matplotlib.pyplot as plt

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, AcadosSimSolver, AcadosSim
import numpy as np
from casadi import SX, vertcat
import os

def main(use_cython=True):
    # (very) simple crane model
    beta = 0.001
    k = 0.9
    a_max = 10
    dt_max = 2.0

    # states
    p1 = SX.sym('p1')
    v1 = SX.sym('v1')
    p2 = SX.sym('p2')
    v2 = SX.sym('v2')

    x = vertcat(p1, v1, p2, v2)

    # controls
    a = SX.sym('a')
    dt = SX.sym('dt')

    u = vertcat(a, dt)

    f_expl = dt*vertcat(v1, a, v2, -beta*v2-k*(p2 - p1))

    model = AcadosModel()

    model.f_expl_expr = f_expl
    model.x = x
    model.u = u
    model.name = 'crane_time_opt'

    # create ocp object to formulate the OCP

    x0 = np.array([2.0, 0.0, 2.0, 0.0])
    xf = np.array([0.0, 0.0, 0.0, 0.0])

    ocp = AcadosOcp()
    ocp.model = model

    # N - maximum number of bangs
    N = 7
    Tf = N
    nx = model.x.size()[0]
    nu = model.u.size()[0]

    # set dimensions
    ocp.dims.N = N

    # set cost
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    ocp.model.cost_expr_ext_cost = dt
    ocp.model.cost_expr_ext_cost_e = 0

    ocp.constraints.lbu = np.array([-a_max, 0.0])
    ocp.constraints.ubu = np.array([+a_max, dt_max])
    ocp.constraints.idxbu = np.array([0, 1])

    ocp.constraints.x0 = x0
    ocp.constraints.lbx_e = xf
    ocp.constraints.ubx_e = xf
    ocp.constraints.idxbx_e = np.array([0, 1, 2, 3])

    # set prediction horizon
    ocp.solver_options.tf = Tf

    # set options
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'#'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 3
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    ocp.solver_options.nlp_solver_max_iter = 5000
    ocp.solver_options.nlp_solver_tol_stat = 1e-6
    ocp.solver_options.levenberg_marquardt = 0.1
    ocp.solver_options.sim_method_num_steps = 15
    ocp.solver_options.qp_solver_iter_max = 100
    ocp.code_export_directory = 'c_generated_code'
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.exact_hess_constr = 0
    ocp.solver_options.exact_hess_dyn = 0

    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file='acados_ocp.json')
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        ocp_solver = AcadosOcpSolver.create_cython_solver('acados_ocp.json')
    else: # ctypes
        ## Note: skip generate and build assuming this is done before (in cython run)
        ocp_solver = AcadosOcpSolver(ocp, json_file='acados_ocp.json', build=False, generate=False)

    ocp_solver.reset()

    for i, tau in enumerate(np.linspace(0, 1, N)):
        ocp_solver.set(i, 'x', (1-tau)*x0 + tau*xf)
        ocp_solver.set(i, 'u', np.array([0.1, 0.5]))

    simX = np.zeros((N+1, nx))
    simU = np.zeros((N, nu))

    status = ocp_solver.solve()

    if status != 0:
        ocp_solver.print_statistics()
        raise Exception(f'acados returned status {status}.')

    # get solution
    for i in range(N):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N,:] = ocp_solver.get(N, "x")

    dts = simU[:, 1]

    print("acados solved OCP successfully, creating integrator to simulate the solution")
    print(f"optimal time: {sum(dts)} s")

    # simulate on finer grid
    sim = AcadosSim()

    # set model
    sim.model = model

    # set options
    sim.solver_options.integrator_type = 'ERK'
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 3
    sim.solver_options.T = 1.0 # dummy value

    dt_approx = 0.0005

    dts_fine = np.zeros((N,))
    Ns_fine = np.zeros((N,), dtype='int16')

    # compute number of simulation steps for bang interval + dt_fine
    for i in range(N):
        N_approx = max(int(dts[i]/dt_approx), 1)
        dts_fine[i] = dts[i]/N_approx
        Ns_fine[i] = int(round(dts[i]/dts_fine[i]))

    N_fine = int(np.sum(Ns_fine))

    simU_fine = np.zeros((N_fine, nu))
    ts_fine = np.zeros((N_fine+1, ))
    simX_fine = np.zeros((N_fine+1, nx))
    simX_fine[0, :] = x0

    acados_integrator = AcadosSimSolver(sim)

    k = 0
    for i in range(N):
        u = simU[i, 0]
        acados_integrator.set("u", np.hstack((u, np.ones(1, ))))

        # set simulation time
        acados_integrator.set("T", dts_fine[i])

        for j in range(Ns_fine[i]):
            acados_integrator.set("x", simX_fine[k,:])
            status = acados_integrator.solve()
            if status != 0:
                raise Exception(f'acados returned status {status}.')

            simX_fine[k+1,:] = acados_integrator.get("x")
            simU_fine[k, :] = u
            ts_fine[k+1] = ts_fine[k] + dts_fine[i]

            k += 1

    # visualize
    if not os.environ.get('ACADOS_ON_CI'):
        plt.figure()

        state_labels = ['p1', 'v1', 'p2', 'v2']

        for i,l in enumerate(state_labels):
            plt.subplot(5, 1, i+1)

            plt.plot(ts_fine, simX_fine[:, i], label='time optimal solution')
            plt.grid(True)
            plt.ylabel(l)
            if i ==0:
                plt.legend(loc=1)

        plt.subplot(5, 1, 5)
        plt.step(ts_fine, np.hstack((simU_fine[:, 0], simU_fine[-1, 0])), '-', where='post')
        plt.grid(True)
        plt.ylabel('a')
        plt.xlabel('t')

        plt.show()

if __name__ == "__main__":
    for use_cython in [True, False]:
        main(use_cython=use_cython)

