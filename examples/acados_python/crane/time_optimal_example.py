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
from casadi import SX, vertcat, diag
from typing import Tuple

def plot_crane_trajectories(ts, simX, simU):

    state_labels = ['$p_1$', '$v_1$', '$p_2$', '$v_2$']
    fig, axes = plt.subplots(nrows=len(state_labels)+1, ncols=1, sharex=True)

    for i, l in enumerate(state_labels):
        axes[i].plot(ts, simX[:, i])
        axes[i].grid(True)
        axes[i].set_ylabel(l)

    axes[0].set_title('time optimal solution')

    axes[-1].step(ts, np.hstack((simU[:, 0], simU[-1, 0])), '-', where='post')
    axes[-1].grid(True)
    axes[-1].set_ylabel('a')
    axes[-1].set_xlabel('t')
    axes[-1].set_xlim(ts[0], ts[-1])

    fig.align_ylabels()

    plt.show()

def setup_solver_and_integrator(x0: np.ndarray, xf: np.ndarray, N: int, use_cython: bool = True) -> Tuple[AcadosOcpSolver, AcadosSimSolver]:

    # (very) simple crane model
    beta = 0.001
    k = 0.9
    a_max = 10
    dt_max = 2.0
    dt_min = 1e-3

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

    ocp = AcadosOcp()
    ocp.model = model
    ocp.dims.N = N
    ocp.solver_options.tf = N

    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    ocp.model.cost_expr_ext_cost = dt
    ocp.model.cost_expr_ext_cost_e = 0

    ocp.model.cost_expr_ext_cost_custom_hess = diag(vertcat(SX.zeros(1, 1), 1./(dt), SX.zeros(model.x.rows(), 1)))

    ocp.constraints.lbu = np.array([-a_max, dt_min])
    ocp.constraints.ubu = np.array([+a_max, dt_max])
    ocp.constraints.idxbu = np.array([0, 1])

    ocp.constraints.x0 = x0
    ocp.constraints.lbx_e = xf
    ocp.constraints.ubx_e = xf
    ocp.constraints.idxbx_e = np.array([0, 1, 2, 3])

    # set options
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'#'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 3
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    ocp.solver_options.nlp_solver_max_iter = 5000
    ocp.solver_options.nlp_solver_tol_stat = 1e-6
    ocp.solver_options.sim_method_num_steps = 15
    ocp.solver_options.qp_solver_iter_max = 100
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

    # integrator
    sim = AcadosSim()

    sim.model = model
    sim.solver_options.integrator_type = 'ERK'
    sim.solver_options.num_stages = 4
    sim.solver_options.num_steps = 3
    sim.solver_options.T = 1.0 # dummy value
    integrator = AcadosSimSolver(sim)

    return ocp_solver, integrator

def main(use_cython=True):

    nu = 2
    nx = 4

    N = 7 # N - maximum number of bangs

    x0 = np.array([2.0, 0.0, 2.0, 0.0])
    xf = np.zeros((nx,))

    ocp_solver, integrator = setup_solver_and_integrator(x0, xf, N, use_cython)

    # initialization
    for i, tau in enumerate(np.linspace(0, 1, N)):
        ocp_solver.set(i, 'x', (1-tau)*x0 + tau*xf)
        ocp_solver.set(i, 'u', np.array([0.1, 0.5]))

    status = ocp_solver.solve()

    if status != 0:
        ocp_solver.print_statistics()
        raise Exception(f'acados returned status {status}.')

    # get solution
    simX = np.zeros((N+1, nx))
    simU = np.zeros((N, nu))
    for i in range(N):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N,:] = ocp_solver.get(N, "x")

    dts = simU[:, 1]

    print("acados solved OCP successfully, creating integrator to simulate the solution")
    print(f"optimal time: {sum(dts)} s")

    # simulate on finer grid
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

    k = 0
    for i in range(N):
        u = simU[i, 0]
        integrator.set("u", np.hstack((u, np.ones(1, ))))

        # set simulation time
        integrator.set("T", dts_fine[i])

        for _ in range(Ns_fine[i]):
            integrator.set("x", simX_fine[k,:])
            status = integrator.solve()
            if status != 0:
                raise Exception(f'acados returned status {status}.')

            simX_fine[k+1,:] = integrator.get("x")
            simU_fine[k, :] = u
            ts_fine[k+1] = ts_fine[k] + dts_fine[i]

            k += 1

    plot_crane_trajectories(ts_fine, simX_fine, simU_fine)


if __name__ == "__main__":
    for use_cython in [True, False]:
        main(use_cython=use_cython)

