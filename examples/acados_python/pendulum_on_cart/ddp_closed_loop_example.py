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
sys.path.insert(0, 'common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from pendulum_model import export_pendulum_ode_model
from utils import plot_pendulum
import numpy as np
import scipy.linalg


def main():
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx
    N_horizon = 20


    # set cost module
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)

    ocp.cost.W_e = Q_mat

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
    x0 = np.array([0.0, .3*np.pi, 0.0, 0.0])
    ocp.constraints.x0 = x0
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.idxbx_e = np.arange(nx)
    ocp.constraints.lbx_e = 0*x0
    ocp.constraints.ubx_e = 0*x0

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'DDP'

    ocp.solver_options.qp_solver_cond_N = N_horizon

    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = Tf

    # create ocp solver
    ocp.solver_options.globalization = "MERIT_BACKTRACKING"
    # ocp.solver_options.globalization_alpha_min = 1e-3
    ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.with_adaptive_levenberg_marquardt = True
    ocp.solver_options.print_level = 1
    ocp.translate_to_feasibility_problem(keep_x0=True, keep_cost=False)
    acados_ocp_solver = AcadosOcpSolver(ocp)

    sim = AcadosSim()
    sim.model = model
    sim.solver_options.integrator_type = ocp.solver_options.integrator_type
    sim.solver_options.T = Tf/N_horizon

    acados_integrator = AcadosSimSolver(sim)

    # initialize solver with linspace
    x_init = np.linspace(x0, np.zeros((nx,)), N_horizon+1)
    for i in range(N_horizon):
        acados_ocp_solver.set(i, "x", x_init[i])

    Nsim = 100
    simX = np.zeros((Nsim+1, nx))
    simU = np.zeros((Nsim, nu))

    simX[0,:] = x0

    # closed loop
    for i in range(Nsim):

        # TODO: this should not be necessary!
        acados_ocp_solver.set(0, "x", simX[i, :])
        # solve ocp and get next control input
        simU[i,:] = acados_ocp_solver.solve_for_x0(simX[i, :], fail_on_nonzero_status=False)
        acados_ocp_solver.print_statistics()
        status = acados_ocp_solver.get_status()
        if status != 0:
            raise Exception(f'acados returned status {status}.')

        # simulate system
        simX[i+1, :] = acados_integrator.simulate(x=simX[i, :], u=simU[i,:])

    # assert np.allclose(simX[-1, :], np.zeros((nx,)), atol=1e-5)

    # plot results
    plot_pendulum(np.linspace(0, Tf/N_horizon*Nsim, Nsim+1), Fmax, simU, simX)

if __name__ == '__main__':
    main()