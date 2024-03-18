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

from acados_template import AcadosOcp, AcadosOcpSolver
from utils import plot_furuta_pendulum
from furuta_model import get_furuta_model
from integrator_experiment import setup_acados_integrator
import numpy as np
import scipy.linalg
from casadi import vertcat

def setup(x0, umax, dt_0, N_horizon, Tf, RTI=False):
    ocp = AcadosOcp()

    model = get_furuta_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    ocp.dims.N = N_horizon

    # set cost module
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    Q_mat = np.diag([50., 500., 1., 1.])
    R_mat = np.diag([1e3])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    ocp.model.cost_y_expr = vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    ocp.constraints.lbu = np.array([-umax])
    ocp.constraints.ubu = np.array([+umax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = x0

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'

    # NOTE we use a nonuniform grid!
    ocp.solver_options.time_steps = np.array([dt_0] + [(Tf-dt_0)/(N_horizon-1)]*(N_horizon-1))
    ocp.solver_options.sim_method_num_steps = np.array([1] + [2]*(N_horizon-1))
    ocp.solver_options.levenberg_marquardt = 1e-6
    ocp.solver_options.nlp_solver_max_iter = 10

    ocp.solver_options.nlp_solver_type = 'SQP_RTI' if RTI else 'SQP'
    ocp.solver_options.qp_solver_cond_N = N_horizon

    ocp.solver_options.tf = Tf

    solver_json = 'acados_ocp_' + model.name + '.json'
    ocp_solver = AcadosOcpSolver(ocp, json_file = solver_json)

    # setup plant simulator
    integrator = setup_acados_integrator(model, dt_0, num_stages=2, num_steps=2, integrator_type="IRK",
                            newton_iter=20, newton_tol=1e-10)

    return ocp_solver, integrator


def main(use_RTI=False):

    x0 = np.array([0.0, np.pi, 0.0, 0.0])
    umax = .45

    Tf = .350       # total prediction time
    N_horizon = 8   # number of shooting intervals
    dt_0 = 0.025    # sampling time = length of first shooting interval

    ocp_solver, integrator = setup(x0, umax, dt_0, N_horizon, Tf, use_RTI)

    nx = ocp_solver.acados_ocp.dims.nx
    nu = ocp_solver.acados_ocp.dims.nu

    Nsim = 50
    simX = np.zeros((Nsim+1, nx))
    simU = np.zeros((Nsim, nu))

    simX[0,:] = x0

    if use_RTI:
        t_preparation = np.zeros((Nsim))
        t_feedback = np.zeros((Nsim))

    else:
        t = np.zeros((Nsim))

    # do some initial iterations to start with a good initial guess
    num_iter_initial = 10
    for _ in range(num_iter_initial):
        ocp_solver.solve_for_x0(x0_bar = x0, fail_on_nonzero_status=False)

    # closed loop
    for i in range(Nsim):

        if use_RTI:
            # preparation phase
            ocp_solver.options_set('rti_phase', 1)
            status = ocp_solver.solve()
            t_preparation[i] = ocp_solver.get_stats('time_tot')

            # set initial state
            ocp_solver.set(0, "lbx", simX[i, :])
            ocp_solver.set(0, "ubx", simX[i, :])

            # feedback phase
            ocp_solver.options_set('rti_phase', 2)
            status = ocp_solver.solve()
            t_feedback[i] = ocp_solver.get_stats('time_tot')

            simU[i, :] = ocp_solver.get(0, "u")

        else:
            # solve ocp and get next control input
            simU[i,:] = ocp_solver.solve_for_x0(x0_bar = simX[i, :], fail_on_nonzero_status=False)

            t[i] = ocp_solver.get_stats('time_tot')

        # simulate system
        simX[i+1, :] = integrator.simulate(x=simX[i, :], u=simU[i,:])

    # evaluate timings
    if use_RTI:
        # scale to milliseconds
        t_preparation *= 1000
        t_feedback *= 1000
        print(f'Computation time in preparation phase in ms: \
                min {np.min(t_preparation):.3f} median {np.median(t_preparation):.3f} max {np.max(t_preparation):.3f}')
        print(f'Computation time in feedback phase in ms:    \
                min {np.min(t_feedback):.3f} median {np.median(t_feedback):.3f} max {np.max(t_feedback):.3f}')
    else:
        # scale to milliseconds
        t *= 1000
        print(f'Computation time in ms: min {np.min(t):.3f} median {np.median(t):.3f} max {np.max(t):.3f}')

    # plot results
    plot_furuta_pendulum(np.linspace(0, (Tf/N_horizon)*Nsim, Nsim+1),simX, simU, umax, plt_show=True)

    ocp_solver = None


if __name__ == '__main__':
    main(use_RTI=False)
    main(use_RTI=True)

