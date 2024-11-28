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

from acados_template import AcadosOcp, AcadosOcpSolver, latexify_plot
from pendulum_model import export_pendulum_ode_model
import numpy as np
import matplotlib.pyplot as plt
import casadi as ca
import scipy

latexify_plot()

def setup_solver(N: int, dt: float, u_max: float = 60):
    ocp = AcadosOcp()

    model = export_pendulum_ode_model()
    ocp.model = model

    Tf = N*dt
    ocp.solver_options.N_horizon = N

    nx = model.x.rows()
    nu = model.u.rows()

    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-1])

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref = np.zeros((nx+nu, ))
    ocp.cost.yref_e = np.zeros((nx, ))
    ocp.cost.W_e = Q_mat
    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)

    # set constraints
    ocp.constraints.lbu = np.array([-u_max])
    ocp.constraints.ubu = np.array([+u_max])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.nlp_solver_max_iter = 600

    ocp.solver_options.tf = Tf

    return AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')


def main():

    dt = 0.05
    N = 25

    ocp_solver = setup_solver(N, dt)

    nx = ocp_solver.acados_ocp.dims.nx

    num_grid = 50
    thetas = np.linspace(0, 0.2*np.pi, num_grid)
    optimal_value_fun = np.zeros((num_grid))
    optimal_value_grad = np.zeros((num_grid))

    x0 = np.zeros((nx,))

    for n, tau in enumerate(np.linspace(0, 1, N+1)):
        x0[1] = tau*thetas[-1]
        ocp_solver.set(n, 'x', x0)

    # state value function gradient
    for k, theta in enumerate(thetas):
        print(f'Solving OCP for {theta=}')

        x0[1] = theta
        _ = ocp_solver.solve_for_x0(x0)
        optimal_value_fun[k] = ocp_solver.get_cost()
        optimal_value_grad[k] = ocp_solver.eval_and_get_optimal_value_gradient(with_respect_to='initial_state')[1]

    cd_optimal_value_grad = (optimal_value_fun[2:]-optimal_value_fun[:-2])/(thetas[2:]-thetas[:-2])

    assert np.allclose(optimal_value_grad[1:-1], cd_optimal_value_grad, rtol=1e-2, atol=1e-2)


    # state-action value function gradient (aka Q-function)
    u = ocp_solver.get(0, 'u').item()
    us = np.linspace(0.9*u, 1.1*u, num_grid)

    Q_fun = np.zeros((num_grid,))
    Q_grad = np.zeros((num_grid,))

    for k, u0 in enumerate(us):
        print(f'Solving OCP for {u0=}')

        ocp_solver.constraints_set(0, 'lbu', u0)
        ocp_solver.constraints_set(0, 'ubu', u0)
        status = ocp_solver.solve()
        Q_fun[k] = ocp_solver.get_cost()
        Q_grad[k] = ocp_solver.eval_and_get_optimal_value_gradient(with_respect_to='initial_control')[0]

    cd_Q_grad = (Q_fun[2:]-Q_fun[:-2])/(us[2:]-us[:-2])

    assert np.allclose(Q_grad[1:-1], cd_Q_grad, rtol=1e-2, atol=1e-2)

    _, axes = plt.subplots(nrows=4, ncols=1, figsize=(3, 8))

    axes[0].plot(thetas, optimal_value_fun)
    axes[1].plot(thetas, optimal_value_grad, label='exact')

    axes[1].plot(thetas[1:-1], cd_optimal_value_grad, label='central differences')
    axes[0].set_ylabel(r'optimal value $V^*(\bar{x}_0(\theta))$')
    axes[1].set_ylabel(r'$\nabla_{\theta} V^*(\bar{x}_0(\theta))$')
    axes[0].set_xlabel(r'$\theta$')
    axes[1].set_xlabel(r'$\theta$')

    axes[2].plot(us, Q_fun)
    axes[3].plot(us, Q_grad, label='exact')

    axes[3].plot(us[1:-1], cd_Q_grad, label='central differences')
    axes[2].set_ylabel(r'optimal state-action value $Q(\bar{x}_0, \bar{u}_0)$')
    axes[3].set_ylabel(r'$\nabla_{\theta} Q(\bar{x}_0, \bar{u}_0)$')
    axes[2].set_xlabel(r'$u$')
    axes[3].set_xlabel(r'$u$')

    for i in range(4):
        axes[i].grid()
    axes[1].legend()
    axes[3].legend()

    axes[0].set_xlim(thetas[0], thetas[-1])
    axes[1].set_xlim(thetas[0], thetas[-1])
    axes[2].set_xlim(us[0], us[-1])
    axes[3].set_xlim(us[0], us[-1])

    plt.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()