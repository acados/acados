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
import casadi as ca
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, plot_trajectories, plot_convergence
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
import matplotlib.pyplot as plt

VARIANTS = ['exact Hess', 'inexact Hess', 'zero-order']

def export_pendulum_ode_model_with_discrete_dyn(dT, with_custom_jacobian=True, with_custom_hessian=False) -> AcadosModel:

    model: AcadosModel = export_pendulum_ode_model()

    x = model.x
    u = model.u

    ode = ca.Function('ode', [x, u], [model.f_expl_expr])
    # set up RK4
    integrator_type = 'RK4'
    if integrator_type == 'Euler':
        xf = x + dT * ode(x, u)
    elif integrator_type == 'RK4':
        k1 = ode(x, u)
        k2 = ode(x+dT/2*k1,u)
        k3 = ode(x+dT/2*k2,u)
        k4 = ode(x+dT*k3,  u)
        xf = x + dT/6 * (k1 + 2*k2 + 2*k3 + k4)

    model.disc_dyn_expr = xf

    xss = np.zeros((model.x.size()[0],))
    uss = np.zeros((model.u.size()[0],))
    if with_custom_jacobian:
        disc_dyn_jac_ux_fun = ca.Function('disc_dyn_jac_ux_fun', [x, u], [ca.jacobian(xf, ca.vertcat(u, x))])
        disc_dyn_jac_ux_fun_evaluated = disc_dyn_jac_ux_fun(xss, uss)
        model.disc_dyn_custom_jac_ux_expr = disc_dyn_jac_ux_fun_evaluated
    if with_custom_hessian:
        model.pi = ca.SX.sym('pi', model.x.size()[0])
        ux = ca.vertcat(u, x)
        adj_ux = ca.jtimes(model.disc_dyn_expr, ux, model.pi, True)
        model.disc_dyn_custom_hess_ux_expr = 0.7 * ca.jacobian(adj_ux, ux) # make this inexact

    return model

def create_ocp_solver(variant: str) -> AcadosOcpSolver:
    if variant not in VARIANTS:
        raise ValueError("Unknown variant")

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    Tf = 1.0
    N_horizon = 20

    # set model
    integrator_type = 'DISCRETE'
    if integrator_type == 'DISCRETE':
        with_custom_jacobian = (variant == 'zero-order')
        with_custom_hessian = (variant != 'exact Hess')
        model = export_pendulum_ode_model_with_discrete_dyn(Tf/N_horizon,
                                            with_custom_jacobian=with_custom_jacobian,
                                            with_custom_hessian=with_custom_hessian)
    else:
        raise NotImplementedError("This example only supports DISCRETE integrator type.")

    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # set dimensions
    ocp.solver_options.N_horizon = N_horizon

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
    ocp.constraints.x0 = np.array([0.0, 0.1 * np.pi, 0.0, 0.0])
    ocp.constraints.idxbu = np.array([0])

    ocp.solver_options.integrator_type = integrator_type
    ocp.solver_options.print_level = 1
    # ocp.solver_options.with_anderson_acceleration = True
    # ocp.solver_options.anderson_activation_threshold = anderson_activation_threshold

    if variant in ['exact Hess', 'inexact Hess']:
        ocp.solver_options.hessian_approx = 'EXACT'
        # ocp.solver_options.regularize_method = 'CONVEXIFY'
    else:
        ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.nlp_solver_type = 'SQP'

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    return ocp_solver

def main():
    solutions = []
    labels = []
    kkt_norm_list = []

    for variant in VARIANTS:
        print(f'Solving variant: {variant}')

        labels.append(variant)

        ocp_solver = create_ocp_solver(variant)
        ocp = ocp_solver.acados_ocp

        status = ocp_solver.solve()

        if status != 0:
            raise Exception(f'acados returned status {status}.')

        # get solution
        sol = ocp_solver.get_iterate()
        solutions.append(sol)

        res_all = ocp_solver.get_stats("res_all")
        kkt_norms = np.linalg.norm(res_all, axis=1)
        kkt_norm_list.append(kkt_norms)

        del ocp_solver

    idx_plot_traj = [0, 1, 2]
    idx_plot_conv = [0, 1, 2]

    traj_fig_filename = None
    conv_fig_filename = None

    # tests:
    assert len(kkt_norm_list[0]) < len(kkt_norm_list[1]), "exact Hess should converge faster than inexact Hess"
    assert solutions[0].allclose(solutions[1], atol=1e-6), "exact Hess and inexact Hess solutions should be close"
    assert not solutions[0].allclose(solutions[2], atol=1e-3), "zero-order and exact Hess solutions should not be close"

    plot_trajectories(
        x_traj_list=[np.array(solutions[i].x) for i in idx_plot_traj],
        u_traj_list=[np.array(solutions[i].u) for i in idx_plot_traj],
        linestyle_list=['-', '--', ':', '-.'],
        labels_list=[labels[i] for i in idx_plot_traj],
        idxpx=[1],
        x_labels=ocp.model.x_labels,
        u_labels=ocp.model.u_labels,
        time_traj_list=[ocp.solver_options.shooting_nodes for _ in idx_plot_traj],
        idxbu=ocp.constraints.idxbu,
        lbu=ocp.constraints.lbu,
        ubu=ocp.constraints.ubu,
        show_plot=False,
        single_column=True,
        bbox_to_anchor=(.7, 0.),
        # figsize=(3., 2.5),
        ncol_legend=1,
        color_list=['C2', 'C3', 'C4'],
        fig_filename=traj_fig_filename,
    )


    plot_convergence(
        [kkt_norm_list[i] for i in idx_plot_conv],
        [labels[i] for i in idx_plot_conv],
        show_plot=False,
        # figsize=(2.8, 3.5),
        fig_filename=conv_fig_filename,
        )
    plt.show()


if __name__ == "__main__":
    main()
