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

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, latexify_plot
from pendulum_model import export_pendulum_ode_model, export_linearized_pendulum, export_pendulum_ode_model_with_discrete_rk4, export_linearized_pendulum_ode_model_with_discrete_rk4

from utils import plot_pendulum
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import casadi as ca

X0 = np.array([0.5, 0.0001, 0.00001, 0.00001])
FMAX = 80
T_HORIZON = 1.0
N = 20

# TODO:
# - compare correctness via finite differences
# - exact Hessian CasADi check
# - check definiteness

def create_solver_and_integrator(hessian_approx, linearized_dynamics=False, discrete=False):
    ocp = AcadosOcp()

    ocp.solver_options.integrator_type = 'DISCRETE'

    if linearized_dynamics:
        if discrete:
            model = export_linearized_pendulum_ode_model_with_discrete_rk4(T_HORIZON/N, X0, np.array([1.]))
        else:
            model = export_linearized_pendulum(X0, np.array([1.]))
            ocp.solver_options.integrator_type = 'ERK'

    else:
        if discrete:
            model = export_pendulum_ode_model_with_discrete_rk4(T_HORIZON/N)
        else:
            model = export_pendulum_ode_model()
            ocp.solver_options.integrator_type = 'ERK'


    model.name += f'{hessian_approx}'

    ocp.model = model

    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = nx + nu
    ny_e = nx
    ocp.dims.N = N

    # set cost module
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 5*np.diag([1e-4])
    # R = 2*np.diag([1e0])

    ocp.cost.W = scipy.linalg.block_diag(Q, R)
    ocp.cost.W_e = Q

    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u**2)

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    # ocp.constraints.lbu = np.array([-FMAX])
    # ocp.constraints.ubu = np.array([+FMAX])
    # ocp.constraints.idxbu = np.array([0])
    ocp.constraints.x0 = X0

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = hessian_approx
    # ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI
    ocp.solver_options.sim_method_num_steps = 5
    ocp.solver_options.tol = 1e-6
    # ocp.solver_options.sim_method_newton_tol = 1e-5

    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.qp_solver_warm_start = 0

    ocp.solver_options.qp_solver_iter_max = 500
    ocp.solver_options.nlp_solver_max_iter = 1000
    # ocp.solver_options.globalization = 'MERIT_BACKTRACKING'
    if hessian_approx == 'EXACT':
        ocp.solver_options.nlp_solver_step_length = 0.0
        ocp.solver_options.nlp_solver_max_iter = 1
        ocp.solver_options.tol = 1e-10
    # set prediction horizon
    ocp.solver_options.tf = T_HORIZON

    acados_ocp_solver = AcadosOcpSolver(ocp, json_file=f'solver_f{hessian_approx}.json')
    # acados_integrator = AcadosSimSolver(ocp)
    return acados_ocp_solver



def sensitivity_experiment(linearized_dynamics=False, discrete=False):
    acados_ocp_solver_exact = create_solver_and_integrator(hessian_approx='EXACT', linearized_dynamics=linearized_dynamics, discrete=discrete)
    acados_ocp_solver_gn = create_solver_and_integrator(hessian_approx='GAUSS_NEWTON', linearized_dynamics=linearized_dynamics, discrete=discrete)

    nval = 250
    # idxp = 0
    # p_max = 1.05

    idxp = 1
    p_max = 0.3
    p_vals = np.linspace(0, p_max, nval)
    u0_values = np.zeros(nval)
    du0_dp_values = np.zeros(nval)
    x0 = X0.copy()

    latexify_plot()

    for i, p0 in enumerate(p_vals):
        x0[idxp] = p0
        u0 = acados_ocp_solver_gn.solve_for_x0(x0)
        u0_values[i] = u0

        print(f'solve sens for {p0=}')
        acados_ocp_solver_gn.store_iterate(filename='iterate.json', overwrite=True)
        acados_ocp_solver_exact.load_iterate(filename='iterate.json')
        acados_ocp_solver_exact.set(0, 'u', u0+1e-7)
        acados_ocp_solver_exact.solve_for_x0(x0, fail_on_nonzero_status=False, print_stats_on_failure=False)
        acados_ocp_solver_exact.eval_param_sens(index=idxp)

        residuals = acados_ocp_solver_exact.get_stats("residuals")
        print(f"residuals sensitivity_solver {residuals}")

        du0_dp_values[i] = acados_ocp_solver_exact.get(0, "sens_u")

    # plot_tangents(p_vals, u0_values, du0_dp_values)

    # Finite difference comparison
    du0_dp_finite_diff = np.gradient(u0_values, p_vals[1]-p_vals[0])

    plt.figure()
    plt.plot(p_vals, du0_dp_values, '--', label='acados')
    plt.plot(p_vals, du0_dp_finite_diff, ':', label='finite differences')
    plt.xlabel('p')
    plt.grid()
    plt.legend()
    plt.ylabel(r'$\partial_p u_0$')
    plt.show()

    err = np.max(np.abs(du0_dp_finite_diff - du0_dp_values))
    rel_err = np.max(np.abs(du0_dp_finite_diff - du0_dp_values) /
                     np.abs(du0_dp_finite_diff))
    print(f"max error acados vs finite differences: abs: {err:.2e} rel {rel_err:.2e}")


def plot_tangents(p_vals, u0_values, du0_dp_values):
    plt.figure()
    for i, p0 in enumerate(p_vals):
        if i % 50 == 25:
            u0 = u0_values[i]
            taylor_0 = u0 + du0_dp_values[i] * (p_vals[0] - p0)
            taylor_1 = u0 + du0_dp_values[i] * (p_vals[-1] - p0)
            plt.scatter(p0, u0, marker='*', color="C1")
            plt.plot([p_vals[0], p_vals[-1]], [taylor_0, taylor_1], color="C1", alpha=0.3)

    plt.plot(p_vals, u0_values, ':', color="C0")
    plt.xlabel('$p$')
    plt.ylabel('$u_0$')
    plt.grid()

if __name__ == "__main__":
    sensitivity_experiment(linearized_dynamics=False, discrete=True)
