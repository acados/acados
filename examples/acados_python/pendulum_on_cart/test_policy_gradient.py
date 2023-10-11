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

from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, AcadosSimSolver, latexify_plot, casadi_length
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
# - exact Hessian CasADi check
# - compare correctness via finite differences -> done
# - check definiteness -> a lot of indefiniteness

def create_ocp_description(hessian_approx, linearized_dynamics=False, discrete=False) -> AcadosOcp:
    ocp = AcadosOcp()

    if discrete:
        ocp.solver_options.integrator_type = 'DISCRETE'
        if linearized_dynamics:
            model: AcadosModel = export_linearized_pendulum_ode_model_with_discrete_rk4(T_HORIZON/N, X0, np.array([1.]))
        else:
            model: AcadosModel = export_pendulum_ode_model_with_discrete_rk4(T_HORIZON/N)

    else:
        ocp.solver_options.integrator_type = 'ERK'
        if linearized_dynamics:
            model: AcadosModel = export_linearized_pendulum(X0, np.array([1.]))
        else:
            model: AcadosModel = export_pendulum_ode_model()


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

    ocp.cost.yref = np.zeros((ny, ))
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
        ocp.solver_options.nlp_solver_step_length = 0.5
        ocp.solver_options.nlp_solver_max_iter = 1
        ocp.solver_options.tol = 1e-14
    # set prediction horizon
    ocp.solver_options.tf = T_HORIZON
    return ocp


def create_solver(hessian_approx, linearized_dynamics=False, discrete=False) -> AcadosOcpSolver:
    ocp = create_ocp_description(hessian_approx, linearized_dynamics=linearized_dynamics, discrete=discrete)

    acados_ocp_solver = AcadosOcpSolver(ocp, json_file=f'solver_f{hessian_approx}.json')

    return acados_ocp_solver



def create_casadi_solver(linearized_dynamics=False, discrete=False):
    ocp = create_ocp_description(hessian_approx='EXACT', linearized_dynamics=linearized_dynamics, discrete=discrete)

    model = ocp.model
    if discrete:
        f_discrete_fun = ca.Function('f_discrete_fun', [model.x, model.u], [model.disc_dyn_expr])
    else:
        raise NotImplementedError()

    if ocp.cost.cost_type == "NONLINEAR_LS":
        cost_res = ocp.model.cost_y_expr - ocp.cost.yref
        cost_term = .5 * (T_HORIZON/N) * (cost_res.T @ ocp.cost.W @ cost_res)
        path_cost_fun = ca.Function('path_cost_fun', [model.x, model.u], [cost_term])
    else:
        raise NotImplementedError()

    if ocp.cost.cost_type_e == "LINEAR_LS":
        cost_res = ocp.cost.yref_e - ocp.cost.Vx_e @ model.x
        cost_term = .5 * cost_res.T @ ocp.cost.W_e @ cost_res
        terminal_cost_fun = ca.Function('terminal_cost_fun', [model.x], [cost_term])
    else:
        raise NotImplementedError()

    # create CasADi NLP
    nlp_vars = []
    nlp_constr = []
    nx = casadi_length(model.x)
    nu = casadi_length(model.u)

    nlp_params = []
    x0 = ca.SX.sym('x0', nx)
    nlp_params += [x0]
    cost = 0.0

    for k in range(N):
        xk = ca.SX.sym('xk', nx)
        uk = ca.SX.sym('uk', nu)
        nlp_vars += [uk, xk]
        # dynamics
        if k == 0:
            nlp_constr += [xk - x0]
        else:
            nlp_constr += [xnext - xk]
        xnext = f_discrete_fun(xk, uk)
        cost += path_cost_fun(xk, uk)

    # terminal
    xk = ca.SX.sym('xk', nx)
    nlp_constr += [xnext - xk]
    nlp_vars += [xk]
    cost += terminal_cost_fun(xk)

    constr_vec = ca.vertcat(*nlp_constr)
    nlp_x = ca.vertcat(*nlp_vars)
    nlp = {
            'x': nlp_x,
            'f': cost,
            'g': constr_vec,
            'p': ca.vertcat(*nlp_params),
           }
    casadi_solver = ca.nlpsol('casadi_solver', 'ipopt', nlp)

    # exact hessian
    ng = casadi_length(constr_vec)
    # lam_g = ca.SX.sym('lam_g', ng)
    # lagrangian = cost + lam_g.T @ constr_vec
    # lag_hess_expr, _ = ca.hessian(lagrangian, nlp_x)
    # lag_hess_fun = ca.Function('lag_hess_fun', [nlp_x, lam_g], [lag_hess_expr])
    lag_hess_fun = casadi_solver.get_function('nlp_hess_l')

    lbg = np.zeros(ng)
    ubg = np.zeros(ng)
    return casadi_solver, lag_hess_fun, lbg, ubg



def sensitivity_experiment(linearized_dynamics=False, discrete=False, show=True):
    acados_ocp_solver_exact = create_solver(hessian_approx='EXACT', linearized_dynamics=linearized_dynamics, discrete=discrete)
    acados_ocp_solver_gn = create_solver(hessian_approx='GAUSS_NEWTON', linearized_dynamics=linearized_dynamics, discrete=discrete)

    nval = 250
    # idxp = 0
    # p_max = 1.05

    idxp = 1
    p_max = 0.4
    p_vals = np.linspace(0, p_max, nval)
    u0_values = np.zeros(nval)
    du0_dp_values = np.zeros(nval)
    x0 = X0.copy()
    exact_hessian_status = 0.0*p_vals

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
        exact_hessian_status[i] = acados_ocp_solver_exact.get_stats('qp_stat')[-1]

        # check_hessian_last_qp(acados_ocp_solver_exact)

        residuals = acados_ocp_solver_exact.get_stats("residuals")
        print(f"residuals sensitivity_solver {residuals}")

        du0_dp_values[i] = acados_ocp_solver_exact.get(0, "sens_u")

    # plot_tangents(p_vals, u0_values, du0_dp_values)

    # Finite difference comparison
    du0_dp_finite_diff = np.gradient(u0_values, p_vals[1]-p_vals[0])

    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
    axes[0].plot(p_vals, u0_values)
    axes[0].grid()
    axes[0].set_ylabel("$u_0$")
    axes[1].plot(p_vals, du0_dp_values, '--', label='acados')
    axes[1].plot(p_vals, du0_dp_finite_diff, ':', label='finite differences')
    axes[2].set_xlabel('p')
    axes[1].grid()
    axes[1].legend()
    axes[1].set_ylabel(r'$\partial_p u_0$')

    axes[2].plot(p_vals, exact_hessian_status)
    axes[2].set_ylabel("QP status")
    axes[2].grid()

    if show:
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

def check_hessian_last_qp(solver: AcadosOcpSolver):
    for i in range(N+1):
        hess_block = get_hessian_block(solver, i)
        eig_vals, _ = np.linalg.eig(hess_block)
        if any(eig_vals < 0.0):
            print(f"found eig < 0 at node {i} {eig_vals=}")

def get_hessian_block(solver, i) -> np.ndarray:
    Q_mat = solver.get_from_qp_in(i, 'Q')
    R_mat = solver.get_from_qp_in(i, 'R')
    S_mat = solver.get_from_qp_in(i, 'S')
    hess_block = scipy.linalg.block_diag(
        R_mat, Q_mat
    )
    nu = R_mat.shape[0]
    hess_block[nu:, :nu] = S_mat
    hess_block[:nu, nu:] = S_mat.T
    return hess_block

def compare_acados_casadi_hessians(linearized_dynamics=False, discrete=False):
    casadi_solver, lag_hess_fun, lbg, ubg = create_casadi_solver(linearized_dynamics=linearized_dynamics, discrete=discrete)

    acados_ocp_solver_gn = create_solver(hessian_approx='GAUSS_NEWTON', linearized_dynamics=linearized_dynamics, discrete=discrete)
    acados_ocp_solver_exact = create_solver(hessian_approx='EXACT', linearized_dynamics=linearized_dynamics, discrete=discrete)

    x0 = X0.copy()
    x0[1] = 0.1

    nlp_sol = casadi_solver(p=x0, lbg=lbg, ubg=ubg)
    # print(f"{nlp_sol=}")
    casadi_hess = lag_hess_fun(x=nlp_sol['x'], p=x0, lam_f=1.0, lam_g=nlp_sol['lam_g'])['triu_hess_gamma_x_x'].full()

    acados_ocp_solver_gn.solve_for_x0(x0)
    acados_ocp_solver_gn.store_iterate(filename='iterate.json', overwrite=True)
    acados_ocp_solver_exact.load_iterate(filename='iterate.json')
    acados_ocp_solver_exact.solve_for_x0(x0, fail_on_nonzero_status=False, print_stats_on_failure=False)

    offset = 0
    for i in range(N+1):
        hess_block_acados = get_hessian_block(acados_ocp_solver_exact, i)
        nv = hess_block_acados.shape[0]
        hess_block_casadi = casadi_hess[offset:offset+nv, offset:offset+nv]
        hess_error_norm = np.max(np.abs(hess_block_acados - hess_block_casadi))
        print(f"hess block {i} error {hess_error_norm:.2e}")
        if hess_error_norm > 0:
            print(f"diff\n{hess_block_acados - hess_block_casadi}")
            print(f"\nhess acados\n{hess_block_acados}")
            breakpoint()

if __name__ == "__main__":
    # sensitivity_experiment(linearized_dynamics=False, discrete=True)
    compare_acados_casadi_hessians(linearized_dynamics=True, discrete=True)
