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
sys.path.insert(0, '../common')

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
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    # R = 2*np.diag([1e0])
    R = 5*np.diag([1e-4])
    W = scipy.linalg.block_diag(Q, R)
    y_expr = ca.vertcat(model.x, model.u**2)
    if hessian_approx == 'GAUSS_NEWTON':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'

        ocp.cost.W = W
        ocp.cost.W_e = Q

        ocp.model.cost_y_expr = y_expr

        ocp.cost.Vx_e = np.eye(nx)

        ocp.cost.yref  = np.zeros((ny, ))
        ocp.cost.yref_e = np.zeros((ny_e, ))
    else:
        ocp.cost.cost_type = 'EXTERNAL'
        ocp.cost.cost_type_e = 'EXTERNAL'
        ocp.model.cost_expr_ext_cost = 0.5 * (y_expr.T @ W @ y_expr)
        ocp.model.cost_expr_ext_cost_e = 0.5* model.x.T @ Q @ model.x

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
    ocp.solver_options.qp_solver_ric_alg = 0
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
    elif ocp.cost.cost_type == "EXTERNAL":
        path_cost_fun = ca.Function('path_cost_fun', [model.x, model.u], [(T_HORIZON/N) * model.cost_expr_ext_cost])
    else:
        raise NotImplementedError()

    if ocp.cost.cost_type_e == "LINEAR_LS":
        cost_res = ocp.cost.yref_e - ocp.cost.Vx_e @ model.x
        cost_term = .5 * cost_res.T @ ocp.cost.W_e @ cost_res
        terminal_cost_fun = ca.Function('terminal_cost_fun', [model.x], [cost_term])
    elif ocp.cost.cost_type_e == "EXTERNAL":
        terminal_cost_fun = ca.Function('terminal_cost_fun', [model.x], [model.cost_expr_ext_cost_e])
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
    casadi_solver = ca.nlpsol('casadi_solver', 'ipopt', nlp, {'ipopt': {'print_level': 0}, 'print_time': False})

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
    casadi_solver, lag_hess_fun, lbg, ubg = create_casadi_solver(linearized_dynamics=linearized_dynamics, discrete=discrete)

    nval = 100
    # idxp = 0
    # p_max = 1.05

    idxp = 1
    p_max = 0.4
    p_vals = np.linspace(0, p_max, nval)

    du0_dp_values = np.zeros(nval)
    x0 = X0.copy()
    exact_hessian_status = np.zeros(nval)

    # for hessian check
    min_eigv_vals = np.zeros(nval)
    min_abs_eigv_vals = np.zeros(nval)
    min_abs_eigv_proj = np.zeros(nval)
    min_eigv_proj = np.zeros(nval)
    cond_proj_hess_vals = np.zeros(nval)
    min_eigP_vals = np.zeros(nval)
    min_abs_eigP_vals = np.zeros(nval)
    u0_values = np.zeros(nval)
    hess_errors = np.zeros(nval)

    latexify_plot()
    for i, p0 in enumerate(p_vals):
        x0[idxp] = p0
        u0 = acados_ocp_solver_gn.solve_for_x0(x0)
        u0_values[i] = u0

    du0_dp_finite_diff = np.gradient(u0_values, p_vals[1]-p_vals[0])
    for i, p0 in enumerate(p_vals):
        x0[idxp] = p0
        u0 = acados_ocp_solver_gn.solve_for_x0(x0)
        acados_ocp_solver_gn.store_iterate(filename='iterate.json', overwrite=True, verbose=False)
        acados_ocp_solver_exact.load_iterate(filename='iterate.json', verbose=False)
        acados_ocp_solver_exact.set(0, 'u', u0+1e-7)
        acados_ocp_solver_exact.solve_for_x0(x0, fail_on_nonzero_status=False, print_stats_on_failure=False)
        acados_ocp_solver_exact.eval_param_sens(index=idxp)

        exact_hessian_status[i] = acados_ocp_solver_exact.get_stats('qp_stat')[-1]

        residuals = acados_ocp_solver_exact.get_stats("residuals")
        print(f"residuals sensitivity_solver {residuals}")

        du0_dp_values[i] = acados_ocp_solver_exact.get(0, "sens_u")

        # solve with casadi and compare hessians
        nlp_sol = casadi_solver(p=x0, lbg=lbg, ubg=ubg)
        casadi_hess_l = lag_hess_fun(x=nlp_sol['x'], p=x0, lam_f=1.0, lam_g=nlp_sol['lam_g'])['triu_hess_gamma_x_x']
        casadi_hess = ca.triu2symm(ca.triu(casadi_hess_l)).full()
        min_eigv_vals[i], min_abs_eigv_vals[i], hess_errors[i], min_abs_eigv_proj[i], min_eigv_proj[i], cond_proj_hess_vals[i], min_eigP_vals[i], min_abs_eigP_vals[i] = compare_hessian(casadi_hess, acados_ocp_solver_exact)

        K_mat, K_regularized = compute_K(acados_ocp_solver_exact)

        # du0_dp_values[i] = K_regularized[0][idxp]

        if np.abs(du0_dp_values[i] - K_mat[0][idxp]) > 1e-5:
            print(f"K and du0_dp differ too much")
            print(f"{du0_dp_values[i]=} {K_mat[0][idxp]=}")

    max_hess_error = np.max(hess_errors)
    if max_hess_error > 1e-4:
        raise Exception(f"Hessian error {max_hess_error} > 1e-4 when comparing to casadi.")

    solution_sens_mean_diff = np.mean(np.abs(du0_dp_values -du0_dp_finite_diff))
    if solution_sens_mean_diff > 1.0:
        raise Exception(f"Mean of solution sensitivity difference wrt finite differences {solution_sens_mean_diff} > 1.0.")

    # plot_tangents(p_vals, u0_values, du0_dp_values)
    # Finite difference comparison
    fig, axes = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(9.,9.))
    isub = 0
    axes[isub].plot(p_vals, u0_values)
    axes[isub].grid()
    axes[isub].set_ylabel("$u_0$")
    axes[isub].set_xlim([p_vals[0], p_vals[-1]])

    isub += 1
    axes[isub].plot(p_vals, du0_dp_values, '--', label='acados')
    axes[isub].plot(p_vals, du0_dp_finite_diff, ':', label='finite differences')
    axes[isub].grid()
    axes[isub].legend()
    axes[isub].set_ylabel(r'$\partial_p u_0$')

    # isub += 1
    # axes[isub].plot(p_vals, exact_hessian_status)
    # axes[isub].set_ylabel("QP status")
    # axes[isub].grid()

    isub += 1
    axes[isub].plot(p_vals, min_eigv_vals, label='full hess')
    axes[isub].plot(p_vals, min_eigv_proj, label='proj hess')
    axes[isub].plot(p_vals, min_eigP_vals, label='$P$ Riccati')
    axes[isub].set_ylabel("min eigval")
    axes[isub].legend()
    axes[isub].grid()

    isub += 1
    axes[isub].plot(p_vals, cond_proj_hess_vals, label='proj hess')
    axes[isub].legend()
    axes[isub].grid()
    axes[isub].set_yscale('log')
    axes[isub].set_ylabel("cond")

    isub += 1
    axes[isub].plot(p_vals, min_abs_eigv_vals, '--', label='full hess')
    axes[isub].plot(p_vals, min_abs_eigv_proj, label='proj hess')
    axes[isub].plot(p_vals, min_abs_eigP_vals, label='$P$ Riccati')
    axes[isub].set_ylabel("abs eigval")
    axes[isub].set_yscale('log')
    axes[isub].legend()
    axes[isub].grid()
    axes[-1].set_xlabel('p')

    filename = f'sensitivity_experiment_{"lindyn" if linearized_dynamics else "nonlindyn"}.pdf'
    plt.savefig(filename)
    print(f"stored figure in {filename}")

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


def compute_K(solver: AcadosOcpSolver):
    R_mat = solver.get_from_qp_in(0, 'R')
    S_mat = solver.get_from_qp_in(0, 'S')
    A_mat = solver.get_from_qp_in(0, 'A')
    B_mat = solver.get_from_qp_in(0, 'B')

    P_mat = solver.get_from_qp_in(1, 'P')

    M_mat = R_mat + B_mat.T @ P_mat @ B_mat

    eigvals, eigvecs = np.linalg.eigh(M_mat)

    # print("eigenvalues\n", eigvals)

    eigvals_regularized = [ev if np.abs(ev) > 1e-4 else np.sign(ev)*1e-4 for ev in eigvals]

    K_mat = -np.linalg.solve(M_mat, S_mat + B_mat.T @ P_mat @ A_mat)

    P_mat_regularized = eigvecs @ np.diag(eigvals_regularized) @ eigvecs.T
    K_regularized = -np.linalg.solve(P_mat_regularized, S_mat + B_mat.T @ P_mat @ A_mat)
    return K_mat, K_regularized


def run_hessian_comparison(linearized_dynamics=False, discrete=False):
    casadi_solver, lag_hess_fun, lbg, ubg = create_casadi_solver(linearized_dynamics=linearized_dynamics, discrete=discrete)

    acados_ocp_solver_gn = create_solver(hessian_approx='GAUSS_NEWTON', linearized_dynamics=linearized_dynamics, discrete=discrete)
    acados_ocp_solver_exact = create_solver(hessian_approx='EXACT', linearized_dynamics=linearized_dynamics, discrete=discrete)

    x0 = X0.copy()
    x0[1] = 0.1

    nlp_sol = casadi_solver(p=x0, lbg=lbg, ubg=ubg)
    # print(f"{nlp_sol=}")
    casadi_hess_l = lag_hess_fun(x=nlp_sol['x'], p=x0, lam_f=1.0, lam_g=nlp_sol['lam_g'])['triu_hess_gamma_x_x']
    casadi_hess = ca.triu2symm(ca.triu(casadi_hess_l)).full()
    acados_ocp_solver_gn.solve_for_x0(x0)
    acados_ocp_solver_gn.store_iterate(filename='iterate.json', overwrite=True)
    acados_ocp_solver_exact.load_iterate(filename='iterate.json')
    acados_ocp_solver_exact.solve_for_x0(x0, fail_on_nonzero_status=False, print_stats_on_failure=False)

    _, _, _ = compare_hessian(casadi_hess, acados_ocp_solver_exact)

def compare_hessian(casadi_hess, acados_solver: AcadosOcpSolver):
    offset = 0
    min_eigv_total = 1e12
    min_abs_eigv = 1e12
    hess_error_norm_total = 0.0

    for i in range(N+1):
        hess_block_acados = acados_solver.get_hessian_block(i)
        nv = hess_block_acados.shape[0]
        hess_block_casadi = casadi_hess[offset:offset+nv, offset:offset+nv]
        hess_error_norm = np.max(np.abs(hess_block_acados - hess_block_casadi))
        offset += nv
        hess_error_norm_total = max(hess_error_norm, hess_error_norm_total)
        # print(f"hess block {i} error {hess_error_norm:.2e}")

        eigv = np.linalg.eigvals(hess_block_acados)
        min_eigv = np.min(eigv)
        min_eigv_total = min(min_eigv, min_eigv_total)
        min_abs_eigv = min(min_abs_eigv, np.min(np.abs(eigv)))

    # check projected Hessian & P matrices
    min_abs_eig_proj_hess = 1e12
    min_eig_proj_hess = 1e12
    min_eig_P = 1e12
    min_abs_eig_P = 1e12
    cond_proj_hess = 0.
    for i in range(1, N):
        P_mat = acados_solver.get_from_qp_in(i, 'P')
        B_mat = acados_solver.get_from_qp_in(i-1, 'B')
        # Lr: lower triangular decomposition of R within Riccati != R in qp_in!
        Lr = acados_solver.get_from_qp_in(i-1, 'Lr')
        R_ric = Lr @ Lr.T
        proj_hess_block = R_ric + B_mat.T @ P_mat @ B_mat
        eigv = np.linalg.eigvals(proj_hess_block)
        min_eigv = np.min(eigv)
        min_eig_proj_hess = min(min_eigv, min_eig_proj_hess)
        min_abs_eig_proj_hess = min(min_abs_eig_proj_hess, np.min(np.abs(eigv)))
        cond_proj_hess = max(cond_proj_hess, np.linalg.cond(proj_hess_block))
        # P
        eigv = np.linalg.eigvals(P_mat)
        min_eig_P = min(min_eig_P, np.min(eigv))
        min_abs_eig_P = min(min_abs_eig_P, np.min(np.abs(eigv)))

    return min_eigv_total, min_abs_eigv, hess_error_norm_total, min_abs_eig_proj_hess, min_eig_proj_hess, cond_proj_hess, min_eig_P, min_abs_eig_P

if __name__ == "__main__":
    sensitivity_experiment(linearized_dynamics=False, discrete=True)
    # run_hessian_comparison(linearized_dynamics=False, discrete=True)
