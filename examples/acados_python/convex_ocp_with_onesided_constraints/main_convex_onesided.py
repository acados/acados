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

from typing import List
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, latexify_plot, ACADOS_INFTY, AcadosOcpFlattenedIterate, AcadosOcpIterates
import numpy as np
from matplotlib import pyplot as plt
import casadi as cs
from scipy.linalg import block_diag

## This example demonstrates the convergence of the solver for a convex double integrator OCP.
# The dynamics are linear, the cost is quadratic and the constraints are convex.
# The constraints can be formulated as nonlinear constraints (BGH) and as convex-over-nonlinear constraints (BGP).
# The BGP formulation results in an exact Hessian for the constraints in this example.
# The BGH formulation can be used with an exact Hessian or a Gauss-Newton Hessian approximation.
# Using the BGH formulation with a Gauss-Newton Hessian approximation requires globalization and does not converge well in this example.

def export_double_integrator_model():

    model_name = 'double_integrator_2d'
    nx = 4

    # define states
    pos = cs.SX.sym('pos', 2)
    vel = cs.SX.sym('vel', 2)
    states = cs.vertcat(pos, vel)

    # define controls
    acc = cs.SX.sym('acc', 2)
    controls = cs.vertcat(acc)

    # define dynamics expression
    states_dot = cs.SX.sym('xdot', nx)

    ## dynamics
    f_expl = cs.vertcat(vel, acc)
    f_impl = states_dot - f_expl

    model = AcadosModel()
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = states
    model.xdot = states_dot
    model.u = controls
    model.name = model_name

    return model


def solve_ocp(modification=1, constraint_formulation="BGH", hessian_approx="EXACT", globalization = 'FIXED_STEP', plot=False) -> tuple[AcadosOcpIterates, float]:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_double_integrator_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    Tf = 5.0
    N = 50

    Xi_0 = [0,0] # Initial position
    Vi_0 = [-4,6] #  Initial velocity
    u_init = np.array([-0.0, -0.0])

    qp_solver_iter_max = 100
    qp_solver = 'PARTIAL_CONDENSING_HPIPM'

    if modification == 1:
        # Change initial position
        Xi_0 = [-10, 0]
    elif modification == 2:
        # Change initial velocity
        Vi_0 = [-3,6]
    elif modification == 3:
        u_init = np.array([-1.0, -1.0])
    elif modification == 4:
        Tf = 10.0
    # other QP solver settings
    elif modification == 5:
        # fails with infeasible QP
        qp_solver_iter_max = 200
        qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    elif modification == 6:
        # works fine
        qp_solver_iter_max = 200
        qp_solver = 'FULL_CONDENSING_DAQP'
    elif modification == 7:
        # works fine
        qp_solver_iter_max = 200
        qp_solver = 'PARTIAL_CONDENSING_OSQP'
    elif modification == 8:
        # works fine
        qp_solver_iter_max = 200
        qp_solver = 'PARTIAL_CONDENSING_QPOASES'
    elif modification == 9:
        # works fine
        qp_solver_iter_max = 200
        qp_solver = 'FULL_CONDENSING_HPIPM'

    x0 = np.array(Xi_0 + Vi_0)

    ###########################################################################
    # Define cost
    ###########################################################################

    cost_type = "EXTERNAL"
    # cost_type = "LINEAR_LS"
    cost_type = "LLS_SMALL"
    W_x = cs.diag(cs.vertcat(100, 100, 0, 0))
    W_u = 1.* cs.diag(cs.vertcat(1, 1))

    P_des = cs.vertcat(100, -50) # Desired Position
    V_des = cs.vertcat(0, 0) # Desired velocity, but no weight is applied to it
    if cost_type == "EXTERNAL":
        X_des = cs.vertcat(P_des, V_des)

        cost_x = (X_des - model.x).T @ W_x @ (X_des - model.x)

        ocp.cost.cost_type = 'EXTERNAL'
        ocp.model.cost_expr_ext_cost = cost_x + model.u.T @ W_u @ model.u
        ocp.cost.cost_type_e = 'EXTERNAL'
        ocp.model.cost_expr_ext_cost_e = cost_x
    elif cost_type == "LINEAR_LS":
        ocp.cost.cost_type = 'LINEAR_LS'
        ny = nx + nu
        ocp.cost.Vx = np.zeros((ny, nx))
        ocp.cost.Vx[:nx, :] = np.eye(nx)
        ocp.cost.Vu = np.zeros((ny, nu))
        ocp.cost.Vu[nx:, :] = np.eye(nu)
        ocp.cost.W = 2 * block_diag(W_x.full(), W_u.full())
        ocp.cost.Vx_e = np.eye(nx)
        ocp.cost.W_e = 2 * W_x.full()
        ocp.cost.yref = np.concatenate((P_des, V_des, np.zeros((nu, 1)))).flatten()
        ocp.cost.yref_e = np.concatenate((P_des, V_des)).flatten()
    elif cost_type == "LLS_SMALL":
        ocp.cost.cost_type = 'LINEAR_LS'
        ny = 4
        ocp.cost.Vx = np.zeros((ny, nx))
        ocp.cost.Vx[:2, :2] = np.eye(2)
        ocp.cost.Vu = np.zeros((ny, nu))
        ocp.cost.Vu[2:, :] = np.eye(nu)
        ocp.cost.W = 2 * block_diag(cs.diag(cs.vertcat(100, 100)), W_u.full())
        ocp.cost.yref = np.concatenate((P_des, np.zeros((nu, 1)))).flatten()
        ocp.cost.Vx_e = ocp.cost.Vx[:2, :]
        ocp.cost.W_e = 2 * cs.diag(cs.vertcat(100, 100)).full()
        ocp.cost.yref_e = P_des.full().flatten()
    else:
        raise Exception(f'Unknown cost type {cost_type}.')

    ###########################################################################
    # Define constraints
    ###########################################################################

    # Initial conditions
    ocp.constraints.x0 = x0

    # Nonlinear Constraints
    max_velocity_xy = 10
    max_acc_xy = 3
    eps = 1e-16

    lh = np.array([-ACADOS_INFTY, -ACADOS_INFTY])  # ACADOS_INFTY means corresponding constraints are ignored in acados.
    uh = np.array([max_velocity_xy, max_acc_xy])
    if constraint_formulation.startswith("BGH"):
        velocity_and_acceleration_norms = cs.vertcat(cs.sqrt(cs.sumsqr(ocp.model.x[2:])+eps),
                                                     cs.sqrt(cs.sumsqr(ocp.model.u)+eps))
        ocp.model.con_h_expr = velocity_and_acceleration_norms
        ocp.constraints.uh = uh
        ocp.constraints.lh = lh

        ocp.model.con_h_expr_0 = velocity_and_acceleration_norms
        ocp.constraints.uh_0 = uh
        ocp.constraints.lh_0 = lh

    elif constraint_formulation == "BGH_reverse":
        ocp.model.con_h_expr = -velocity_and_acceleration_norms
        ocp.constraints.uh = -lh
        ocp.constraints.lh = -uh

        ocp.model.con_h_expr_0 = velocity_and_acceleration_norms
        ocp.constraints.uh_0 = -lh
        ocp.constraints.lh_0 = -uh

    elif constraint_formulation == "BGP":
        ocp.model.con_r_in_phi = cs.SX.sym('con_r', 4, 1)
        ocp.model.con_phi_expr = cs.vertcat(cs.sumsqr(ocp.model.con_r_in_phi[:2]),
                                                    cs.sumsqr(ocp.model.con_r_in_phi[2:]))
        ocp.model.con_phi_expr = cs.vertcat(cs.sqrt(cs.sumsqr(ocp.model.con_r_in_phi[:2]) + eps),
                                            cs.sqrt(cs.sumsqr(ocp.model.con_r_in_phi[2:]) + eps))
        ocp.model.con_r_expr = cs.vertcat(ocp.model.x[2:], ocp.model.u)
        ocp.constraints.uphi = uh
        ocp.constraints.lphi = lh

        ocp.model.con_r_in_phi_0 = ocp.model.con_r_in_phi
        ocp.model.con_phi_expr_0 = ocp.model.con_phi_expr
        ocp.model.con_r_expr_0 = ocp.model.con_r_expr
        ocp.constraints.uphi_0 = uh
        ocp.constraints.lphi_0 = lh
        ocp.solver_options.exact_hess_constr = 0
    else:
        raise Exception(f'Unknown constraint formulation {constraint_formulation}.')

    # Define bounds on u around the circle, helps with convergence.
    # Without this, exact hessian solver needs line search.
    ocp.constraints.lbu = -(max_acc_xy + 0.1) * np.ones((nu,))
    ocp.constraints.ubu = (max_acc_xy + 0.1) * np.ones((nu,))
    ocp.constraints.idxbu = np.arange(nu)

    ###########################################################################
    # set solver options
    ###########################################################################
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.hessian_approx = hessian_approx
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 1
    ocp.solver_options.qp_solver_iter_max = qp_solver_iter_max
    ocp.solver_options.nlp_solver_max_iter = 1000
    ocp.solver_options.nlp_solver_type = 'SQP'
    ocp.solver_options.globalization = globalization
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.qp_solver_mu0 = 1e3
    ocp.solver_options.store_iterates = True
    ocp.solver_options.eval_residual_at_max_iter = True
    ocp.solver_options.nlp_solver_ext_qp_res = 1

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N

    ocp_solver = AcadosOcpSolver(ocp)

    sol_X = np.zeros((N+1, nx))
    sol_U = np.zeros((N, nu))

    for i in range(N):
        ocp_solver.set(i, "x", x0)
        ocp_solver.set(i, "u", u_init)
    ocp_solver.set(N, "x", x0)

    # Solve the problem
    ocp_solver.solve()
    ocp_solver.print_statistics()

    qp_res_ineq = ocp_solver.get_stats("qp_res_ineq")
    if qp_res_ineq[-1] > 1e-6:
        raise ValueError(f"qp_res_ineq at last iteration is {qp_res_ineq[-1]}, which is larger than 1e-6.")

    # get solution
    for i in range(N):
        sol_X[i,:] = ocp_solver.get(i, "x")
        sol_U[i,:] = ocp_solver.get(i, "u")
    sol_X[N,:] = ocp_solver.get(N, "x")

    print("Initial state: ", sol_X[0,:])
    print("Initial control: ", sol_U[0,:])
    print("Terminal state: ", sol_X[N,:])
    if plot:
        plot_ocp_solution(ocp, sol_X, sol_U, N)

    iterates = ocp_solver.get_iterates()
    residual = np.max(ocp_solver.get_residuals())
    del ocp_solver
    return iterates, residual

def convergence_plot(ll_iterates: List[AcadosOcpIterates], labels: List[str], solution: AcadosOcpFlattenedIterate):
    latexify_plot()
    plt.figure()
    for iterates, label in zip(ll_iterates, labels):
        flat_iterates = [it.flatten() for it in iterates.iterate_list]
        n_iter = len(flat_iterates)
        error_norm = np.zeros(n_iter)
        for i in range(n_iter):
            iterate = flat_iterates[i]
            error_norm[i] = max(np.linalg.norm(solution.x - iterate.x, np.inf),
                                np.linalg.norm(solution.u - iterate.u, np.inf),)
        plt.plot(error_norm[0:-1], error_norm[1:], '-o', label=label)
    plt.xlabel('$||w^n - w^*||_{\infty}$')
    plt.ylabel('$||w^{n+1} - w^*||_{\infty}$')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.grid()
    plt.show()


def plot_ocp_solution(ocp: AcadosOcp, sol_X, sol_U, N):
    latexify_plot()
    plt.figure()
    plt.plot(sol_X[:,0], sol_X[:,1], '-o')
    plt.title('Position trajectory')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()

    plt.figure()
    plt.title('Acceleration and velocity norm')
    acc_norm = [np.linalg.norm(sol_U[i]) for i in range(N)] + [None]
    vel_norm = [np.linalg.norm(sol_X[i, 2:]) for i in range(N+1)]
    plt.plot(ocp.solver_options.shooting_nodes, acc_norm, '-o', label='acc norm')
    plt.plot(ocp.solver_options.shooting_nodes, vel_norm, '-o', label='vel norm')
    plt.grid()
    plt.legend()

    plt.show()


def main(modification=1):
    settings = [("BGP", "EXACT", "FIXED_STEP"), ("BGH", "EXACT", "FIXED_STEP"), ("BGH_reverse", "EXACT", "FIXED_STEP"), ("BGH", "GAUSS_NEWTON", "MERIT_BACKTRACKING")]
    ll_iterates = []
    residuals = []
    for (contraint_formulation, hessian_approx, globalization) in settings:
        iterates, residual = solve_ocp(modification=modification, constraint_formulation=contraint_formulation, hessian_approx=hessian_approx, globalization=globalization, plot=False)
        ll_iterates.append(iterates)
        residuals.append(residual)
        # checks
        n_iter = len(iterates.iterate_list)
        print(f"Settings: {contraint_formulation}, {hessian_approx}, {globalization}, got {n_iter} iterations.")
        if hessian_approx == "EXACT":
            assert residual < 1e-6, f"Residual {residual} too large."
            assert n_iter < 9, f"Number of iterations {n_iter} too large."
            if contraint_formulation == "BGP":
                sol_bgp = iterates.iterate_list[-1].flatten()
            elif contraint_formulation == "BGH":
                sol_bgh = iterates.iterate_list[-1].flatten()
            elif contraint_formulation == "BGH_reverse":
                sol_bgh_reverse = iterates.iterate_list[-1].flatten()
        else:
            assert n_iter > 1e3, f"Number of iterations {n_iter} too small, expected no convergence."

    assert np.allclose(sol_bgp.x, sol_bgh.x, atol=1e-6), f"Solution BGP and BGH differ."
    assert np.allclose(sol_bgp.u, sol_bgh.u, atol=1e-6), f"Solution BGP and BGH differ."
    assert np.allclose(sol_bgp.x, sol_bgh_reverse.x, atol=1e-6), f"Solution BGP and BGH_reverse differ."
    assert np.allclose(sol_bgp.u, sol_bgh_reverse.u, atol=1e-6), f"Solution BGP and BGH_reverse differ."

    min_residual = min(residuals)
    min_residual_idx = residuals.index(min_residual)
    solution = ll_iterates[min_residual_idx].iterate_list[-1].flatten()

    convergence_plot(ll_iterates, labels=[f"{s[0]}, {s[1]}" for s in settings], solution=solution)


if __name__ == '__main__':
    main()
