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


from acados_template import AcadosModel, AcadosOcp, AcadosMultiphaseOcp, AcadosOcpSolver, casadi_length
import numpy as np
import casadi as ca
from scipy.linalg import block_diag
from polynom_utils import plot_open_loop_trajectory_pwpol_u, export_pendulum_ode_model

def create_ocp_formulation_without_opts(cost_type, degree_u_polynom, explicit_symmetric_penalties=True, penalty_type='L2') -> (AcadosOcp, ca.Function):
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    nu_original = casadi_length(model.u)

    ocp.model: AcadosModel = model
    nx = casadi_length(model.x)

    # set cost
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-1])
    W_mat = block_diag(Q_mat, R_mat)

    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref = np.zeros(nx+nu_original)
    ocp.cost.yref_e = np.zeros(nx)
    if cost_type == 'NONLINEAR_LS':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'
        ocp.cost.W = W_mat
        ocp.cost.W_e = Q_mat
    elif cost_type == 'CONVEX_OVER_NONLINEAR':
        ocp.cost.cost_type = 'CONVEX_OVER_NONLINEAR'
        ocp.cost.cost_type_e = 'CONVEX_OVER_NONLINEAR'
        conl_res = ca.SX.sym('residual_conl', ocp.model.cost_y_expr.shape)
        conl_res_e = ca.SX.sym('residual_conl', ocp.model.cost_y_expr_e.shape)
        ocp.model.cost_r_in_psi_expr = conl_res
        ocp.model.cost_r_in_psi_expr_e = conl_res_e
        ocp.model.cost_psi_expr = .5 * conl_res.T @ W_mat @ conl_res
        ocp.model.cost_psi_expr_e = .5 * conl_res_e.T @ Q_mat @ conl_res_e
    elif cost_type == 'NLS_TO_CONL':
        ocp.cost.cost_type = 'NONLINEAR_LS'
        ocp.cost.cost_type_e = 'NONLINEAR_LS'
        ocp.cost.W = W_mat
        ocp.cost.W_e = Q_mat
        ocp.translate_nls_cost_to_conl()
    else:
        raise ValueError(f'cost_type {cost_type} not supported.')

    # set constraints as penalties
    Fmax = 80
    if penalty_type == 'L2':
        weight = 1e5
        if explicit_symmetric_penalties:
            ocp.formulate_constraint_as_L2_penalty(ocp.model.u, weight, Fmax, -Fmax)
        else:
            ocp.formulate_constraint_as_L2_penalty(ocp.model.u, weight, Fmax, None)
            ocp.formulate_constraint_as_L2_penalty(ocp.model.u, weight, None, -Fmax)
    elif penalty_type == 'Huber':
        # combined Huber and L2 penalty
        ocp.solver_options.nlp_solver_max_iter = 200
        weight = 1e0
        if explicit_symmetric_penalties:
            ocp.formulate_constraint_as_Huber_penalty(ocp.model.u, weight, Fmax, -Fmax, huber_delta=1e0, use_xgn=True)
            ocp.formulate_constraint_as_L2_penalty(ocp.model.u, 1e5, Fmax, -Fmax)
        else:
            raise NotImplementedError('Huber penalty not implemented for non-explicit symmetric penalties.')

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    evaluate_polynomial_u_fun = model.augment_model_with_polynomial_control(degree_u_polynom)

    return ocp, evaluate_polynomial_u_fun



def create_ocp_solver(cost_type, N_horizon, degree_u_polynom,
                      explicit_symmetric_penalties=True,
                      penalty_type='L2', nlp_solver_max_iter=None) -> (AcadosOcpSolver, ca.Function):
    ocp: AcadosOcp
    ocp, evaluate_polynomial_u_fun = create_ocp_formulation_without_opts(cost_type, degree_u_polynom, explicit_symmetric_penalties=explicit_symmetric_penalties, penalty_type=penalty_type)

    ocp.dims.N = N_horizon

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP, FULL_CONDENSING_DAQP
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'IRK'
    # ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.cost_discretization = 'INTEGRATOR'

    if nlp_solver_max_iter is not None:
        ocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter

    # set prediction horizon
    T_horizon = 1.0
    n_short = 1
    dt_short = 0.02
    n_long = N_horizon-n_short
    ocp.solver_options.time_steps = np.array(n_short * [dt_short] + n_long * [(T_horizon - dt_short * n_short) / n_long])
    ocp.solver_options.tf = T_horizon

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    return ocp_solver, evaluate_polynomial_u_fun


def create_mocp_solver(cost_type, N_list, degrees_u_polynom, explicit_symmetric_penalties=True, penalty_type='L2', nlp_solver_max_iter=None) -> (AcadosOcpSolver, ca.Function):

    N_horizon = sum(N_list)

    mocp = AcadosMultiphaseOcp(N_list)
    polynomial_u_funs = []

    for i, degree_u_polynom in enumerate(degrees_u_polynom):
        ocp, evaluate_polynomial_u_fun = create_ocp_formulation_without_opts(cost_type, degree_u_polynom, explicit_symmetric_penalties=explicit_symmetric_penalties, penalty_type=penalty_type)

        polynomial_u_funs.append(evaluate_polynomial_u_fun)
        mocp.set_phase(ocp, i)

    # set options
    mocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP, FULL_CONDENSING_DAQP
    mocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    mocp.solver_options.integrator_type = 'IRK'
    # mocp.solver_options.print_level = 1
    mocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    mocp.solver_options.cost_discretization = 'INTEGRATOR'
    if nlp_solver_max_iter is not None:
        mocp.solver_options.nlp_solver_max_iter = nlp_solver_max_iter

    # set prediction horizon
    T_horizon = 1.0
    n_short = 1
    dt_short = 0.02
    n_long = N_horizon-n_short
    mocp.solver_options.time_steps = np.array(n_short * [dt_short] + n_long * [(T_horizon - dt_short * n_short) / n_long])
    mocp.solver_options.tf = T_horizon

    mocp.solver_options.sim_method_num_stages = np.array(n_short * [2] + n_long * [4])

    ocp_solver = AcadosOcpSolver(mocp, verbose=True)

    return ocp_solver, polynomial_u_funs



def main_mocp(cost_type='NONLINEAR_LS', explicit_symmetric_penalties=True, penalty_type='L2'):
    N_list = [1, 5]
    N_horizon = sum(N_list)
    degrees_u_polynom = [0, 4]

    # create solver and extract
    ocp_solver, polynomial_u_funs = create_mocp_solver(cost_type, N_list, degrees_u_polynom, explicit_symmetric_penalties=explicit_symmetric_penalties, penalty_type=penalty_type)
    ocp = ocp_solver.acados_ocp
    model = ocp.model[0]
    nu_original = model.nu_original
    nx = casadi_length(model.x)

    # solve OCP
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    cost_val = ocp_solver.get_cost()
    print(f"found optimal cost: {cost_val:.4e}")

    # get solution
    simX = np.ndarray((N_horizon+1, nx))
    simU = []
    for i in range(N_horizon):
        simX[i,:] = ocp_solver.get(i, "x")
        simU.append(ocp_solver.get(i, "u"))
    simX[N_horizon,:] = ocp_solver.get(N_horizon, "x")

    # compute U fine:
    n_pol_eval = 20
    U_fine = []
    i_ocp = 0
    for i_phase in range(len(N_list)):
        for _ in range(N_list[i_phase]):
            dt = ocp.solver_options.time_steps[i_ocp]
            u_coeff = simU[i_ocp]
            U_interval = np.zeros((n_pol_eval+1, nu_original))
            for j in range(n_pol_eval+1):
                t = (j/n_pol_eval) * dt
                U_interval[j, :] = polynomial_u_funs[i_phase](u_coeff, t)
            U_fine.append(U_interval)
            i_ocp += 1
    ocp_solver.store_iterate('mocp_sol.json', overwrite=True)

    # plot
    plot_open_loop_trajectory_pwpol_u(ocp.solver_options.shooting_nodes, simX, U_fine, title=f'controls: piecewise polynomials with degree {degrees_u_polynom}')


def main_ocp(cost_type='NONLINEAR_LS', explicit_symmetric_penalties=True, penalty_type='L2'):
    N_horizon = 100
    degree_u_polynom = 0

    # create solver and extract
    ocp_solver, evaluate_polynomial_u_fun = create_ocp_solver(cost_type, N_horizon, degree_u_polynom, explicit_symmetric_penalties=explicit_symmetric_penalties, penalty_type=penalty_type)
    ocp = ocp_solver.acados_ocp
    model = ocp.model
    nu_original = model.nu_original
    nu = casadi_length(model.u)
    nx = casadi_length(model.x)

    # solve OCP
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    cost_val = ocp_solver.get_cost()
    print(f"found optimal cost: {cost_val:.4e}")

    # get solution
    simX = np.ndarray((N_horizon+1, nx))
    simU = np.ndarray((N_horizon, nu))
    for i in range(N_horizon):
        simX[i,:] = ocp_solver.get(i, "x")
        simU[i,:] = ocp_solver.get(i, "u")
    simX[N_horizon,:] = ocp_solver.get(N_horizon, "x")

    # compute U fine:
    n_pol_eval = 20
    U_fine = []
    for i, dt in enumerate(ocp.solver_options.time_steps):
        u_coeff = simU[i, :]
        U_interval = np.zeros((n_pol_eval+1, nu_original))
        for j in range(n_pol_eval+1):
            t = (j/n_pol_eval) * dt
            U_interval[j, :] = evaluate_polynomial_u_fun(u_coeff, t)
        U_fine.append(U_interval)

    ocp_solver.store_iterate('ocp_sol.json', overwrite=True)

    # plot
    plot_open_loop_trajectory_pwpol_u(ocp.solver_options.shooting_nodes, simX, U_fine, title=f'controls: piecewise polynomials of degree {degree_u_polynom}')


if __name__ == '__main__':
    # main(cost_type="CONVEX_OVER_NONLINEAR", explicit_symmetric_penalties=True)
    # main_ocp(cost_type="CONVEX_OVER_NONLINEAR", explicit_symmetric_penalties=True, penalty_type='Huber')
    main_ocp(cost_type="CONVEX_OVER_NONLINEAR", explicit_symmetric_penalties=True, penalty_type='L2')
    main_mocp(cost_type="CONVEX_OVER_NONLINEAR", explicit_symmetric_penalties=True, penalty_type='L2')
