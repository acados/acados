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

sys.path.insert(0, '../getting_started')

from matplotlib import pyplot as plt

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim, ACADOS_INFTY
from acados_template.mpc_utils import AcadosCostConstraintEvaluator
from pendulum_model import export_pendulum_ode_model
from utils_eval import plot_pendulum_eval
import numpy as np
import scipy.linalg
import math
from casadi import vertcat, SX

# Define constraints to make evaluation more challenging
constraint_par = {'omega_dot_min_1': -4,
                  'omega_dot_min_2': -6,
                  'omega_dot_max_1': 1.,
                  'omega_dot_max_2': 2.,
                  'iter_omega_change': 6,
                  'v_max': 5,
                  'F_max': 80}


def setup(x0, N_horizon, Tf, td, RTI=False, parametric_constraints=False):
    # create ocp object to formulate the OCP
    global constraint_par
    ocp = AcadosOcp()

    # set model
    ocp.model = model = export_pendulum_ode_model()

    if parametric_constraints:
        omega_dot_min = SX.sym('omega_dot_min')
        omega_dot_max = SX.sym('omega_dot_max')

        ocp.model.p = omega_dot_min
        ocp.model.p_global = omega_dot_max
        ocp.model.con_h_expr = vertcat(omega_dot_min - model.x[3],  # assuming ub=0, lb=-inf
                                       model.x[3] - omega_dot_max)  # assuming ub=0, lb=-inf
        ocp.parameter_values = np.array([constraint_par['omega_dot_min_1']])
        ocp.p_global_values = np.array([constraint_par['omega_dot_max_1']])

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    # set cost module
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    Q_mat = 2 * np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2 * np.diag([1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat * 10

    ocp.model.cost_y_expr = vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    # set constraints
    ocp.constraints.lbu = np.array([-constraint_par['F_max']])
    ocp.constraints.ubu = np.array([+constraint_par['F_max']])

    ocp.constraints.idxbx = np.array([2])
    ocp.constraints.lbx = np.array([-constraint_par['v_max']])
    ocp.constraints.ubx = np.array([constraint_par['v_max']])

    ocp.constraints.idxbx_e = np.array([2])
    ocp.constraints.lbx_e = np.array([-constraint_par['v_max']])
    ocp.constraints.ubx_e = np.array([constraint_par['v_max']])

    ocp.constraints.idxsbx = np.array([0])
    ocp.constraints.lsbx = np.zeros((1,))
    ocp.constraints.usbx = np.zeros((1,))

    if parametric_constraints:
        ocp.constraints.uh = np.array([0, 0])
        ocp.constraints.lh = np.array([-ACADOS_INFTY, -ACADOS_INFTY])
        ocp.constraints.idxsh = np.array([0, 1])
        ocp.cost.zu = ocp.cost.zl = 2e3 * np.ones((3,))
        ocp.cost.Zu = ocp.cost.Zl = 5e3 * np.ones((3,))
    else:
        ocp.cost.zu = ocp.cost.zl = 2e3 * np.ones((1,))
        ocp.cost.Zu = ocp.cost.Zl = 5e3 * np.ones((1,))

    ocp.constraints.x0 = x0
    ocp.constraints.idxbu = np.array([0])

    # set prediction horizon
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = Tf

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'  # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.sim_method_newton_iter = 10

    if RTI:
        ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    else:
        ocp.solver_options.nlp_solver_type = 'SQP'
        ocp.solver_options.globalization = 'MERIT_BACKTRACKING'  # turns on globalization
        ocp.solver_options.nlp_solver_max_iter = 150

    ocp.solver_options.qp_solver_cond_N = N_horizon

    solver_json = 'acados_ocp_' + model.name + '.json'
    acados_ocp_solver = AcadosOcpSolver(ocp, json_file=solver_json, save_p_global=True)

    # create an integrator with the same settings as used in the OCP solver.
    # attention, with p_global, integrator cannot be built from OCP
    sim = AcadosSim()
    sim.model = export_pendulum_ode_model()

    # set simulation time
    sim.solver_options.T = td
    # set options
    sim.solver_options.integrator_type = 'IRK'
    sim.solver_options.num_stages = 3
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 3  # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    acados_integrator = AcadosSimSolver(sim)

    acados_evaluator = AcadosCostConstraintEvaluator(ocp, with_parametric_bounds=False)

    return acados_ocp_solver, acados_integrator, acados_evaluator


def update_comprehensive_dict(comprehensive_dict, input_dict):
    """
    Updates the global comprehensive_dict with the values from the input dictionary.
    Values for each key are stored in lists, which are appended to with each call.

    :param input_dict: A dictionary with keys and values to add to the comprehensive_dict.
    """
    for key, value in input_dict.items():
        if key not in comprehensive_dict:
            comprehensive_dict[key] = []  # Initialize the list if the key doesn't exist
        comprehensive_dict[key].append(value)  # Append the value to the list

    return comprehensive_dict


def main(use_RTI: bool = False, parametric_constraints: bool = True, plot_results: bool = False):
    global constraint_par
    x0 = np.array([0.0, np.pi, 0.0, 0.0])

    Tf = .8
    N_horizon = 40
    td = Tf / N_horizon

    ocp_solver, integrator, evaluator = setup(
        x0, N_horizon, Tf,
        td=Tf / N_horizon,
        RTI=use_RTI,
        parametric_constraints=parametric_constraints)

    nx = ocp_solver.acados_ocp.dims.nx
    nu = ocp_solver.acados_ocp.dims.nu

    Nsim = 100
    simX = np.zeros((Nsim + 1, nx))
    simU = np.zeros((Nsim, nu))
    eval_dict = {}

    simX[0, :] = x0

    # do some initial iterations to start with a good initial guess
    num_iter_initial = 5
    for _ in range(num_iter_initial):
        ocp_solver.solve_for_x0(x0_bar=x0)

    evaluator.update_all(ocp_solver)

    # closed loop
    for i in range(Nsim):

        # change constraint parameter in the middle of the simulation
        if i == constraint_par['iter_omega_change'] and parametric_constraints:
            new_omega_dot_min = np.array([constraint_par['omega_dot_min_2']])
            for j in range(ocp_solver.acados_ocp.dims.N):
                ocp_solver.set(j, "p", new_omega_dot_min)
            ocp_solver.set_p_global_and_precompute_dependencies(np.array([constraint_par['omega_dot_max_2']]))
            evaluator.update_all(ocp_solver)

        if use_RTI:
            # preparation phase
            ocp_solver.options_set('rti_phase', 1)
            ocp_solver.solve()

            # set initial state
            ocp_solver.set(0, "lbx", simX[i, :])
            ocp_solver.set(0, "ubx", simX[i, :])

            # feedback phase
            ocp_solver.options_set('rti_phase', 2)
            ocp_solver.solve()

            simU[i, :] = ocp_solver.get(0, "u")

        else:
            # solve ocp and get next control input
            simU[i, :] = ocp_solver.solve_for_x0(x0_bar=simX[i, :])

        # evaluate the cost of the full trajectory
        solution_obj = ocp_solver.store_iterate_to_obj()
        cost_ext_eval = evaluator.evaluate_ocp_cost(solution_obj)
        cost_int_eval = ocp_solver.get_cost()
        abs_error = np.abs(cost_ext_eval - cost_int_eval)

        # formatted print relative error up to 3 decimal places
        print(f'cost_err_abs: {abs_error:.9f}')
        assert math.isclose(cost_ext_eval, cost_int_eval, abs_tol=1e-8, rel_tol=1e-8)

        # simulate system
        simX[i + 1, :] = integrator.simulate(x=simX[i, :], u=simU[i, :])
        eval_iter = evaluator.evaluate(x=simX[i, :], u=simU[i, :])
        eval_dict = update_comprehensive_dict(eval_dict, eval_iter)

    # plot results
    if not plot_results:
        return
    model = ocp_solver.acados_ocp.model

    fix, axes = plot_pendulum_eval(
        np.linspace(0, td * Nsim, Nsim + 1),
        simU,
        simX,
        eval_dict,
        latexify=False,
        time_label=model.t_label,
        x_labels=model.x_labels,
        u_labels=model.u_labels)

    axes[2].axhline(constraint_par['v_max'], alpha=0.7, color='tab:red')
    axes[2].axhline(-constraint_par['v_max'], alpha=0.7, color='tab:red')

    if parametric_constraints:
        constr_omega_dot_min = np.empty(Nsim)
        constr_omega_dot_min[:constraint_par['iter_omega_change'] + 1] = constraint_par['omega_dot_min_1']
        constr_omega_dot_min[constraint_par['iter_omega_change'] + 1:] = constraint_par['omega_dot_min_2']
        constr_omega_dot_max = np.empty(Nsim)
        constr_omega_dot_max[:constraint_par['iter_omega_change'] + 1] = constraint_par['omega_dot_max_1']
        constr_omega_dot_max[constraint_par['iter_omega_change'] + 1:] = constraint_par['omega_dot_max_2']
        axes[3].plot(np.linspace(0, td * Nsim, Nsim), constr_omega_dot_min, alpha=0.7, color='tab:red')
        axes[3].plot(np.linspace(0, td * Nsim, Nsim), constr_omega_dot_max, alpha=0.7, color='tab:red')


    axes[-1].set_ylim([-1.2 * constraint_par['F_max'], 1.2 * constraint_par['F_max']])
    axes[-1].axhline(constraint_par['F_max'], alpha=0.7, color='tab:red')
    axes[-1].axhline(-constraint_par['F_max'], alpha=0.7, color='tab:red')
    plt.show()


if __name__ == '__main__':

    for parametric_constraints in [True, False]:
        print(f'Parametric_constraints: {parametric_constraints}\n')
        main(use_RTI=True, parametric_constraints=parametric_constraints, plot_results=False)
