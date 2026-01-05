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

from acados_template import AcadosParamManager, AcadosParam, AcadosOcp, AcadosOcpSolver
import numpy as np
import casadi as ca

def setup_qp():

    # Define dimensions
    nx = 2  # number of states
    nu = 1  # number of inputs
    N_horizon = 10  # horizon length

    # Define QP matrices with parameters
    Q_param = AcadosParam("Q", np.ones((nx, nx)))
    R_param = AcadosParam("R", np.ones((nu, nu)))
    A_param = AcadosParam("A", np.ones((nx, nx)))
    B_param = AcadosParam("B", np.ones((nx, nu)))

    param_manager = AcadosParamManager([Q_param, R_param, A_param, B_param], N_horizon)

    Q_expr = param_manager.get_expression("Q")
    R_expr = param_manager.get_expression("R")
    A_expr = param_manager.get_expression("A")
    B_expr = param_manager.get_expression("B")

    x = ca.SX.sym('x', nx)
    u = ca.SX.sym('u', nu)

    # Discrete dynamics expression
    discrete_dyn = A_expr @ x + B_expr @ u

    # External cost expressions
    cost_expr = ca.vertcat(x, u).T @ ca.vertcat(
        ca.horzcat(Q_expr, ca.DM.zeros((nx, nu))),
        ca.horzcat(ca.DM.zeros((nu, nx)), R_expr)
    ) @ ca.vertcat(x, u)

    cost_expr_e = x.T @ Q_expr @ x

    # Create OCP model
    ocp = AcadosOcp()
    ocp.model.name = 'parametric_qp_model'

    # Set symbolic variables
    ocp.model.x = x
    ocp.model.u = u
    ocp.model.p = param_manager.get_p_stagewise_expression()

    # Set discrete dynamics
    ocp.model.disc_dyn_expr = discrete_dyn

    # Set external cost
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost = cost_expr
    ocp.model.cost_expr_ext_cost_e = cost_expr_e

    # Set initial condition
    x0 = np.array([1.0, 0.0])
    ocp.constraints.x0 = x0

    # Set control bounds
    ocp.constraints.lbu = np.array([-1.0])
    ocp.constraints.ubu = np.array([1.0])
    ocp.constraints.idxbu = np.array([0])

    # Horizon configuration
    ocp.solver_options.N_horizon = N_horizon
    ocp.solver_options.tf = N_horizon

    # Solver options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    ocp.solver_options.hessian_approx = 'EXACT'
    ocp.solver_options.integrator_type = 'DISCRETE'

    ocp.parameter_values = param_manager.get_p_stagewise_values(stage=0)

    return ocp, param_manager

def main():

    ocp, param_manager = setup_qp()

    # Create solver
    ocp_solver = AcadosOcpSolver(ocp, json_file="acados_ocp_qp.json")

    # Set parameter values
    param_manager.set_value("Q", np.array([[1.0, 0.0], [0.0, 1.0]]), stage=0)
    param_manager.set_value("R", np.array([[0.1]]), stage=0)

    # Update solver with new parameter values
    ocp_solver.set(0, "p", param_manager.get_p_stagewise_values(stage=0))

    # Solve the QP
    status = ocp_solver.solve()

    if status != 0:
        print(f"Solver failed with status {status}")
    else:
        print("Success")


if __name__ == "__main__":
    main()
