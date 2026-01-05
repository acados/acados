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

from acados_template import AcadosParamManager, AcadosParam, AcadosOcp, AcadosModel, AcadosOcpSolver
import numpy as np

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

    param_manager = AcadosParamManager()

    Q_expr = param_manager.get_parameter_expression("Q")
    R_expr = param_manager.get_parameter_expression("R")
    A_expr = param_manager.get_parameter_expression("A")
    B_expr = param_manager.get_parameter_expression("B")

    ocp = AcadosOcp([Q_param, R_param, A_param, B_param], N_horizon)


    return ocp, param_manager

def main():

    ocp, param_manager = setup_qp()

    # Create solver
    ocp_solver = AcadosOcpSolver(ocp, json_file="acados_ocp_qp.json")

    # Set parameter values
    param_manager.set_param_value("Q", np.array([[1.0, 0.0], [0.0, 1.0]]), stage=0)
    param_manager.set_param_value("R", np.array([[0.1]]), stage=0)

    # Update solver with new parameter values
    ocp_solver.set(0, "p", param_manager.get_p_stagewise_values(stage=0))

    # Solve the QP
    status = ocp_solver.solve()

    if status != 0:
        print(f"Solver failed with status {status}")

