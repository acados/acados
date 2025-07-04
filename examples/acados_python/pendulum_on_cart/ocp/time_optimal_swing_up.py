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

from acados_template import AcadosOcp, AcadosOcpSolver, ACADOS_INFTY, AcadosOcpOptions
from pendulum_model import export_free_time_pendulum_ode_model
import numpy as np
from utils import plot_pendulum
import casadi as ca
from casadi.tools import entry, struct_symSX

def formulate_ocp(opts: AcadosOcpOptions) -> AcadosOcp:
    # create ocp object to formulate the OCP
    N = 100
    nx = 5
    nu = 1
    Tf = 1.0

    # Parameters
    max_f = 5.
    max_x1 = 1.0
    max_v = 2.0

    x1_0 = 0.0
    theta_0 = np.pi
    dx1_0 = 0.0
    dtheta_0 = 0.0

    theta_f = 0.0
    dx1_f = 0.0
    dtheta_f = 0.0

    ocp = AcadosOcp()

    model = export_free_time_pendulum_ode_model()

    if opts.qpscaling_scale_objective:
        model.name += "scaled_objective"
    else:
        model.name += "no_scaling"

    # set model
    ocp.model = model

    # Initial conditions
    ocp.constraints.lbx_0 = np.array([0.0, x1_0, theta_0, dx1_0, dtheta_0])
    ocp.constraints.ubx_0 = np.array([ACADOS_INFTY, x1_0, theta_0, dx1_0, dtheta_0])
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3, 4])

    # Actuator constraints
    ocp.constraints.lbu = np.array([-max_f])
    ocp.constraints.ubu = np.array([+max_f])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.lbx = np.array([0.0, -max_x1, -max_v])
    ocp.constraints.ubx = np.array([ACADOS_INFTY, max_x1, max_v])
    ocp.constraints.idxbx = np.array([0, 1, 3])

    # Terminal constraints
    ocp.constraints.lbx_e = np.array([0.0, theta_f, dx1_f, dtheta_f])
    ocp.constraints.ubx_e = np.array([ACADOS_INFTY, theta_f, dx1_f, dtheta_f])
    ocp.constraints.idxbx_e = np.array([0, 2, 3, 4])

    # Define objective function
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = model.x[0]

    # set solver options
    ocp.solver_options = opts
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    return ocp


def main(scale_objective: bool):

    opts = AcadosOcpOptions()
    opts.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
    opts.qp_solver_mu0 = 1e3
    opts.hessian_approx = 'EXACT'
    opts.regularize_method = 'MIRROR'
    opts.integrator_type = 'ERK'
    opts.print_level = 1
    opts.nlp_solver_max_iter = 1000
    opts.qp_solver_iter_max = 1000

    opts.nlp_solver_type = 'SQP_WITH_FEASIBLE_QP'#'SQP'#

    # Globalization
    opts.globalization = 'FUNNEL_L1PEN_LINESEARCH'#'MERIT_BACKTRACKING'#
    opts.globalization_full_step_dual = True
    opts.globalization_funnel_use_merit_fun_only = False

    # Scaling
    if scale_objective:
        opts.qpscaling_scale_objective = "OBJECTIVE_GERSHGORIN"
    else:
        opts.qpscaling_scale_objective = "NO_OBJECTIVE_SCALING"

    ocp = formulate_ocp(opts)

    N = ocp.solver_options.N_horizon

    # create solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)

    # intialize
    T0 = 1.0
    for i in range(N):
        ocp_solver.set(i, "x", np.array([T0, 0.0, np.pi, 0.0, 0.0]))
        ocp_solver.set(i, "u", np.array([0.0]))
    ocp_solver.set(N, "x", np.array([T0, 0.0, np.pi, 0.0, 0.0]))

    # solve
    status = ocp_solver.solve()

    if scale_objective:
        if status != 0:
            raise ValueError("Objective scaling should make the problem solvable!")
        print(f"Objective scaling worked, status: {status}.")
    else:
        if status == 0:
            raise ValueError("In Byrd-Omojokun equations, the multipliers should get too large which makes HPIPM fail!")
        print("Problem not solvable without objective scaling, as expected.")

if __name__ == "__main__":
    main(scale_objective=False)
    main(scale_objective=True)
