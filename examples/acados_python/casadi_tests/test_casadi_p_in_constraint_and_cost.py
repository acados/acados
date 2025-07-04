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
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosCasadiOcpSolver
from pendulum_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from utils import plot_pendulum
from casadi import SX, vertcat

PLOT = False

def formulate_ocp(p_in_constraint=True, p_in_cost=True):
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model
    x = model.x
    u = model.u

    Tf = 1.0
    nx = x.rows()
    nu = u.rows()
    N = 20

    # set dimensions
    ocp.solver_options.N_horizon = N

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])
    cost_W = scipy.linalg.block_diag(Q, R)

    n_param = 6
    p = SX.sym('p', n_param)
    constraint_quotient = p[0]
    y_param = p[1:nx+nu+1]
    ocp.model.p = p

    # define cost with parametric reference
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.cost.cost_type_e = 'EXTERNAL'

    residual = y_param - vertcat(x, u)
    ocp.model.cost_expr_ext_cost = residual.T @ cost_W @ residual
    res_e = y_param[0:nx] - x
    ocp.model.cost_expr_ext_cost_e = res_e.T @ Q @ res_e

    # set constraints
    Fmax = 80
    ocp.constraints.lh = np.array([-Fmax])
    ocp.constraints.uh = np.array([+Fmax])
    ocp.model.con_h_expr = model.u / constraint_quotient

    ocp.constraints.lh_0 = np.array([-Fmax])
    ocp.constraints.uh_0 = np.array([+Fmax])
    ocp.model.con_h_expr_0 = model.u / constraint_quotient

    p_0 = np.zeros(n_param)
    p_0[0] = 0.5 if p_in_constraint else 1
    p_0[1] = 1 if p_in_cost else 0
    ocp.parameter_values = p_0

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'EXACT' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.regularize_method = 'CONVEXIFY' # GAUSS_NEWTON, EXACT
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization

    # set prediction horizon
    ocp.solver_options.tf = Tf
    return ocp


def main(p_in_constraint=True, p_in_cost=True):
    
    print(f"\n\nRunning example with p_in_constraint={p_in_constraint}, p_in_cost={p_in_cost}")
    ocp = formulate_ocp(p_in_constraint, p_in_cost)

    N = ocp.solver_options.N_horizon
    Tf = ocp.solver_options.tf

    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    status = ocp_solver.solve()

    if status != 0:
        raise Exception(f'acados returned status {status}.')
    result = ocp_solver.store_iterate_to_obj()

    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False)
    casadi_ocp_solver.load_iterate_from_obj(result)
    casadi_ocp_solver.solve()
    result_casadi = casadi_ocp_solver.store_iterate_to_obj()

    result.flatten().allclose(other=result_casadi.flatten())
    
    Fmax = 80
    if PLOT: 
        acados_u = np.array([ocp_solver.get(i, "u") for i in range(N)])
        acados_x = np.array([ocp_solver.get(i, "x") for i in range(N+1)])
        casadi_u = np.array([casadi_ocp_solver.get(i, "u") for i in range(N)])
        casadi_x = np.array([casadi_ocp_solver.get(i, "x") for i in range(N+1)])
        plot_pendulum(np.linspace(0, Tf, N+1), Fmax, acados_u, acados_x, latexify=False)
        plot_pendulum(np.linspace(0, Tf, N+1), Fmax, casadi_u, casadi_x, latexify=False)

if __name__ == "__main__":
    main(p_in_constraint=False, p_in_cost=False)
    main(p_in_constraint=True, p_in_cost=False)
    main(p_in_constraint=False, p_in_cost=True)
    main(p_in_constraint=True, p_in_cost=True)