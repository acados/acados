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

import sys
sys.path.insert(0, '../getting_started')

import numpy as np
import casadi as ca

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSim, AcadosSimSolver, AcadosCasadiOcpSolver
from pendulum_model import export_pendulum_ode_model

def formulate_ocp(Tf: float = 1.0, N: int = 20)-> AcadosOcp:
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    # add dummy control
    model.u = ca.vertcat(model.u, ca.SX.sym('dummy_u'))
    ocp.model = model

    # set h constraints
    ocp.model.con_h_expr_0 = ca.norm_2(model.x)
    ocp.constraints.lh_0 = np.array([0])
    ocp.constraints.uh_0 = np.array([3.16])

    nx = model.x.rows()
    nu = model.u.rows()

    # set prediction horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost matrices
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-2, 1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0, np.pi, 0, 0])  # initial state
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3])

    # set partial bounds for state
    ocp.constraints.lbx = np.array([-10,-10])
    ocp.constraints.ubx = np.array([10,10])
    ocp.constraints.idxbx = np.array([0,3])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization

    return ocp

def compare_integrators(ocp, ocp_solver, acasadi_solver):
    sim = AcadosSim.from_ocp(ocp)
    sim_solver = AcadosSimSolver(sim, verbose=False)

    for i in range(ocp.solver_options.N_horizon):
        x_i = ocp_solver.get(i, 'x')
        u_i = ocp_solver.get(i, 'u')
        dt_i = ocp.solver_options.time_steps[i]

        sim_solver.set("T", dt_i)
        sim_solver.set("x", x_i)
        sim_solver.set("u", u_i)
        sim_solver.solve()
        x_sim = sim_solver.get('x')

        x_casados = acasadi_solver.f_discr_fun(x_i, u_i, dt_i)

        diff = np.linalg.norm(x_sim - x_casados)
        if diff < 1e-6:
            print(f"Step {i}: ||x_sim - x_casados|| = {diff:.6e} < 1e-6, PASSED")
        else:
            raise Exception(f"Step {i}: ||x_sim - x_casados|| = {diff:.6e} >= 1e-6, FAILED")

def main(casadi_solver_name="ipopt", itype="ERK", casados=True):
    ocp = formulate_ocp()
    ocp.model.name = f"model_{itype}"
    ocp.code_gen_options.code_export_directory = f"c_generated_code_{itype}"
    ocp.solver_options.integrator_type = itype
    if itype != "GNSF":
        tau = np.linspace(0, 1, ocp.solver_options.N_horizon+1)
        nodes = ocp.solver_options.tf * tau**3
        ocp.solver_options.time_steps = np.diff(nodes)

    initial_iterate = ocp.create_default_initial_iterate()

    ## solve using acados
    # create acados solver
    ocp_solver = AcadosOcpSolver(ocp, verbose=False)
    # initialize solver
    ocp_solver.set_iterate(initial_iterate)
    # solve with acados
    status = ocp_solver.solve()
    ocp_solver.print_statistics()
    # get solution
    result_acados = ocp_solver.get_iterate()

    print(f"Using Casadi solver {casadi_solver_name}.")
    acasadi_solver = AcadosCasadiOcpSolver(ocp,
                                           verbose=False,
                                           solver=casadi_solver_name,
                                           with_casados = casados)

    compare_integrators(ocp, ocp_solver, acasadi_solver)

    acasadi_solver.set_iterate(result_acados)
    status = acasadi_solver.solve()
    print(f"Casadi solver returned status {status}.")
    nlp_iter_ca = acasadi_solver.get_stats("nlp_iter")
    print(f"Casadi solver finished after {nlp_iter_ca} iterations.")

    x_acados = ocp_solver.get_flat('x')
    u_acados = ocp_solver.get_flat('u')
    lambda_acados = ocp_solver.get_flat('lam')
    pi_acados = ocp_solver.get_flat('pi')

    x_casadi = acasadi_solver.get_flat('x')
    u_casadi = acasadi_solver.get_flat('u')
    lambda_casadi = acasadi_solver.get_flat('lam')
    pi_casadi = acasadi_solver.get_flat('pi')

    diff_x = np.linalg.norm(x_acados - x_casadi)
    print(f"||x_acados - x_casadi|| = {diff_x:.6e}")
    diff_u = np.linalg.norm(u_acados - u_casadi)
    print(f"||u_acados - u_casadi|| = {diff_u:.6e}")
    diff_lambda = np.linalg.norm(lambda_acados - lambda_casadi)
    print(f"||lambda_acados - lambda_casadi|| = {diff_lambda:.6e}")
    diff_pi = np.linalg.norm(pi_acados - pi_casadi)
    print(f"||pi_acados - pi_casadi|| = {diff_pi:.6e}")

    max_diff = max(diff_x, diff_u)
    if max_diff > 5e-4:
        raise Exception(f"Max difference between acados and casadi solver is {max_diff:.6e} >= 5e-4, FAILED")

if __name__ == "__main__":
    intergrator_types = ["ERK", "IRK", "GNSF"]
    for itype in intergrator_types:
        print(f"Testing integrator type {itype}.")
        main(itype=itype, casados=True)
        if itype == "ERK":
            print(f"Testing integrator type {itype} without casados.")
            main(itype=itype, casados=False)