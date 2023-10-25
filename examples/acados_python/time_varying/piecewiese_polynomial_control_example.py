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


from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, casadi_length
import numpy as np
import casadi as ca
from scipy.linalg import block_diag
from utils import plot_open_loop_trajectory_pwpol_u, export_pendulum_ode_model

def augment_model_with_polynomial_control(model: AcadosModel, d: int = 1, delta_T=None):
    # add time to model
    if model.t == []:
        model.t = ca.SX.sym('t')
    t = model.t

    if delta_T is None:
        delta_T = ca.SX.sym('delta_T')

    u_old = model.u
    nu_original = casadi_length(model.u)

    u_coeff = ca.SX.sym('u_coeff', (d+1) * nu_original)
    u_new = ca.SX.zeros(nu_original, 1)
    for i in range(d+1):
        u_new += (t / delta_T) ** i * u_coeff[i*nu_original:(i+1)*nu_original]

    evaluate_polynomial_u_fun = ca.Function("evaluate_polynomial_u", [u_coeff, t, delta_T], [u_new])

    model.f_impl_expr = ca.substitute(model.f_impl_expr, u_old, u_new)
    model.cost_y_expr = ca.substitute(model.cost_y_expr, u_old, u_new)

    model.u = u_coeff

    model.p = ca.vertcat(model.p, delta_T)

    return model, evaluate_polynomial_u_fun


def main():
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    nu_original = model.u.size()[0]

    ocp.model = model

    T_horizon = 1.0
    nx = model.x.size()[0]
    N_horizon = 20
    d = 3

    # set dimensions
    ocp.dims.N = N_horizon

    # set cost
    Q_mat = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R_mat = 2*np.diag([1e-1])

    # the 'EXTERNAL' cost type can be used to define general cost terms
    # NOTE: This leads to additional (exact) hessian contributions when using GAUSS_NEWTON hessian.
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x

    ocp.cost.W = block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat
    ocp.cost.yref = np.zeros(nx+nu_original)
    ocp.cost.yref_e = np.zeros(nx)

    # set constraints
    Fmax = 80
    # ocp.constraints.lbu = np.array([-Fmax])
    # ocp.constraints.ubu = np.array([+Fmax])
    # ocp.constraints.idxbu = np.array([0])
    # formulate u bounds as barrier
    ocp.cost.W = block_diag(ocp.cost.W, 1e5)
    ocp.cost.yref = np.concatenate((ocp.cost.yref, np.zeros(nu_original)))
    ocp.model.cost_y_expr = ca.vertcat(ocp.model.cost_y_expr,
                                       ca.fmax(ca.fmax(0.0, (ocp.model.u - Fmax)), -Fmax - ocp.model.u))

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    # PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_QPOASES, FULL_CONDENSING_HPIPM,
    # PARTIAL_CONDENSING_QPDUNES, PARTIAL_CONDENSING_OSQP, FULL_CONDENSING_DAQP
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'IRK'
    # ocp.solver_options.print_level = 1
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.cost_discretization = 'INTEGRATOR'

    # set prediction horizon
    n_short = 3
    dt_short = 0.02
    n_long = N_horizon-n_short
    ocp.solver_options.time_steps = np.array(n_short * [dt_short] + n_long * [(T_horizon - dt_short * n_short) / n_long])
    ocp.solver_options.tf = T_horizon

    ocp.parameter_values = np.array([T_horizon / N_horizon])
    model, evaluate_polynomial_u_fun = augment_model_with_polynomial_control(model, d=d)
    ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')
    for i in range(N_horizon):
        ocp_solver.set(i, 'p', ocp.solver_options.time_steps[i])

    nu = model.u.size()[0]
    simX = np.ndarray((N_horizon+1, nx))
    simU = np.ndarray((N_horizon, nu))


    status = ocp_solver.solve()
    ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    cost_val = ocp_solver.get_cost()
    print(f"found optimal cost: {cost_val:.4e}")

    # get solution
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
            U_interval[j, :] = evaluate_polynomial_u_fun(u_coeff, t, dt)
        U_fine.append(U_interval)


    shooting_nodes = np.concatenate((np.array([0.]), np.cumsum(ocp.solver_options.time_steps)))
    plot_open_loop_trajectory_pwpol_u(shooting_nodes, simX, U_fine, title=f'controls: piecewise polynomials of degree {d}')


if __name__ == '__main__':
    main()