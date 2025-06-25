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
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel, AcadosCasadiOcpSolver
import numpy as np
import sys
sys.path.insert(0, '../getting_started')
from utils import plot_pendulum

from casadi import MX, vertcat, sin, cos
import casadi as ca
import time

PLOT = False

def export_pendulum_ode_model() -> AcadosModel:
    # constants
    m_cart = 1. # mass of the cart [kg]

    # parameters
    g = MX.sym("g")
    p = g
    m = MX.sym("m")
    l = MX.sym("l")
    p_global = ca.vertcat(m,l)

    # set up states & controls
    x1      = MX.sym('x1')
    theta   = MX.sym('theta')
    v1      = MX.sym('v1')
    dtheta  = MX.sym('dtheta')

    x = vertcat(x1, theta, v1, dtheta)

    F = MX.sym('F')
    u = vertcat(F)

    # xdot
    nx = x.shape[0]
    xdot = MX.sym('xdot', nx)

    # dynamics
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    denominator = m_cart + m - m*cos_theta**2
    f_expl = vertcat(v1,
                     dtheta,
                     (-m*l*sin_theta*dtheta**2 + m*g*cos_theta*sin_theta+F)/denominator,
                     (-m*l*cos_theta*sin_theta*dtheta**2 + F*cos_theta+(m_cart+m)*g*sin_theta)/(l*denominator)
                     )

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.p = p
    model.p_global = p_global
    model.name = 'pendulum_ode'

    # store meta information
    model.x_labels = ['$x$ [m]', r'$\theta$ [rad]', '$v$ [m]', r'$\dot{\theta}$ [rad/s]']
    model.u_labels = ['$F$']
    model.t_label = '$t$ [s]'

    return model


def ocp_formulation() -> AcadosOcp:

    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model()
    ocp.model = model

    # dimensions
    nx = model.x.rows()
    nu = model.u.rows()

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    # path cost
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q, R).full()

    # terminal cost
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    ocp.parameter_values = np.array([9.81])
    ocp.p_global_values = np.array([0.1, 0.8])
    
    return ocp


def main(stage_varying=True):

    print(f"\n\nRunning example with stage_varying_p={stage_varying}")
    # create ocp
    ocp = ocp_formulation()

    Tf = 1.0
    N_horizon = 20

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON' # 'GAUSS_NEWTON', 'EXACT'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP' # SQP_RTI, SQP
    ocp.solver_options.globalization = 'MERIT_BACKTRACKING' # turns on globalization

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N_horizon

    # create acados solver
    print(f"Creating ocp solver with p_global = {ocp.model.p_global}, p = {ocp.model.p}")
    ocp_solver = AcadosOcpSolver(ocp, generate=True, build=True, verbose=False, save_p_global=True)
    if stage_varying:
        for i in range(0, N_horizon+1):
            ocp_solver.set(stage_= i, field_= 'p', value_=np.array([9.81+i*0.3]))
    status = ocp_solver.solve()

    if status != 0:
        raise Exception(f'acados returned status {status}.')
    result = ocp_solver.store_iterate_to_obj()

    casadi_ocp_solver = AcadosCasadiOcpSolver(ocp, verbose=False)
    if stage_varying:
        for i in range(0, N_horizon+1):
            casadi_ocp_solver.set(stage= i, field= 'p', value_=np.array([9.81+i*0.3]))
    casadi_ocp_solver.load_iterate_from_obj(result)
    casadi_ocp_solver.solve()
    result_casadi = casadi_ocp_solver.store_iterate_to_obj()

    result.flatten().allclose(other=result_casadi.flatten())

    if PLOT:
        acados_u = np.array([ocp_solver.get(i, "u") for i in range(N_horizon)])
        acados_x = np.array([ocp_solver.get(i, "x") for i in range(N_horizon+1)])
        casadi_u = np.array([casadi_ocp_solver.get(i, "u") for i in range(N_horizon)])
        casadi_x = np.array([casadi_ocp_solver.get(i, "x") for i in range(N_horizon+1)])
        plot_pendulum(ocp.solver_options.shooting_nodes, ocp.constraints.ubu[0], acados_u, acados_x, x_labels=ocp.model.x_labels, u_labels=ocp.model.u_labels)
        plot_pendulum(ocp.solver_options.shooting_nodes, ocp.constraints.ubu[0], casadi_u, casadi_x, x_labels=ocp.model.x_labels, u_labels=ocp.model.u_labels)

if __name__ == "__main__":
    main(stage_varying=False)
    main(stage_varying=True)
