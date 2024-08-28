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

print("use CasADi branch se")

import sys
sys.path.insert(0, '../common')

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
import numpy as np
import scipy.linalg
from utils import plot_pendulum

from casadi import MX, vertcat, sin, cos, Function
import casadi as ca

knots = [[0,0,0,0,0.2,0.5,0.8,1,1,1,1],[0,0,0,0.1,0.5,0.9,1,1,1]]
np.random.seed(1)
data = np.random.random((7,6,2)).ravel(order='F')

def export_pendulum_ode_model(lut=True) -> AcadosModel:

    model_name = 'pendulum'

    # constants
    m_cart = 1. # mass of the cart [kg]

    # parameters
    g = MX.sym("g")
    p = g

    m = MX.sym("m")
    l = MX.sym("l")
    p_slow = [m, l]

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

    # parameters
    if lut:
        # Coefficient of B-spline
        C = MX.sym("C",data.shape[0],1)
        p_slow += [C]

    # dynamics
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    denominator = m_cart + m - m*cos_theta*cos_theta
    f_expl = vertcat(v1,
                     dtheta,
                     (-m*l*sin_theta*dtheta*dtheta + m*g*cos_theta*sin_theta+F)/denominator,
                     (-m*l*cos_theta*sin_theta*dtheta*dtheta + F*cos_theta+(m_cart+m)*g*sin_theta)/(l*denominator)
                     )

    if lut:
        x_in = ca.vertcat(u/100+0.5,theta/np.pi+0.5)

        # Disturb the dynamics by a sprinkle of bspline
        f_expl[2:4] += 0.01*ca.bspline(x_in,C,knots,[3,2],2)

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    # model.z = z
    model.p = p
    model.p_slow = ca.vcat(p_slow)
    model.name = model_name

    return model


def main(use_cython=False, lut=True, use_p_slow=True):

    print(f"\n\nRunning example with lut={lut}, use_p_slow={use_p_slow}")
    # create ocp object to formulate the OCP
    ocp = AcadosOcp()

    # set model
    model = export_pendulum_ode_model(lut=lut)
    ocp.model = model

    # parameter values
    p_values = np.array([9.81])
    p_slow_values = np.array([0.1, 0.8])
    if lut:
        p_slow_values = np.concatenate([p_slow_values, data])
    if not use_p_slow:
        model.p = vertcat(model.p, model.p_slow)
        model.p_slow = None
        p_values = np.concatenate([p_values, p_slow_values])
        p_slow_values = np.array([])

    # tmp hack
    # model.p = vertcat(model.p, model.p_slow)
    # p_values = np.concatenate([p_values, p_slow_values])

    ocp.parameter_values = p_values

    Tf = 1.0
    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx
    N_horizon = 20

    # set cost
    Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
    R = 2*np.diag([1e-2])

    ocp.cost.W_e = Q
    ocp.cost.W = scipy.linalg.block_diag(Q, R)

    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    # set constraints
    Fmax = 80
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([+Fmax])
    ocp.constraints.idxbu = np.array([0])

    ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])

    # set options
    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.print_level = 0
    ocp.solver_options.nlp_solver_type = 'SQP_RTI' # SQP_RTI, SQP

    ocp.solver_options.custom_update_filename = 'custom_update_function.c'
    ocp.solver_options.custom_update_header_filename = 'custom_update_function.h'

    # set prediction horizon
    ocp.solver_options.tf = Tf
    ocp.solver_options.N_horizon = N_horizon

    ocp.solver_options.custom_update_copy = False
    ocp.solver_options.custom_templates = [
        ('custom_update_p_slow_template.in.c', 'custom_update_function.c'),
        ('custom_update_function_zoro_template.in.h', 'custom_update_function.h'),
    ]

    # create ocp solver
    print(f"Creating ocp solver with p_slow = {model.p_slow}, p = {model.p}")

    solver_json = 'acados_ocp_' + model.name + '.json'
    if use_cython:
        AcadosOcpSolver.generate(ocp, json_file=solver_json)
        AcadosOcpSolver.build(ocp.code_export_directory, with_cython=True)
        ocp_solver = AcadosOcpSolver.create_cython_solver(solver_json)
    else:
        ocp_solver = AcadosOcpSolver(ocp, json_file = solver_json)

    # call SQP_RTI solver in the loop:
    residuals = []

    if use_p_slow:
        ocp_solver.custom_update(p_slow_values)
    for i in range(20):
        status = ocp_solver.solve()
        # ocp_solver.print_statistics() # encapsulates: stat = ocp_solver.get_stats("statistics")
        residuals+= list(ocp_solver.get_residuals())

    print(residuals)

    return residuals


if __name__ == "__main__":
    ref_nolut = main(use_cython=False, use_p_slow=False, lut=False)
    res_nolut = main(use_cython=False, use_p_slow=True, lut=False)
    np.testing.assert_almost_equal(ref_nolut,res_nolut)
    ref_lut = main(use_cython=False, use_p_slow=False, lut=True)
    res_lut = main(use_cython=False, use_p_slow=True, lut=True)
    np.testing.assert_almost_equal(ref_lut, res_lut)

    with np.testing.assert_raises(Exception):
        np.testing.assert_almost_equal(ref_lut, ref_nolut)

    # main(use_cython=True)
