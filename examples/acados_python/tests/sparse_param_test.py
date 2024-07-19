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
sys.path.insert(0, '../pendulum_on_cart/common')

from acados_template import AcadosSim, AcadosSimSolver, AcadosModel, sim_get_default_cmake_builder
from utils import plot_pendulum

import casadi as ca
import numpy as np


def export_pendulum_ode_model() -> AcadosModel:
    model_name = 'pendulum'

    # constants
    m_cart = 1. # mass of the cart [kg]
    m = 0.1 # mass of the ball [kg]
    g = 9.81 # gravity constant [m/s^2]
    l = 0.8 # length of the rod [m]

    # set up states & controls
    x1      = ca.SX.sym('x1')
    theta   = ca.SX.sym('theta')
    v1      = ca.SX.sym('v1')
    dtheta  = ca.SX.sym('dtheta')

    x = ca.vertcat(x1, theta, v1, dtheta)

    F = ca.SX.sym('F')
    u = ca.vertcat(F)

    # xdot
    x1_dot      = ca.SX.sym('x1_dot')
    theta_dot   = ca.SX.sym('theta_dot')
    v1_dot      = ca.SX.sym('v1_dot')
    dtheta_dot  = ca.SX.sym('dtheta_dot')

    xdot = ca.vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot)

    # parameters
    # NOTE: such sparse parameter formulations are not recommended to use in acados!
    p = ca.SX.zeros(1000, 1)
    p[-1] = ca.SX.sym('p')
    m_cart = p[-1]
    p = ca.sparsify(p)
    # dynamics
    cos_theta = ca.cos(theta)
    sin_theta = ca.sin(theta)
    denominator = m_cart + m - m*cos_theta*cos_theta
    f_expl = ca.vertcat(v1,
                     dtheta,
                     (-m*l*sin_theta*dtheta*dtheta + m*g*cos_theta*sin_theta+F)/denominator,
                     (-m*l*cos_theta*sin_theta*dtheta*dtheta + F*cos_theta+(m_cart+m)*g*sin_theta)/(l*denominator)
                     )

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    # model.z = z
    model.p = p
    model.name = model_name

    return model


def sparse_param_test():
    sim = AcadosSim()

    # model
    model = export_pendulum_ode_model()
    sim.model = model

    Tf = 0.1
    nx = model.x.rows()
    nu = model.u.rows()
    N = 200

    # set simulation time
    sim.solver_options.T = Tf
    # set options
    sim.solver_options.num_stages = 7
    sim.solver_options.num_steps = 3
    sim.solver_options.newton_iter = 10 # for implicit integrator
    sim.solver_options.collocation_type = "GAUSS_RADAU_IIA"
    sim.solver_options.integrator_type = "IRK" # ERK, IRK, GNSF
    sim.solver_options.sens_forw = True
    sim.solver_options.sens_adj = True
    sim.solver_options.sens_hess = False
    sim.solver_options.sens_algebraic = False
    sim.solver_options.output_z = False
    sim.solver_options.sim_method_jac_reuse = False
    sim.parameter_values = np.zeros(model.p.shape).flatten()
    sim.parameter_values[-1] = 1.0


    # create
    cmake_builder = sim_get_default_cmake_builder()
    acados_integrator = AcadosSimSolver(sim)

    simX = np.zeros((N+1, nx))
    x0 = np.array([0.0, np.pi+1, 0.0, 0.0])
    u0 = np.array([0.0])
    acados_integrator.set("u", u0)

    simX[0,:] = x0

    ## Single test call
    import time

    t0 = time.time()
    acados_integrator.set("seed_adj", np.ones((nx, 1)))
    acados_integrator.set("x", x0)
    acados_integrator.set("u", u0)
    status = acados_integrator.solve()
    time_external = time.time() - t0

    S_forw = acados_integrator.get("S_forw")
    Sx = acados_integrator.get("Sx")
    Su = acados_integrator.get("Su")
    S_hess = acados_integrator.get("S_hess")
    S_adj = acados_integrator.get("S_adj")
    print(f"\ntimings of last call to acados_integrator: with Python interface, set and get {time_external*1e3:.4f}ms")


    # get timings (of last call)
    CPUtime = acados_integrator.get("CPUtime")
    LAtime = acados_integrator.get("LAtime")
    ADtime = acados_integrator.get("ADtime")
    print(f"\ntimings of last call to acados_integrator: overall CPU: {CPUtime*1e3:.4f} ms, linear algebra {LAtime*1e3:.4f} ms, external functions {ADtime*1e3:.4f} ms")

    print("S_forw, sensitivities of simulation result wrt x,u:\n", S_forw)
    print("Sx, sensitivities of simulation result wrt x:\n", Sx)
    print("Su, sensitivities of simulation result wrt u:\n", Su)
    print("S_adj, adjoint sensitivities:\n", S_adj)
    print("S_hess, second order sensitivities:\n", S_hess)

    # turn off sensitivity propagation when not needed
    acados_integrator.options_set('sens_forw', False)
    acados_integrator.options_set('sens_adj', False)
    acados_integrator.options_set('sens_hess', False)

    # call in loop:
    for i in range(N):
        # set initial state
        acados_integrator.set("x", simX[i,:])
        # solve
        status = acados_integrator.solve()
        # get solution
        simX[i+1,:] = acados_integrator.get("x")

    if status != 0:
        raise Exception(f'acados returned status {status}.')

    if not np.allclose(S_adj, np.array([1., 0.32978441, 1.1, 1.0644604, 0.03091267])):
        raise Exception("adjoint sensitivities should match reference solution.")

    # plot results
    plot_pendulum(np.linspace(0, N*Tf, N+1), 10, np.zeros((N, nu)), simX, latexify=False)


if __name__ == '__main__':
    sparse_param_test()