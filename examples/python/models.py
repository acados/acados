#
#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren, Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor, Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan, Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
from casadi import cos, Function, sin, SX, vertcat

def chen_model():
    """ The following ODE model comes from Chen1998. """
    nx, nu = (2, 1)
    x = SX.sym('x', nx)
    u = SX.sym('u', nu)
    mu = 0.5
    rhs = vertcat(x[1] + u*(mu + (1.-mu)*x[0]), x[0] + u*(mu - 4.*(1.-mu)*x[1]))
    return Function('chen', [x, u], [rhs], ['x', 'u'], ['xdot']), nx, nu

def pendulum_model():
    """ Nonlinear inverse pendulum model. """
    M = 1    # mass of the cart [kg]
    m = 0.1  # mass of the ball [kg]
    g = 9.81 # gravity constant [m/s^2]
    l = 0.8  # length of the rod [m]

    p = SX.sym('p')         # horizontal displacement [m]
    theta = SX.sym('theta') # angle with the vertical [rad]
    v = SX.sym('v')         # horizontal velocity [m/s]
    omega = SX.sym('omega') # angular velocity [rad/s]
    F = SX.sym('F')         # horizontal force [N]

    ode_rhs = vertcat(v,
                      omega,
                      (- l*m*sin(theta)*omega**2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta)**2),
                      (- l*m*cos(theta)*sin(theta)*omega**2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta)**2)))

    nx = 4
    # for IRK
    xdot = SX.sym('xdot', nx, 1)
    z = SX.sym('z',0,1)
    return (Function('pendulum', [vertcat(p, theta, v, omega), F], [ode_rhs], ['x', 'u'], ['xdot']),
            nx, # number of states
            1,  # number of controls
            Function('impl_pendulum', [vertcat(p, theta, v, omega), F, xdot, z], [ode_rhs-xdot],
                                    ['x', 'u','xdot','z'], ['rhs']))
