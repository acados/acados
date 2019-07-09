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
import matplotlib.pyplot as plt
from numpy import array, diag, eye, zeros
from scipy.linalg import block_diag

from acados import ocp_nlp
from casadi import SX, Function, vertcat

nx, nu = 2, 1

x = SX.sym('x', nx)
u = SX.sym('u', nu)

ode_fun = Function('ode_fun', [x, u], [vertcat(x[1], u)], ['x', 'u'], ['rhs'])

N = 15

nlp = ocp_nlp(N, nx, nu)
nlp.set_dynamics(ode_fun, {'integrator': 'rk4', 'step': 0.1})

q, r = 1, 1
P = eye(nx)

nlp.set_stage_cost(eye(nx+nu), zeros(nx+nu), diag([q, q, r]))
nlp.set_terminal_cost(eye(nx), zeros(nx), P)

x0 = array([1, 1])
nlp.set_field("lbx", 0, x0)
nlp.set_field("ubx", 0, x0)

nlp.initialize_solver("sqp", {"qp_solver": "qpoases"})

output = nlp.solve()

print("states:")
print(output.states())
print()

print("controls:")
print(output.controls())
print()

# from models import chen_model

# N = 10
# ode_fun, nx, nu = chen_model()
# nlp = ocp_nlp({'N': N, 'nx': nx, 'nu': nu, 'ng': N*[1] + [0]})

# # ODE Model
# step = 0.1
# nlp.model_set(ode_fun, step)

# # Cost function
# Q = diag([1.0, 1.0])
# R = 1e-2
# x = SX.sym('x', nx)
# u = SX.sym('u', nu)
# u_N = SX.sym('u', 0)
# f = ocp_nlp_function(Function('ls_cost', [x, u], [vertcat(x, u)]))
# f_N = ocp_nlp_function(Function('ls_cost_N', [x, u_N], [x]))
# ls_cost = ocp_nlp_ls_cost(N, N*[f]+[f_N])
# ls_cost.ls_cost_matrix = N*[block_diag(Q, R)] + [Q]
# nlp.set_cost(ls_cost)

# # Constraints
# g = ocp_nlp_function(Function('path_constraint', [x, u], [u]))
# g_N = ocp_nlp_function(Function('path_constraintN', [x, u], [SX([])]))
# nlp.set_path_constraints(N*[g] + [g_N])
# for i in range(N):
#     nlp.lg[i] = -0.5
#     nlp.ug[i] = +0.5

# solver = ocp_nlp_solver('sqp', nlp, {'integrator_steps': 2, 'qp_solver':'condensing_qpoases', 'sensitivity_method': 'gauss-newton'})

# # Simulation
# STATES = [array([0.1, 0.1])]
# CONTROLS = []
# for i in range(20):
#     state_traj, control_traj = solver.evaluate(STATES[-1])
#     STATES += [state_traj[1]]
#     CONTROLS += [control_traj[0]]

# plt.ion()
# plt.plot([x[0] for x in STATES], [x[1] for x in STATES])
# plt.axis('equal')
