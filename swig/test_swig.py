from acados import *
from numpy import array, ones, zeros, diag

# qp_in = ocp_qp_in({'N':5, 'nx':2, 'nu':1})
# x0 = array([1,1])

# # specify OCP
# qp_in.A = array([[0,1],[0,0]])
# qp_in.B = array([[0],[1]])
# qp_in.Q = diag([1,1])
# qp_in.R = diag([1])

# # specify initial condition
# x0 = array([1,1])
# qp_in.lb[0] = x0
# qp_in.ub[0] = x0

# solver = ocp_qp_solver("condensing_qpoases", qp_in)
# result = solver.solve()
# print(result)

from casadi import *

# This model comes from Chen1998

nx = 2
nu = 1

x = SX.sym('x', nx)
u = SX.sym('u', nu)

mu = 0.5
rhs = vertcat(x[1] + u*(mu + (1.-mu)*x[0]), x[0] + u*(mu - 4.*(1.-mu)*x[1]))


print(rhs)

# ode = Function('ode', [x, u], [rhs])
# jac_ode = Function('jac_ode', [x, u], [rhs, jacobian(rhs, x)])

# Sx = SX.sym('Sx', nx, nx)
# Su = SX.sym('Su', nx, nu)

# vde_x = jtimes(rhs, x, Sx)
# vde_u = jacobian(rhs, u) + jtimes(rhs, x, Su)

# vde = Function('vde', [x, Sx, Su, u], [rhs, vde_x, vde_u])

# ode.generate('ode.c')
# jac_ode.generate('jac_ode.c')
# vde.generate('vde.c')

# nlp = ocp_nlp_in()
# a = nlp.set_model('ode', 5)
# print(a)
