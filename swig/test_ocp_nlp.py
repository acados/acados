from acados import *
from casadi import *

nx = 2
nu = 1

nlp = ocp_nlp_in({'N':5, 'nx':nx, 'nu':nu})

# This model comes from Chen1998

x = SX.sym('x', nx)
u = SX.sym('u', nu)

mu = 0.5
rhs = vertcat(x[1] + u*(mu + (1.-mu)*x[0]), x[0] + u*(mu - 4.*(1.-mu)*x[1]))

Sx = SX.sym('Sx', nx, nx)
Su = SX.sym('Su', nx, nu)

vde_x = jtimes(rhs, x, Sx)
vde_u = jacobian(rhs, u) + jtimes(rhs, x, Su)
vde = Function('vde', [x, Sx, Su, u], [rhs, vde_x, vde_u])
vde.generate('vde.c')
nlp.set_model('vde')

solver = ocp_nlp_solver("gauss-newton-sqp", nlp)
print(solver.solve())
