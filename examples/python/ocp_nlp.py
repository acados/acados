from acados import *
from casadi import *
from numpy import array, diag
from scipy.linalg import block_diag
import matplotlib.pyplot as plt

N = 10
nx = 2
nu = 1

nlp = ocp_nlp({'N':N, 'nx':nx, 'nu':nu})
# Specify initial condition
current_state = array([0.1, 0.1])
nlp.lb[0] = current_state
nlp.ub[0] = current_state
# Weighting matrix
Q = diag([1.0, 1.0])
R = 1e-2
nlp.ls_cost_matrix = N*[block_diag(Q, R)] + [Q]

# The following ODE model comes from Chen1998
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

solver = ocp_nlp_solver('gauss-newton-sqp', nlp)

STATES = [current_state]

for i in range(10):
    current_state = solver.solve(current_state)[1]
    print(current_state)
    STATES.append(current_state)

plt.ion()
plt.plot([x[0] for x in STATES], [x[1] for x in STATES])
plt.axis('equal')

