import matplotlib.pyplot as plt
from numpy import array, diag
from scipy.linalg import block_diag

from acados import ocp_nlp, ocp_nlp_solver
from casadi import SX, Function, vertcat

N = 10
nx = 2
nu = 1

nlp = ocp_nlp({'N': N, 'nx': nx, 'nu': nu})
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
ode_fun = Function('ode_fun', [x, u], [rhs])
step = 0.1
nlp.set_model(ode_fun, step)

solver = ocp_nlp_solver('gauss-newton-sqp', nlp, {'integrator_steps': 2})

STATES = [current_state]
for i in range(50):
    state_traj, control_traj = solver.solve(current_state)
    current_state = state_traj[1]
    STATES += [current_state]
    print(STATES)

plt.ion()
plt.plot([x[0] for x in STATES], [x[1] for x in STATES])
plt.axis('equal')
