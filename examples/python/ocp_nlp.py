import matplotlib.pyplot as plt
from numpy import array, diag
from scipy.linalg import block_diag

from acados import ocp_nlp, ocp_nlp_solver
from models import chen_model

N = 10
ode_fun, nx, nu = chen_model()
nlp = ocp_nlp({'N': N, 'nx': nx, 'nu': nu})

# ODE Model
step = 0.1
nlp.set_model(ode_fun, step)

# Cost function
Q = diag([1.0, 1.0])
R = 1e-2
nlp.ls_cost_matrix = N*[block_diag(Q, R)] + [Q]

solver = ocp_nlp_solver('gauss-newton-sqp', nlp, {'integrator_steps': 2})

# Simulation
STATES = [array([0.1, 0.1])]
for i in range(50):
    state_traj, control_traj = solver.evaluate(STATES[-1])
    STATES += [state_traj[1]]

plt.ion()
plt.plot([x[0] for x in STATES], [x[1] for x in STATES])
plt.axis('equal')
