import matplotlib.pyplot as plt
from numpy import array, diag
from scipy.linalg import block_diag

from casadi import SX, Function, vertcat
from acados import ocp_nlp_function, ocp_nlp_ls_cost, ocp_nlp, ocp_nlp_solver
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
x = SX.sym('x',nx)
u = SX.sym('u',nu)
uN = SX.sym('u',0)
F = ocp_nlp_function(Function('ls_cost', [x,u], [vertcat(x,u)]))
FN = ocp_nlp_function(Function('ls_costN', [x,uN], [x]))
ls_cost = ocp_nlp_ls_cost(N, N*[F]+[FN])
ls_cost.ls_cost_matrix = N*[block_diag(Q, R)] + [Q]
nlp.set_cost(ls_cost)

# Constraints
g = SX([]);
G = ocp_nlp_function(Function('path_constraint', [x,u], [g]))
GN = ocp_nlp_function(Function('path_constraintN', [x,uN], [g]))
path_constraints = N*[G] + [GN]
nlp.set_path_constraints(path_constraints)

solver = ocp_nlp_solver('sqp', nlp, {'integrator_steps': 2, 'qp_solver': 'condensing_qpoases', 'sensitivity_method': 'gauss-newton'})

# Simulation
STATES = [array([0.1, 0.1])]
for i in range(50):
    state_traj, control_traj = solver.evaluate(STATES[-1])
    STATES += [state_traj[1]]

plt.ion()
plt.plot([x[0] for x in STATES], [x[1] for x in STATES])
plt.axis('equal')
