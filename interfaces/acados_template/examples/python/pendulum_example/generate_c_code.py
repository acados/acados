from acados_template import *
import acados_template as at
from export_ode_model import *
import numpy as np
import scipy.linalg
from ctypes import *

# create render arguments
ra = acados_ocp_nlp()

# export model 
model = export_ode_model()

# set model_name 
ra.model_name = model.name

Tf = 1.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
nyN = nx
N = 10

# set ocp_nlp_dimensions
nlp_dims     = ra.dims
nlp_dims.nx  = nx 
nlp_dims.ny  = ny 
nlp_dims.nyN = nyN 
nlp_dims.nbx = 0
nlp_dims.nbu = nu 
nlp_dims.nu  = model.u.size()[0]
nlp_dims.N   = N

# set weighting matrices
nlp_cost = ra.cost
Q = np.eye(4)
Q[0,0] = 1e3
Q[1,1] = 1e-2
Q[2,2] = 1e3
Q[3,3] = 1e-2

R = np.eye(1)
R[0,0] = 1e-2

nlp_cost.W = scipy.linalg.block_diag(Q, R) 

Vx = np.zeros((ny, nx))
Vx[0,0] = 1.0
Vx[1,1] = 1.0
Vx[2,2] = 1.0
Vx[3,3] = 1.0

nlp_cost.Vx = Vx

Vu = np.zeros((ny, nu))
Vu[4,0] = 1.0
nlp_cost.Vu = Vu

nlp_cost.WN = Q 

VxN = np.zeros((nyN, nx))
VxN[0,0] = 1.0
VxN[1,1] = 1.0
VxN[2,2] = 1.0
VxN[3,3] = 1.0

nlp_cost.VxN = VxN

nlp_cost.yref  = np.zeros((ny, ))
nlp_cost.yrefN = np.zeros((nyN, ))

# setting bounds
Fmax = 80.0
nlp_con = ra.constraints
nlp_con.lbu = np.array([-Fmax])
nlp_con.ubu = np.array([+Fmax])
nlp_con.x0 = np.array([0.0, 0.0, 3.14, 0.0])
nlp_con.idxbu = np.array([1])

# set constants
ra.constants['PI'] = 3.1415926535897932

# set QP solver
ra.solver_config.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
ra.solver_config.hessian_approx = 'GAUSS_NEWTON'
ra.solver_config.integrator_type = 'ERK'

# set prediction horizon
ra.solver_config.tf = Tf
ra.solver_config.nlp_solver_type = 'SQP'

# set header path
ra.acados_include_path  = '/usr/local/include'
ra.acados_lib_path      = '/usr/local/lib'

# json_layout = acados_ocp2json_layout(ra)
# with open('acados_layout.json', 'w') as f:
#     json.dump(json_layout, f, default=np_array_to_list)

# import pdb; pdb.set_trace()
acados_solver = generate_solver(model, ra, json_file = 'acados_ocp.json')

Nsim = 100

simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))

for i in range(Nsim):
    status = acados_solver.solve()

    # get solution
    x0 = acados_solver.get(0, "x")
    u0 = acados_solver.get(0, "u")
    
    for j in range(nx):
        simX[i,j] = x0[j]

    for j in range(nu):
        simU[i,j] = u0[j]
    
    # update initial condition
    x0 = acados_solver.get(1, "x")

    acados_solver.set(0, "lbx", x0)
    acados_solver.set(0, "ubx", x0)

# plot results
import matplotlib
import matplotlib.pyplot as plt
t = np.linspace(0.0, Tf/N, Nsim)
plt.subplot(2, 1, 1)
plt.step(t, simU, 'r')
plt.title('closed-loop simulation')
plt.ylabel('u')
plt.xlabel('t')
plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(t, simX[:,2])
plt.ylabel('theta')
plt.xlabel('t')
plt.grid(True)
plt.show()

