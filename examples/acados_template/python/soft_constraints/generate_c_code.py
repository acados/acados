from acados_template import *
import acados_template as at
from export_ode_model import *
import numpy as np
import scipy.linalg
from ctypes import *

FORMULATION = 1 # 0 for linear soft bounds, 
            # 1 for equivalent nonlinear soft constraint

def export_nonlinear_constraint():

    con_name = 'nl_con'

    # set up states & controls
    x1      = SX.sym('x1')
    theta   = SX.sym('theta')
    v1      = SX.sym('v1')
    dtheta  = SX.sym('dtheta')
    
    x = vertcat(x1, v1, theta, dtheta)

    # controls
    F = SX.sym('F')
    u = vertcat(F)

    # voltage sphere
    constraint = acados_constraint()

    constraint.expr = u
    constraint.x = x
    constraint.u = u
    constraint.nc = 1
    constraint.name = con_name

    return constraint

# create render arguments
ocp = acados_ocp_nlp()

# export model 
model = export_ode_model()

# set model_name 
ocp.model_name = model.name

Tf = 1.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
ny_e = nx
N = 100

# set ocp_nlp_dimensions
nlp_dims     = ocp.dims
nlp_dims.nx  = nx 
nlp_dims.ny  = ny 
nlp_dims.ny_e = ny_e 
nlp_dims.nbx = 0
nlp_dims.ns = nu 
nlp_dims.nu  = model.u.size()[0]
nlp_dims.N   = N

if FORMULATION == 1:
    nlp_dims.nh   = model.u.size()[0]
    nlp_dims.nsh  = model.u.size()[0] 
    nlp_dims.nbu  = 0 
    nlp_dims.nsbu = 0 
else:
    nlp_dims.nh   = 0
    nlp_dims.nsh  = 0
    nlp_dims.nbu  = model.u.size()[0]
    nlp_dims.nsbu = model.u.size()[0]

# set weighting matrices
nlp_cost = ocp.cost
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

nlp_cost.W_e = Q 

Vx_e = np.zeros((ny_e, nx))
Vx_e[0,0] = 1.0
Vx_e[1,1] = 1.0
Vx_e[2,2] = 1.0
Vx_e[3,3] = 1.0

nlp_cost.Vx_e = Vx_e

nlp_cost.yref  = np.zeros((ny, ))
nlp_cost.yref_e = np.zeros((ny_e, ))

nlp_cost.zl = 50*np.ones((1, ))
nlp_cost.Zl = 0*np.ones((1, 1))
nlp_cost.zu = 50*np.ones((1, ))
nlp_cost.Zu = 0*np.ones((1, 1))

# setting bounds
Fmax = 80.0
nlp_con = ocp.constraints
nlp_con.x0 = np.array([0.0, 0.0, 3.14, 0.0])


constraint = export_nonlinear_constraint()
if FORMULATION == 1:
    nlp_con.lh = np.array([-Fmax])
    nlp_con.uh = np.array([+Fmax])
    nlp_con.lsh = 0*np.array([-Fmax])
    nlp_con.ush = 0*np.array([+Fmax])
    nlp_con.idxsh = np.array([0])
else:
    nlp_con.lbu = np.array([-Fmax])
    nlp_con.ubu = np.array([+Fmax])
    nlp_con.lsbu = 0*np.array([-Fmax])
    nlp_con.usbu = 0*np.array([+Fmax])
    nlp_con.idxbu = np.array([0])
    nlp_con.idxsbu = np.array([0])

# set QP solver
ocp.solver_config.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
ocp.solver_config.hessian_approx = 'GAUSS_NEWTON'
ocp.solver_config.integrator_type = 'ERK'

# set prediction horizon
ocp.solver_config.tf = Tf
ocp.solver_config.nlp_solver_type = 'SQP'

# set header path
ocp.acados_include_path  = '/usr/local/include'
ocp.acados_lib_path      = '/usr/local/lib'

# json_layout = acados_ocp2json_layout(ocp)
# with open('acados_layout.json', 'w') as f:
#     json.dump(json_layout, f, default=np_array_to_list)
# exit()

if FORMULATION == 1:
    ocp.con_h_name = 'nl_con'
    acados_solver = generate_solver(model, ocp, con_h=constraint, json_file = 'acados_ocp.json')
else:
    acados_solver = generate_solver(model, ocp, json_file = 'acados_ocp.json')

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
plt.step(t, simU, color='r')
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

