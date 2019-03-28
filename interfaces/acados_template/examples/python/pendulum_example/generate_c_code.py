from acados_template import *
import acados_template as at
import numpy as np
from ctypes import *
import matplotlib
import matplotlib.pyplot as plt
import scipy.linalg

USE_JSON_DUMP = 1

def export_ode_model():

    model_name = 'pendulum_ode'

    # constants
    M = 1.
    m = 0.1
    g = 9.81
    l = 0.8

    # set up states & controls
    x1      = SX.sym('x1')
    theta   = SX.sym('theta')
    v1      = SX.sym('v1')
    dtheta  = SX.sym('dtheta')
    
    x = vertcat(x1, v1, theta, dtheta)

    # controls
    F = SX.sym('F')
    u = vertcat(F)
    
    # xdot
    x1_dot      = SX.sym('x1_dot')
    theta_dot   = SX.sym('theta_dot')
    v1_dot      = SX.sym('v1_dot')
    dtheta_dot  = SX.sym('dtheta_dot')

    xdot = vertcat(x1_dot, theta_dot, v1_dot, dtheta_dot)
    
    # algebraic variables
    z = []

    # parameters
    p = []
    
    # dynamics     
    denominator = M + m - m*cos(theta)*cos(theta)
    f_expl = vertcat(v1, (-m*l*sin(theta)*dtheta*dtheta + m*g*cos(theta)*sin(theta)+F)/denominator, dtheta, (-m*l*cos(theta)*sin(theta)*dtheta*dtheta + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator))
    
    f_impl = xdot - f_expl
   
    model = acados_dae()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = z
    model.p = p
    model.name = model_name

    return model 

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

VxN = np.zeros((ny, nx))
VxN[0,0] = 1.0
VxN[1,1] = 1.0
VxN[2,2] = 1.0
VxN[3,3] = 1.0

nlp_cost.VxN = VxN

nlp_cost.yref  = np.zeros((ny, 1))
nlp_cost.yrefN = np.zeros((nyN, 1))

# setting bounds
Fmax = 80.0
nlp_con = ra.constraints
nlp_con.lbu = np.array([-Fmax])
nlp_con.ubu = np.array([+Fmax])
nlp_con.x0 = np.array([0.0, 0.0, 3.14, 0.0])
nlp_con.idxbu = np.array([1])

# set constants
const1 = ocp_nlp_constant()
const1.name  = 'PI'
const1.value = 3.1415926535897932
ra.constants = [const1]

# set QP solver
ra.solver_config.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
# ra.solver_config.qp_solver = 'FULL_CONDENSING_QPOASES'
ra.solver_config.hessian_approx = 'GAUSS_NEWTON'
# ra.solver_config.hessian_approx = 'EXACT'
ra.solver_config.integrator_type = 'ERK'
# ra.solver_config.integrator_type = 'IRK'

# set prediction horizon
ra.solver_config.tf = Tf
ra.solver_config.nlp_solver_type = 'SQP'

# set header path
ra.acados_include_path = '/usr/local/include'
ra.acados_lib_path = '/usr/local/lib'

if USE_JSON_DUMP == 1: 
    name_file = 'acados_ocp'
    ocp_nlp = ra
    ocp_nlp.cost = ra.cost.__dict__
    ocp_nlp.constraints = ra.constraints.__dict__
    ocp_nlp.solver_config = ra.solver_config.__dict__
    ocp_nlp.dims = ra.dims.__dict__
    ocp_nlp = ocp_nlp.__dict__
    ocp_nlp = rename_keys(ocp_nlp)
    with open(name_file, 'w') as f:
        json.dump(ocp_nlp, f, default=np_array_to_list)

    ra = ocp_nlp_as_object(ocp_nlp)
    ra.cost = ocp_nlp_as_object(ra.cost)
    ra.constraints = ocp_nlp_as_object(ra.constraints)
    ra.solver_config = ocp_nlp_as_object(ra.solver_config)
    ra.dims = ocp_nlp_as_object(ra.dims)

generate_solver(model, ra)

# make 
os.chdir('c_generated_code')
os.system('make')
os.system('make shared_lib')
os.chdir('..')

acados   = CDLL('c_generated_code/acados_solver_pendulum_ode.so')

acados.acados_create()

nlp_opts = acados.acados_get_nlp_opts()
nlp_dims = acados.acados_get_nlp_dims()
nlp_config = acados.acados_get_nlp_config()
nlp_out = acados.acados_get_nlp_out()
nlp_in = acados.acados_get_nlp_in()

# closed loop simulation TODO(add proper simulation)
Nsim = 100

lb0 = np.ascontiguousarray(np.zeros((5,1)), dtype=np.float64)
ub0 = np.ascontiguousarray(np.zeros((5,1)), dtype=np.float64)
lb0 = cast(lb0.ctypes.data, POINTER(c_double))
ub0 = cast(ub0.ctypes.data, POINTER(c_double))

x0 = np.ascontiguousarray(np.zeros((4,1)), dtype=np.float64)
x0 = cast(x0.ctypes.data, POINTER(c_double))
u0 = np.ascontiguousarray(np.zeros((1,1)), dtype=np.float64)
u0 = cast(u0.ctypes.data, POINTER(c_double))

simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))

for i in range(Nsim):
    acados.acados_solve()

    # get solution
    acados.ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "x", x0);
    acados.ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "u", u0);
    
    for j in range(nx):
        simX[i,j] = x0[j]
    for j in range(nu):
        simU[i,j] = u0[j]
    
    # update initial condition
    acados.ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 1, "x", x0);

    # ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lb", lb0);

    acados.ocp_nlp_constraints_model_set.argtypes = [c_void_p, c_void_p, c_void_p, c_int, c_char_p, POINTER(c_double)]
    field_name = "lbx"
    arg = field_name.encode('utf-8')
    # acados.ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, 0, arg, x0);
    acados.ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, arg, x0);
    field_name = "ubx"
    arg = field_name.encode('utf-8')
    # acados.ocp_nlp_constraints_bounds_set(nlp_config, nlp_dims, nlp_in, 0, arg, x0);
    acados.ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, arg, x0);

# plot results
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

# free memory
acados.acados_free()


