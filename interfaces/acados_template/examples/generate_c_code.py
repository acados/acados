from jinja2 import Environment, FileSystemLoader
from acados_template import *
import acados_template as at
import numpy as np

def export_ode_model():

    model_name = 'pendulum_ode'

    ## constants
    M = 1.
    m = 0.1
    g = 9.81
    l = 0.8

    ## set up States & Controls
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
    
    ## algebraic variables
    z = []

    # parameters
    p = []
    
    ## dynamics     
    denominator = M + m - m*cos(theta)*cos(theta)
    f_expl = vertcat(v1, (-m*l*sin(theta)*dtheta*dtheta + m*g*cos(theta)*sin(theta)+F)/denominator, dtheta, (-m*l*cos(theta)*sin(theta)*dtheta*dtheta + F*cos(theta)+(M+m)*g*sin(theta))/(l*denominator))
    
    f_impl = xdot - f_expl
   
    model = ode_model()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = z
    model.p = p
    model.name = model_name

    return model 

# setting up loader and environment
acados_path = os.path.dirname(os.path.abspath(at.__file__))
file_loader = FileSystemLoader(acados_path + '/c_templates')
env = Environment(loader = file_loader)
template = env.get_template('main.in.c')

# create render arguments
ra = ocp_nlp_render_arguments()

# export model 
model = export_ode_model()

# set model_name 
ra.model_name = model.name

# set ocp_nlp_dimensions
nlp_dims = ra.dims
nlp_dims.nx = model.x.size()[0]
nlp_dims.nbx = 0
nlp_dims.nbu = model.u.size()[0]
nlp_dims.nu = model.u.size()[0]
nlp_dims.N  = 100

# set weighting matrices
nlp_cost = ra.cost
Q = np.eye(4)
Q[0,0] = 1e3
Q[1,1] = 1e-2
Q[2,2] = 1e3
Q[3,3] = 1e-2

R = np.eye(1)
R[0,0] = 1e-2

nlp_cost.Q = Q 
nlp_cost.R = R

# setting bounds
Fmax = 80.0
nlp_con = ra.constraints
nlp_con.lbu = np.array([-Fmax])
nlp_con.ubu = np.array([+Fmax])
nlp_con.x0 = np.array([0.0, 0.0, 3.14, 0.0])

# set constants
const1 = ocp_nlp_constant()
const1.name  = 'PI'
const1.value = 3.1415926535897932
ra.constants = [const1]

# set QP solver
# ra.solver_config.qp_solver = 'PARTIAL_CONDENSING_HPIPM'
ra.solver_config.qp_solver = 'FULL_CONDENSING_QPOASES'
ra.solver_config.hessian_approx = 'GAUSS_NEWTON'
# ra.solver_config.hessian_approx = 'EXACT'
ra.solver_config.integrator_type = 'ERK'
# ra.solver_config.integrator_type = 'IRK'

# set prediction horizon
ra.solver_config.tf = 1.0

# explicit model -- generate C code
generate_c_code_explicit_ode(model);

# implicit model -- generate C code
opts = dict(generate_hess=1)
generate_c_code_implicit_ode(model, opts);

# set header path
ra.acados_include_path = '~/.local/include'
ra.acados_lib_path = '~/.local/lib'

# check render arguments
check_ra(ra)

# render source template
output = template.render(ra=ra)
# output file
out_file = open('./c_generated_code/main_' + model.name + '.c', 'w+')
out_file.write(output)

# render header template
template = env.get_template('model.in.h')
output = template.render(ra=ra)
# output file
out_file = open('./c_generated_code/' + model.name + '_model/' + model.name + '_model.h', 'w+')
out_file.write(output)

# render Makefile template
template = env.get_template('Makefile.in')
output = template.render(ra=ra)
# output file
out_file = open('./c_generated_code/Makefile', 'w+')
out_file.write(output)
