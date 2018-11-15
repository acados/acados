from casadi import *
from generate_c_code_explicit_ode import *
class ode_model():
    def __init__(self):
        self.f_impl_expr = None
        self.f_expl_expr = None
        self.x = None
        self.xdot = None
        self.u = None
        self.z = None
        self.name = None

def export_ode_model():
    # This function generates an implicit ODE / index-1 DAE model,
    # which consists of a CasADi expression f_impl_expr, f_expl_expr
    # that depends on the symbolic CasADi variables x, xdot, u, z,
    # and a model name, which will be used as a prefix for generated C
    # functions later on

    # This model is based on the explicit pendulum model
    # but formulated implicitly to test implicit integrators with it.

    model_name = 'pendulum_ode'

    ## Constants
    M = 1
    m = 0.1
    g = 9.81
    l = 0.8

    ## Set up States & Controls
    x1      = SX.sym('x1')
    theta   = SX.sym('theta')
    v1      = SX.sym('v1')
    dtheta  = SX.sym('dtheta')
    
    x = vertcat(x1, v1, theta, dtheta)

    # Controls
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
    
    ## Dynamics     
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
    model.name = model_name

    # Explicit Model -- Generate C Code
    generate_c_code_explicit_ode( model );

    return 
