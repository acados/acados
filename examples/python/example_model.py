from casadi import SX, Function, vertcat

def example_model():
    """ The following ODE model comes from Chen1998. """
    nx = 2
    nu = 1
    x = SX.sym('x', nx)
    u = SX.sym('u', nu)
    mu = 0.5
    rhs = vertcat(x[1] + u*(mu + (1.-mu)*x[0]), x[0] + u*(mu - 4.*(1.-mu)*x[1]))
    return nx, nu, Function('ode_fun', [x, u], [rhs])
