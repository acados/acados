from casadi import cos, Function, sin, SX, vertcat

def chen_model():
    """ The following ODE model comes from Chen1998. """
    nx, nu = (2, 1)
    x = SX.sym('x', nx)
    u = SX.sym('u', nu)
    mu = 0.5
    rhs = vertcat(x[1] + u*(mu + (1.-mu)*x[0]), x[0] + u*(mu - 4.*(1.-mu)*x[1]))
    return Function('chen', [x, u], [rhs], ['x', 'u'], ['xdot']), nx, nu

def pendulum_model():
    """ Nonlinear inverse pendulum model. """
    M = 1    # mass of the cart [kg]
    m = 0.1  # mass of the ball [kg]
    g = 9.81 # gravity constant [m/s^2]
    l = 0.8  # length of the rod [m]

    p = SX.sym('p')         # horizontal displacement [m]
    theta = SX.sym('theta') # angle with the vertical [rad]
    v = SX.sym('v')         # horizontal velocity [m/s]
    omega = SX.sym('omega') # angular velocity [rad/s]
    F = SX.sym('F')         # horizontal force [N]

    ode_rhs = vertcat(v,
                      omega,
                      (- l*m*sin(theta)*omega**2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta)**2),
                      (- l*m*cos(theta)*sin(theta)*omega**2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta)**2)))

    nx = 4
    # for IRK
    xdot = SX.sym('xdot', nx, 1)
    z = SX.sym('z',0,1)
    return (Function('pendulum', [vertcat(p, theta, v, omega), F], [ode_rhs], ['x', 'u'], ['xdot']),
            nx, # number of states
            1,  # number of controls
            Function('impl_pendulum', [vertcat(p, theta, v, omega), F, xdot, z], [ode_rhs-xdot],
                                    ['x', 'u','xdot','z'], ['rhs']))
