from acados_template import AcadosModel
from casadi import SX, vertcat, sin, cos

# Reference for model equations:
# http://users.isr.ist.utl.pt/~jag/publications/08-JETC-RCarona-vcontrol.pdf

def export_robot_model() -> AcadosModel:
    model_name = "unicycle_ode"

    # set up states & controls
    x = SX.sym("x")
    y = SX.sym("y")
    v = SX.sym("x_d")
    theta = SX.sym("theta")
    theta_d = SX.sym("theta_d")

    x = vertcat(x, y, v, theta, theta_d)

    F = SX.sym("F")
    T = SX.sym("T")
    u = vertcat(F, T)

    # xdot
    x_dot = SX.sym("x_dot")
    y_dot = SX.sym("y_dot")
    v_dot = SX.sym("v_dot")
    theta_dot = SX.sym("theta_dot")
    theta_ddot = SX.sym("theta_ddot")

    xdot = vertcat(x_dot, y_dot, v_dot, theta_dot, theta_ddot)

    # algebraic variables
    # z = None

    # parameters
    p = []

    # dynamics
    f_expl = vertcat(v * cos(theta), v * sin(theta), F, theta_d, T)

    f_impl = xdot - f_expl

    model = AcadosModel()

    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    # model.z = z
    model.p = p
    model.name = model_name

    return model
